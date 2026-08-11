#include "swiftig/alignment.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace swiftig {
namespace {

bool is_canonical_base(char base) {
    return base == 'A' || base == 'C' || base == 'G' || base == 'T';
}

bool is_canonical_match(char query, char reference) {
    return query == reference && is_canonical_base(query);
}

void append_cigar_op(std::string& cigar, std::size_t count, char op) {
    if (count == 0) return;
    cigar += std::to_string(count);
    cigar.push_back(op);
}

std::string make_cigar(
    const std::vector<char>& operations,
    std::size_t query_start,
    std::size_t query_end,
    std::size_t query_size,
    std::size_t reference_start,
    std::size_t reference_end,
    std::size_t reference_size) {
    std::string cigar;
    append_cigar_op(cigar, query_start, 'S');
    append_cigar_op(cigar, reference_start, 'N');
    std::size_t run = 0;
    char previous = '\0';
    for (const char op : operations) {
        if (op == previous) {
            ++run;
        } else {
            append_cigar_op(cigar, run, previous);
            previous = op;
            run = 1;
        }
    }
    append_cigar_op(cigar, run, previous);
    append_cigar_op(cigar, query_size - query_end, 'S');
    append_cigar_op(cigar, reference_size - reference_end, 'N');
    return cigar;
}

}  // namespace

Alignment local_align_affine(
    const std::string& query,
    const std::string& reference,
    const Scoring& scoring,
    int estimated_diagonal,
    int band_width) {
    Alignment result;
    if (query.empty() || reference.empty()) return result;
    const bool banded = band_width >= 0;
    const std::size_t trace_columns = banded
        ? static_cast<std::size_t>(2 * band_width + 1) : reference.size() + 1;
    if (trace_columns == 0 || query.size() >
        (std::numeric_limits<std::size_t>::max() / trace_columns) - 1) {
        return result;
    }
    const std::size_t rows = query.size() + 1;
    const std::size_t columns = reference.size() + 1;
    std::vector<std::uint8_t> trace(rows * trace_columns, 0);
    std::vector<int> previous(columns, 0), current(columns, 0);
    constexpr int neg_inf = std::numeric_limits<int>::min() / 4;
    std::vector<int> insertion_previous(columns, neg_inf), insertion_current(columns, neg_inf);

    const auto trace_index = [&](std::size_t i, std::size_t j) -> std::optional<std::size_t> {
        if (!banded) return i * trace_columns + j;
        const long theoretical_low = static_cast<long>(i) -
            static_cast<long>(estimated_diagonal) - band_width;
        const long offset = static_cast<long>(j) - theoretical_low;
        if (offset < 0 || offset >= static_cast<long>(trace_columns)) return std::nullopt;
        return i * trace_columns + static_cast<std::size_t>(offset);
    };

    int best_score = 0;
    std::size_t best_i = 0;
    std::size_t best_j = 0;
    std::size_t previous_low = 1;
    std::size_t previous_high = reference.size();

    for (std::size_t i = 1; i < rows; ++i) {
        std::size_t j_low = 1;
        std::size_t j_high = reference.size();
        if (band_width >= 0) {
            const long center = static_cast<long>(i) - static_cast<long>(estimated_diagonal);
            const long low = center - band_width;
            const long high = center + band_width;
            if (high < 1 || low > static_cast<long>(reference.size())) {
                previous_low = 1;
                previous_high = 0;
                continue;
            }
            j_low = static_cast<std::size_t>(std::max<long>(1, low));
            j_high = static_cast<std::size_t>(std::min<long>(static_cast<long>(reference.size()), high));
            if (j_low > j_high) {
                previous.swap(current);
                insertion_previous.swap(insertion_current);
                continue;
            }
        }
        int deletion = neg_inf;
        for (std::size_t j = j_low; j <= j_high; ++j) {
            const bool previous_has_j = j >= previous_low && j <= previous_high;
            const bool previous_has_diagonal = j > 0 && j - 1 >= previous_low && j - 1 <= previous_high;
            const int previous_score = previous_has_j ? previous[j] : 0;
            const int previous_insertion = previous_has_j ? insertion_previous[j] : neg_inf;
            const int diagonal_score = previous_has_diagonal ? previous[j - 1] : 0;
            const int insertion_open = previous_score + scoring.gap_open;
            const int insertion_extend = previous_insertion + scoring.gap_extend;
            const bool insertion_is_extension = insertion_extend > insertion_open;
            const int insertion = insertion_is_extension ? insertion_extend : insertion_open;
            insertion_current[j] = insertion;

            const int deletion_open = (j == j_low ? 0 : current[j - 1]) + scoring.gap_open;
            const int deletion_extend = deletion + scoring.gap_extend;
            const bool deletion_is_extension = deletion_extend > deletion_open;
            deletion = deletion_is_extension ? deletion_extend : deletion_open;

            const char query_base = query[i - 1];
            const char reference_base = reference[j - 1];
            const int substitution_delta = is_canonical_match(query_base, reference_base)
                ? scoring.match
                : (!is_canonical_base(query_base) || !is_canonical_base(reference_base)
                    ? 0 : scoring.mismatch);
            const int substitution = diagonal_score + substitution_delta;
            int score = 0;
            std::uint8_t direction = 0;
            if (substitution > score) {
                score = substitution;
                direction = 1;
            }
            if (insertion > score) {
                score = insertion;
                direction = 2;
            }
            if (deletion > score) {
                score = deletion;
                direction = 3;
            }
            current[j] = score;
            if (const auto index = trace_index(i, j)) {
                trace[*index] = static_cast<std::uint8_t>(
                    direction | (insertion_is_extension ? 4U : 0U) |
                    (deletion_is_extension ? 8U : 0U));
            }
            if (score > best_score) {
                best_score = score;
                best_i = i;
                best_j = j;
            }
        }
        previous.swap(current);
        insertion_previous.swap(insertion_current);
        previous_low = j_low;
        previous_high = j_high;
    }

    if (best_score <= 0) return result;
    enum class State { Match, Insertion, Deletion };
    State state = State::Match;
    std::size_t i = best_i;
    std::size_t j = best_j;
    std::string aligned_query;
    std::string aligned_reference;
    std::vector<char> operations;
    aligned_query.reserve(query.size());
    aligned_reference.reserve(reference.size());

    while (i > 0 || j > 0) {
        const auto index = trace_index(i, j);
        if (!index) break;
        const auto code = trace[*index];
        if (state == State::Match) {
            const auto direction = static_cast<std::uint8_t>(code & 3U);
            if (direction == 0) break;
            if (direction == 2) {
                state = State::Insertion;
                continue;
            }
            if (direction == 3) {
                state = State::Deletion;
                continue;
            }
            aligned_query.push_back(query[i - 1]);
            aligned_reference.push_back(reference[j - 1]);
            operations.push_back('M');
            if (is_canonical_match(query[i - 1], reference[j - 1])) ++result.matches;
            else ++result.mismatches;
            --i;
            --j;
        } else if (state == State::Insertion) {
            if (i == 0) break;
            const bool extends = (code & 4U) != 0;
            aligned_query.push_back(query[i - 1]);
            aligned_reference.push_back('-');
            operations.push_back('I');
            ++result.insertions;
            --i;
            if (!extends) state = State::Match;
        } else {
            if (j == 0) break;
            const bool extends = (code & 8U) != 0;
            aligned_query.push_back('-');
            aligned_reference.push_back(reference[j - 1]);
            operations.push_back('D');
            ++result.deletions;
            --j;
            if (!extends) state = State::Match;
        }
    }

    std::reverse(aligned_query.begin(), aligned_query.end());
    std::reverse(aligned_reference.begin(), aligned_reference.end());
    std::reverse(operations.begin(), operations.end());
    result.score = best_score;
    result.query_start = i;
    result.query_end = best_i;
    result.reference_start = j;
    result.reference_end = best_j;
    result.aligned_query = std::move(aligned_query);
    result.aligned_reference = std::move(aligned_reference);
    result.cigar = make_cigar(
        operations, result.query_start, result.query_end, query.size(),
        result.reference_start, result.reference_end, reference.size());
    return result;
}

void refresh_airr_cigar(
    Alignment& alignment,
    std::size_t query_size,
    std::size_t reference_size) {
    std::vector<char> operations;
    operations.reserve(alignment.aligned_query.size());
    for (std::size_t i = 0; i < alignment.aligned_query.size(); ++i) {
        const char q = alignment.aligned_query[i];
        const char r = alignment.aligned_reference[i];
        if (q == '-' && r != '-') operations.push_back('D');
        else if (q != '-' && r == '-') operations.push_back('I');
        else if (q != '-' && r != '-') operations.push_back('M');
    }
    alignment.cigar = make_cigar(
        operations, alignment.query_start, alignment.query_end, query_size,
        alignment.reference_start, alignment.reference_end, reference_size);
}

std::optional<std::size_t> map_reference_to_query(
    const Alignment& alignment,
    std::size_t reference_position) {
    if (!alignment.valid() || reference_position < alignment.reference_start ||
        reference_position >= alignment.reference_end) return std::nullopt;
    std::size_t query_position = alignment.query_start;
    std::size_t current_reference = alignment.reference_start;
    for (std::size_t column = 0; column < alignment.aligned_query.size(); ++column) {
        const char q = alignment.aligned_query[column];
        const char r = alignment.aligned_reference[column];
        if (r != '-' && current_reference == reference_position) {
            if (q == '-') return std::nullopt;
            return query_position;
        }
        if (q != '-') ++query_position;
        if (r != '-') ++current_reference;
    }
    return std::nullopt;
}

}  // namespace swiftig
