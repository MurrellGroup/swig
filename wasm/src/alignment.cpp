#include "swiftig/alignment.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_set>
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

std::string extended_geometry_key(
    const std::string& query,
    const std::string& reference,
    const Alignment& alignment) {
    const auto left = std::min(alignment.query_start, alignment.reference_start);
    const auto right = std::min(
        query.size() - alignment.query_end,
        reference.size() - alignment.reference_end);
    std::string key;
    key.reserve(
        2 * (left + alignment.aligned_query.size() + right) + 1);
    key.append(query, alignment.query_start - left, left);
    key += alignment.aligned_query;
    key.append(query, alignment.query_end, right);
    key.push_back('\n');
    key.append(reference, alignment.reference_start - left, left);
    key += alignment.aligned_reference;
    key.append(reference, alignment.reference_end, right);
    return key;
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

std::vector<Alignment> local_align_affine_paths(
    const std::string& query,
    const std::string& reference,
    const Scoring& scoring,
    int estimated_diagonal,
    int band_width,
    std::size_t maximum_paths,
    int score_tolerance,
    std::size_t maximum_trace_states,
    TracebackSearchStats* search_stats) {
    std::vector<Alignment> results;
    if (search_stats) *search_stats = {};
    if (maximum_paths == 0 || query.empty() || reference.empty()) return results;
    if (maximum_paths == 1) {
        auto best = local_align_affine(
            query, reference, scoring, estimated_diagonal, band_width);
        if (best.valid()) results.push_back(std::move(best));
        if (search_stats) search_stats->distinct_paths = results.size();
        return results;
    }

    const bool banded = band_width >= 0;
    const std::size_t trace_columns = banded
        ? static_cast<std::size_t>(2 * band_width + 1) : reference.size() + 1;
    if (trace_columns == 0 || query.size() >
        (std::numeric_limits<std::size_t>::max() / trace_columns) - 1) {
        return results;
    }
    const std::size_t rows = query.size() + 1;
    const auto cells = rows * trace_columns;
    constexpr int neg_inf = std::numeric_limits<int>::min() / 4;
    std::vector<int> match(cells, 0);
    std::vector<int> insertion(cells, neg_inf);
    std::vector<int> deletion(cells, neg_inf);
    std::vector<std::uint8_t> trace(cells, 0);

    const auto trace_index = [&](std::size_t i, std::size_t j) -> std::optional<std::size_t> {
        if (!banded) return i * trace_columns + j;
        const long theoretical_low = static_cast<long>(i) -
            static_cast<long>(estimated_diagonal) - band_width;
        const long offset = static_cast<long>(j) - theoretical_low;
        if (offset < 0 || offset >= static_cast<long>(trace_columns)) return std::nullopt;
        return i * trace_columns + static_cast<std::size_t>(offset);
    };
    const auto match_score = [&](std::size_t i, std::size_t j) {
        if (const auto index = trace_index(i, j)) return match[*index];
        return 0;
    };
    const auto insertion_score = [&](std::size_t i, std::size_t j) {
        if (const auto index = trace_index(i, j)) return insertion[*index];
        return neg_inf;
    };
    const auto deletion_score = [&](std::size_t i, std::size_t j) {
        if (const auto index = trace_index(i, j)) return deletion[*index];
        return neg_inf;
    };
    const auto substitution_delta = [&](std::size_t i, std::size_t j) {
        const char query_base = query[i - 1];
        const char reference_base = reference[j - 1];
        return is_canonical_match(query_base, reference_base)
            ? scoring.match
            : (!is_canonical_base(query_base) || !is_canonical_base(reference_base)
                ? 0 : scoring.mismatch);
    };

    int best_score = 0;
    std::size_t best_i = 0;
    std::size_t best_j = 0;
    for (std::size_t i = 1; i < rows; ++i) {
        std::size_t j_low = 1;
        std::size_t j_high = reference.size();
        if (banded) {
            const long center = static_cast<long>(i) - static_cast<long>(estimated_diagonal);
            const long low = center - band_width;
            const long high = center + band_width;
            if (high < 1 || low > static_cast<long>(reference.size())) continue;
            j_low = static_cast<std::size_t>(std::max<long>(1, low));
            j_high = static_cast<std::size_t>(std::min<long>(
                static_cast<long>(reference.size()), high));
            if (j_low > j_high) continue;
        }
        for (std::size_t j = j_low; j <= j_high; ++j) {
            const auto index = *trace_index(i, j);
            const int insertion_open = match_score(i - 1, j) + scoring.gap_open;
            const int insertion_extend = insertion_score(i - 1, j) + scoring.gap_extend;
            const bool insertion_is_extension = insertion_extend > insertion_open;
            insertion[index] = insertion_is_extension ? insertion_extend : insertion_open;

            const int deletion_open = match_score(i, j - 1) + scoring.gap_open;
            const int deletion_extend = deletion_score(i, j - 1) + scoring.gap_extend;
            const bool deletion_is_extension = deletion_extend > deletion_open;
            deletion[index] = deletion_is_extension ? deletion_extend : deletion_open;

            const int substitution = match_score(i - 1, j - 1) + substitution_delta(i, j);
            int score = 0;
            std::uint8_t direction = 0;
            if (substitution > score) {
                score = substitution;
                direction = 1;
            }
            if (insertion[index] > score) {
                score = insertion[index];
                direction = 2;
            }
            if (deletion[index] > score) {
                score = deletion[index];
                direction = 3;
            }
            match[index] = score;
            trace[index] = static_cast<std::uint8_t>(
                direction | (insertion_is_extension ? 4U : 0U) |
                (deletion_is_extension ? 8U : 0U));
            if (score > best_score) {
                best_score = score;
                best_i = i;
                best_j = j;
            }
        }
    }
    if (best_score <= 0) return results;

    enum class TraceState : std::uint8_t { Match, Insertion, Deletion };
    const auto deterministic_best = [&]() {
        Alignment result;
        TraceState state = TraceState::Match;
        std::size_t i = best_i;
        std::size_t j = best_j;
        std::string aligned_query;
        std::string aligned_reference;
        std::vector<char> operations;
        while (i > 0 || j > 0) {
            const auto index = trace_index(i, j);
            if (!index) break;
            const auto code = trace[*index];
            if (state == TraceState::Match) {
                const auto direction = static_cast<std::uint8_t>(code & 3U);
                if (direction == 0) break;
                if (direction == 2) {
                    state = TraceState::Insertion;
                    continue;
                }
                if (direction == 3) {
                    state = TraceState::Deletion;
                    continue;
                }
                aligned_query.push_back(query[i - 1]);
                aligned_reference.push_back(reference[j - 1]);
                operations.push_back('M');
                if (is_canonical_match(query[i - 1], reference[j - 1])) ++result.matches;
                else ++result.mismatches;
                --i;
                --j;
            } else if (state == TraceState::Insertion) {
                if (i == 0) break;
                const bool extends = (code & 4U) != 0;
                aligned_query.push_back(query[i - 1]);
                aligned_reference.push_back('-');
                operations.push_back('I');
                ++result.insertions;
                --i;
                if (!extends) state = TraceState::Match;
            } else {
                if (j == 0) break;
                const bool extends = (code & 8U) != 0;
                aligned_query.push_back('-');
                aligned_reference.push_back(reference[j - 1]);
                operations.push_back('D');
                ++result.deletions;
                --j;
                if (!extends) state = TraceState::Match;
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
    }();
    results.push_back(deterministic_best);
    std::unordered_set<std::string> geometries;
    geometries.insert(extended_geometry_key(query, reference, deterministic_best));

    struct PartialTrace {
        TraceState state = TraceState::Match;
        std::size_t i = 0;
        std::size_t j = 0;
        int loss = 0;
        std::uint64_t serial = 0;
        std::string aligned_query_reverse;
        std::string aligned_reference_reverse;
        std::vector<char> operations_reverse;
    };
    struct WorseTrace {
        bool operator()(const PartialTrace& left, const PartialTrace& right) const noexcept {
            if (left.loss != right.loss) return left.loss > right.loss;
            return left.serial > right.serial;
        }
    };
    std::priority_queue<PartialTrace, std::vector<PartialTrace>, WorseTrace> frontier;
    std::uint64_t serial = 0;
    frontier.push(PartialTrace{TraceState::Match, best_i, best_j, 0, serial++, {}, {}, {}});
    std::size_t expanded = 0;
    score_tolerance = std::max(score_tolerance, 0);

    auto push_if_allowed = [&](PartialTrace candidate, int delta) {
        if (delta < 0) delta = 0;
        candidate.loss += delta;
        if (candidate.loss > score_tolerance) return;
        candidate.serial = serial++;
        frontier.push(std::move(candidate));
    };

    while (!frontier.empty() && results.size() < maximum_paths &&
           expanded < maximum_trace_states) {
        auto current = frontier.top();
        frontier.pop();
        ++expanded;

        if (current.state == TraceState::Match) {
            const int current_score = match_score(current.i, current.j);
            if (current_score <= 0 || current.i == 0 || current.j == 0) {
                if (current.aligned_query_reverse.empty()) continue;
                Alignment alignment;
                alignment.score = best_score - current.loss;
                alignment.query_start = current.i;
                alignment.query_end = best_i;
                alignment.reference_start = current.j;
                alignment.reference_end = best_j;
                alignment.aligned_query = std::move(current.aligned_query_reverse);
                alignment.aligned_reference = std::move(current.aligned_reference_reverse);
                std::reverse(alignment.aligned_query.begin(), alignment.aligned_query.end());
                std::reverse(alignment.aligned_reference.begin(), alignment.aligned_reference.end());
                auto operations = std::move(current.operations_reverse);
                std::reverse(operations.begin(), operations.end());
                for (std::size_t column = 0; column < alignment.aligned_query.size(); ++column) {
                    const char q = alignment.aligned_query[column];
                    const char r = alignment.aligned_reference[column];
                    if (q == '-' && r != '-') ++alignment.deletions;
                    else if (q != '-' && r == '-') ++alignment.insertions;
                    else if (is_canonical_match(q, r)) ++alignment.matches;
                    else ++alignment.mismatches;
                }
                alignment.cigar = make_cigar(
                    operations, alignment.query_start, alignment.query_end, query.size(),
                    alignment.reference_start, alignment.reference_end, reference.size());
                const auto key = extended_geometry_key(query, reference, alignment);
                if (geometries.insert(key).second) results.push_back(std::move(alignment));
                continue;
            }

            const int diagonal = match_score(current.i - 1, current.j - 1) +
                substitution_delta(current.i, current.j);
            if (diagonal > 0) {
                auto next = current;
                next.state = TraceState::Match;
                next.aligned_query_reverse.push_back(query[current.i - 1]);
                next.aligned_reference_reverse.push_back(reference[current.j - 1]);
                next.operations_reverse.push_back('M');
                --next.i;
                --next.j;
                push_if_allowed(std::move(next), current_score - diagonal);
            }
            const int insert = insertion_score(current.i, current.j);
            if (insert > 0) {
                auto next = current;
                next.state = TraceState::Insertion;
                push_if_allowed(std::move(next), current_score - insert);
            }
            const int remove = deletion_score(current.i, current.j);
            if (remove > 0) {
                auto next = current;
                next.state = TraceState::Deletion;
                push_if_allowed(std::move(next), current_score - remove);
            }
        } else if (current.state == TraceState::Insertion) {
            if (current.i == 0) continue;
            const int current_score = insertion_score(current.i, current.j);
            const int opened = match_score(current.i - 1, current.j) + scoring.gap_open;
            const int extended = insertion_score(current.i - 1, current.j) + scoring.gap_extend;
            current.aligned_query_reverse.push_back(query[current.i - 1]);
            current.aligned_reference_reverse.push_back('-');
            current.operations_reverse.push_back('I');
            --current.i;
            if (opened > neg_inf / 2) {
                auto next = current;
                next.state = TraceState::Match;
                push_if_allowed(std::move(next), current_score - opened);
            }
            if (extended > neg_inf / 2) {
                auto next = current;
                next.state = TraceState::Insertion;
                push_if_allowed(std::move(next), current_score - extended);
            }
        } else {
            if (current.j == 0) continue;
            const int current_score = deletion_score(current.i, current.j);
            const int opened = match_score(current.i, current.j - 1) + scoring.gap_open;
            const int extended = deletion_score(current.i, current.j - 1) + scoring.gap_extend;
            current.aligned_query_reverse.push_back('-');
            current.aligned_reference_reverse.push_back(reference[current.j - 1]);
            current.operations_reverse.push_back('D');
            --current.j;
            if (opened > neg_inf / 2) {
                auto next = current;
                next.state = TraceState::Match;
                push_if_allowed(std::move(next), current_score - opened);
            }
            if (extended > neg_inf / 2) {
                auto next = current;
                next.state = TraceState::Deletion;
                push_if_allowed(std::move(next), current_score - extended);
            }
        }
    }
    if (search_stats) {
        search_stats->expanded_states = expanded;
        search_stats->distinct_paths = results.size();
        search_stats->state_limit_hit = expanded >= maximum_trace_states &&
            results.size() < maximum_paths && !frontier.empty();
    }
    return results;
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
