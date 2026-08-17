#include "swiftig/engine.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <optional>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "swiftig/alignment.hpp"
#include "swiftig/allele_tree.hpp"

#ifndef SWIG_V_TREE_ROOT_ALIGNMENTS
#define SWIG_V_TREE_ROOT_ALIGNMENTS 3
#endif

#ifndef SWIG_V_TREE_WEAK_ROOT_ALIGNMENTS
#define SWIG_V_TREE_WEAK_ROOT_ALIGNMENTS 8
#endif

#ifndef SWIG_V_TREE_FINAL_REALIGN
#define SWIG_V_TREE_FINAL_REALIGN 0
#endif

#ifndef SWIG_V_TREE_EXHAUSTIVE_LEAVES
#define SWIG_V_TREE_EXHAUSTIVE_LEAVES 0
#endif

#ifndef SWIG_V_TREE_TRACEBACKS
#define SWIG_V_TREE_TRACEBACKS 1
#endif

#ifndef SWIG_V_TREE_TRACEBACK_TOLERANCE
#define SWIG_V_TREE_TRACEBACK_TOLERANCE 0
#endif

#ifndef SWIG_V_TREE_TRACE_STATE_LIMIT
#define SWIG_V_TREE_TRACE_STATE_LIMIT 8192
#endif

namespace swiftig {
namespace {

bool same_or_unknown_locus(const std::string& a, const std::string& b) {
    return a.empty() || b.empty() || a == b;
}

bool locus_has_d(const std::string& locus) {
    return locus.empty() || locus == "IGH" || locus == "TRB" || locus == "TRD";
}

bool supports_rearranged_pair(const SegmentHit& v, const SegmentHit& j) {
    // A rearranged V/J pair must face the recombination junction: the V hit
    // reaches its reference 3' end and the J hit reaches its reference 5' end.
    constexpr std::size_t v_end_slack = 45;
    constexpr std::size_t j_start_slack = 30;
    const auto v_aligned = static_cast<std::size_t>(
        v.alignment.matches + v.alignment.mismatches);
    const auto j_aligned = static_cast<std::size_t>(
        j.alignment.matches + j.alignment.mismatches);
    if (v_aligned + j_aligned < 80) return false;
    const bool faces_junction =
        v.alignment.reference_end + v_end_slack >= v.gene->sequence.size() &&
        j.alignment.reference_start <= j_start_slack;
    // Strong split alignments are retained for large biological/technical
    // deletions even if local scoring chooses the V-side HSP before the anchor.
    return faces_junction || (v_aligned >= 80 && j_aligned >= 30);
}

std::size_t aligned_bases(const Alignment& alignment) {
    return static_cast<std::size_t>(alignment.matches + alignment.mismatches);
}

bool canonical_base(char base) {
    return base == 'A' || base == 'C' || base == 'G' || base == 'T';
}

int substitution_score(char query, char reference, const Scoring& scoring) {
    if (!canonical_base(query) || !canonical_base(reference)) return 0;
    return query == reference ? scoring.match : scoring.mismatch;
}

std::pair<int, int> substitution_counts(char query, char reference) {
    if (query == '-' || reference == '-') return {0, 0};
    return canonical_base(query) && query == reference
        ? std::pair<int, int>{1, 0}
        : std::pair<int, int>{0, 1};
}

/**
 * A local alignment is allowed to reset across a short run of terminal SHM.
 * Once a strong V or J core has been selected, restore co-linear outer bases
 * where both the read and reference still exist.  This changes only endpoint
 * materialization (linear in the omitted prefix/suffix); candidate discovery
 * and the dynamic-programming search remain untouched.
 */
void extend_terminal_substitutions(
    const std::string& query,
    const std::string& reference,
    const Scoring& scoring,
    Alignment& alignment,
    bool extend_left,
    bool extend_right) {
    if (!alignment.valid()) return;
    if (extend_left) {
        const auto length = std::min(alignment.query_start, alignment.reference_start);
        if (length > 0) {
            const auto query_start = alignment.query_start - length;
            const auto reference_start = alignment.reference_start - length;
            for (std::size_t offset = 0; offset < length; ++offset) {
                const char q = query[query_start + offset];
                const char r = reference[reference_start + offset];
                alignment.score += substitution_score(q, r, scoring);
                const auto [matches, mismatches] = substitution_counts(q, r);
                alignment.matches += matches;
                alignment.mismatches += mismatches;
            }
            alignment.aligned_query.insert(0, query, query_start, length);
            alignment.aligned_reference.insert(0, reference, reference_start, length);
            alignment.query_start = query_start;
            alignment.reference_start = reference_start;
        }
    }
    if (extend_right) {
        const auto length = std::min(
            query.size() - alignment.query_end,
            reference.size() - alignment.reference_end);
        if (length > 0) {
            for (std::size_t offset = 0; offset < length; ++offset) {
                const char q = query[alignment.query_end + offset];
                const char r = reference[alignment.reference_end + offset];
                alignment.score += substitution_score(q, r, scoring);
                const auto [matches, mismatches] = substitution_counts(q, r);
                alignment.matches += matches;
                alignment.mismatches += mismatches;
            }
            alignment.aligned_query.append(query, alignment.query_end, length);
            alignment.aligned_reference.append(reference, alignment.reference_end, length);
            alignment.query_end += length;
            alignment.reference_end += length;
        }
    }
    refresh_airr_cigar(alignment, query.size(), reference.size());
}

void rank_segment_hits(std::vector<SegmentHit>& hits) {
    std::sort(hits.begin(), hits.end(), [](const SegmentHit& left, const SegmentHit& right) {
        if (left.alignment.score != right.alignment.score) return left.alignment.score > right.alignment.score;
        if (aligned_bases(left.alignment) != aligned_bases(right.alignment)) return aligned_bases(left.alignment) > aligned_bases(right.alignment);
        if (left.alignment.identity() != right.alignment.identity()) return left.alignment.identity() > right.alignment.identity();
        return left.gene->name < right.gene->name;
    });
}

struct FixedAlignmentPath {
    std::string query;
    std::string root_reference;
    std::size_t query_origin = 0;
    std::size_t reference_origin = 0;
    std::vector<std::uint32_t> query_prefix;
    std::vector<std::uint32_t> reference_prefix;
    std::vector<std::uint32_t> aligned_prefix;
    std::vector<int> reference_column;
};

FixedAlignmentPath extend_fixed_path(
    const std::string& query,
    const Gene& root,
    const Alignment& alignment) {
    FixedAlignmentPath path;
    const auto left = std::min(alignment.query_start, alignment.reference_start);
    const auto right = std::min(
        query.size() - alignment.query_end,
        root.sequence.size() - alignment.reference_end);
    path.query_origin = alignment.query_start - left;
    path.reference_origin = alignment.reference_start - left;
    path.query.reserve(left + alignment.aligned_query.size() + right);
    path.root_reference.reserve(left + alignment.aligned_reference.size() + right);
    path.query.append(query, path.query_origin, left);
    path.query += alignment.aligned_query;
    path.query.append(query, alignment.query_end, right);
    path.root_reference.append(root.sequence, path.reference_origin, left);
    path.root_reference += alignment.aligned_reference;
    path.root_reference.append(root.sequence, alignment.reference_end, right);

    const auto columns = path.query.size();
    path.query_prefix.assign(columns + 1, 0);
    path.reference_prefix.assign(columns + 1, 0);
    path.aligned_prefix.assign(columns + 1, 0);
    path.reference_column.assign(root.sequence.size(), -1);
    std::size_t reference_position = path.reference_origin;
    for (std::size_t column = 0; column < columns; ++column) {
        path.query_prefix[column + 1] = path.query_prefix[column] +
            static_cast<std::uint32_t>(path.query[column] != '-');
        path.reference_prefix[column + 1] = path.reference_prefix[column] +
            static_cast<std::uint32_t>(path.root_reference[column] != '-');
        path.aligned_prefix[column + 1] = path.aligned_prefix[column] +
            static_cast<std::uint32_t>(
                path.query[column] != '-' && path.root_reference[column] != '-');
        if (path.root_reference[column] != '-' && reference_position < root.sequence.size()) {
            path.reference_column[reference_position++] = static_cast<int>(column);
        }
    }
    return path;
}

struct ScoredSpan {
    int score = std::numeric_limits<int>::min() / 4;
    std::size_t start = 0;
    std::size_t end = 0;
};

bool better_span(const ScoredSpan& candidate, const ScoredSpan& current) {
    if (candidate.score != current.score) return candidate.score > current.score;
    if (candidate.end != current.end) return candidate.end < current.end;
    return candidate.start < current.start;
}

struct RangeScore {
    std::size_t left = 0;
    std::size_t right = 0;
    int sum = 0;
    ScoredSpan prefix;
    ScoredSpan suffix;
    ScoredSpan best;
};

RangeScore merge_range_scores(const RangeScore& left, const RangeScore& right) {
    RangeScore result;
    result.left = left.left;
    result.right = right.right;
    result.sum = left.sum + right.sum;

    result.prefix = left.prefix;
    const ScoredSpan joined_prefix{
        left.sum + right.prefix.score, left.left, right.prefix.end};
    if (better_span(joined_prefix, result.prefix)) result.prefix = joined_prefix;

    result.suffix = right.suffix;
    const ScoredSpan joined_suffix{
        left.suffix.score + right.sum, left.suffix.start, right.right};
    if (better_span(joined_suffix, result.suffix)) result.suffix = joined_suffix;

    result.best = left.best;
    if (better_span(right.best, result.best)) result.best = right.best;
    const ScoredSpan crossing{
        left.suffix.score + right.prefix.score,
        left.suffix.start,
        right.prefix.end};
    if (better_span(crossing, result.best)) result.best = crossing;
    return result;
}

class FixedPathScorer {
public:
    FixedPathScorer(const FixedAlignmentPath& path, const Scoring& scoring)
        : path_(path), scoring_(scoring) {
        std::vector<int> values(path.query.size(), 0);
        char previous_gap = '\0';
        for (std::size_t column = 0; column < values.size(); ++column) {
            const char query_base = path.query[column];
            const char reference_base = path.root_reference[column];
            const char gap = query_base == '-' ? 'D' : reference_base == '-' ? 'I' : '\0';
            if (gap != '\0') {
                values[column] = gap == previous_gap ? scoring.gap_extend : scoring.gap_open;
            } else {
                values[column] = substitution_score(query_base, reference_base, scoring);
            }
            previous_gap = gap;
        }
        tree_.resize(std::max<std::size_t>(1, values.size() * 4));
        if (!values.empty()) build(1, 0, values.size(), values);
    }

    [[nodiscard]] ScoredSpan best() const noexcept {
        return path_.query.empty() ? ScoredSpan{} : tree_[1].best;
    }

    void set_reference_base(std::size_t reference_position, char reference_base) {
        if (reference_position >= path_.reference_column.size()) return;
        const int signed_column = path_.reference_column[reference_position];
        if (signed_column < 0) return;
        const auto column = static_cast<std::size_t>(signed_column);
        const int value = substitution_score(path_.query[column], reference_base, scoring_);
        update(1, column, value);
    }

private:
    const FixedAlignmentPath& path_;
    const Scoring& scoring_;
    std::vector<RangeScore> tree_;

    void build(
        std::size_t node,
        std::size_t left,
        std::size_t right,
        const std::vector<int>& values) {
        if (right - left == 1) {
            const ScoredSpan span{values[left], left, right};
            tree_[node] = RangeScore{left, right, values[left], span, span, span};
            return;
        }
        const auto middle = left + (right - left) / 2;
        build(node * 2, left, middle, values);
        build(node * 2 + 1, middle, right, values);
        tree_[node] = merge_range_scores(tree_[node * 2], tree_[node * 2 + 1]);
    }

    void update(std::size_t node, std::size_t position, int value) {
        auto& range = tree_[node];
        if (range.right - range.left == 1) {
            const ScoredSpan span{value, range.left, range.right};
            range = RangeScore{range.left, range.right, value, span, span, span};
            return;
        }
        if (position < tree_[node * 2].right) update(node * 2, position, value);
        else update(node * 2 + 1, position, value);
        range = merge_range_scores(tree_[node * 2], tree_[node * 2 + 1]);
    }
};

Alignment materialize_fixed_tree_alignment(
    const FixedAlignmentPath& path,
    const Gene& allele,
    const ScoredSpan& span,
    std::size_t query_size) {
    Alignment result;
    if (span.score <= 0 || span.start >= span.end || span.end > path.query.size()) return result;
    result.score = span.score;
    result.query_start = path.query_origin + path.query_prefix[span.start];
    result.query_end = path.query_origin + path.query_prefix[span.end];
    result.reference_start = path.reference_origin + path.reference_prefix[span.start];
    result.reference_end = path.reference_origin + path.reference_prefix[span.end];
    result.aligned_query = path.query.substr(span.start, span.end - span.start);
    result.aligned_reference = path.root_reference.substr(span.start, span.end - span.start);
    std::size_t reference_position = result.reference_start;
    for (std::size_t column = 0; column < result.aligned_reference.size(); ++column) {
        const char query_base = result.aligned_query[column];
        char& reference_base = result.aligned_reference[column];
        if (reference_base != '-' && reference_position < allele.sequence.size()) {
            reference_base = allele.sequence[reference_position++];
        }
        if (query_base == '-' && reference_base != '-') {
            ++result.deletions;
        } else if (query_base != '-' && reference_base == '-') {
            ++result.insertions;
        } else if (query_base != '-' && reference_base != '-') {
            if (canonical_base(query_base) && query_base == reference_base) ++result.matches;
            else ++result.mismatches;
        }
    }
    refresh_airr_cigar(result, query_size, allele.sequence.size());
    return result;
}

std::size_t longest_exact_run(const Alignment& alignment) {
    std::size_t longest = 0;
    std::size_t current = 0;
    for (std::size_t i = 0; i < alignment.aligned_query.size(); ++i) {
        if (alignment.aligned_query[i] != '-' &&
            alignment.aligned_query[i] == alignment.aligned_reference[i]) {
            longest = std::max(longest, ++current);
        } else {
            current = 0;
        }
    }
    return longest;
}

bool is_cys_codon(const std::string& sequence, std::size_t position) {
    return position + 2 < sequence.size() && sequence[position] == 'T' &&
        sequence[position + 1] == 'G' &&
        (sequence[position + 2] == 'T' || sequence[position + 2] == 'C');
}

bool is_w_or_f_codon(const std::string& sequence, std::size_t position) {
    if (position + 2 >= sequence.size()) return false;
    if (sequence.compare(position, 3, "TGG") == 0) return true;
    return sequence[position] == 'T' && sequence[position + 1] == 'T' &&
        (sequence[position + 2] == 'T' || sequence[position + 2] == 'C');
}

bool is_gly_codon(const std::string& sequence, std::size_t position) {
    return position + 2 < sequence.size() && sequence[position] == 'G' && sequence[position + 1] == 'G';
}

std::optional<std::size_t> v_anchor_in_reference(const Gene& gene) {
    const auto& sequence = gene.sequence;
    if (sequence.size() < 3) return std::nullopt;
    const std::size_t begin = sequence.size() > 120 ? sequence.size() - 120 : 0;
    std::optional<std::size_t> framed_anchor;
    std::optional<std::size_t> any_anchor;
    const std::size_t frame = gene.coding_frame_start >= 0
        ? static_cast<std::size_t>(gene.coding_frame_start) : 0;
    for (std::size_t i = begin; i + 2 < sequence.size(); ++i) {
        if (!is_cys_codon(sequence, i)) continue;
        any_anchor = i;
        if (i % 3 == frame) framed_anchor = i;
    }
    return framed_anchor ? framed_anchor : any_anchor;
}

std::optional<std::size_t> j_anchor_in_reference(const std::string& sequence) {
    const std::size_t end = std::min<std::size_t>(sequence.size(), 120);
    for (std::size_t i = 0; i + 5 < end; ++i) {
        if (is_w_or_f_codon(sequence, i) && is_gly_codon(sequence, i + 3)) return i;
    }
    return std::nullopt;
}

std::optional<std::size_t> fallback_v_anchor(const std::string& sequence, const Alignment& alignment) {
    const std::size_t begin = alignment.query_end > 120 ? alignment.query_end - 120 : alignment.query_start;
    const std::size_t end = std::min(sequence.size(), alignment.query_end + 18);
    std::optional<std::size_t> anchor;
    for (std::size_t i = begin; i + 2 < end; ++i) {
        if (is_cys_codon(sequence, i)) anchor = i;
    }
    return anchor;
}

std::optional<std::size_t> fallback_j_anchor(
    const std::string& sequence,
    const Alignment& alignment,
    std::optional<std::size_t> cys_anchor) {
    std::size_t begin = alignment.query_start > 12 ? alignment.query_start - 12 : 0;
    if (cys_anchor) begin = std::max(begin, *cys_anchor + 3);
    const std::size_t end = std::min(sequence.size(), alignment.query_start + 120);
    for (std::size_t i = begin; i + 5 < end; ++i) {
        if (cys_anchor && (i - *cys_anchor) % 3 != 0) continue;
        if (is_w_or_f_codon(sequence, i) && is_gly_codon(sequence, i + 3)) return i;
    }
    return std::nullopt;
}

bool net_frameshift(const Alignment& alignment) {
    int shift = (alignment.insertions - alignment.deletions) % 3;
    if (shift < 0) shift += 3;
    return shift != 0;
}

void merge_equivalent_calls(SegmentHit& selected, const std::vector<SegmentHit>& hits) {
    selected.call = selected.gene->name;
    for (const auto& hit : hits) {
        if (hit.gene == selected.gene || hit.alignment.score != selected.alignment.score ||
            aligned_bases(hit.alignment) != aligned_bases(selected.alignment) ||
            hit.alignment.query_start != selected.alignment.query_start ||
            hit.alignment.query_end != selected.alignment.query_end ||
            !same_or_unknown_locus(hit.gene->locus, selected.gene->locus)) continue;
        selected.call += "," + hit.gene->name;
    }
}

std::vector<SegmentHit> uncertain_alternatives(
    const SegmentHit& selected,
    const std::vector<SegmentHit>& hits) {
    std::unordered_set<std::string> selected_calls;
    std::size_t start = 0;
    while (start <= selected.call.size()) {
        const auto end = selected.call.find(',', start);
        selected_calls.insert(selected.call.substr(
            start, end == std::string::npos ? std::string::npos : end - start));
        if (end == std::string::npos) break;
        start = end + 1;
    }
    const int tolerance = std::max(2, std::abs(selected.alignment.score) / 100);
    std::vector<SegmentHit> alternatives;
    for (const auto& hit : hits) {
        if (!hit.gene || selected_calls.contains(hit.gene->name)) continue;
        if (hit.alignment.score + tolerance < selected.alignment.score) continue;
        alternatives.push_back(hit);
    }
    return alternatives;
}

std::optional<std::pair<std::size_t, std::size_t>> map_reference_interval(
    const Alignment& alignment,
    std::size_t start,
    std::size_t end) {
    if (start >= end) return std::nullopt;
    std::optional<std::size_t> query_start;
    for (std::size_t position = start; position < end; ++position) {
        query_start = map_reference_to_query(alignment, position);
        if (query_start) break;
    }
    std::optional<std::size_t> query_last;
    for (std::size_t position = end; position > start; --position) {
        query_last = map_reference_to_query(alignment, position - 1);
        if (query_last) break;
    }
    if (!query_start || !query_last || *query_last < *query_start) return std::nullopt;
    return std::pair<std::size_t, std::size_t>{*query_start, *query_last + 1};
}

}  // namespace

AnnotationEngine::AnnotationEngine(const GermlineDatabase& database, EngineOptions options)
    : database_(database),
      options_(std::move(options)),
      v_statistics_(reference_database_statistics(database.v.genes())),
      d_statistics_(reference_database_statistics(database.d.genes())),
      j_statistics_(reference_database_statistics(database.j.genes())),
      c_statistics_(reference_database_statistics(database.c.genes())) {}

std::vector<SegmentHit> AnnotationEngine::align_candidates(
    const std::string& query,
    const SegmentIndex& index,
    std::size_t top_n,
    const Scoring& scoring,
    std::size_t min_length,
    const std::string& locus_filter,
    const std::vector<Candidate>* candidate_hints) const {
    if (index.empty() || query.empty() || top_n == 0) return {};
    const auto canonical_bases = static_cast<std::size_t>(std::count_if(
        query.begin(), query.end(), [](char base) {
            return base == 'A' || base == 'C' || base == 'G' || base == 'T';
        }));
    const auto minimum_possible_matches = static_cast<std::size_t>(
        std::ceil(static_cast<double>(min_length) * options_.min_identity));
    if (canonical_bases < minimum_possible_matches) return {};
    const std::size_t pool_size = std::max<std::size_t>(top_n * 4, 16);
    const std::size_t candidate_limit = std::min(
        index.genes().size(), pool_size * (locus_filter.empty() ? 1 : 2));
    auto candidates = candidate_hints && !candidate_hints->empty()
        ? *candidate_hints : index.candidates(query, candidate_limit);
    std::vector<Candidate> work;
    work.reserve(candidates.empty() ? index.genes().size() : candidates.size());
    for (const auto& candidate : candidates) {
        const auto& gene = index.genes()[candidate.gene_index];
        if (locus_filter.empty() || same_or_unknown_locus(locus_filter, gene.locus)) work.push_back(candidate);
    }
    if (work.empty()) {
        for (std::uint32_t i = 0; i < index.genes().size(); ++i) {
            const auto& gene = index.genes()[i];
            if (locus_filter.empty() || same_or_unknown_locus(locus_filter, gene.locus)) {
                work.push_back(Candidate{i, std::numeric_limits<int>::max(), 0});
            }
        }
    }
    if (work.empty()) return {};
    // A strong seed mode normally needs only a small refinement set. Retain
    // the full safety pool when the short-seed rescue tier was needed, or for
    // exhaustive seedless fallback, so robustness is paid for only on hard
    // reads rather than on every ordinary repertoire record.
    const bool exhaustive = work.front().diagonal == std::numeric_limits<int>::max();
    if (!exhaustive && !work.front().weak_seed_signal) {
        auto refinement_limit = std::max<std::size_t>(top_n, 3);
        if (options_.assigner_strategy == AssignerStrategy::Aer &&
            index.kmer_size() == 9 && min_length == options_.min_v_length &&
            work.size() > refinement_limit) {
            const auto relative_margin =
                static_cast<std::uint32_t>(work.front().votes) * 5U / 100U;
            const auto vote_margin = std::max<std::uint32_t>(8U, relative_margin);
            const auto cap = std::min<std::size_t>(work.size(), 16);
            while (refinement_limit < cap &&
                   static_cast<std::uint32_t>(work[refinement_limit].votes) + vote_margin >=
                       static_cast<std::uint32_t>(work.front().votes)) {
                ++refinement_limit;
            }
        }
        if (work.size() > refinement_limit) work.resize(refinement_limit);
    }
    std::vector<SegmentHit> hits;
    hits.reserve(work.size());
    for (const auto& candidate : work) {
        const auto& gene = index.genes()[candidate.gene_index];
        const bool seeded = candidate.diagonal != std::numeric_limits<int>::max();
        const int adaptive_band = seeded
            ? std::max(options_.band_width, std::min(
                options_.max_band_width, candidate.diagonal_span + 8)) : -1;
        auto alignment = local_align_affine(
            query, gene.sequence, scoring,
            seeded ? candidate.diagonal : 0,
            adaptive_band);
        if (!alignment.valid() || aligned_bases(alignment) < min_length ||
            alignment.identity() < options_.min_identity) continue;
        hits.push_back(SegmentHit{&gene, std::move(alignment), gene.name, query.size()});
    }
    std::sort(hits.begin(), hits.end(), [](const SegmentHit& a, const SegmentHit& b) {
        if (a.alignment.score != b.alignment.score) return a.alignment.score > b.alignment.score;
        if (aligned_bases(a.alignment) != aligned_bases(b.alignment)) {
            return aligned_bases(a.alignment) > aligned_bases(b.alignment);
        }
        if (a.alignment.identity() != b.alignment.identity()) {
            return a.alignment.identity() > b.alignment.identity();
        }
        return a.gene->name < b.gene->name;
    });
    if (hits.size() > top_n) hits.resize(top_n);
    return hits;
}

std::vector<SegmentHit> AnnotationEngine::align_v_allele_tree(
    const std::string& query,
    std::size_t top_n,
    const Scoring& scoring,
    std::size_t min_length) const {
    std::vector<SegmentHit> hits;
    AlleleTreeSearchStats stats;
    stats.queries = 1;
    const auto& tree = database_.v_tree;
    const auto& root_index = tree.roots();
    if (query.empty() || top_n == 0 || tree.empty() || root_index.empty()) {
        tree.record_search(stats);
        return hits;
    }

    const auto canonical_bases = static_cast<std::size_t>(std::count_if(
        query.begin(), query.end(), canonical_base));
    const auto minimum_possible_matches = static_cast<std::size_t>(
        std::ceil(static_cast<double>(min_length) * options_.min_identity));
    if (canonical_bases < minimum_possible_matches) {
        tree.record_search(stats);
        return hits;
    }

    const auto candidate_limit = std::min<std::size_t>(root_index.genes().size(), 16);
    auto candidates = root_index.candidates(query, candidate_limit);
    if (candidates.empty()) {
        candidates.reserve(root_index.genes().size());
        for (std::uint32_t index = 0; index < root_index.genes().size(); ++index) {
            candidates.push_back(Candidate{index, std::numeric_limits<int>::max(), 0});
        }
    }
    stats.root_candidates = static_cast<std::uint32_t>(candidates.size());
    const bool exhaustive = candidates.front().diagonal == std::numeric_limits<int>::max();
    const bool weak = !exhaustive && candidates.front().weak_seed_signal;
    const auto configured_root_limit = weak
        ? static_cast<std::size_t>(SWIG_V_TREE_WEAK_ROOT_ALIGNMENTS)
        : static_cast<std::size_t>(SWIG_V_TREE_ROOT_ALIGNMENTS);
    const auto root_limit = exhaustive
        ? candidates.size()
        : std::min(candidates.size(), configured_root_limit);

    struct RootAlignment {
        std::size_t cluster = 0;
        Candidate candidate;
        Alignment alignment;
        FixedAlignmentPath path;
    };
    std::vector<RootAlignment> root_alignments;
    root_alignments.reserve(root_limit);
    for (std::size_t candidate_index = 0; candidate_index < root_limit; ++candidate_index) {
        const auto& candidate = candidates[candidate_index];
        if (candidate.gene_index >= tree.clusters().size()) continue;
        const auto& cluster = tree.clusters()[candidate.gene_index];
        const auto& root = database_.v.genes()[cluster.root_gene_index];
        const bool seeded = candidate.diagonal != std::numeric_limits<int>::max();
        const int adaptive_band = seeded
            ? std::max(options_.band_width, std::min(
                options_.max_band_width, candidate.diagonal_span + 8)) : -1;
        auto alignment = local_align_affine(
            query, root.sequence, scoring,
            seeded ? candidate.diagonal : 0,
            adaptive_band);
        ++stats.root_alignments;
        if (!alignment.valid() || aligned_bases(alignment) < min_length) continue;
        ++stats.root_tracebacks;
        root_alignments.push_back(RootAlignment{
            candidate.gene_index,
            candidate,
            std::move(alignment),
            {},
        });
        root_alignments.back().path = extend_fixed_path(
            query, root, root_alignments.back().alignment);
    }

    struct ProxyHit {
        std::size_t root_alignment = 0;
        std::uint32_t gene_index = 0;
        ScoredSpan span;
    };
    std::vector<ProxyHit> proxies;
    const auto score_root_path = [&](std::size_t root_alignment_index) {
        const auto& root = root_alignments[root_alignment_index];
        const auto& cluster = tree.clusters()[root.cluster];
        if (cluster.nodes.empty()) return;
        ++stats.clusters_scored;
        stats.nodes_scored += static_cast<std::uint32_t>(cluster.nodes.size());
        FixedPathScorer scorer(root.path, scoring);
        std::vector<std::vector<std::size_t>> children(cluster.nodes.size());
        for (std::size_t node_index = 1; node_index < cluster.nodes.size(); ++node_index) {
            children[cluster.nodes[node_index].parent].push_back(node_index);
        }
        const auto visit = [&](auto&& self, std::size_t node_index) -> void {
            const auto& node = cluster.nodes[node_index];
            proxies.push_back(ProxyHit{
                root_alignment_index,
                node.gene_index,
                scorer.best(),
            });
            for (const auto child_index : children[node_index]) {
                const auto& child = cluster.nodes[child_index];
                for (const auto& mutation : child.edge_mutations) {
                    ++stats.mutation_updates;
                    scorer.set_reference_base(mutation.position, mutation.child_base);
                }
                self(self, child_index);
                for (const auto& mutation : child.edge_mutations) {
                    scorer.set_reference_base(mutation.position, mutation.parent_base);
                }
            }
        };
        visit(visit, 0);
    };
    for (std::size_t root_alignment_index = 0;
         root_alignment_index < root_alignments.size(); ++root_alignment_index) {
        score_root_path(root_alignment_index);
    }
    const auto proxy_better = [&](const ProxyHit& left, const ProxyHit& right) {
        if (left.span.score != right.span.score) return left.span.score > right.span.score;
        const auto& left_path = root_alignments[left.root_alignment].path;
        const auto& right_path = root_alignments[right.root_alignment].path;
        const auto left_aligned = left_path.aligned_prefix[left.span.end] -
            left_path.aligned_prefix[left.span.start];
        const auto right_aligned = right_path.aligned_prefix[right.span.end] -
            right_path.aligned_prefix[right.span.start];
        if (left_aligned != right_aligned) return left_aligned > right_aligned;
        const auto& left_name = database_.v.genes()[left.gene_index].name;
        const auto& right_name = database_.v.genes()[right.gene_index].name;
        if (left_name != right_name) return left_name < right_name;
        return left.root_alignment < right.root_alignment;
    };
    std::sort(proxies.begin(), proxies.end(), proxy_better);

#if SWIG_V_TREE_TRACEBACKS > 1
    // First score every selected root once. Only if the provisional winning
    // allele actually contains an indel do we reopen that one root DP and
    // enumerate its near-optimal traceback geometries. The selected-path
    // trigger covered 33/38 fixed-path discrepancies in the diagnostic data,
    // while avoiding multipath work on ordinary substitution-only queries.
    std::optional<std::size_t> multipath_root_alignment;
    for (const auto& proxy : proxies) {
        const auto& gene = database_.v.genes()[proxy.gene_index];
        const auto& root = root_alignments[proxy.root_alignment];
        const auto provisional = materialize_fixed_tree_alignment(
            root.path, gene, proxy.span, query.size());
        if (!provisional.valid() || aligned_bases(provisional) < min_length ||
            provisional.identity() < options_.min_identity) continue;
        if (provisional.insertions + provisional.deletions > 0) {
            multipath_root_alignment = proxy.root_alignment;
        }
        break;
    }
    if (multipath_root_alignment) {
        const auto seed_index = *multipath_root_alignment;
        const auto seed = root_alignments[seed_index];
        const auto& cluster = tree.clusters()[seed.cluster];
        const auto& root_gene = database_.v.genes()[cluster.root_gene_index];
        const bool seeded = seed.candidate.diagonal != std::numeric_limits<int>::max();
        const int adaptive_band = seeded
            ? std::max(options_.band_width, std::min(
                options_.max_band_width, seed.candidate.diagonal_span + 8)) : -1;
        TracebackSearchStats traceback_stats;
        auto alternatives = local_align_affine_paths(
            query, root_gene.sequence, scoring,
            seeded ? seed.candidate.diagonal : 0,
            adaptive_band,
            static_cast<std::size_t>(SWIG_V_TREE_TRACEBACKS),
            SWIG_V_TREE_TRACEBACK_TOLERANCE,
            static_cast<std::size_t>(SWIG_V_TREE_TRACE_STATE_LIMIT),
            &traceback_stats);
        ++stats.root_alignments;
        ++stats.multipath_searches;
        stats.trace_states += static_cast<std::uint32_t>(traceback_stats.expanded_states);
        stats.trace_limit_hits += static_cast<std::uint32_t>(traceback_stats.state_limit_hit);
        // alternatives[0] is the already-scored deterministic optimum.
        for (std::size_t path_index = 1; path_index < alternatives.size(); ++path_index) {
            auto& alignment = alternatives[path_index];
            if (!alignment.valid() || aligned_bases(alignment) < min_length) continue;
            root_alignments.push_back(RootAlignment{
                seed.cluster,
                seed.candidate,
                std::move(alignment),
                {},
            });
            root_alignments.back().path = extend_fixed_path(
                query, root_gene, root_alignments.back().alignment);
            ++stats.root_tracebacks;
            score_root_path(root_alignments.size() - 1);
        }
        std::sort(proxies.begin(), proxies.end(), proxy_better);
    }
#endif

    hits.reserve(top_n);
#if SWIG_V_TREE_EXHAUSTIVE_LEAVES
    // Validation build only: align every leaf admitted by the root search.
    // Comparing this ranking with the sparse-delta ranking directly measures
    // the approximation made by holding the root alignment fixed.
    std::vector<SegmentHit> exhaustive_hits;
    exhaustive_hits.reserve(proxies.size());
    for (const auto& proxy : proxies) {
        const auto& gene = database_.v.genes()[proxy.gene_index];
        const auto& root = root_alignments[proxy.root_alignment];
        const bool seeded = root.candidate.diagonal != std::numeric_limits<int>::max();
        const int adaptive_band = seeded
            ? std::max(options_.band_width, std::min(
                options_.max_band_width, root.candidate.diagonal_span + 8)) : -1;
        auto alignment = local_align_affine(
            query, gene.sequence, scoring,
            seeded ? root.candidate.diagonal : 0,
            adaptive_band);
        ++stats.final_realignments;
        if (!alignment.valid() || aligned_bases(alignment) < min_length ||
            alignment.identity() < options_.min_identity) continue;
        exhaustive_hits.push_back(SegmentHit{&gene, std::move(alignment), gene.name, query.size()});
    }
    std::sort(exhaustive_hits.begin(), exhaustive_hits.end(), [](const SegmentHit& left, const SegmentHit& right) {
        if (left.alignment.score != right.alignment.score) {
            return left.alignment.score > right.alignment.score;
        }
        if (aligned_bases(left.alignment) != aligned_bases(right.alignment)) {
            return aligned_bases(left.alignment) > aligned_bases(right.alignment);
        }
        if (left.alignment.identity() != right.alignment.identity()) {
            return left.alignment.identity() > right.alignment.identity();
        }
        return left.gene->name < right.gene->name;
    });
    if (exhaustive_hits.size() > top_n) exhaustive_hits.resize(top_n);
    hits = std::move(exhaustive_hits);
#else
    std::unordered_set<std::uint32_t> emitted_genes;
    for (const auto& proxy : proxies) {
        if (emitted_genes.contains(proxy.gene_index)) continue;
        const auto& gene = database_.v.genes()[proxy.gene_index];
        const auto& root = root_alignments[proxy.root_alignment];
#if SWIG_V_TREE_FINAL_REALIGN
        const bool seeded = root.candidate.diagonal != std::numeric_limits<int>::max();
        const int adaptive_band = seeded
            ? std::max(options_.band_width, std::min(
                options_.max_band_width, root.candidate.diagonal_span + 8)) : -1;
        auto alignment = local_align_affine(
            query, gene.sequence, scoring,
            seeded ? root.candidate.diagonal : 0,
            adaptive_band);
        ++stats.final_realignments;
#else
        auto alignment = materialize_fixed_tree_alignment(
            root.path, gene, proxy.span, query.size());
#endif
        if (!alignment.valid() || aligned_bases(alignment) < min_length ||
            alignment.identity() < options_.min_identity) continue;
        emitted_genes.insert(proxy.gene_index);
        hits.push_back(SegmentHit{&gene, std::move(alignment), gene.name, query.size()});
        if (hits.size() == top_n) break;
    }
#endif
    std::sort(hits.begin(), hits.end(), [](const SegmentHit& left, const SegmentHit& right) {
        if (left.alignment.score != right.alignment.score) {
            return left.alignment.score > right.alignment.score;
        }
        if (aligned_bases(left.alignment) != aligned_bases(right.alignment)) {
            return aligned_bases(left.alignment) > aligned_bases(right.alignment);
        }
        if (left.alignment.identity() != right.alignment.identity()) {
            return left.alignment.identity() > right.alignment.identity();
        }
        return left.gene->name < right.gene->name;
    });
    tree.record_search(stats);
    return hits;
}

AnnotationEngine::OrientationResult AnnotationEngine::annotate_orientation(
    const std::string& sequence,
    const OrientationHints* hints) const {
    OrientationResult result;
    result.oriented_sequence = sequence;
    auto v_hits = options_.assigner_strategy == AssignerStrategy::RiatMp
        ? align_v_allele_tree(
            sequence, options_.top_v, options_.v_scoring, options_.min_v_length)
        : align_candidates(
            sequence, database_.v, options_.top_v, options_.v_scoring,
            options_.min_v_length, "", hints ? &hints->v : nullptr);
    auto j_hits = align_candidates(
        sequence, database_.j, options_.top_j, options_.j_scoring, options_.min_j_length,
        "", hints ? &hints->j : nullptr);
    for (auto& hit : v_hits) extend_terminal_substitutions(
        sequence, hit.gene->sequence, options_.v_scoring, hit.alignment, true, false);
    for (auto& hit : j_hits) extend_terminal_substitutions(
        sequence, hit.gene->sequence, options_.j_scoring, hit.alignment, false, true);
    rank_segment_hits(v_hits);
    rank_segment_hits(j_hits);

    double best_pair_score = -std::numeric_limits<double>::infinity();
    const std::size_t overlap = options_.allow_vdj_overlap
        ? static_cast<std::size_t>(std::max(options_.max_vdj_overlap, 0)) : 0;
    for (std::size_t vi = 0; vi <= v_hits.size(); ++vi) {
        for (std::size_t ji = 0; ji <= j_hits.size(); ++ji) {
            const SegmentHit* v = vi < v_hits.size() ? &v_hits[vi] : nullptr;
            const SegmentHit* j = ji < j_hits.size() ? &j_hits[ji] : nullptr;
            if (!v && !j) continue;
            if (v && j) {
                if (!same_or_unknown_locus(v->gene->locus, j->gene->locus)) continue;
                if (!supports_rearranged_pair(*v, *j)) continue;
                if (v->alignment.query_start > j->alignment.query_start) continue;
                if (v->alignment.query_end > j->alignment.query_start + overlap) continue;
                const auto distance = j->alignment.query_start > v->alignment.query_end
                    ? j->alignment.query_start - v->alignment.query_end : 0;
                if (distance > options_.max_junction_span) continue;
            }
            double score = static_cast<double>((v ? v->alignment.score : 0) + (j ? j->alignment.score : 0));
            if (v && j) score += 24.0;
            if (score > best_pair_score) {
                best_pair_score = score;
                result.v = v ? std::optional<SegmentHit>(*v) : std::nullopt;
                result.j = j ? std::optional<SegmentHit>(*j) : std::nullopt;
            }
        }
    }

    if (result.v && result.j && !database_.d.empty()) {
        const std::string locus = !result.v->gene->locus.empty()
            ? result.v->gene->locus : result.j->gene->locus;
        if (locus_has_d(locus)) {
            const std::size_t window_start = result.v->alignment.query_end > overlap
                ? result.v->alignment.query_end - overlap : 0;
            const std::size_t window_end = std::min(
                sequence.size(), result.j->alignment.query_start + overlap);
            if (window_end > window_start) {
                const std::string window = sequence.substr(window_start, window_end - window_start);
                auto d_hits = align_candidates(
                    window, database_.d, options_.top_d, options_.d_scoring,
                    options_.min_d_match, locus);
                std::vector<SegmentHit> viable_d_hits;
                for (auto hit : d_hits) {
                    if (longest_exact_run(hit.alignment) < options_.min_d_match) continue;
                    hit.alignment.query_start += window_start;
                    hit.alignment.query_end += window_start;
                    refresh_airr_cigar(hit.alignment, sequence.size(), hit.gene->sequence.size());
                    viable_d_hits.push_back(std::move(hit));
                }
                if (!viable_d_hits.empty()) {
                    result.d = viable_d_hits.front();
                    merge_equivalent_calls(*result.d, viable_d_hits);
                    result.d_alternatives = uncertain_alternatives(*result.d, viable_d_hits);
                }
            }
        }
    }

    const std::string locus = result.v ? result.v->gene->locus : (result.j ? result.j->gene->locus : "");
    if (!database_.c.empty() && result.j) {
        auto c_hits = align_candidates(
            sequence, database_.c, options_.top_c, options_.c_scoring, options_.min_c_length, locus);
        std::vector<SegmentHit> viable_c_hits;
        for (const auto& hit : c_hits) {
            // NCBI C references may carry one leading J-derived base to complete a codon.
            if (!result.j || hit.alignment.query_start + 1 >= result.j->alignment.query_end) {
                viable_c_hits.push_back(hit);
            }
        }
        if (!viable_c_hits.empty()) {
            result.c = viable_c_hits.front();
            merge_equivalent_calls(*result.c, viable_c_hits);
            result.c_alternatives = uncertain_alternatives(*result.c, viable_c_hits);
        }
    }

    result.rank_score = best_pair_score;
    if (result.v) {
        merge_equivalent_calls(*result.v, v_hits);
        result.v_alternatives = uncertain_alternatives(*result.v, v_hits);
    }
    if (result.j) {
        merge_equivalent_calls(*result.j, j_hits);
        result.j_alternatives = uncertain_alternatives(*result.j, j_hits);
    }
    if (result.d) result.rank_score += result.d->alignment.score * 0.5;
    if (result.c) result.rank_score += result.c->alignment.score * 0.25;
    return result;
}

Annotation AnnotationEngine::annotate(const SequenceRecord& record, const AnnotationHints* hints) const {
    Annotation annotation;
    annotation.sequence_id = record.id;
    annotation.sequence = record.sequence;
    annotation.quality = record.quality;
    OrientationResult chosen;
    chosen.rank_score = -std::numeric_limits<double>::infinity();
    const std::string reversed = options_.search_reverse ? reverse_complement(record.sequence) : std::string{};
    bool run_forward = options_.search_forward;
    bool run_reverse = options_.search_reverse;
    if (run_forward && run_reverse && hints == nullptr) {
        const auto forward_strength = orientation_seed_strength(record.sequence);
        const auto reverse_strength = orientation_seed_strength(reversed);
        if (forward_strength >= 6 && forward_strength > reverse_strength * 2 + 4) run_reverse = false;
        else if (reverse_strength >= 6 && reverse_strength > forward_strength * 2 + 4) run_forward = false;
    }
    if (run_forward) chosen = annotate_orientation(
        record.sequence, hints ? &hints->forward : nullptr);
    if (run_reverse) {
        auto reverse = annotate_orientation(reversed, hints ? &hints->reverse : nullptr);
        if (!run_forward || reverse.rank_score > chosen.rank_score) {
            chosen = std::move(reverse);
            annotation.rev_comp = true;
        }
    }
    annotation.oriented_sequence = std::move(chosen.oriented_sequence);
    annotation.v = std::move(chosen.v);
    annotation.d = std::move(chosen.d);
    annotation.j = std::move(chosen.j);
    annotation.c = std::move(chosen.c);
    annotation.v_alternatives = std::move(chosen.v_alternatives);
    annotation.d_alternatives = std::move(chosen.d_alternatives);
    annotation.j_alternatives = std::move(chosen.j_alternatives);
    annotation.c_alternatives = std::move(chosen.c_alternatives);
    const auto assign_support = [&](std::optional<SegmentHit>& hit,
                                    const Scoring& scoring,
                                    const ReferenceDatabaseStatistics& statistics) {
        if (!hit) return;
        const auto query_length = hit->search_query_length > 0
            ? hit->search_query_length : annotation.oriented_sequence.size();
        hit->support = calibrated_alignment_evalue(
            hit->alignment.score, query_length, statistics, scoring);
    };
    assign_support(annotation.v, options_.v_scoring, v_statistics_);
    assign_support(annotation.d, options_.d_scoring, d_statistics_);
    assign_support(annotation.j, options_.j_scoring, j_statistics_);
    assign_support(annotation.c, options_.c_scoring, c_statistics_);
    if (annotation.v && !annotation.v->gene->locus.empty()) annotation.locus = annotation.v->gene->locus;
    else if (annotation.j) annotation.locus = annotation.j->gene->locus;
    if (annotation.v) {
        annotation.v_annotation_source = annotation.v->gene->annotation_source;
        if (annotation.v->gene->region_bounds[9] > 0) annotation.region_definition = "IMGT";
    }
    if (annotation.j) annotation.j_annotation_source = annotation.j->gene->annotation_source;

    if (annotation.v && annotation.j) {
        const auto v_end = annotation.v->alignment.query_end;
        const auto j_start = annotation.j->alignment.query_start;
        if (annotation.d) {
            const auto d_start = annotation.d->alignment.query_start;
            const auto d_end = annotation.d->alignment.query_end;
            if (d_start > v_end) annotation.np1 = annotation.oriented_sequence.substr(v_end, d_start - v_end);
            if (j_start > d_end) annotation.np2 = annotation.oriented_sequence.substr(d_end, j_start - d_end);
        } else if (j_start > v_end) {
            annotation.np1 = annotation.oriented_sequence.substr(v_end, j_start - v_end);
        }
    }
    stitch_alignment(annotation);
    annotate_junction(annotation);
    if (annotation.v && annotation.sequence_aa.empty()) {
        const auto coding_frame = annotation.v->gene->coding_frame_start >= 0
            ? static_cast<std::size_t>(annotation.v->gene->coding_frame_start) : 0;
        const auto reference_phase = (annotation.v->alignment.reference_start + 3 - coding_frame) % 3;
        const auto frame = (annotation.v->alignment.query_start + 3 - reference_phase) % 3;
        annotation.sequence_aa = translate_dna(annotation.oriented_sequence, frame);
        annotation.sequence_frame = static_cast<int>(frame + 1);
        annotation.v_frameshift = net_frameshift(annotation.v->alignment);
    }
    annotate_v_regions(annotation);
    return annotation;
}

void AnnotationEngine::annotate_v_regions(Annotation& annotation) const {
    if (!annotation.v) return;
    const auto& bounds = annotation.v->gene->region_bounds;
    std::array<RegionCall*, 5> regions{
        &annotation.fwr1, &annotation.cdr1, &annotation.fwr2, &annotation.cdr2, &annotation.fwr3};
    for (std::size_t i = 0; i < regions.size(); ++i) {
        if (bounds[i * 2] < 0 || bounds[i * 2 + 1] <= bounds[i * 2]) continue;
        const auto interval = map_reference_interval(
            annotation.v->alignment,
            static_cast<std::size_t>(bounds[i * 2]),
            static_cast<std::size_t>(bounds[i * 2 + 1]));
        if (!interval || interval->second > annotation.oriented_sequence.size()) continue;
        auto& region = *regions[i];
        region.sequence = annotation.oriented_sequence.substr(
            interval->first, interval->second - interval->first);
        std::size_t region_frame = 0;
        if (annotation.sequence_frame) {
            const auto sequence_frame = static_cast<std::size_t>(*annotation.sequence_frame - 1);
            region_frame = (sequence_frame + 3 - interval->first % 3) % 3;
        }
        region.sequence_aa = translate_dna(region.sequence, region_frame);
        region.start = interval->first + 1;
        region.end = interval->second;
    }
}

std::uint32_t AnnotationEngine::orientation_seed_strength(const std::string& sequence) const {
    std::uint32_t strength = 0;
    const auto v = options_.assigner_strategy == AssignerStrategy::RiatMp
        ? database_.v_tree.roots().candidates(sequence, 1)
        : database_.v.candidates(sequence, 1);
    const auto j = database_.j.candidates(sequence, 1);
    if (!v.empty()) strength += v.front().votes;
    if (!j.empty()) strength += j.front().votes;
    return strength;
}

void AnnotationEngine::stitch_alignment(Annotation& annotation) const {
    std::vector<const SegmentHit*> hits;
    if (annotation.v) hits.push_back(&*annotation.v);
    if (annotation.d) hits.push_back(&*annotation.d);
    if (annotation.j) hits.push_back(&*annotation.j);
    if (hits.empty()) return;
    std::sort(hits.begin(), hits.end(), [](const SegmentHit* a, const SegmentHit* b) {
        return a->alignment.query_start < b->alignment.query_start;
    });
    std::size_t cursor = hits.front()->alignment.query_start;
    for (const auto* hit : hits) {
        const auto& alignment = hit->alignment;
        if (alignment.query_start < cursor) continue;
        if (alignment.query_start > cursor) {
            const auto length = alignment.query_start - cursor;
            annotation.sequence_alignment += annotation.oriented_sequence.substr(cursor, length);
            annotation.germline_alignment.append(length, 'N');
        }
        annotation.sequence_alignment += alignment.aligned_query;
        annotation.germline_alignment += alignment.aligned_reference;
        cursor = alignment.query_end;
    }
}

void AnnotationEngine::annotate_junction(Annotation& annotation) const {
    if (!annotation.v || !annotation.j) return;
    const auto& sequence = annotation.oriented_sequence;
    std::optional<std::size_t> cys;
    std::optional<std::size_t> v_reference_anchor;
    if (annotation.v->gene->region_bounds[9] >= 3) {
        const auto candidate = static_cast<std::size_t>(annotation.v->gene->region_bounds[9] - 3);
        if (is_cys_codon(annotation.v->gene->sequence, candidate)) v_reference_anchor = candidate;
    }
    if (v_reference_anchor) {
        cys = map_reference_to_query(annotation.v->alignment, *v_reference_anchor);
    }

    std::optional<std::size_t> j_anchor;
    std::optional<std::size_t> j_reference_anchor;
    if (annotation.j->gene->cdr3_stop >= 0 &&
        static_cast<std::size_t>(annotation.j->gene->cdr3_stop) + 1 <
            annotation.j->gene->sequence.size()) {
        // IgBLAST .aux stores the 0-based final CDR3 base, immediately before
        // the conserved J W/F anchor codon.
        const auto candidate = static_cast<std::size_t>(annotation.j->gene->cdr3_stop) + 1;
        if (is_w_or_f_codon(annotation.j->gene->sequence, candidate) &&
            is_gly_codon(annotation.j->gene->sequence, candidate + 3)) {
            j_reference_anchor = candidate;
        }
    }
    if (j_reference_anchor) {
        j_anchor = map_reference_to_query(annotation.j->alignment, *j_reference_anchor);
    }
    if (!cys || !j_anchor || *j_anchor <= *cys || *j_anchor + 3 > sequence.size() ||
        *j_anchor - *cys > options_.max_junction_span + 6) return;

    const bool in_frame = (*j_anchor - *cys) % 3 == 0;
    annotation.vj_in_frame = in_frame;
    annotation.junction = sequence.substr(*cys, *j_anchor + 3 - *cys);
    annotation.junction_aa = translate_dna(annotation.junction);
    if (*j_anchor > *cys + 3) {
        annotation.cdr3 = sequence.substr(*cys + 3, *j_anchor - (*cys + 3));
        annotation.cdr3_aa = translate_dna(annotation.cdr3);
        annotation.cdr3_start = *cys + 4;
        annotation.cdr3_end = *j_anchor;
    }
    if (annotation.j->alignment.query_end > *j_anchor) {
        annotation.fwr4.sequence = sequence.substr(
            *j_anchor, annotation.j->alignment.query_end - *j_anchor);
        annotation.fwr4.sequence_aa = translate_dna(annotation.fwr4.sequence);
        annotation.fwr4.start = *j_anchor + 1;
        annotation.fwr4.end = annotation.j->alignment.query_end;
    }

    const std::size_t frame = *cys % 3;
    annotation.sequence_frame = static_cast<int>(frame + 1);
    annotation.sequence_aa = translate_dna(sequence, frame);
    std::size_t coding_start = annotation.v->alignment.query_start;
    while (coding_start < *j_anchor && coding_start % 3 != frame) ++coding_start;
    const auto coding = sequence.substr(coding_start, *j_anchor + 3 - coding_start);
    const bool has_stop = translate_dna(coding).find('*') != std::string::npos;
    annotation.stop_codon = has_stop;
    annotation.v_frameshift = net_frameshift(annotation.v->alignment);
    annotation.j_frameshift = net_frameshift(annotation.j->alignment);
    annotation.productive = in_frame && !has_stop && !*annotation.v_frameshift && !*annotation.j_frameshift;
    annotation.complete_vdj = annotation.v->alignment.reference_start <= 2 &&
        annotation.j->alignment.reference_end + 3 >= annotation.j->gene->sequence.size();
    if (annotation.d && annotation.d->gene->coding_frame_start >= 0) {
        int phase = static_cast<int>(annotation.d->alignment.reference_start) -
            annotation.d->gene->coding_frame_start -
            static_cast<int>(annotation.d->alignment.query_start) + static_cast<int>(frame);
        phase %= 3;
        if (phase < 0) phase += 3;
        annotation.d_frame = phase + 1;
    }

    if (!annotation.sequence_alignment.empty() && annotation.v) {
        const auto alignment_start = annotation.v->alignment.query_start;
        std::size_t query_position = alignment_start;
        std::optional<std::size_t> cys_column;
        for (std::size_t column = 0; column < annotation.sequence_alignment.size(); ++column) {
            if (annotation.sequence_alignment[column] == '-') continue;
            if (query_position == *cys) {
                cys_column = column;
                break;
            }
            ++query_position;
        }
        if (!cys_column) return;
        const auto alignment_frame = *cys_column % 3;
        annotation.sequence_alignment_aa = translate_dna(
            annotation.sequence_alignment, alignment_frame, true);
        annotation.germline_alignment_aa = translate_dna(
            annotation.germline_alignment, alignment_frame, true);
    }
}

}  // namespace swiftig
