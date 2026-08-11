#include "swiftig/engine.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "swiftig/alignment.hpp"

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
    : database_(database), options_(std::move(options)) {}

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
        const auto refinement_limit = std::max<std::size_t>(top_n, 3);
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
        hits.push_back(SegmentHit{&gene, std::move(alignment), gene.name});
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

AnnotationEngine::OrientationResult AnnotationEngine::annotate_orientation(
    const std::string& sequence,
    const OrientationHints* hints) const {
    OrientationResult result;
    result.oriented_sequence = sequence;
    const auto v_hits = align_candidates(
        sequence, database_.v, options_.top_v, options_.v_scoring, options_.min_v_length,
        "", hints ? &hints->v : nullptr);
    const auto j_hits = align_candidates(
        sequence, database_.j, options_.top_j, options_.j_scoring, options_.min_j_length,
        "", hints ? &hints->j : nullptr);

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
    if (annotation.v && !annotation.v->gene->locus.empty()) annotation.locus = annotation.v->gene->locus;
    else if (annotation.j) annotation.locus = annotation.j->gene->locus;

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
    annotate_v_regions(annotation);
    annotate_junction(annotation);
    if (annotation.v && annotation.sequence_aa.empty()) {
        const auto frame = (annotation.v->alignment.query_start + 3 -
            (annotation.v->alignment.reference_start % 3)) % 3;
        annotation.sequence_aa = translate_dna(annotation.oriented_sequence, frame);
        annotation.v_frameshift = net_frameshift(annotation.v->alignment);
    }
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
        region.sequence_aa = translate_dna(region.sequence);
        region.start = interval->first + 1;
        region.end = interval->second;
    }
}

std::uint32_t AnnotationEngine::orientation_seed_strength(const std::string& sequence) const {
    std::uint32_t strength = 0;
    const auto v = database_.v.candidates(sequence, 1);
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
        v_reference_anchor = static_cast<std::size_t>(annotation.v->gene->region_bounds[9] - 3);
    } else {
        v_reference_anchor = v_anchor_in_reference(*annotation.v->gene);
    }
    if (v_reference_anchor) {
        cys = map_reference_to_query(annotation.v->alignment, *v_reference_anchor);
    }
    if (!cys) cys = fallback_v_anchor(sequence, annotation.v->alignment);

    std::optional<std::size_t> j_anchor;
    std::optional<std::size_t> j_reference_anchor;
    if (annotation.j->gene->cdr3_stop >= 0 &&
        static_cast<std::size_t>(annotation.j->gene->cdr3_stop) + 1 <
            annotation.j->gene->sequence.size()) {
        // IgBLAST .aux stores the 0-based final CDR3 base, immediately before
        // the conserved J W/F anchor codon.
        j_reference_anchor = static_cast<std::size_t>(annotation.j->gene->cdr3_stop) + 1;
    } else {
        j_reference_anchor = j_anchor_in_reference(annotation.j->gene->sequence);
    }
    if (j_reference_anchor) {
        j_anchor = map_reference_to_query(annotation.j->alignment, *j_reference_anchor);
    }
    if (!j_anchor) j_anchor = fallback_j_anchor(sequence, annotation.j->alignment, cys);
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
        const auto alignment_frame = (*cys - alignment_start) % 3;
        annotation.sequence_alignment_aa = translate_dna(
            annotation.sequence_alignment, alignment_frame, true);
        annotation.germline_alignment_aa = translate_dna(
            annotation.germline_alignment, alignment_frame, true);
    }
}

}  // namespace swiftig
