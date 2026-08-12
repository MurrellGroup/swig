#include "swiftig/double_d.hpp"

#include <algorithm>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "swiftig/alignment.hpp"

namespace swiftig {
namespace {

bool canonical(const std::string& sequence) {
    return std::all_of(sequence.begin(), sequence.end(), [](char base) {
        return base == 'A' || base == 'C' || base == 'G' || base == 'T';
    });
}

bool d_locus(const std::string& locus) {
    return locus == "IGH" || locus == "TRB" || locus == "TRD";
}

bool matching_locus(const Gene& gene, const std::string& locus) {
    return gene.locus.empty() || locus.empty() || gene.locus == locus;
}

Alignment extend_seed(
    const std::string& query,
    const Gene& gene,
    std::size_t query_seed,
    std::size_t reference_seed,
    std::size_t seed_length) {
    const auto& reference = gene.sequence;
    std::size_t best_left = 0;
    int best_left_score = 0;
    int running = 0;
    const auto left_limit = std::min(query_seed, reference_seed);
    for (std::size_t step = 1; step <= left_limit; ++step) {
        running += query[query_seed - step] == reference[reference_seed - step] ? 2 : -2;
        if (running > best_left_score) {
            best_left_score = running;
            best_left = step;
        }
    }

    std::size_t best_right = 0;
    int best_right_score = 0;
    running = 0;
    const auto query_after = query_seed + seed_length;
    const auto reference_after = reference_seed + seed_length;
    const auto right_limit = std::min(query.size() - query_after, reference.size() - reference_after);
    for (std::size_t step = 1; step <= right_limit; ++step) {
        running += query[query_after + step - 1] == reference[reference_after + step - 1] ? 2 : -2;
        if (running > best_right_score) {
            best_right_score = running;
            best_right = step;
        }
    }

    Alignment alignment;
    alignment.query_start = query_seed - best_left;
    alignment.query_end = query_after + best_right;
    alignment.reference_start = reference_seed - best_left;
    alignment.reference_end = reference_after + best_right;
    alignment.aligned_query = query.substr(
        alignment.query_start, alignment.query_end - alignment.query_start);
    alignment.aligned_reference = reference.substr(
        alignment.reference_start, alignment.reference_end - alignment.reference_start);
    for (std::size_t index = 0; index < alignment.aligned_query.size(); ++index) {
        if (alignment.aligned_query[index] == alignment.aligned_reference[index]) ++alignment.matches;
        else ++alignment.mismatches;
    }
    alignment.score = alignment.matches * 2 - alignment.mismatches * 2;
    refresh_airr_cigar(alignment, query.size(), reference.size());
    return alignment;
}

int delta_distance(
    const std::string& span,
    const std::vector<Gene>& genes,
    const std::string& locus,
    std::size_t delta,
    std::size_t minimum_length) {
    if (span.size() <= delta) return -1;
    const auto length = span.size() - delta;
    if (length < minimum_length) return -1;
    int best = std::numeric_limits<int>::max();
    for (std::size_t span_start = 0; span_start + length <= span.size(); ++span_start) {
        for (const auto& gene : genes) {
            if (!matching_locus(gene, locus) || gene.sequence.size() < length) continue;
            for (std::size_t gene_start = 0; gene_start + length <= gene.sequence.size(); ++gene_start) {
                int distance = 0;
                for (std::size_t position = 0; position < length; ++position) {
                    distance += span[span_start + position] != gene.sequence[gene_start + position];
                    if (distance >= best) break;
                }
                best = std::min(best, distance);
            }
        }
    }
    return best == std::numeric_limits<int>::max() ? -1 : best;
}

std::string calls_for_hit(const SegmentHit& selected, const std::vector<SegmentHit>& hits) {
    std::vector<std::string> names;
    for (const auto& hit : hits) {
        if (hit.alignment.score != selected.alignment.score ||
            hit.alignment.query_start != selected.alignment.query_start ||
            hit.alignment.query_end != selected.alignment.query_end ||
            hit.alignment.reference_start != selected.alignment.reference_start ||
            hit.alignment.reference_end != selected.alignment.reference_end) continue;
        names.push_back(hit.gene->name);
    }
    std::sort(names.begin(), names.end());
    names.erase(std::unique(names.begin(), names.end()), names.end());
    std::string result;
    for (const auto& name : names) {
        if (!result.empty()) result += ',';
        result += name;
    }
    return result.empty() ? selected.gene->name : result;
}

std::string number(double value) {
    std::ostringstream stream;
    stream.setf(std::ios::fixed);
    stream.precision(6);
    stream << value;
    auto result = stream.str();
    while (!result.empty() && result.back() == '0') result.pop_back();
    if (!result.empty() && result.back() == '.') result.pop_back();
    return result;
}

void field(std::ostream& output, const std::string& value, bool& first) {
    if (!first) output.put('\t');
    first = false;
    output << value;
}

std::string alternatives(const std::vector<DoubleDAlternative>& values) {
    std::ostringstream output;
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index) output << ';';
        output << values[index].d_call << '+' << values[index].d2_call << '|'
               << values[index].pair_score << '|' << values[index].score_gain << '|'
               << values[index].pseudo_distance;
    }
    return output.str();
}

}  // namespace

DoubleDScreener::DoubleDScreener(
    const GermlineDatabase& database,
    DoubleDOptions options)
    : database_(database), options_(std::move(options)) {
    if (options_.seed_length < 6 || options_.seed_length > 24) {
        options_.seed_length = 11;
    }
    for (std::size_t gene_index = 0; gene_index < database_.d.genes().size(); ++gene_index) {
        const auto& sequence = database_.d.genes()[gene_index].sequence;
        if (sequence.size() < options_.seed_length) continue;
        for (std::size_t position = 0; position + options_.seed_length <= sequence.size(); ++position) {
            const auto seed = sequence.substr(position, options_.seed_length);
            if (canonical(seed)) seeds_[seed].push_back(SeedHit{gene_index, position});
        }
    }
}

std::optional<DoubleDCall> DoubleDScreener::screen(const Annotation& annotation) const {
    if (options_.mode == DoubleDMode::Off || !annotation.v || !annotation.j ||
        database_.d.empty() || !d_locus(annotation.locus)) return std::nullopt;
    const auto window_start = annotation.v->alignment.query_end;
    const auto window_end = annotation.j->alignment.query_start;
    if (window_end <= window_start || window_end > annotation.oriented_sequence.size()) {
        return std::nullopt;
    }
    const auto vj_span = window_end - window_start;
    if (options_.mode == DoubleDMode::LongSpan && vj_span < options_.minimum_vj_span) {
        return std::nullopt;
    }
    const auto query = annotation.oriented_sequence.substr(window_start, vj_span);
    if (query.size() < options_.seed_length * 2) return std::nullopt;

    std::vector<SegmentHit> hits;
    std::unordered_set<std::string> diagonals;
    for (std::size_t query_position = 0;
         query_position + options_.seed_length <= query.size(); ++query_position) {
        const auto seed = query.substr(query_position, options_.seed_length);
        if (!canonical(seed)) continue;
        const auto found = seeds_.find(seed);
        if (found == seeds_.end()) continue;
        for (const auto& seed_hit : found->second) {
            const auto& gene = database_.d.genes()[seed_hit.gene];
            if (!matching_locus(gene, annotation.locus)) continue;
            const auto diagonal = static_cast<long long>(query_position) -
                static_cast<long long>(seed_hit.position);
            const auto key = std::to_string(seed_hit.gene) + ':' + std::to_string(diagonal);
            if (!diagonals.insert(key).second) continue;
            auto alignment = extend_seed(
                query, gene, query_position, seed_hit.position, options_.seed_length);
            alignment.query_start += window_start;
            alignment.query_end += window_start;
            refresh_airr_cigar(
                alignment, annotation.oriented_sequence.size(), gene.sequence.size());
            hits.push_back(SegmentHit{&gene, std::move(alignment), gene.name});
        }
    }
    if (hits.size() < 2) return std::nullopt;
    std::sort(hits.begin(), hits.end(), [](const SegmentHit& left, const SegmentHit& right) {
        if (left.alignment.query_start != right.alignment.query_start) {
            return left.alignment.query_start < right.alignment.query_start;
        }
        if (left.alignment.query_end != right.alignment.query_end) {
            return left.alignment.query_end > right.alignment.query_end;
        }
        if (left.alignment.score != right.alignment.score) {
            return left.alignment.score > right.alignment.score;
        }
        return left.gene->name < right.gene->name;
    });

    int best_single_score = 0;
    for (const auto& hit : hits) best_single_score = std::max(best_single_score, hit.alignment.score);
    struct Pair {
        std::size_t first = 0;
        std::size_t second = 0;
        int score = 0;
        int gain = 0;
        int pseudo = -1;
    };
    std::vector<Pair> candidates;
    for (std::size_t first = 0; first < hits.size(); ++first) {
        for (std::size_t second = first + 1; second < hits.size(); ++second) {
            const auto& d = hits[first].alignment;
            const auto& d2 = hits[second].alignment;
            if (d.query_end > d2.query_start) continue;
            const auto pair_score = d.score + d2.score;
            const auto gain = pair_score - best_single_score;
            if (gain < options_.minimum_score_gain) continue;
            candidates.push_back(Pair{first, second, pair_score, gain, -1});
        }
    }
    if (candidates.empty()) return std::nullopt;
    const auto pair_order = [&](const Pair& left, const Pair& right) {
        if (left.score != right.score) return left.score > right.score;
        if (left.gain != right.gain) return left.gain > right.gain;
        const auto left_bases = hits[left.first].alignment.matches + hits[left.second].alignment.matches;
        const auto right_bases = hits[right.first].alignment.matches + hits[right.second].alignment.matches;
        if (left_bases != right_bases) return left_bases > right_bases;
        const auto left_gap = hits[left.second].alignment.query_start - hits[left.first].alignment.query_end;
        const auto right_gap = hits[right.second].alignment.query_start - hits[right.first].alignment.query_end;
        if (left_gap != right_gap) return left_gap < right_gap;
        return hits[left.first].gene->name + hits[left.second].gene->name <
            hits[right.first].gene->name + hits[right.second].gene->name;
    };
    std::sort(candidates.begin(), candidates.end(), pair_order);

    // The pseudo-tandem test is substantially more expensive than indexed seed
    // lookup. Evaluate candidates in evidence order and stop after the range
    // that can be reported as a near-tied alternative. Many seed combinations
    // share the same outer span, so cache their delta-distance as well.
    std::vector<Pair> pairs;
    std::unordered_map<std::string, int> pseudo_cache;
    for (auto candidate : candidates) {
        if (!pairs.empty() && candidate.score + 2 < pairs.front().score) break;
        const auto& d = hits[candidate.first].alignment;
        const auto& d2 = hits[candidate.second].alignment;
        const auto span_key = std::to_string(d.query_start) + ':' + std::to_string(d2.query_end);
        const auto cached = pseudo_cache.find(span_key);
        if (cached != pseudo_cache.end()) candidate.pseudo = cached->second;
        else {
            const auto span = annotation.oriented_sequence.substr(
                d.query_start, d2.query_end - d.query_start);
            candidate.pseudo = delta_distance(
                span, database_.d.genes(), annotation.locus,
                options_.pseudo_trim, options_.seed_length);
            pseudo_cache.emplace(span_key, candidate.pseudo);
        }
        if (candidate.pseudo >= 0 && candidate.pseudo <= options_.maximum_pseudo_mismatches) continue;
        pairs.push_back(candidate);
    }
    if (pairs.empty()) return std::nullopt;
    std::sort(pairs.begin(), pairs.end(), pair_order);

    const auto& selected = pairs.front();
    DoubleDCall result;
    result.d = hits[selected.first];
    result.d2 = hits[selected.second];
    result.d.call = calls_for_hit(result.d, hits);
    result.d2.call = calls_for_hit(result.d2, hits);
    result.baseline_d_call = annotation.d
        ? (annotation.d->call.empty() ? annotation.d->gene->name : annotation.d->call) : "";
    result.vj_span = vj_span;
    result.inter_d_length = result.d2.alignment.query_start - result.d.alignment.query_end;
    result.pair_score = selected.score;
    result.best_single_score = best_single_score;
    result.score_gain = selected.gain;
    result.pseudo_distance = selected.pseudo;
    for (std::size_t index = 1; index < pairs.size() && result.alternatives.size() < 12; ++index) {
        const auto& pair = pairs[index];
        if (pair.score + 2 < selected.score) break;
        result.alternatives.push_back(DoubleDAlternative{
            calls_for_hit(hits[pair.first], hits),
            calls_for_hit(hits[pair.second], hits),
            pair.score, pair.gain, pair.pseudo});
    }
    return result;
}

void write_double_d_header(std::ostream& output) {
    output << "swig_batch_record_index\tsequence_id\tlocus\tv_call\tstandard_d_call\td_call\td2_call\tj_call"
        "\td_score\td_identity\td_cigar\td2_score\td2_identity\td2_cigar"
        "\td_sequence_start\td_sequence_end\td_germline_start\td_germline_end"
        "\td2_sequence_start\td2_sequence_end\td2_germline_start\td2_germline_end"
        "\td_sequence_alignment\td_germline_alignment\td2_sequence_alignment\td2_germline_alignment"
        "\tnp2\tnp3\tnp2_length\tnp3_length"
        "\tswig_double_d_vj_span\tswig_double_d_seed_length\tswig_double_d_pair_score"
        "\tswig_double_d_best_single_score\tswig_double_d_score_gain"
        "\tswig_double_d_pseudo_distance\tswig_double_d_mode\tswig_double_d_alternatives\n";
}

void write_double_d_record(
    std::ostream& output,
    const Annotation& annotation,
    const DoubleDCall& call,
    const DoubleDOptions& options,
    std::size_t batch_record_index) {
    bool first = true;
    const auto hit_call = [](const SegmentHit& hit) {
        return hit.call.empty() ? hit.gene->name : hit.call;
    };
    const auto optional_call = [](const std::optional<SegmentHit>& hit) {
        return hit ? (hit->call.empty() ? hit->gene->name : hit->call) : std::string{};
    };
    const auto& d = call.d.alignment;
    const auto& d2 = call.d2.alignment;
    field(output, std::to_string(batch_record_index), first);
    field(output, annotation.sequence_id, first);
    field(output, annotation.locus, first);
    field(output, optional_call(annotation.v), first);
    field(output, call.baseline_d_call, first);
    field(output, hit_call(call.d), first);
    field(output, hit_call(call.d2), first);
    field(output, optional_call(annotation.j), first);
    field(output, std::to_string(d.score), first);
    field(output, number(d.identity()), first);
    field(output, d.cigar, first);
    field(output, std::to_string(d2.score), first);
    field(output, number(d2.identity()), first);
    field(output, d2.cigar, first);
    field(output, std::to_string(d.query_start + 1), first);
    field(output, std::to_string(d.query_end), first);
    field(output, std::to_string(d.reference_start + 1), first);
    field(output, std::to_string(d.reference_end), first);
    field(output, std::to_string(d2.query_start + 1), first);
    field(output, std::to_string(d2.query_end), first);
    field(output, std::to_string(d2.reference_start + 1), first);
    field(output, std::to_string(d2.reference_end), first);
    field(output, d.aligned_query, first);
    field(output, d.aligned_reference, first);
    field(output, d2.aligned_query, first);
    field(output, d2.aligned_reference, first);
    const auto np2 = annotation.oriented_sequence.substr(
        d.query_end, d2.query_start - d.query_end);
    const auto j_start = annotation.j ? annotation.j->alignment.query_start : d2.query_end;
    const auto np3 = j_start > d2.query_end
        ? annotation.oriented_sequence.substr(d2.query_end, j_start - d2.query_end) : std::string{};
    field(output, np2, first);
    field(output, np3, first);
    field(output, std::to_string(np2.size()), first);
    field(output, std::to_string(np3.size()), first);
    field(output, std::to_string(call.vj_span), first);
    field(output, std::to_string(options.seed_length), first);
    field(output, std::to_string(call.pair_score), first);
    field(output, std::to_string(call.best_single_score), first);
    field(output, std::to_string(call.score_gain), first);
    field(output, call.pseudo_distance < 0 ? "" : std::to_string(call.pseudo_distance), first);
    field(output, options.mode == DoubleDMode::All ? "all" : "long_span", first);
    field(output, alternatives(call.alternatives), first);
    output.put('\n');
}

}  // namespace swiftig
