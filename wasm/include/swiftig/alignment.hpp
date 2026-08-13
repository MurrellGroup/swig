#pragma once

#include <string>
#include <vector>

#include "swiftig/types.hpp"

namespace swiftig {

struct TracebackSearchStats {
    std::size_t expanded_states = 0;
    std::size_t distinct_paths = 0;
    bool state_limit_hit = false;
};

Alignment local_align_affine(
    const std::string& query,
    const std::string& reference,
    const Scoring& scoring,
    int estimated_diagonal = 0,
    int band_width = -1);

// Enumerate distinct, near-optimal traceback geometries from one affine local
// alignment matrix. `score_tolerance` is the maximum raw-score loss from the
// optimal root alignment. This does not align any descendant reference.
std::vector<Alignment> local_align_affine_paths(
    const std::string& query,
    const std::string& reference,
    const Scoring& scoring,
    int estimated_diagonal,
    int band_width,
    std::size_t maximum_paths,
    int score_tolerance,
    std::size_t maximum_trace_states = 8192,
    TracebackSearchStats* search_stats = nullptr);

void refresh_airr_cigar(
    Alignment& alignment,
    std::size_t query_size,
    std::size_t reference_size);

std::optional<std::size_t> map_reference_to_query(
    const Alignment& alignment,
    std::size_t reference_position);

}  // namespace swiftig
