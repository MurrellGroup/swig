#pragma once

#include <string>

#include "swiftig/types.hpp"

namespace swiftig {

Alignment local_align_affine(
    const std::string& query,
    const std::string& reference,
    const Scoring& scoring,
    int estimated_diagonal = 0,
    int band_width = -1);

void refresh_airr_cigar(
    Alignment& alignment,
    std::size_t query_size,
    std::size_t reference_size);

std::optional<std::size_t> map_reference_to_query(
    const Alignment& alignment,
    std::size_t reference_position);

}  // namespace swiftig
