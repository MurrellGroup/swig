#pragma once

#include <cstddef>
#include <optional>
#include <vector>

#include "swiftig/types.hpp"

namespace swiftig {

struct ReferenceDatabaseStatistics {
    std::size_t sequence_count = 0;
    std::size_t total_length = 0;
};

[[nodiscard]] ReferenceDatabaseStatistics reference_database_statistics(
    const std::vector<Gene>& genes) noexcept;

// Returns the expected number of chance local alignments scoring at least
// `score` under a uniform independent-nucleotide null. Constants are fitted
// offline for each scoring tuple shipped by SwiftIG. An uncalibrated tuple
// returns nullopt rather than emitting a misleading number.
[[nodiscard]] std::optional<double> calibrated_alignment_evalue(
    int score,
    std::size_t query_length,
    const ReferenceDatabaseStatistics& database,
    const Scoring& scoring) noexcept;

}  // namespace swiftig
