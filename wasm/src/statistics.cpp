#include "swiftig/statistics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace swiftig {
namespace {

struct Calibration {
    Scoring scoring;
    double lambda;
    double k;
    double alpha;
    double beta;
};

// Deterministic uniform-DNA simulations of the exact affine recurrence in
// alignment.cpp, fitted to P(max_score < S) = exp(-K A exp(-lambda S)). See
// scripts/calibrate-airr-support.cpp and scripts/fit-airr-support.py.
constexpr std::array<Calibration, 5> calibrations{{
    {{2, -3, -5, -1}, 0.615091626817, 0.657543879520, 2.29229884023, -9.32267423552},
    {{2, -3, -13, -1}, 0.633731436761, 0.382178144270, 0.383131230945, 0.317757738848},
    {{2, -3, -17, -2}, 0.633731436761, 0.381864932099, 0.379080052598, 0.341396551998},
    {{2, -4, -11, -1}, 0.645352027689, 0.280867618578, 0.230455684630, -3.33906776426},
    {{2, -4, -13, -1}, 0.645457954698, 0.281099152429, 0.231233105738, -3.35008024475},
}};

bool same_scoring(const Scoring& left, const Scoring& right) noexcept {
    return left.match == right.match && left.mismatch == right.mismatch &&
        left.gap_open == right.gap_open && left.gap_extend == right.gap_extend;
}

const Calibration* find_calibration(const Scoring& scoring) noexcept {
    for (const auto& calibration : calibrations) {
        if (same_scoring(scoring, calibration.scoring)) return &calibration;
    }
    return nullptr;
}

double finite_length_adjustment(
    double query_length,
    double database_length,
    double database_sequences,
    const Calibration& calibration) noexcept {
    if (query_length <= 1.0 || database_length <= database_sequences) return 0.0;
    const double upper = std::max(0.0, std::min(
        query_length - 1.0,
        (database_length - database_sequences) / database_sequences));
    double adjustment = 0.0;
    // Damped fixed-point iteration is stable for the calibrated parameter
    // range and mirrors the standard BLAST finite-size search-space equation.
    for (int iteration = 0; iteration < 64; ++iteration) {
        const double area = std::max(
            1.0,
            (query_length - adjustment) *
                (database_length - database_sequences * adjustment));
        double candidate = (calibration.alpha / calibration.lambda) *
            (std::log(calibration.k) + std::log(area)) + calibration.beta;
        candidate = std::clamp(candidate, 0.0, upper);
        if (std::abs(candidate - adjustment) < 1e-10) return candidate;
        adjustment = 0.5 * (adjustment + candidate);
    }
    return adjustment;
}

}  // namespace

ReferenceDatabaseStatistics reference_database_statistics(
    const std::vector<Gene>& genes) noexcept {
    ReferenceDatabaseStatistics result;
    result.sequence_count = genes.size();
    for (const auto& gene : genes) result.total_length += gene.sequence.size();
    return result;
}

std::optional<double> calibrated_alignment_evalue(
    int score,
    std::size_t query_length,
    const ReferenceDatabaseStatistics& database,
    const Scoring& scoring) noexcept {
    const auto* calibration = find_calibration(scoring);
    if (!calibration || score <= 0 || query_length == 0 ||
        database.sequence_count == 0 || database.total_length == 0) return std::nullopt;

    const double m = static_cast<double>(query_length);
    const double n = static_cast<double>(database.total_length);
    const double sequences = static_cast<double>(database.sequence_count);
    const double adjustment = finite_length_adjustment(m, n, sequences, *calibration);
    const double effective_query = std::max(1.0, m - adjustment);
    const double effective_database = std::max(1.0, n - sequences * adjustment);
    const double log_evalue = std::log(calibration->k) +
        std::log(effective_query) + std::log(effective_database) -
        calibration->lambda * static_cast<double>(score);
    if (log_evalue <= std::log(std::numeric_limits<double>::denorm_min())) return 0.0;
    if (log_evalue >= std::log(std::numeric_limits<double>::max())) {
        return std::numeric_limits<double>::max();
    }
    return std::exp(log_evalue);
}

}  // namespace swiftig
