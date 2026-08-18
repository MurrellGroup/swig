#pragma once

#include <array>
#include <string>
#include <vector>

#include "swiftig/types.hpp"

namespace swiftig {

struct TracebackSearchStats {
    std::size_t expanded_states = 0;
    std::size_t distinct_paths = 0;
    bool state_limit_hit = false;
};

// Scratch storage for the optimized affine kernel. One workspace is reused by
// an AnnotationEngine for an entire batch, avoiding hundreds of thousands of
// trace/score-vector allocations. It contains no result state.
struct AlignmentWorkspace {
    std::vector<std::uint8_t> trace;
    std::vector<int> previous;
    std::vector<int> current;
    std::vector<int> insertion_previous;
    std::vector<int> insertion_current;
    std::vector<char> operations;
    using ScoreVector4 = int __attribute__((vector_size(16)));
    std::vector<ScoreVector4> batch_previous;
    std::vector<ScoreVector4> batch_current;
    std::vector<ScoreVector4> batch_insertion_previous;
    std::vector<ScoreVector4> batch_insertion_current;
    std::vector<std::uint32_t> batch_trace;
};

struct AlignmentScoreRequest {
    const std::string* reference = nullptr;
    int estimated_diagonal = 0;
    int band_width = -1;
};

Alignment local_align_affine(
    const std::string& query,
    const std::string& reference,
    const Scoring& scoring,
    int estimated_diagonal = 0,
    int band_width = -1);

// Result-equivalent allocation-light kernel. The original function above is
// retained as the reference implementation and remains the Standard path.
Alignment local_align_affine_fast(
    const std::string& query,
    const std::string& reference,
    const Scoring& scoring,
    int estimated_diagonal,
    int band_width,
    AlignmentWorkspace& workspace);

// Compute four independent affine Smith-Waterman scores in SIMD lanes. This
// is used only to prove which candidates can reach the output cutoff; selected
// candidates still receive a full reference-equivalent traceback.
std::array<int, 4> local_align_affine_scores4(
    const std::string& query,
    const std::array<AlignmentScoreRequest, 4>& requests,
    std::size_t request_count,
    const Scoring& scoring,
    AlignmentWorkspace& workspace);

// Full SIMD DP plus independent deterministic tracebacks for up to four
// candidates. Every lane is byte-for-byte equivalent to the scalar kernel.
std::array<Alignment, 4> local_align_affine_fast4(
    const std::string& query,
    const std::array<AlignmentScoreRequest, 4>& requests,
    std::size_t request_count,
    const Scoring& scoring,
    AlignmentWorkspace& workspace);

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
