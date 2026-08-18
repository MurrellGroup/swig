#pragma once

#include <optional>
#include <string>
#include <vector>

#include "swiftig/alignment.hpp"
#include "swiftig/index.hpp"
#include "swiftig/statistics.hpp"
#include "swiftig/types.hpp"

namespace swiftig {

struct OrientationHints {
    std::vector<Candidate> v;
    std::vector<Candidate> j;
};

struct AnnotationHints {
    OrientationHints forward;
    OrientationHints reverse;
};

class AnnotationEngine {
public:
    AnnotationEngine(const GermlineDatabase& database, EngineOptions options = {});
    [[nodiscard]] Annotation annotate(
        const SequenceRecord& record,
        const AnnotationHints* hints = nullptr) const;

private:
    struct OrientationResult {
        std::string oriented_sequence;
        std::optional<SegmentHit> v;
        std::optional<SegmentHit> d;
        std::optional<SegmentHit> j;
        std::optional<SegmentHit> c;
        std::vector<SegmentHit> v_alternatives;
        std::vector<SegmentHit> d_alternatives;
        std::vector<SegmentHit> j_alternatives;
        std::vector<SegmentHit> c_alternatives;
        double rank_score = 0.0;
    };

    const GermlineDatabase& database_;
    EngineOptions options_;
    ReferenceDatabaseStatistics v_statistics_;
    ReferenceDatabaseStatistics d_statistics_;
    ReferenceDatabaseStatistics j_statistics_;
    ReferenceDatabaseStatistics c_statistics_;
    mutable AlignmentWorkspace alignment_workspace_;

    [[nodiscard]] OrientationResult annotate_orientation(
        const std::string& sequence,
        const OrientationHints* hints = nullptr) const;
    [[nodiscard]] std::vector<SegmentHit> align_candidates(
        const std::string& query,
        const SegmentIndex& index,
        std::size_t top_n,
        const Scoring& scoring,
        std::size_t min_length,
        const std::string& locus_filter = "",
        const std::vector<Candidate>* candidate_hints = nullptr) const;
    [[nodiscard]] std::vector<SegmentHit> align_v_allele_tree(
        const std::string& query,
        std::size_t top_n,
        const Scoring& scoring,
        std::size_t min_length,
        const std::vector<Candidate>* candidate_hints = nullptr) const;
    void annotate_junction(Annotation& annotation) const;
    void annotate_v_regions(Annotation& annotation) const;
    void stitch_alignment(Annotation& annotation) const;
    [[nodiscard]] OrientationHints orientation_hints(const std::string& sequence) const;
    [[nodiscard]] static std::uint32_t orientation_seed_strength(
        const OrientationHints& hints) noexcept;
};

}  // namespace swiftig
