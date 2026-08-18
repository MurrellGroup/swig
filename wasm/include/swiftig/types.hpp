#pragma once

#include <cstddef>
#include <cstdint>
#include <array>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace swiftig {

inline constexpr const char* kVersion = "0.6.0";

enum class Segment : std::uint8_t { V = 0, D = 1, J = 2, C = 3 };

enum class AssignerStrategy : std::uint8_t {
    Standard = 0,
    RiatMp = 1,
    Aer = 2,
};

struct SequenceRecord {
    std::string id;
    std::string description;
    std::string sequence;
    std::string quality;
};

struct Gene {
    std::string name;
    std::string sequence;
    std::string locus;
    int coding_frame_start = -1;
    int cdr3_stop = -1;
    std::string chain_type;
    std::array<int, 10> region_bounds{{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};
    std::string annotation_source;
};

struct RegionCall {
    std::string sequence;
    std::string sequence_aa;
    std::optional<std::size_t> start;
    std::optional<std::size_t> end;
};

struct Scoring {
    int match = 2;
    int mismatch = -3;
    int gap_open = -5;
    int gap_extend = -1;
};

struct Alignment {
    int score = 0;
    std::size_t query_start = 0;
    std::size_t query_end = 0;
    std::size_t reference_start = 0;
    std::size_t reference_end = 0;
    int matches = 0;
    int mismatches = 0;
    int insertions = 0;
    int deletions = 0;
    std::string aligned_query;
    std::string aligned_reference;
    std::string cigar;

    [[nodiscard]] bool valid() const noexcept {
        return score > 0 && query_end > query_start && reference_end > reference_start;
    }

    [[nodiscard]] std::size_t query_length() const noexcept {
        return query_end - query_start;
    }

    [[nodiscard]] std::size_t reference_length() const noexcept {
        return reference_end - reference_start;
    }

    [[nodiscard]] std::size_t alignment_columns() const noexcept {
        return static_cast<std::size_t>(matches + mismatches + insertions + deletions);
    }

    [[nodiscard]] double identity() const noexcept {
        const auto n = alignment_columns();
        return n == 0 ? 0.0 : static_cast<double>(matches) / static_cast<double>(n);
    }
};

struct SegmentHit {
    const Gene* gene = nullptr;
    Alignment alignment;
    std::string call;
    // Length of the query span searched to produce this hit. D is searched in
    // the V/J-bounded junction; V, J, and C use the complete oriented query.
    std::size_t search_query_length = 0;
    // Calibrated BLAST-form expectation value for this SwiftIG score. Missing
    // when a caller supplies a scoring tuple without an offline calibration.
    std::optional<double> support;

    [[nodiscard]] bool valid() const noexcept { return gene != nullptr && alignment.valid(); }
};

struct Annotation {
    std::string sequence_id;
    std::string sequence;
    std::string quality;
    std::string oriented_sequence;
    std::string sequence_aa;
    bool rev_comp = false;
    std::optional<bool> productive;
    std::optional<bool> vj_in_frame;
    std::optional<bool> stop_codon;
    std::optional<bool> complete_vdj;
    std::optional<bool> v_frameshift;
    std::optional<bool> j_frameshift;
    std::string locus;
    std::optional<SegmentHit> v;
    std::optional<SegmentHit> d;
    std::optional<SegmentHit> j;
    std::optional<SegmentHit> c;
    std::vector<SegmentHit> v_alternatives;
    std::vector<SegmentHit> d_alternatives;
    std::vector<SegmentHit> j_alternatives;
    std::vector<SegmentHit> c_alternatives;
    std::string sequence_alignment;
    std::string sequence_alignment_aa;
    std::string germline_alignment;
    std::string germline_alignment_aa;
    std::string junction;
    std::string junction_aa;
    std::string cdr3;
    std::string cdr3_aa;
    RegionCall fwr1;
    RegionCall cdr1;
    RegionCall fwr2;
    RegionCall cdr2;
    RegionCall fwr3;
    RegionCall fwr4;
    std::string np1;
    std::string np2;
    std::optional<std::size_t> cdr3_start;
    std::optional<std::size_t> cdr3_end;
    std::optional<int> d_frame;
    std::optional<int> sequence_frame;
    std::string region_definition;
    std::string v_annotation_source;
    std::string j_annotation_source;
};

struct EngineOptions {
    AssignerStrategy assigner_strategy = AssignerStrategy::Standard;
    // AER and RIAT-MP use result-equivalent allocation-light kernels by
    // default. The reference kernels remain selectable for validation.
    bool optimized_kernels = true;
    Scoring v_scoring{2, -3, -5, -1};
    // Jointly calibrated on the supplied low-SHM and IgG simulations. D uses
    // a six-base exact-support floor; the strong gap-open penalty avoids
    // explaining N-additions as internal D gaps.
    Scoring d_scoring{2, -3, -13, -1};
    // J tolerates distributed SHM but strongly penalizes a new gap and mildly
    // penalizes its extension, which favored the correct J allele across both
    // supplied mutation regimes.
    Scoring j_scoring{2, -3, -17, -2};
    Scoring c_scoring{2, -3, -5, -1};
    std::size_t top_v = 3;
    std::size_t top_d = 2;
    std::size_t top_j = 2;
    std::size_t top_c = 3;
    std::size_t min_v_length = 24;
    std::size_t min_d_match = 6;
    std::size_t min_j_length = 10;
    std::size_t min_c_length = 30;
    double min_identity = 0.60;
    int band_width = 28;
    int max_band_width = 256;
    int max_vdj_overlap = 0;
    std::size_t max_junction_span = 225;
    bool search_forward = true;
    bool search_reverse = true;
    bool allow_vdj_overlap = false;
};

struct Candidate {
    std::uint32_t gene_index = 0;
    int diagonal = std::numeric_limits<int>::max();
    std::uint16_t votes = 0;
    int diagonal_span = 0;
    bool weak_seed_signal = false;
};

std::string reverse_complement(const std::string& sequence);
std::string normalize_dna(const std::string& sequence, bool strip_gaps = false);
std::string infer_locus(const std::string& gene_name);
std::string translate_dna(const std::string& sequence, std::size_t frame = 0, bool keep_gaps = false);

}  // namespace swiftig
