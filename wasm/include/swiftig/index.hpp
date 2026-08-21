#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "swiftig/types.hpp"
#include "swiftig/allele_tree.hpp"

namespace swiftig {

class SegmentIndex {
public:
    SegmentIndex() = default;
    SegmentIndex(std::vector<Gene> genes, std::uint32_t kmer_size);

    void reset(std::vector<Gene> genes, std::uint32_t kmer_size);
    [[nodiscard]] std::vector<Candidate> candidates(
        const std::string& query,
        std::size_t limit,
        std::size_t max_seed_frequency = 96,
        std::size_t stride = 1,
        bool force_fallback = false) const;
    // Allocation-light, result-equivalent candidate search used by the AER
    // and RIAT-MP production kernels. `candidates` above is deliberately kept
    // as the reference implementation for equivalence testing.
    [[nodiscard]] std::vector<Candidate> candidates_fast(
        const std::string& query,
        std::size_t limit,
        std::size_t max_seed_frequency = 96,
        std::size_t stride = 1,
        bool force_fallback = false) const;
    [[nodiscard]] const std::vector<Gene>& genes() const noexcept { return genes_; }
    [[nodiscard]] std::vector<Gene>& mutable_genes() noexcept { return genes_; }
    [[nodiscard]] std::uint32_t kmer_size() const noexcept { return kmer_size_; }
    [[nodiscard]] bool empty() const noexcept { return genes_.empty(); }

private:
    struct SeedHit {
        std::uint32_t gene = 0;
        std::uint32_t position = 0;
    };

    std::vector<Gene> genes_;
    std::uint32_t kmer_size_ = 0;
    std::uint32_t fallback_kmer_size_ = 0;
    std::size_t maximum_gene_length_ = 0;
    // The unordered maps are retained for the reference candidate kernel.
    // Production AER/RIAT lookups use compact direct-address CSR tables when
    // k is small enough (all bundled indexes currently use k <= 9).
    std::unordered_map<std::uint64_t, std::vector<SeedHit>> seeds_;
    std::unordered_map<std::uint64_t, std::vector<SeedHit>> fallback_seeds_;
    std::vector<std::uint32_t> fast_seed_offsets_;
    std::vector<SeedHit> fast_seed_hits_;
    std::vector<std::uint32_t> fast_fallback_seed_offsets_;
    std::vector<SeedHit> fast_fallback_seed_hits_;
    void rebuild();
};

class GermlineDatabase {
public:
    SegmentIndex v;
    AlleleTreeIndex v_tree;
    SegmentIndex d;
    SegmentIndex j;
    SegmentIndex c;

    static GermlineDatabase from_fastas(
        const std::string& v_path,
        const std::string& d_path,
        const std::string& j_path,
        const std::string& c_path = "");
    static GermlineDatabase load(const std::string& path);
    void save(const std::string& path) const;
    std::size_t apply_auxiliary(const std::string& path);
    std::size_t apply_internal_data(const std::string& path);
    std::size_t apply_d_frame_data(const std::string& path);
    [[nodiscard]] std::size_t gene_count() const noexcept;
};

}  // namespace swiftig
