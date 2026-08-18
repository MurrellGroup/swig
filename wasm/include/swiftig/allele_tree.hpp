#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <limits>
#include <vector>

#include "swiftig/types.hpp"

namespace swiftig {

class SegmentIndex;

struct AlleleTreeMutation {
    std::uint16_t position = 0;
    char parent_base = 'N';
    char child_base = 'N';
};

struct AlleleTreeNode {
    std::uint32_t gene_index = 0;
    std::uint32_t parent = 0;
    std::vector<AlleleTreeMutation> edge_mutations;
    std::uint32_t first_child = std::numeric_limits<std::uint32_t>::max();
    std::uint32_t next_sibling = std::numeric_limits<std::uint32_t>::max();
};

struct AlleleTreeCluster {
    std::uint32_t root_gene_index = 0;
    std::vector<AlleleTreeNode> nodes;
};

struct AlleleTreeSearchStats {
    std::uint32_t queries = 0;
    std::uint32_t root_candidates = 0;
    std::uint32_t root_alignments = 0;
    std::uint32_t root_tracebacks = 0;
    std::uint32_t multipath_searches = 0;
    std::uint32_t trace_states = 0;
    std::uint32_t trace_limit_hits = 0;
    std::uint32_t clusters_scored = 0;
    std::uint32_t nodes_scored = 0;
    std::uint32_t mutation_updates = 0;
    std::uint32_t final_realignments = 0;
};

class AlleleTreeIndex {
public:
    AlleleTreeIndex();
    ~AlleleTreeIndex();
    AlleleTreeIndex(AlleleTreeIndex&&) noexcept;
    AlleleTreeIndex& operator=(AlleleTreeIndex&&) noexcept;
    AlleleTreeIndex(const AlleleTreeIndex&) = delete;
    AlleleTreeIndex& operator=(const AlleleTreeIndex&) = delete;

    void reset(const std::vector<Gene>& genes, std::size_t cluster_radius = 12);

    [[nodiscard]] bool empty() const noexcept;
    [[nodiscard]] const SegmentIndex& roots() const noexcept;
    [[nodiscard]] const std::vector<AlleleTreeCluster>& clusters() const noexcept {
        return clusters_;
    }
    [[nodiscard]] std::size_t cluster_radius() const noexcept { return cluster_radius_; }
    [[nodiscard]] std::size_t edge_count() const noexcept { return edge_count_; }
    [[nodiscard]] std::size_t mutation_count() const noexcept { return mutation_count_; }
    [[nodiscard]] std::size_t maximum_edge_mutations() const noexcept {
        return maximum_edge_mutations_;
    }

    void record_search(const AlleleTreeSearchStats& update) const noexcept;
    void reset_search_stats() const noexcept { search_stats_ = {}; }
    [[nodiscard]] AlleleTreeSearchStats search_stats() const noexcept { return search_stats_; }

private:
    std::unique_ptr<SegmentIndex> roots_;
    std::vector<AlleleTreeCluster> clusters_;
    std::size_t cluster_radius_ = 0;
    std::size_t edge_count_ = 0;
    std::size_t mutation_count_ = 0;
    std::size_t maximum_edge_mutations_ = 0;
    mutable AlleleTreeSearchStats search_stats_;
};

}  // namespace swiftig
