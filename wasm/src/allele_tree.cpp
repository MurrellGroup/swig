#include "swiftig/allele_tree.hpp"

#include "swiftig/index.hpp"

#include <algorithm>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace swiftig {
namespace {

std::size_t hamming_distance(const std::string& left, const std::string& right) {
    if (left.size() != right.size()) return std::numeric_limits<std::size_t>::max();
    std::size_t distance = 0;
    for (std::size_t position = 0; position < left.size(); ++position) {
        distance += left[position] != right[position];
    }
    return distance;
}

void saturating_add(std::uint32_t& value, std::uint32_t update) noexcept {
    value = update > std::numeric_limits<std::uint32_t>::max() - value
        ? std::numeric_limits<std::uint32_t>::max()
        : value + update;
}

}  // namespace

AlleleTreeIndex::AlleleTreeIndex() : roots_(std::make_unique<SegmentIndex>()) {}
AlleleTreeIndex::~AlleleTreeIndex() = default;
AlleleTreeIndex::AlleleTreeIndex(AlleleTreeIndex&&) noexcept = default;
AlleleTreeIndex& AlleleTreeIndex::operator=(AlleleTreeIndex&&) noexcept = default;

bool AlleleTreeIndex::empty() const noexcept { return clusters_.empty(); }
const SegmentIndex& AlleleTreeIndex::roots() const noexcept { return *roots_; }

void AlleleTreeIndex::reset(const std::vector<Gene>& genes, std::size_t cluster_radius) {
    roots_->reset({}, 0);
    clusters_.clear();
    cluster_radius_ = cluster_radius;
    edge_count_ = 0;
    mutation_count_ = 0;
    maximum_edge_mutations_ = 0;
    search_stats_ = {};
    if (genes.empty()) return;

    const std::size_t gene_count = genes.size();
    std::vector<std::uint16_t> distances(gene_count * gene_count,
        std::numeric_limits<std::uint16_t>::max());
    for (std::size_t left = 0; left < gene_count; ++left) {
        distances[left * gene_count + left] = 0;
        for (std::size_t right = left + 1; right < gene_count; ++right) {
            const auto distance = hamming_distance(genes[left].sequence, genes[right].sequence);
            const auto stored = static_cast<std::uint16_t>(std::min<std::size_t>(
                distance, std::numeric_limits<std::uint16_t>::max()));
            distances[left * gene_count + right] = stored;
            distances[right * gene_count + left] = stored;
        }
    }
    const auto distance = [&](std::size_t left, std::size_t right) {
        return static_cast<std::size_t>(distances[left * gene_count + right]);
    };

    std::map<std::size_t, std::vector<std::uint32_t>> length_groups;
    for (std::uint32_t gene_index = 0; gene_index < genes.size(); ++gene_index) {
        length_groups[genes[gene_index].sequence.size()].push_back(gene_index);
    }

    std::vector<Gene> root_genes;
    for (const auto& [unused_length, group] : length_groups) {
        static_cast<void>(unused_length);
        std::vector<bool> assigned(group.size(), false);
        std::size_t remaining = group.size();
        while (remaining != 0) {
            std::size_t best_group_position = group.size();
            std::size_t best_neighbours = 0;
            std::size_t best_distance_sum = std::numeric_limits<std::size_t>::max();
            for (std::size_t candidate_position = 0; candidate_position < group.size(); ++candidate_position) {
                if (assigned[candidate_position]) continue;
                std::size_t neighbours = 0;
                std::size_t distance_sum = 0;
                for (std::size_t other_position = 0; other_position < group.size(); ++other_position) {
                    if (assigned[other_position]) continue;
                    const auto value = distance(group[candidate_position], group[other_position]);
                    if (value > cluster_radius) continue;
                    ++neighbours;
                    distance_sum += value;
                }
                const bool better = neighbours > best_neighbours ||
                    (neighbours == best_neighbours && distance_sum < best_distance_sum) ||
                    (neighbours == best_neighbours && distance_sum == best_distance_sum &&
                     (best_group_position == group.size() ||
                      genes[group[candidate_position]].name < genes[group[best_group_position]].name));
                if (better) {
                    best_group_position = candidate_position;
                    best_neighbours = neighbours;
                    best_distance_sum = distance_sum;
                }
            }

            const auto root_gene_index = group[best_group_position];
            std::vector<std::uint32_t> members;
            for (std::size_t position = 0; position < group.size(); ++position) {
                if (assigned[position] || distance(root_gene_index, group[position]) > cluster_radius) continue;
                assigned[position] = true;
                --remaining;
                members.push_back(group[position]);
            }

            AlleleTreeCluster cluster;
            cluster.root_gene_index = root_gene_index;
            cluster.nodes.push_back(AlleleTreeNode{root_gene_index, 0, {}});
            std::vector<bool> in_tree(members.size(), false);
            const auto root_member = std::find(members.begin(), members.end(), root_gene_index);
            in_tree[static_cast<std::size_t>(root_member - members.begin())] = true;
            std::size_t tree_size = 1;
            while (tree_size < members.size()) {
                std::size_t best_parent_node = 0;
                std::size_t best_child_position = members.size();
                std::size_t best_edge = std::numeric_limits<std::size_t>::max();
                for (std::size_t parent_node = 0; parent_node < cluster.nodes.size(); ++parent_node) {
                    const auto parent_gene = cluster.nodes[parent_node].gene_index;
                    for (std::size_t child_position = 0; child_position < members.size(); ++child_position) {
                        if (in_tree[child_position]) continue;
                        const auto edge = distance(parent_gene, members[child_position]);
                        const bool better = edge < best_edge ||
                            (edge == best_edge &&
                             (best_child_position == members.size() ||
                              genes[members[child_position]].name <
                                  genes[members[best_child_position]].name));
                        if (better) {
                            best_parent_node = parent_node;
                            best_child_position = child_position;
                            best_edge = edge;
                        }
                    }
                }

                const auto parent_gene = cluster.nodes[best_parent_node].gene_index;
                const auto child_gene = members[best_child_position];
                std::vector<AlleleTreeMutation> mutations;
                mutations.reserve(best_edge);
                for (std::size_t position = 0; position < genes[parent_gene].sequence.size(); ++position) {
                    const char parent_base = genes[parent_gene].sequence[position];
                    const char child_base = genes[child_gene].sequence[position];
                    if (parent_base == child_base) continue;
                    mutations.push_back(AlleleTreeMutation{
                        static_cast<std::uint16_t>(position), parent_base, child_base});
                }
                maximum_edge_mutations_ = std::max(maximum_edge_mutations_, mutations.size());
                mutation_count_ += mutations.size();
                ++edge_count_;
                cluster.nodes.push_back(AlleleTreeNode{
                    child_gene, static_cast<std::uint32_t>(best_parent_node), std::move(mutations)});
                in_tree[best_child_position] = true;
                ++tree_size;
            }

            // Persist the DFS adjacency once. The hot RIAT-MP search used to
            // allocate a vector-of-vectors and rebuild this identical topology
            // for every read and every selected root alignment.
            for (std::size_t node_index = cluster.nodes.size(); node_index-- > 1;) {
                auto& node = cluster.nodes[node_index];
                auto& parent = cluster.nodes[node.parent];
                node.next_sibling = parent.first_child;
                parent.first_child = static_cast<std::uint32_t>(node_index);
            }

            root_genes.push_back(genes[root_gene_index]);
            clusters_.push_back(std::move(cluster));
        }
    }
    roots_->reset(std::move(root_genes), 9);
}

void AlleleTreeIndex::record_search(const AlleleTreeSearchStats& update) const noexcept {
    saturating_add(search_stats_.queries, update.queries);
    saturating_add(search_stats_.root_candidates, update.root_candidates);
    saturating_add(search_stats_.root_alignments, update.root_alignments);
    saturating_add(search_stats_.root_tracebacks, update.root_tracebacks);
    saturating_add(search_stats_.multipath_searches, update.multipath_searches);
    saturating_add(search_stats_.trace_states, update.trace_states);
    saturating_add(search_stats_.trace_limit_hits, update.trace_limit_hits);
    saturating_add(search_stats_.clusters_scored, update.clusters_scored);
    saturating_add(search_stats_.nodes_scored, update.nodes_scored);
    saturating_add(search_stats_.mutation_updates, update.mutation_updates);
    saturating_add(search_stats_.final_realignments, update.final_realignments);
}

}  // namespace swiftig
