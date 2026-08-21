#include "swiftig/index.hpp"

#include <algorithm>
#include <array>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

#include "swiftig/fasta.hpp"

namespace swiftig {
namespace {

int encode_base(char base) {
    switch (base) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

template <typename Fn>
void for_each_kmer(const std::string& sequence, std::uint32_t k, std::size_t stride, Fn&& fn) {
    if (k == 0 || k > 31 || sequence.size() < k) return;
    const std::uint64_t mask = (std::uint64_t{1} << (2 * k)) - 1;
    std::uint64_t code = 0;
    std::uint32_t valid = 0;
    for (std::size_t i = 0; i < sequence.size(); ++i) {
        const int base = encode_base(sequence[i]);
        if (base < 0) {
            code = 0;
            valid = 0;
            continue;
        }
        code = ((code << 2U) | static_cast<std::uint64_t>(base)) & mask;
        if (valid < k) ++valid;
        if (valid == k) {
            const std::size_t position = i + 1 - k;
            if (position % std::max<std::size_t>(stride, 1) == 0) fn(code, position);
        }
    }
}

int floor_div(int value, int divisor) {
    return value >= 0 ? value / divisor : -((-value + divisor - 1) / divisor);
}

// A query touches thousands of (gene, diagonal-bin) pairs, but the old
// node-based unordered_map allocated and chased pointers for every one of
// them. This reusable open-addressed table keeps keys, votes, and occupancy
// stamps in contiguous arrays. Generation stamps make reset O(touched), with
// no clearing or allocation after the largest query seen by a worker.
class FlatVoteTable {
public:
    void reset(
        std::size_t expected_entries,
        std::size_t gene_count,
        int minimum_bin,
        int maximum_bin) {
        const auto bin_count = maximum_bin >= minimum_bin
            ? static_cast<std::size_t>(maximum_bin - minimum_bin + 1) : 0;
        const bool dense_allowed = gene_count != 0 && bin_count != 0 &&
            gene_count <= 2'000'000 / bin_count;
        dense_ = dense_allowed;
        if (dense_) {
            dense_minimum_bin_ = minimum_bin;
            dense_bin_count_ = bin_count;
            const auto cells = gene_count * bin_count;
            if (cells > dense_votes_.size()) {
                dense_votes_.resize(cells);
                dense_stamps_.assign(cells, 0);
                generation_ = 1;
            } else if (++generation_ == 0) {
                std::fill(dense_stamps_.begin(), dense_stamps_.end(), 0);
                generation_ = 1;
            }
            touched_.clear();
            return;
        }
        std::size_t requested = 64;
        const auto target = std::max<std::size_t>(expected_entries * 2, 64);
        while (requested < target) requested <<= 1U;
        if (requested > keys_.size()) {
            keys_.resize(requested);
            votes_.resize(requested);
            stamps_.assign(requested, 0);
            generation_ = 1;
        } else if (++generation_ == 0) {
            std::fill(stamps_.begin(), stamps_.end(), 0);
            generation_ = 1;
        }
        touched_.clear();
        mask_ = keys_.size() - 1;
    }

    void add(std::uint64_t key, std::uint32_t weight) {
        if (dense_) {
            const auto gene = static_cast<std::uint32_t>(key >> 32U);
            const auto bin = static_cast<std::int32_t>(key & 0xffffffffU);
            const auto relative_bin = static_cast<long long>(bin) - dense_minimum_bin_;
            if (relative_bin < 0 ||
                relative_bin >= static_cast<long long>(dense_bin_count_)) return;
            const auto slot = static_cast<std::size_t>(gene) * dense_bin_count_ +
                static_cast<std::size_t>(relative_bin);
            if (slot >= dense_votes_.size()) return;
            if (dense_stamps_[slot] != generation_) {
                dense_stamps_[slot] = generation_;
                dense_votes_[slot] = static_cast<std::uint16_t>(weight);
                touched_.push_back(static_cast<std::uint32_t>(slot));
            } else {
                dense_votes_[slot] = static_cast<std::uint16_t>(
                    std::min<std::uint32_t>(
                        std::numeric_limits<std::uint16_t>::max(),
                        static_cast<std::uint32_t>(dense_votes_[slot]) + weight));
            }
            return;
        }
        // SplitMix64's finalizer gives good low bits for the power-of-two mask
        // even though gene ids occupy the high half of these compound keys.
        std::uint64_t mixed = key;
        mixed ^= mixed >> 30U;
        mixed *= 0xbf58476d1ce4e5b9ULL;
        mixed ^= mixed >> 27U;
        mixed *= 0x94d049bb133111ebULL;
        mixed ^= mixed >> 31U;
        std::size_t slot = static_cast<std::size_t>(mixed) & mask_;
        while (stamps_[slot] == generation_ && keys_[slot] != key) {
            slot = (slot + 1) & mask_;
        }
        if (stamps_[slot] != generation_) {
            if ((touched_.size() + 1) * 2 >= keys_.size()) {
                grow();
                add(key, weight);
                return;
            }
            stamps_[slot] = generation_;
            keys_[slot] = key;
            votes_[slot] = weight;
            touched_.push_back(static_cast<std::uint32_t>(slot));
        } else {
            votes_[slot] = std::min<std::uint32_t>(
                std::numeric_limits<std::uint16_t>::max(), votes_[slot] + weight);
        }
    }

    template <typename Fn>
    void for_each(Fn&& fn) const {
        if (dense_) {
            for (const auto slot : touched_) {
                const auto gene = static_cast<std::uint32_t>(slot / dense_bin_count_);
                const auto bin = static_cast<std::int32_t>(
                    slot % dense_bin_count_) + dense_minimum_bin_;
                const std::uint64_t key = (static_cast<std::uint64_t>(gene) << 32U) |
                    static_cast<std::uint32_t>(bin);
                fn(key, static_cast<std::uint32_t>(dense_votes_[slot]));
            }
            return;
        }
        for (const auto slot : touched_) fn(keys_[slot], votes_[slot]);
    }

private:
    void grow() {
        std::vector<std::pair<std::uint64_t, std::uint32_t>> entries;
        entries.reserve(touched_.size());
        for (const auto slot : touched_) entries.emplace_back(keys_[slot], votes_[slot]);
        const auto new_size = std::max<std::size_t>(64, keys_.size() * 2);
        keys_.resize(new_size);
        votes_.resize(new_size);
        stamps_.assign(new_size, 0);
        touched_.clear();
        generation_ = 1;
        mask_ = new_size - 1;
        for (const auto& [key, votes] : entries) add(key, votes);
    }

    std::vector<std::uint64_t> keys_;
    std::vector<std::uint32_t> votes_;
    std::vector<std::uint32_t> stamps_;
    std::vector<std::uint32_t> touched_;
    std::vector<std::uint16_t> dense_votes_;
    std::vector<std::uint32_t> dense_stamps_;
    std::size_t dense_bin_count_ = 0;
    int dense_minimum_bin_ = 0;
    bool dense_ = false;
    std::uint32_t generation_ = 0;
    std::size_t mask_ = 0;
};

FlatVoteTable& vote_workspace() {
#if defined(__wasm__)
    // Each web/CLI worker owns an independent single-threaded WASM instance.
    // Avoid pulling in unsupported TLS destructors for that target.
    static FlatVoteTable table;
#else
    static thread_local FlatVoteTable table;
#endif
    return table;
}

// Reuse per-gene reduction storage just as the vote table is reused.  A
// generation stamp means each query initializes only genes that actually
// received a seed vote; the common path performs no heap allocation or
// O(database-size) clearing after warm-up.
class CandidateWorkspace {
public:
    void reset(std::size_t gene_count) {
        if (gene_count > stamps_.size()) {
            stamps_.assign(gene_count, 0);
            best_votes_.resize(gene_count);
            best_diagonal_.resize(gene_count);
            diagonal_span_.resize(gene_count);
            generation_ = 1;
        } else if (++generation_ == 0) {
            std::fill(stamps_.begin(), stamps_.end(), 0);
            generation_ = 1;
        }
        touched_.clear();
    }

    void consider(std::uint32_t gene, std::uint32_t votes, int diagonal) {
        if (gene >= stamps_.size()) return;
        if (stamps_[gene] != generation_) {
            stamps_[gene] = generation_;
            best_votes_[gene] = votes;
            best_diagonal_[gene] = diagonal;
            diagonal_span_[gene] = 0;
            touched_.push_back(gene);
            return;
        }
        if (votes > best_votes_[gene] ||
            (votes == best_votes_[gene] && diagonal < best_diagonal_[gene])) {
            best_votes_[gene] = votes;
            best_diagonal_[gene] = diagonal;
        }
    }

    [[nodiscard]] bool active(std::uint32_t gene) const {
        return gene < stamps_.size() && stamps_[gene] == generation_;
    }
    [[nodiscard]] std::uint32_t votes(std::uint32_t gene) const {
        return best_votes_[gene];
    }
    [[nodiscard]] int diagonal(std::uint32_t gene) const {
        return best_diagonal_[gene];
    }
    void include_diagonal(std::uint32_t gene, int diagonal) {
        const long long difference =
            static_cast<long long>(diagonal) - best_diagonal_[gene];
        const auto distance = static_cast<int>(std::min<long long>(
            std::abs(difference), std::numeric_limits<int>::max()));
        diagonal_span_[gene] = std::max(diagonal_span_[gene], distance);
    }
    [[nodiscard]] int span(std::uint32_t gene) const { return diagonal_span_[gene]; }
    [[nodiscard]] const std::vector<std::uint32_t>& touched() const { return touched_; }

private:
    std::vector<std::uint32_t> stamps_;
    std::vector<std::uint32_t> best_votes_;
    std::vector<int> best_diagonal_;
    std::vector<int> diagonal_span_;
    std::vector<std::uint32_t> touched_;
    std::uint32_t generation_ = 0;
};

CandidateWorkspace& candidate_workspace() {
#if defined(__wasm__)
    static CandidateWorkspace workspace;
#else
    static thread_local CandidateWorkspace workspace;
#endif
    return workspace;
}

#ifndef SWIG_WEB

template <typename T>
void write_pod(std::ostream& output, const T& value) {
    output.write(reinterpret_cast<const char*>(&value), sizeof(T));
    if (!output) throw std::runtime_error("failed while writing index");
}

template <typename T>
T read_pod(std::istream& input) {
    T value{};
    input.read(reinterpret_cast<char*>(&value), sizeof(T));
    if (!input) throw std::runtime_error("truncated SwiftIG index");
    return value;
}

void write_string(std::ostream& output, const std::string& value) {
    if (value.size() > std::numeric_limits<std::uint32_t>::max()) {
        throw std::runtime_error("index string is too large");
    }
    const auto size = static_cast<std::uint32_t>(value.size());
    write_pod(output, size);
    output.write(value.data(), static_cast<std::streamsize>(value.size()));
    if (!output) throw std::runtime_error("failed while writing index string");
}

std::string read_string(std::istream& input) {
    const auto size = read_pod<std::uint32_t>(input);
    if (size > 16U * 1024U * 1024U) throw std::runtime_error("invalid string size in SwiftIG index");
    std::string value(size, '\0');
    input.read(value.data(), static_cast<std::streamsize>(size));
    if (!input) throw std::runtime_error("truncated SwiftIG index string");
    return value;
}

void write_segment(std::ostream& output, const SegmentIndex& segment) {
    write_pod(output, segment.kmer_size());
    const auto count = static_cast<std::uint32_t>(segment.genes().size());
    write_pod(output, count);
    for (const auto& gene : segment.genes()) {
        write_string(output, gene.name);
        write_string(output, gene.sequence);
        write_string(output, gene.locus);
        write_pod<std::int32_t>(output, static_cast<std::int32_t>(gene.coding_frame_start));
        write_pod<std::int32_t>(output, static_cast<std::int32_t>(gene.cdr3_stop));
        write_string(output, gene.chain_type);
        for (const int bound : gene.region_bounds) {
            write_pod<std::int32_t>(output, static_cast<std::int32_t>(bound));
        }
    }
}

SegmentIndex read_segment(std::istream& input, std::uint32_t version) {
    const auto k = read_pod<std::uint32_t>(input);
    const auto count = read_pod<std::uint32_t>(input);
    if (k > 31 || count > 10'000'000U) throw std::runtime_error("invalid SwiftIG index metadata");
    std::vector<Gene> genes;
    genes.reserve(count);
    for (std::uint32_t i = 0; i < count; ++i) {
        Gene gene{read_string(input), read_string(input), read_string(input), -1, -1, {}};
        if (version >= 2U) {
            gene.coding_frame_start = read_pod<std::int32_t>(input);
            gene.cdr3_stop = read_pod<std::int32_t>(input);
            gene.chain_type = read_string(input);
        }
        if (version >= 3U) {
            for (auto& bound : gene.region_bounds) bound = read_pod<std::int32_t>(input);
        }
        genes.push_back(std::move(gene));
    }
    return SegmentIndex(std::move(genes), k);
}

#endif

}  // namespace

SegmentIndex::SegmentIndex(std::vector<Gene> genes, std::uint32_t kmer_size) {
    reset(std::move(genes), kmer_size);
}

void SegmentIndex::reset(std::vector<Gene> genes, std::uint32_t kmer_size) {
    if (kmer_size > 31) {
        genes.clear();
        kmer_size = 0;
    }
    genes_ = std::move(genes);
    kmer_size_ = genes_.empty() ? 0 : kmer_size;
    rebuild();
}

void SegmentIndex::rebuild() {
    seeds_.clear();
    fallback_seeds_.clear();
    fast_seed_offsets_.clear();
    fast_seed_hits_.clear();
    fast_fallback_seed_offsets_.clear();
    fast_fallback_seed_hits_.clear();
    if (kmer_size_ == 0) return;
    fallback_kmer_size_ = kmer_size_ > 4 ? kmer_size_ - 2 : 0;
    std::size_t total_bases = 0;
    maximum_gene_length_ = 0;
    for (const auto& gene : genes_) {
        total_bases += gene.sequence.size();
        maximum_gene_length_ = std::max(maximum_gene_length_, gene.sequence.size());
    }
    seeds_.reserve(total_bases / 3 + 1);
    fallback_seeds_.reserve(total_bases / 3 + 1);
    for (std::uint32_t gene_index = 0; gene_index < genes_.size(); ++gene_index) {
        for_each_kmer(genes_[gene_index].sequence, kmer_size_, 1,
            [&](std::uint64_t code, std::size_t position) {
                seeds_[code].push_back(SeedHit{gene_index, static_cast<std::uint32_t>(position)});
            });
        if (fallback_kmer_size_ != 0) {
            for_each_kmer(genes_[gene_index].sequence, fallback_kmer_size_, 1,
                [&](std::uint64_t code, std::size_t position) {
                    fallback_seeds_[code].push_back(
                        SeedHit{gene_index, static_cast<std::uint32_t>(position)});
                });
        }
    }

    // Every valid DNA k-mer is already an integer in [0, 4^k).  For the
    // production k values (<= 9), a direct-address offset table is both much
    // smaller than a dense table of hit vectors and avoids one hash-table
    // lookup, a pointer chase, and a branch for every query k-mer.  Hits live
    // contiguously in CSR order.  Keep the maps above intact for the reference
    // implementation and fall back to them if a caller constructs a very
    // large-k index.
    const auto build_fast_seed_table = [](
        const auto& source,
        std::uint32_t k,
        std::vector<std::uint32_t>& offsets,
        std::vector<SeedHit>& hits) {
        constexpr std::uint64_t maximum_codes = 1ULL << 20U;
        if (k == 0 || k > 31) return;
        const auto code_count = 1ULL << (2U * k);
        if (code_count > maximum_codes) return;
        std::size_t total_hits = 0;
        for (const auto& [unused, bucket] : source) {
            static_cast<void>(unused);
            total_hits += bucket.size();
        }
        if (total_hits > std::numeric_limits<std::uint32_t>::max()) return;
        offsets.assign(static_cast<std::size_t>(code_count) + 1, 0);
        for (const auto& [code, bucket] : source) {
            if (code >= code_count) continue;
            offsets[static_cast<std::size_t>(code) + 1] =
                static_cast<std::uint32_t>(bucket.size());
        }
        for (std::size_t i = 1; i < offsets.size(); ++i) offsets[i] += offsets[i - 1];
        hits.resize(offsets.back());
        for (const auto& [code, bucket] : source) {
            if (code >= code_count) continue;
            std::copy(bucket.begin(), bucket.end(),
                hits.begin() + offsets[static_cast<std::size_t>(code)]);
        }
    };
    build_fast_seed_table(seeds_, kmer_size_, fast_seed_offsets_, fast_seed_hits_);
    build_fast_seed_table(
        fallback_seeds_, fallback_kmer_size_,
        fast_fallback_seed_offsets_, fast_fallback_seed_hits_);
}

std::vector<Candidate> SegmentIndex::candidates(
    const std::string& query,
    std::size_t limit,
    std::size_t max_seed_frequency,
    std::size_t stride,
    bool force_fallback) const {
    if (limit == 0 || kmer_size_ == 0 || query.size() < kmer_size_) return {};
    constexpr int bin_width = 4;
    std::unordered_map<std::uint64_t, std::uint32_t> vote_bins;
    vote_bins.reserve(query.size() * 2);
    const auto collect = [&](const auto& seed_index, std::uint32_t k, std::uint32_t weight,
                             std::size_t frequency_limit) {
        if (k == 0) return;
        for_each_kmer(query, k, stride, [&](std::uint64_t code, std::size_t query_position) {
            const auto found = seed_index.find(code);
            if (found == seed_index.end() || found->second.size() > frequency_limit) return;
            for (const auto& hit : found->second) {
                const int diagonal = static_cast<int>(query_position) - static_cast<int>(hit.position);
                const int bin = floor_div(diagonal, bin_width);
                const std::uint64_t key = (static_cast<std::uint64_t>(hit.gene) << 32U) |
                    static_cast<std::uint32_t>(bin);
                auto& count = vote_bins[key];
                count = std::min<std::uint32_t>(
                    std::numeric_limits<std::uint16_t>::max(), count + weight);
            }
        });
    };
    collect(seeds_, kmer_size_, 4, max_seed_frequency);
    std::uint32_t strongest_primary_bin = 0;
    for (const auto& [unused, votes] : vote_bins) {
        static_cast<void>(unused);
        strongest_primary_bin = std::max(strongest_primary_bin, votes);
    }
    // Shorter seeds rescue high-SHM and deliberately primary-seedless reads,
    // but traversing their much larger hit lists on an ordinary exact/near-
    // exact query is pure overhead. Primary votes are weighted by four, so
    // 16*k is still well below the normal signal from a 120 nt V amplicon,
    // while catching adversarial reads whose primary seeds match only decoys.
    const auto weak_primary_threshold = std::max<std::uint32_t>(
        32U, kmer_size_ * 16U);
    const bool weak_primary_signal = strongest_primary_bin < weak_primary_threshold;
    const bool use_fallback = force_fallback || weak_primary_signal;
    if (use_fallback) {
        collect(fallback_seeds_, fallback_kmer_size_, 1, max_seed_frequency * 4);
    }

    std::vector<Candidate> ranked;
    ranked.reserve(vote_bins.size());
    for (const auto& [key, votes] : vote_bins) {
        const auto gene = static_cast<std::uint32_t>(key >> 32U);
        const auto bin = static_cast<std::int32_t>(key & 0xffffffffU);
        ranked.push_back(Candidate{
            gene, static_cast<int>(bin) * bin_width,
            static_cast<std::uint16_t>(std::min<std::uint32_t>(
                votes, std::numeric_limits<std::uint16_t>::max())), 0,
            use_fallback});
    }
    std::sort(ranked.begin(), ranked.end(), [](const Candidate& a, const Candidate& b) {
        if (a.votes != b.votes) return a.votes > b.votes;
        if (a.gene_index != b.gene_index) return a.gene_index < b.gene_index;
        return a.diagonal < b.diagonal;
    });

    std::vector<Candidate> result;
    result.reserve(std::min(limit, genes_.size()));
    std::unordered_set<std::uint32_t> used_genes;
    for (auto candidate : ranked) {
        if (!used_genes.insert(candidate.gene_index).second) continue;
        const auto strong_threshold = std::max<std::uint16_t>(4, candidate.votes / 3);
        for (const auto& other : ranked) {
            if (other.gene_index != candidate.gene_index || other.votes < strong_threshold) continue;
            const long long difference = static_cast<long long>(other.diagonal) - candidate.diagonal;
            const auto distance = static_cast<int>(std::min<long long>(
                std::abs(difference), std::numeric_limits<int>::max()));
            candidate.diagonal_span = std::max(candidate.diagonal_span, distance);
        }
        result.push_back(std::move(candidate));
        if (result.size() == limit) break;
    }
    return result;
}

std::vector<Candidate> SegmentIndex::candidates_fast(
    const std::string& query,
    std::size_t limit,
    std::size_t max_seed_frequency,
    std::size_t stride,
    bool force_fallback) const {
    if (limit == 0 || kmer_size_ == 0 || query.size() < kmer_size_) return {};
    constexpr int bin_width = 4;
    auto& vote_bins = vote_workspace();
    // Empirically, 8 bins per query base is already conservative for the full
    // human V collection. The 0.5 load factor leaves room for fallback seeds.
    const int minimum_bin = floor_div(
        -static_cast<int>(std::min<std::size_t>(
            maximum_gene_length_, static_cast<std::size_t>(std::numeric_limits<int>::max()))),
        bin_width) - 1;
    const int maximum_bin = static_cast<int>(std::min<std::size_t>(
        query.size() / bin_width + 1,
        static_cast<std::size_t>(std::numeric_limits<int>::max())));
    vote_bins.reset(
        std::max<std::size_t>(query.size() * 8, genes_.size() * 2),
        genes_.size(), minimum_bin, maximum_bin);
    const auto collect = [&vote_bins, &query, stride](
        const auto& seed_index,
        const std::vector<std::uint32_t>& offsets,
        const std::vector<SeedHit>& dense_hits,
        std::uint32_t k,
        std::uint32_t weight,
        std::size_t frequency_limit) {
        if (k == 0) return;
        for_each_kmer(query, k, stride, [&](std::uint64_t code, std::size_t query_position) {
            if (!offsets.empty() && code + 1 < offsets.size()) {
                const auto begin = offsets[static_cast<std::size_t>(code)];
                const auto end = offsets[static_cast<std::size_t>(code) + 1];
                if (static_cast<std::size_t>(end - begin) > frequency_limit) return;
                for (std::uint32_t i = begin; i < end; ++i) {
                    const auto& hit = dense_hits[i];
                    const int diagonal =
                        static_cast<int>(query_position) - static_cast<int>(hit.position);
                    const int bin = floor_div(diagonal, bin_width);
                    const std::uint64_t key =
                        (static_cast<std::uint64_t>(hit.gene) << 32U) |
                        static_cast<std::uint32_t>(bin);
                    vote_bins.add(key, weight);
                }
                return;
            }
            const auto found = seed_index.find(code);
            if (found == seed_index.end() || found->second.size() > frequency_limit) return;
            for (const auto& hit : found->second) {
                const int diagonal = static_cast<int>(query_position) - static_cast<int>(hit.position);
                const int bin = floor_div(diagonal, bin_width);
                const std::uint64_t key = (static_cast<std::uint64_t>(hit.gene) << 32U) |
                    static_cast<std::uint32_t>(bin);
                vote_bins.add(key, weight);
            }
        });
    };
    collect(seeds_, fast_seed_offsets_, fast_seed_hits_,
        kmer_size_, 4, max_seed_frequency);
    std::uint32_t strongest_primary_bin = 0;
    vote_bins.for_each([&](std::uint64_t, std::uint32_t votes) {
        strongest_primary_bin = std::max(strongest_primary_bin, votes);
    });
    const auto weak_primary_threshold = std::max<std::uint32_t>(
        32U, kmer_size_ * 16U);
    const bool weak_primary_signal = strongest_primary_bin < weak_primary_threshold;
    const bool use_fallback = force_fallback || weak_primary_signal;
    if (use_fallback) {
        collect(fallback_seeds_, fast_fallback_seed_offsets_, fast_fallback_seed_hits_,
            fallback_kmer_size_, 1, max_seed_frequency * 4);
    }

    // The reference implementation sorts every occupied diagonal bin, then
    // repeatedly scans all bins to collapse them by gene. Aggregate the same
    // sufficient statistics in two linear passes and sort only one entry per
    // gene. This is exactly the same ranking and diagonal-span calculation.
    auto& candidate_reduction = candidate_workspace();
    candidate_reduction.reset(genes_.size());
    vote_bins.for_each([&](std::uint64_t key, std::uint32_t votes) {
        const auto gene = static_cast<std::uint32_t>(key >> 32U);
        if (gene >= genes_.size()) return;
        const auto bin = static_cast<std::int32_t>(key & 0xffffffffU);
        const int diagonal = static_cast<int>(bin) * bin_width;
        candidate_reduction.consider(gene, votes, diagonal);
    });
    vote_bins.for_each([&](std::uint64_t key, std::uint32_t votes) {
        const auto gene = static_cast<std::uint32_t>(key >> 32U);
        if (!candidate_reduction.active(gene)) return;
        const auto threshold = std::max<std::uint32_t>(
            4, candidate_reduction.votes(gene) / 3);
        if (votes < threshold) return;
        const auto bin = static_cast<std::int32_t>(key & 0xffffffffU);
        candidate_reduction.include_diagonal(gene, static_cast<int>(bin) * bin_width);
    });

    std::vector<Candidate> ranked;
    ranked.reserve(candidate_reduction.touched().size());
    for (const auto gene : candidate_reduction.touched()) {
        ranked.push_back(Candidate{
            gene,
            candidate_reduction.diagonal(gene),
            static_cast<std::uint16_t>(std::min<std::uint32_t>(
                candidate_reduction.votes(gene), std::numeric_limits<std::uint16_t>::max())),
            candidate_reduction.span(gene),
            use_fallback,
        });
    }
    std::sort(ranked.begin(), ranked.end(), [](const Candidate& a, const Candidate& b) {
        if (a.votes != b.votes) return a.votes > b.votes;
        if (a.gene_index != b.gene_index) return a.gene_index < b.gene_index;
        return a.diagonal < b.diagonal;
    });
    if (ranked.size() > std::min(limit, genes_.size())) {
        ranked.resize(std::min(limit, genes_.size()));
    }
    return ranked;
}

#ifndef SWIG_WEB

GermlineDatabase GermlineDatabase::from_fastas(
    const std::string& v_path,
    const std::string& d_path,
    const std::string& j_path,
    const std::string& c_path) {
    GermlineDatabase database;
    database.v.reset(read_germline_fasta(v_path), 0);
    database.v_tree.reset(database.v.genes(), 12);
    database.d.reset(read_germline_fasta(d_path), 5);
    database.j.reset(read_germline_fasta(j_path), 7);
    database.c.reset(read_germline_fasta(c_path), 9);
    if (database.v.empty()) throw std::runtime_error("a V germline FASTA is required");
    if (database.j.empty()) throw std::runtime_error("a J germline FASTA is required");
    return database;
}

void GermlineDatabase::save(const std::string& path) const {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output) throw std::runtime_error("cannot create index: " + path);
    static constexpr std::array<char, 8> magic{'S', 'W', 'I', 'F', 'T', 'I', 'G', '\0'};
    output.write(magic.data(), static_cast<std::streamsize>(magic.size()));
    write_pod<std::uint32_t>(output, 3U);
    write_segment(output, v);
    write_segment(output, d);
    write_segment(output, j);
    write_segment(output, c);
    output.flush();
    if (!output) throw std::runtime_error("failed to finish index: " + path);
}

GermlineDatabase GermlineDatabase::load(const std::string& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("cannot open index: " + path);
    std::array<char, 8> magic{};
    input.read(magic.data(), static_cast<std::streamsize>(magic.size()));
    static constexpr std::array<char, 8> expected{'S', 'W', 'I', 'F', 'T', 'I', 'G', '\0'};
    if (!input || magic != expected) throw std::runtime_error("not a SwiftIG index: " + path);
    const auto version = read_pod<std::uint32_t>(input);
    if (version != 1U && version != 2U && version != 3U) {
        throw std::runtime_error("unsupported SwiftIG index version: " + std::to_string(version));
    }
    GermlineDatabase database;
    database.v = read_segment(input, version);
    database.d = read_segment(input, version);
    database.j = read_segment(input, version);
    database.c = read_segment(input, version);
    if (database.v.empty() || database.j.empty()) throw std::runtime_error("index has no V or J genes");
    return database;
}

std::size_t GermlineDatabase::apply_auxiliary(const std::string& path) {
    if (path.empty()) return 0;
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open IgBLAST auxiliary file: " + path);
    struct AuxiliaryEntry {
        int frame = -1;
        std::string chain;
        int cdr3_stop = -1;
    };
    std::unordered_map<std::string, AuxiliaryEntry> entries;
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const auto first = line.find_first_not_of(" \t");
        if (first == std::string::npos || line[first] == '#') continue;
        std::istringstream stream(line);
        std::string name;
        AuxiliaryEntry entry;
        if (!(stream >> name >> entry.frame >> entry.chain)) {
            throw std::runtime_error("invalid auxiliary record at " + path + ":" + std::to_string(line_number));
        }
        if (!(stream >> entry.cdr3_stop)) entry.cdr3_stop = -1;
        if (entry.frame < -1 || entry.frame > 2 || entry.cdr3_stop < -1) {
            throw std::runtime_error("invalid auxiliary coordinate at " + path + ":" + std::to_string(line_number));
        }
        entries[name] = std::move(entry);
    }
    std::size_t applied = 0;
    for (auto& gene : j.mutable_genes()) {
        const auto found = entries.find(gene.name);
        if (found == entries.end()) continue;
        gene.coding_frame_start = found->second.frame;
        gene.cdr3_stop = found->second.cdr3_stop;
        gene.chain_type = found->second.chain;
        if (gene.locus.empty()) {
            static const std::unordered_map<std::string, std::string> chain_loci{
                {"JH", "IGH"}, {"JK", "IGK"}, {"JL", "IGL"}, {"JA", "TRA"},
                {"JB", "TRB"}, {"JD", "TRD"}, {"JG", "TRG"}};
            if (const auto locus = chain_loci.find(gene.chain_type); locus != chain_loci.end()) {
                gene.locus = locus->second;
            }
        }
        ++applied;
    }
    if (applied == 0) {
        throw std::runtime_error("auxiliary file has no identifiers matching the J germline database: " + path);
    }
    return applied;
}

std::size_t GermlineDatabase::apply_internal_data(const std::string& path) {
    if (path.empty()) return 0;
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open IgBLAST internal-data file: " + path);
    struct InternalEntry {
        std::array<int, 10> bounds{{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};
        std::string chain;
        int frame = -1;
    };
    std::unordered_map<std::string, InternalEntry> entries;
    std::unordered_map<std::string, InternalEntry> genes;
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const auto first = line.find_first_not_of(" \t");
        if (first == std::string::npos || line[first] == '#') continue;
        std::istringstream stream(line);
        std::string name;
        InternalEntry entry;
        if (!(stream >> name)) {
            throw std::runtime_error("invalid internal-data record at " + path + ":" + std::to_string(line_number));
        }
        for (std::size_t i = 0; i < entry.bounds.size(); i += 2) {
            int start = 0;
            int stop = 0;
            if (!(stream >> start >> stop) || start < 1 || stop < start) {
                throw std::runtime_error("invalid FWR/CDR coordinates at " + path + ":" +
                    std::to_string(line_number));
            }
            entry.bounds[i] = start - 1;
            entry.bounds[i + 1] = stop;
        }
        if (!(stream >> entry.chain >> entry.frame) || entry.frame < -1 || entry.frame > 2) {
            throw std::runtime_error("invalid chain/frame at " + path + ":" + std::to_string(line_number));
        }
        entries[name] = entry;
        const auto star = name.find('*');
        genes.emplace(name.substr(0, star), entry);
    }

    std::size_t applied = 0;
    for (auto& gene : v.mutable_genes()) {
        const InternalEntry* entry = nullptr;
        if (const auto exact = entries.find(gene.name); exact != entries.end()) {
            entry = &exact->second;
        } else {
            const auto star = gene.name.find('*');
            if (const auto by_gene = genes.find(gene.name.substr(0, star)); by_gene != genes.end()) {
                entry = &by_gene->second;
            }
        }
        if (!entry) continue;
        gene.region_bounds = entry->bounds;
        gene.chain_type = entry->chain;
        gene.coding_frame_start = entry->frame;
        if (gene.locus.empty()) {
            static const std::unordered_map<std::string, std::string> chain_loci{
                {"VH", "IGH"}, {"VK", "IGK"}, {"VL", "IGL"}, {"VA", "TRA"},
                {"VB", "TRB"}, {"VD", "TRD"}, {"VG", "TRG"}};
            if (const auto locus = chain_loci.find(gene.chain_type); locus != chain_loci.end()) {
                gene.locus = locus->second;
            }
        }
        ++applied;
    }
    if (applied == 0) {
        throw std::runtime_error("internal-data file has no identifiers matching the V germline database: " + path);
    }
    return applied;
}

std::size_t GermlineDatabase::apply_d_frame_data(const std::string& path) {
    if (path.empty()) return 0;
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open IgBLAST D-frame file: " + path);
    std::unordered_map<std::string, int> frames;
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const auto first = line.find_first_not_of(" \t");
        if (first == std::string::npos || line[first] == '#') continue;
        std::istringstream stream(line);
        std::string name;
        int start = -1;
        if (!(stream >> name >> start) || start < 0) {
            throw std::runtime_error("invalid D-frame record at " + path + ":" + std::to_string(line_number));
        }
        frames[name] = start;
    }
    std::size_t applied = 0;
    for (auto& gene : d.mutable_genes()) {
        if (const auto found = frames.find(gene.name); found != frames.end()) {
            gene.coding_frame_start = found->second;
            ++applied;
        }
    }
    if (applied == 0) {
        throw std::runtime_error("D-frame file has no identifiers matching the D germline database: " + path);
    }
    return applied;
}

#endif

std::size_t GermlineDatabase::gene_count() const noexcept {
    return v.genes().size() + d.genes().size() + j.genes().size() + c.genes().size();
}

}  // namespace swiftig
