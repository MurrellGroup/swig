#pragma once

#include <cstddef>
#include <optional>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "swiftig/index.hpp"
#include "swiftig/types.hpp"

namespace swiftig {

enum class DoubleDMode : int { Off = 0, All = 1, LongSpan = 2 };

struct DoubleDOptions {
    DoubleDMode mode = DoubleDMode::Off;
    std::size_t minimum_vj_span = 40;
    std::size_t seed_length = 11;
    std::size_t pseudo_trim = 5;
    int maximum_pseudo_mismatches = 3;
    int minimum_score_gain = 8;
    // AER-R-only rescue: admit mutation-tolerant seed starts, while requiring
    // the resulting extended hit to cover at least the configured seed length.
    bool robust_seed_rescue = false;
};

struct DoubleDAlternative {
    std::string d_call;
    std::string d2_call;
    int pair_score = 0;
    int score_gain = 0;
    int pseudo_distance = -1;
};

struct DoubleDCall {
    SegmentHit d;
    SegmentHit d2;
    std::string baseline_d_call;
    std::size_t vj_span = 0;
    std::size_t inter_d_length = 0;
    int pair_score = 0;
    int best_single_score = 0;
    int score_gain = 0;
    int pseudo_distance = -1;
    std::vector<DoubleDAlternative> alternatives;
};

class DoubleDScreener {
public:
    DoubleDScreener(const GermlineDatabase& database, DoubleDOptions options);
    [[nodiscard]] std::optional<DoubleDCall> screen(const Annotation& annotation) const;

private:
    struct SeedHit {
        std::size_t gene = 0;
        std::size_t position = 0;
    };

    const GermlineDatabase& database_;
    DoubleDOptions options_;
    std::unordered_map<std::string, std::vector<SeedHit>> seeds_;
    std::unordered_map<std::string, std::vector<SeedHit>> fallback_seeds_;
    std::size_t fallback_seed_length_ = 0;
};

void write_double_d_header(std::ostream& output);
void write_double_d_record(
    std::ostream& output,
    const Annotation& annotation,
    const DoubleDCall& call,
    const DoubleDOptions& options,
    std::size_t batch_record_index);

}  // namespace swiftig
