// Deterministic null simulation used to calibrate SwiftIG's AIRR support values.
//
// Build and run from the repository root:
//   c++ -O3 -std=c++20 scripts/calibrate-airr-support.cpp -o /tmp/calibrate-airr-support
//   /tmp/calibrate-airr-support > /tmp/swiftig-null-scores.tsv
//
// The score-only recurrence below intentionally mirrors
// wasm/src/alignment.cpp::local_align_affine without traceback allocation.

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <string>
#include <vector>

namespace {

struct Scoring {
    const char* name;
    int match;
    int mismatch;
    int gap_open;
    int gap_extend;
};

struct Experiment {
    std::size_t query_length;
    std::size_t reference_length;
    std::size_t replicates;
};

int local_score(
    const std::string& query,
    const std::string& reference,
    const Scoring& scoring) {
    constexpr int negative_infinity = std::numeric_limits<int>::min() / 4;
    std::vector<int> previous(reference.size() + 1, 0);
    std::vector<int> current(reference.size() + 1, 0);
    std::vector<int> insertion_previous(reference.size() + 1, negative_infinity);
    std::vector<int> insertion_current(reference.size() + 1, negative_infinity);
    int best = 0;
    for (std::size_t i = 1; i <= query.size(); ++i) {
        int deletion = negative_infinity;
        current[0] = 0;
        insertion_current[0] = negative_infinity;
        for (std::size_t j = 1; j <= reference.size(); ++j) {
            const int insertion = std::max(
                previous[j] + scoring.gap_open,
                insertion_previous[j] + scoring.gap_extend);
            insertion_current[j] = insertion;
            deletion = std::max(
                current[j - 1] + scoring.gap_open,
                deletion + scoring.gap_extend);
            const int substitution = previous[j - 1] +
                (query[i - 1] == reference[j - 1] ? scoring.match : scoring.mismatch);
            current[j] = std::max({0, substitution, insertion, deletion});
            best = std::max(best, current[j]);
        }
        previous.swap(current);
        insertion_previous.swap(insertion_current);
    }
    return best;
}

std::string random_dna(std::mt19937_64& generator, std::size_t length) {
    static constexpr std::array<char, 4> bases{'A', 'C', 'G', 'T'};
    std::string sequence(length, 'A');
    for (char& base : sequence) base = bases[generator() & 3U];
    return sequence;
}

}  // namespace

int main() {
    // Fixed seed and an explicit experiment table make the embedded constants
    // independently reproducible. The lengths span short D/J searches through
    // ordinary V/C refinement without making calibration part of end-user runs.
    std::mt19937_64 generator(0x535749474556414cULL);
    constexpr std::array<Scoring, 5> profiles{{
        {"v_c_default", 2, -3, -5, -1},
        {"d_truth", 2, -3, -13, -1},
        {"j_truth", 2, -3, -17, -2},
        {"d_igblast_agreement", 2, -4, -11, -1},
        {"j_igblast_agreement", 2, -4, -13, -1},
    }};
    constexpr std::array<Experiment, 15> experiments{{
        {20, 20, 30000},
        {24, 40, 25000},
        {40, 24, 25000},
        {40, 64, 12000},
        {64, 40, 12000},
        {64, 96, 5000},
        {96, 64, 5000},
        {96, 160, 1800},
        {160, 96, 1800},
        {160, 256, 1200},
        {256, 160, 1200},
        {256, 320, 500},
        {320, 256, 500},
        {320, 512, 240},
        {512, 320, 240},
    }};

    std::cout << "profile\tmatch\tmismatch\tgap_open\tgap_extend\tquery_length"
                 "\treference_length\treplicates\tscore\tcount\n";
    for (const auto& experiment : experiments) {
        std::array<std::map<int, std::size_t>, profiles.size()> histograms;
        for (std::size_t replicate = 0; replicate < experiment.replicates; ++replicate) {
            const auto query = random_dna(generator, experiment.query_length);
            const auto reference = random_dna(generator, experiment.reference_length);
            for (std::size_t profile = 0; profile < profiles.size(); ++profile) {
                ++histograms[profile][local_score(query, reference, profiles[profile])];
            }
        }
        for (std::size_t profile = 0; profile < profiles.size(); ++profile) {
            const auto& scoring = profiles[profile];
            for (const auto& [score, count] : histograms[profile]) {
                std::cout << scoring.name << '\t' << scoring.match << '\t'
                          << scoring.mismatch << '\t' << scoring.gap_open << '\t'
                          << scoring.gap_extend << '\t' << experiment.query_length << '\t'
                          << experiment.reference_length << '\t' << experiment.replicates << '\t'
                          << score << '\t' << count << '\n';
            }
        }
    }
}
