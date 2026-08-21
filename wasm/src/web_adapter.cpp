#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

#include "swiftig/airr.hpp"
#include "swiftig/double_d.hpp"
#include "swiftig/engine.hpp"
#include "swiftig/index.hpp"
#include "swiftig/types.hpp"

namespace {

using swiftig::AnnotationEngine;
using swiftig::AssignerStrategy;
using swiftig::EngineOptions;
using swiftig::Gene;
using swiftig::GermlineDatabase;
using swiftig::SequenceRecord;

std::unique_ptr<GermlineDatabase> g_database;
std::optional<EngineOptions> g_engine_options_override;
int g_calling_profile = 0;
AssignerStrategy g_assigner_strategy = AssignerStrategy::Standard;
bool g_optimized_kernels = true;
bool g_optimized_output = true;
std::string g_result;
std::string g_double_d_result;
int g_double_d_count = 0;
std::string g_error;

bool fail(std::string message) {
    g_error = std::move(message);
    return false;
}

bool memory_view(const char* data, std::size_t size, std::string_view& output) {
    if (data == nullptr && size != 0) return fail("invalid browser memory range");
    output = {data, size};
    return true;
}

void trim_cr(std::string& value) {
    if (!value.empty() && value.back() == '\r') value.pop_back();
}

void parse_header(std::string_view header, SequenceRecord& record) {
    const auto first = header.find_first_not_of(" \t");
    if (first == std::string_view::npos) return;
    const auto split = header.find_first_of(" \t", first);
    record.id = std::string(header.substr(first, split == std::string_view::npos
        ? std::string_view::npos : split - first));
    if (split != std::string_view::npos) {
        const auto description = header.find_first_not_of(" \t", split);
        if (description != std::string_view::npos) {
            record.description = std::string(header.substr(description));
        }
    }
}

std::string germline_name(const SequenceRecord& record) {
    std::vector<std::string> fields;
    std::stringstream stream(record.id);
    std::string field;
    while (std::getline(stream, field, '|')) fields.push_back(field);
    for (const auto& candidate : fields) {
        if (!swiftig::infer_locus(candidate).empty() && candidate.find('*') != std::string::npos) {
            return candidate;
        }
    }
    for (const auto& candidate : fields) {
        if (!swiftig::infer_locus(candidate).empty()) return candidate;
    }
    return record.id;
}

bool parse_fasta(
    std::string_view input,
    bool references,
    std::vector<SequenceRecord>& records) {
    std::istringstream stream{std::string(input)};
    SequenceRecord current;
    std::string line;
    bool seen_header = false;
    while (std::getline(stream, line)) {
        trim_cr(line);
        if (line.empty()) continue;
        if (line.front() == '>') {
            if (seen_header) {
                current.sequence = swiftig::normalize_dna(current.sequence, references);
                if (current.sequence.empty()) return fail("empty FASTA record: " + current.id);
                records.push_back(std::move(current));
                current = {};
            }
            seen_header = true;
            parse_header(std::string_view(line).substr(1), current);
            if (current.id.empty()) current.id = "sequence_" + std::to_string(records.size() + 1);
        } else {
            if (!seen_header) return fail("expected FASTA header beginning with '>'");
            current.sequence += line;
        }
    }
    if (seen_header) {
        current.sequence = swiftig::normalize_dna(current.sequence, references);
        if (current.sequence.empty()) return fail("empty FASTA record: " + current.id);
        records.push_back(std::move(current));
    }
    return !records.empty() || fail("no FASTA records found");
}

bool parse_fastq(std::string_view input, std::vector<SequenceRecord>& records) {
    std::istringstream stream{std::string(input)};
    std::string line;
    while (std::getline(stream, line)) {
        trim_cr(line);
        if (line.empty()) continue;
        if (line.front() != '@') return fail("expected FASTQ header beginning with '@'");
        SequenceRecord record;
        parse_header(std::string_view(line).substr(1), record);
        if (record.id.empty()) record.id = "sequence_" + std::to_string(records.size() + 1);
        bool found_plus = false;
        while (std::getline(stream, line)) {
            trim_cr(line);
            if (!line.empty() && line.front() == '+') {
                found_plus = true;
                break;
            }
            record.sequence += line;
        }
        if (!found_plus) return fail("truncated FASTQ sequence: " + record.id);
        record.sequence = swiftig::normalize_dna(record.sequence);
        while (record.quality.size() < record.sequence.size() && std::getline(stream, line)) {
            trim_cr(line);
            record.quality += line;
        }
        if (record.quality.size() != record.sequence.size()) {
            return fail("FASTQ sequence/quality length mismatch: " + record.id);
        }
        records.push_back(std::move(record));
    }
    return !records.empty() || fail("no FASTQ records found");
}

std::vector<std::string> split_delimited(const std::string& line, char delimiter) {
    std::vector<std::string> values;
    std::size_t start = 0;
    while (start <= line.size()) {
        const auto end = line.find(delimiter, start);
        values.push_back(line.substr(start, end == std::string::npos ? std::string::npos : end - start));
        if (end == std::string::npos) break;
        start = end + 1;
    }
    return values;
}

bool parse_airr(std::string_view input, std::vector<SequenceRecord>& records) {
    std::istringstream stream{std::string(input)};
    std::string line;
    while (std::getline(stream, line)) {
        trim_cr(line);
        if (!line.empty()) break;
    }
    if (line.empty()) return fail("AIRR table is empty");
    const char delimiter = line.find('\t') != std::string::npos ? '\t' : ',';
    const auto header = split_delimited(line, delimiter);
    const auto find_column = [&](const std::string& name) {
        const auto found = std::find(header.begin(), header.end(), name);
        return found == header.end() ? header.size() : static_cast<std::size_t>(found - header.begin());
    };
    const auto sequence_column = find_column("sequence");
    const auto id_column = find_column("sequence_id");
    const auto quality_column = find_column("quality");
    if (sequence_column == header.size()) return fail("AIRR input requires a 'sequence' column");
    while (std::getline(stream, line)) {
        trim_cr(line);
        if (line.empty()) continue;
        const auto values = split_delimited(line, delimiter);
        if (sequence_column >= values.size() || values[sequence_column].empty()) continue;
        SequenceRecord record;
        record.id = id_column < values.size() && !values[id_column].empty()
            ? values[id_column] : "sequence_" + std::to_string(records.size() + 1);
        record.sequence = swiftig::normalize_dna(values[sequence_column]);
        if (quality_column < values.size()) record.quality = values[quality_column];
        if (!record.sequence.empty()) records.push_back(std::move(record));
    }
    return !records.empty() || fail("AIRR table has no non-empty sequence rows");
}

bool parse_queries(std::string_view input, int format, std::vector<SequenceRecord>& records) {
    const auto first = input.find_first_not_of(" \t\r\n");
    if (first == std::string_view::npos) return fail("query input is empty");
    if (format == 1 || (format == 0 && input[first] == '>')) {
        return parse_fasta(input, false, records);
    }
    if (format == 2 || (format == 0 && input[first] == '@')) return parse_fastq(input, records);
    if (format == 3 || format == 0) return parse_airr(input, records);
    return fail("unsupported query format");
}

bool parse_germline(
    std::string_view input,
    const char* segment,
    bool required,
    std::vector<Gene>& genes) {
    if (input.empty()) {
        if (required) return fail(std::string("a ") + segment + " germline FASTA is required");
        return true;
    }
    std::vector<SequenceRecord> records;
    if (!parse_fasta(input, true, records)) return false;
    genes.reserve(records.size());
    std::unordered_set<std::string> seen;
    for (const auto& record : records) {
        auto name = germline_name(record);
        if (!seen.insert(name).second) {
            return fail(std::string("duplicate ") + segment + " germline identifier: " + name);
        }
        auto locus = swiftig::infer_locus(record.id + " " + record.description);
        if (locus.empty()) locus = swiftig::infer_locus(name);
        Gene gene;
        gene.name = std::move(name);
        gene.sequence = record.sequence;
        gene.locus = std::move(locus);
        const auto marker = record.description.find("SWIGMETA=");
        if (marker != std::string::npos) {
            std::array<int, 13> values{};
            std::size_t cursor = marker + 9;
            bool valid = true;
            for (std::size_t index = 0; index < values.size(); ++index) {
                bool negative = cursor < record.description.size() && record.description[cursor] == '-';
                if (negative) ++cursor;
                if (cursor >= record.description.size() ||
                    record.description[cursor] < '0' || record.description[cursor] > '9') {
                    valid = false;
                    break;
                }
                int value = 0;
                while (cursor < record.description.size() &&
                    record.description[cursor] >= '0' && record.description[cursor] <= '9') {
                    value = value * 10 + (record.description[cursor] - '0');
                    ++cursor;
                }
                values[index] = negative ? -value : value;
                if (index + 1 < values.size()) {
                    if (cursor >= record.description.size() || record.description[cursor] != ',') {
                        valid = false;
                        break;
                    }
                    ++cursor;
                }
            }
            if (!valid || values[0] < -1 || values[0] > 2 || values[1] < -1) {
                return fail(std::string("invalid SWIGMETA annotation for ") + gene.name);
            }
            gene.coding_frame_start = values[0];
            gene.cdr3_stop = values[1];
            for (std::size_t index = 0; index < gene.region_bounds.size(); ++index) {
                gene.region_bounds[index] = values[index + 2];
            }
            switch (values[12]) {
                case 1: gene.annotation_source = "IMGT-gapped"; break;
                case 2: gene.annotation_source = "AIRR-C"; break;
                case 3: gene.annotation_source = "IMGT-boundary-transfer"; break;
                case 4: gene.annotation_source = "validated-J-motif"; break;
                case 5: gene.annotation_source = "provided"; break;
                case 6: gene.annotation_source = "J-anchor-transfer"; break;
                default: gene.annotation_source = "unclassified"; break;
            }
        }
        genes.push_back(std::move(gene));
    }
    return true;
}

std::string error_text() {
    return g_error.empty() ? "SwiftIG could not complete the request" : g_error;
}

EngineOptions configured_options(int minimum_identity_per_mille, int strand) {
    EngineOptions options;
    if (g_engine_options_override) {
        options = *g_engine_options_override;
    } else if (g_calling_profile == 1) {
        // Calibrated solely against the supplied simulated human-IGH IgBLAST
        // calls. This is an explicit compatibility profile, not the default.
        options.d_scoring = {2, -4, -11, -1};
        options.top_d = 3;
        options.min_d_match = 5;
        options.j_scoring = {2, -4, -13, -1};
        options.top_j = 2;
    } else if (g_calling_profile == 2) {
        // Experimental AER-R-only profile. These values are calibrated with
        // ambiguity-aware scores and boundary losses; the production caller
        // rejects this profile for other assignment strategies.
        options.v_scoring = {2, -4, -13, -1};
        options.aer_r_optimized = true;
        options.aer_r_d_presence_penalty = 10;
        options.aer_r_evidence_conditioned_d_penalty = true;
        options.d_scoring = {2, -4, -13, -1};
        options.top_d = 2;
        options.min_d_match = 5;
        options.j_scoring = {2, -4, -17, -2};
        options.top_j = 2;
    }
    options.assigner_strategy = g_assigner_strategy;
    options.optimized_kernels = g_optimized_kernels;
    options.min_identity = std::clamp(minimum_identity_per_mille, 0, 1000) / 1000.0;
    options.search_forward = strand != 2;
    options.search_reverse = strand != 1;
    return options;
}

}  // namespace

extern "C" {

__attribute__((export_name("swig_alloc")))
void* swig_alloc(std::size_t size) { return std::malloc(size); }

__attribute__((export_name("swig_free")))
void swig_free(void* pointer) { std::free(pointer); }

__attribute__((export_name("swig_set_calling_profile")))
int swig_set_calling_profile(int profile) noexcept {
    if (profile < 0 || profile > 2) return -1;
    g_calling_profile = profile;
    return 0;
}

__attribute__((export_name("swig_set_assigner_strategy")))
int swig_set_assigner_strategy(int strategy) noexcept {
    if (strategy < static_cast<int>(AssignerStrategy::Standard) ||
        strategy > static_cast<int>(AssignerStrategy::AerRobust)) return -1;
    g_assigner_strategy = static_cast<AssignerStrategy>(strategy);
    return 0;
}

// Diagnostic switch used by equivalence tests and native profiling. Normal
// web/CLI callers never need to set it; optimized AER/RIAT kernels are default.
__attribute__((export_name("swig_set_optimized_kernels")))
int swig_set_optimized_kernels(int enabled) noexcept {
    if (enabled != 0 && enabled != 1) return -1;
    g_optimized_kernels = enabled != 0;
    return 0;
}

__attribute__((export_name("swig_set_optimized_output")))
int swig_set_optimized_output(int enabled) noexcept {
    if (enabled != 0 && enabled != 1) return -1;
    g_optimized_output = enabled != 0;
    return 0;
}

__attribute__((export_name("swig_init_database")))
int swig_init_database(
    const char* v_data, std::size_t v_size,
    const char* d_data, std::size_t d_size,
    const char* j_data, std::size_t j_size,
    const char* c_data, std::size_t c_size) noexcept {
    g_error.clear();
    std::string_view v_input, d_input, j_input, c_input;
    if (!memory_view(v_data, v_size, v_input) || !memory_view(d_data, d_size, d_input) ||
        !memory_view(j_data, j_size, j_input) || !memory_view(c_data, c_size, c_input)) return -1;
    std::vector<Gene> v_genes, d_genes, j_genes, c_genes;
    if (!parse_germline(v_input, "V", true, v_genes) ||
        !parse_germline(d_input, "D", false, d_genes) ||
        !parse_germline(j_input, "J", true, j_genes) ||
        !parse_germline(c_input, "C", false, c_genes)) return -1;
    auto database = std::make_unique<GermlineDatabase>();
    if (g_assigner_strategy == AssignerStrategy::RiatMp) {
        // RIAT-MP indexes only one representative root per close-allele tree.
        // Descendant alleles are retained as sparse tree leaves and are never
        // independently aligned by the production RIAT-MP path.
        database->v.reset(std::move(v_genes), 0);
        database->v_tree.reset(database->v.genes(), 12);
    } else {
        database->v.reset(std::move(v_genes), 9);
    }
    database->d.reset(std::move(d_genes), 5);
    database->j.reset(std::move(j_genes), 7);
    database->c.reset(std::move(c_genes), 9);
    const auto count = static_cast<int>(database->gene_count());
    g_database = std::move(database);
    g_double_d_result.clear();
    g_double_d_count = 0;
    return count;
}

__attribute__((export_name("swig_annotate")))
int swig_annotate(
    const char* query_data,
    std::size_t query_size,
    int format,
    int minimum_identity_per_mille,
    int strand) noexcept {
    g_error.clear();
    if (!g_database) {
        fail("load a reference database before annotation");
        return -1;
    }
    std::string_view query_input;
    if (!memory_view(query_data, query_size, query_input)) return -1;
    std::vector<SequenceRecord> records;
    if (!parse_queries(query_input, format, records)) return -1;
    const auto options = configured_options(minimum_identity_per_mille, strand);
    AnnotationEngine engine(*g_database, options);
    if (g_optimized_output) {
        g_result.clear();
        if (query_size <=
            (std::numeric_limits<std::size_t>::max() - 4096) / 6) {
            g_result.reserve(query_size * 6 + 4096);
        }
        swiftig::append_airr_header(g_result);
        for (const auto& record : records) {
            swiftig::append_airr_record(g_result, engine.annotate(record));
        }
    } else {
        std::ostringstream output;
        swiftig::write_airr_header(output);
        for (const auto& record : records) {
            swiftig::write_airr_record(output, engine.annotate(record));
        }
        g_result = std::move(output).str();
    }
    g_double_d_result.clear();
    g_double_d_count = 0;
    return static_cast<int>(records.size());
}

__attribute__((export_name("swig_annotate_double_d")))
int swig_annotate_double_d(
    const char* query_data,
    std::size_t query_size,
    int format,
    int minimum_identity_per_mille,
    int strand,
    int mode,
    int minimum_vj_span,
    int seed_length,
    int pseudo_trim,
    int maximum_pseudo_mismatches,
    int minimum_score_gain) noexcept {
    g_error.clear();
    g_double_d_count = 0;
    if (!g_database) {
        fail("load a reference database before annotation");
        return -1;
    }
    std::string_view query_input;
    if (!memory_view(query_data, query_size, query_input)) return -1;
    std::vector<SequenceRecord> records;
    if (!parse_queries(query_input, format, records)) return -1;
    const auto engine_options = configured_options(minimum_identity_per_mille, strand);
    AnnotationEngine engine(*g_database, engine_options);
    swiftig::DoubleDOptions double_d_options;
    double_d_options.mode = mode == 1
        ? swiftig::DoubleDMode::All : mode == 2
            ? swiftig::DoubleDMode::LongSpan : swiftig::DoubleDMode::Off;
    double_d_options.minimum_vj_span = static_cast<std::size_t>(
        std::clamp(minimum_vj_span, 0, 10000));
    double_d_options.seed_length = static_cast<std::size_t>(
        std::clamp(seed_length, 6, 24));
    double_d_options.pseudo_trim = static_cast<std::size_t>(
        std::clamp(pseudo_trim, 0, 24));
    double_d_options.maximum_pseudo_mismatches = std::clamp(
        maximum_pseudo_mismatches, 0, 24);
    double_d_options.minimum_score_gain = std::clamp(minimum_score_gain, 0, 1000);
    double_d_options.robust_seed_rescue =
        engine_options.assigner_strategy == AssignerStrategy::AerRobust;
    swiftig::DoubleDScreener screener(*g_database, double_d_options);
    std::ostringstream reference_airr_output;
    std::ostringstream double_d_output;
    g_result.clear();
    if (g_optimized_output) {
        if (query_size <=
            (std::numeric_limits<std::size_t>::max() - 4096) / 6) {
            g_result.reserve(query_size * 6 + 4096);
        }
        swiftig::append_airr_header(g_result);
    } else {
        swiftig::write_airr_header(reference_airr_output);
    }
    swiftig::write_double_d_header(double_d_output);
    for (std::size_t record_index = 0; record_index < records.size(); ++record_index) {
        const auto annotation = engine.annotate(records[record_index]);
        if (g_optimized_output) swiftig::append_airr_record(g_result, annotation);
        else swiftig::write_airr_record(reference_airr_output, annotation);
        if (const auto call = screener.screen(annotation)) {
            swiftig::write_double_d_record(
                double_d_output, annotation, *call, double_d_options, record_index);
            ++g_double_d_count;
        }
    }
    if (!g_optimized_output) g_result = std::move(reference_airr_output).str();
    g_double_d_result = std::move(double_d_output).str();
    return static_cast<int>(records.size());
}

__attribute__((export_name("swig_result_ptr")))
const char* swig_result_ptr() noexcept { return g_result.data(); }

__attribute__((export_name("swig_result_len")))
std::size_t swig_result_len() noexcept { return g_result.size(); }

__attribute__((export_name("swig_double_d_result_ptr")))
const char* swig_double_d_result_ptr() noexcept { return g_double_d_result.data(); }

__attribute__((export_name("swig_double_d_result_len")))
std::size_t swig_double_d_result_len() noexcept { return g_double_d_result.size(); }

__attribute__((export_name("swig_double_d_count")))
int swig_double_d_count() noexcept { return g_double_d_count; }

__attribute__((export_name("swig_error_ptr")))
const char* swig_error_ptr() noexcept { return g_error.data(); }

__attribute__((export_name("swig_error_len")))
std::size_t swig_error_len() noexcept { return g_error.size(); }

__attribute__((export_name("swig_gene_count")))
std::size_t swig_gene_count(int segment) noexcept {
    if (!g_database) return 0;
    switch (segment) {
        case 0: return g_database->v.genes().size();
        case 1: return g_database->d.genes().size();
        case 2: return g_database->j.genes().size();
        case 3: return g_database->c.genes().size();
        default: return 0;
    }
}

// Experimental root-index diagnostics. These exports let the benchmark prove
// that the complete V collection is not k-mer indexed or independently
// aligned, and account for the sparse work done below the selected roots.
__attribute__((export_name("swig_v_tree_cluster_count")))
std::uint32_t swig_v_tree_cluster_count() noexcept {
    return g_database
        ? static_cast<std::uint32_t>(g_database->v_tree.clusters().size()) : 0;
}

__attribute__((export_name("swig_v_tree_edge_count")))
std::uint32_t swig_v_tree_edge_count() noexcept {
    return g_database ? static_cast<std::uint32_t>(g_database->v_tree.edge_count()) : 0;
}

__attribute__((export_name("swig_v_tree_mutation_count")))
std::uint32_t swig_v_tree_mutation_count() noexcept {
    return g_database ? static_cast<std::uint32_t>(g_database->v_tree.mutation_count()) : 0;
}

__attribute__((export_name("swig_v_tree_max_edge_mutations")))
std::uint32_t swig_v_tree_max_edge_mutations() noexcept {
    return g_database
        ? static_cast<std::uint32_t>(g_database->v_tree.maximum_edge_mutations()) : 0;
}

__attribute__((export_name("swig_v_tree_full_v_kmer_size")))
std::uint32_t swig_v_tree_full_v_kmer_size() noexcept {
    return g_database ? g_database->v.kmer_size() : 0;
}

__attribute__((export_name("swig_v_tree_root_kmer_size")))
std::uint32_t swig_v_tree_root_kmer_size() noexcept {
    return g_database ? g_database->v_tree.roots().kmer_size() : 0;
}

__attribute__((export_name("swig_v_tree_reset_stats")))
void swig_v_tree_reset_stats() noexcept {
    if (g_database) g_database->v_tree.reset_search_stats();
}

#define SWIG_V_TREE_STAT_EXPORT(function_name, field) \
    __attribute__((export_name(#function_name))) \
    std::uint32_t function_name() noexcept { \
        return g_database ? g_database->v_tree.search_stats().field : 0; \
    }

SWIG_V_TREE_STAT_EXPORT(swig_v_tree_stats_queries, queries)
SWIG_V_TREE_STAT_EXPORT(swig_v_tree_stats_root_candidates, root_candidates)
SWIG_V_TREE_STAT_EXPORT(swig_v_tree_stats_root_alignments, root_alignments)
SWIG_V_TREE_STAT_EXPORT(swig_v_tree_stats_root_tracebacks, root_tracebacks)
SWIG_V_TREE_STAT_EXPORT(swig_v_tree_stats_multipath_searches, multipath_searches)
SWIG_V_TREE_STAT_EXPORT(swig_v_tree_stats_trace_states, trace_states)
SWIG_V_TREE_STAT_EXPORT(swig_v_tree_stats_trace_limit_hits, trace_limit_hits)
SWIG_V_TREE_STAT_EXPORT(swig_v_tree_stats_clusters_scored, clusters_scored)
SWIG_V_TREE_STAT_EXPORT(swig_v_tree_stats_nodes_scored, nodes_scored)
SWIG_V_TREE_STAT_EXPORT(swig_v_tree_stats_mutation_updates, mutation_updates)
SWIG_V_TREE_STAT_EXPORT(swig_v_tree_stats_final_realignments, final_realignments)

#undef SWIG_V_TREE_STAT_EXPORT

// Benchmark-only configuration hook. The browser UI never invokes this export;
// normal annotation therefore continues to use EngineOptions defaults. Keeping
// it in the WASM ABI makes simulated-reference tuning exactly reproducible
// without rebuilding one binary for every scoring combination.
__attribute__((export_name("swig_set_tuning_options")))
int swig_set_tuning_options(
    int d_match,
    int d_mismatch,
    int d_gap_open,
    int d_gap_extend,
    int top_d,
    int min_d_match,
    int j_match,
    int j_mismatch,
    int j_gap_open,
    int j_gap_extend,
    int top_j,
    int min_j_length) noexcept {
    EngineOptions options;
    options.d_scoring = {
        std::clamp(d_match, 1, 20),
        std::clamp(d_mismatch, -50, 0),
        std::clamp(d_gap_open, -100, 0),
        std::clamp(d_gap_extend, -50, 0)};
    options.top_d = static_cast<std::size_t>(std::clamp(top_d, 1, 1000));
    options.min_d_match = static_cast<std::size_t>(std::clamp(min_d_match, 1, 100));
    options.j_scoring = {
        std::clamp(j_match, 1, 20),
        std::clamp(j_mismatch, -50, 0),
        std::clamp(j_gap_open, -100, 0),
        std::clamp(j_gap_extend, -50, 0)};
    options.top_j = static_cast<std::size_t>(std::clamp(top_j, 1, 1000));
    options.min_j_length = static_cast<std::size_t>(std::clamp(min_j_length, 1, 500));
    g_engine_options_override = options;
    return 0;
}

__attribute__((export_name("swig_clear_tuning_options")))
void swig_clear_tuning_options() noexcept { g_engine_options_override.reset(); }

// Companion to swig_set_tuning_options. Benchmarks use it for calibration;
// the direct CLI also preserves R-optimized V decisions when a user overrides
// one of the exposed D/J compatibility controls.
__attribute__((export_name("swig_set_v_tuning_options")))
int swig_set_v_tuning_options(
    int v_match,
    int v_mismatch,
    int v_gap_open,
    int v_gap_extend,
    int aer_r_optimized) noexcept {
    if (!g_engine_options_override) g_engine_options_override = EngineOptions{};
    g_engine_options_override->v_scoring = {
        std::clamp(v_match, 1, 20),
        std::clamp(v_mismatch, -50, 0),
        std::clamp(v_gap_open, -100, 0),
        std::clamp(v_gap_extend, -50, 0)};
    g_engine_options_override->aer_r_optimized = aer_r_optimized != 0;
    g_engine_options_override->aer_r_evidence_conditioned_d_penalty =
        aer_r_optimized != 0;
    return 0;
}

__attribute__((export_name("swig_set_aer_r_decision_tuning")))
int swig_set_aer_r_decision_tuning(int d_presence_penalty) noexcept {
    if (!g_engine_options_override) g_engine_options_override = EngineOptions{};
    g_engine_options_override->aer_r_optimized = true;
    // Calibration sweeps request an exact fixed cost. Production profile 2
    // enables the evidence-conditioned two-point relaxation separately.
    g_engine_options_override->aer_r_evidence_conditioned_d_penalty = false;
    g_engine_options_override->aer_r_d_presence_penalty =
        std::clamp(d_presence_penalty, 0, 100);
    return 0;
}

__attribute__((export_name("swig_version")))
const char* swig_version() noexcept { return swiftig::kVersion; }

}  // extern "C"
