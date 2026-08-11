#include "swiftig/fasta.hpp"

#include <algorithm>
#include <cctype>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

namespace swiftig {
namespace {

void trim_cr(std::string& line) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
}

void parse_header(const std::string& header, SequenceRecord& record) {
    const auto first = header.find_first_not_of(" \t");
    if (first == std::string::npos) return;
    const auto split = header.find_first_of(" \t", first);
    record.id = header.substr(first, split == std::string::npos ? std::string::npos : split - first);
    if (split != std::string::npos) {
        const auto desc = header.find_first_not_of(" \t", split);
        if (desc != std::string::npos) record.description = header.substr(desc);
    }
}

std::string germline_name(const SequenceRecord& record) {
    std::vector<std::string> fields;
    std::stringstream stream(record.id);
    std::string field;
    while (std::getline(stream, field, '|')) fields.push_back(field);
    for (const auto& candidate : fields) {
        if (!infer_locus(candidate).empty() && candidate.find('*') != std::string::npos) return candidate;
    }
    for (const auto& candidate : fields) {
        if (!infer_locus(candidate).empty()) return candidate;
    }
    return record.id;
}

}  // namespace

SequenceReader::SequenceReader(const std::string& path, bool references) : references_(references) {
    if (path == "-") {
        input_ = &std::cin;
    } else {
        owned_ = std::make_unique<std::ifstream>(path);
        if (!*owned_) throw std::runtime_error("cannot open sequence file: " + path);
        input_ = owned_.get();
    }
    while (std::getline(*input_, pending_)) {
        trim_cr(pending_);
        if (!pending_.empty()) break;
    }
    if (pending_.empty()) return;
    if (pending_.front() == '>') format_ = Format::Fasta;
    else if (pending_.front() == '@') format_ = Format::Fastq;
    else throw std::runtime_error("expected FASTA ('>') or FASTQ ('@') input");
}

bool SequenceReader::next(SequenceRecord& record) {
    if (format_ == Format::Fasta) return next_fasta(record);
    if (format_ == Format::Fastq) return next_fastq(record);
    return false;
}

bool SequenceReader::next_fasta(SequenceRecord& record) {
    if (pending_.empty()) return false;
    record = {};
    parse_header(pending_.substr(1), record);
    std::string line;
    pending_.clear();
    while (std::getline(*input_, line)) {
        trim_cr(line);
        if (!line.empty() && line.front() == '>') {
            pending_ = line;
            break;
        }
        record.sequence += line;
    }
    ++record_number_;
    if (record.id.empty()) record.id = "sequence_" + std::to_string(record_number_);
    record.sequence = normalize_dna(record.sequence, references_);
    if (record.sequence.empty()) throw std::runtime_error("empty sequence for record: " + record.id);
    return true;
}

bool SequenceReader::next_fastq(SequenceRecord& record) {
    if (pending_.empty()) return false;
    record = {};
    parse_header(pending_.substr(1), record);
    pending_.clear();
    std::string line;
    bool found_plus = false;
    while (std::getline(*input_, line)) {
        trim_cr(line);
        if (!line.empty() && line.front() == '+') {
            found_plus = true;
            break;
        }
        record.sequence += line;
    }
    if (!found_plus) throw std::runtime_error("truncated FASTQ sequence for: " + record.id);
    while (record.quality.size() < record.sequence.size() && std::getline(*input_, line)) {
        trim_cr(line);
        record.quality += line;
    }
    if (record.quality.size() != record.sequence.size()) {
        throw std::runtime_error("FASTQ sequence/quality length mismatch for: " + record.id);
    }
    while (std::getline(*input_, pending_)) {
        trim_cr(pending_);
        if (!pending_.empty()) break;
    }
    if (!pending_.empty() && pending_.front() != '@') {
        throw std::runtime_error("expected '@' at next FASTQ record");
    }
    ++record_number_;
    if (record.id.empty()) record.id = "sequence_" + std::to_string(record_number_);
    record.sequence = normalize_dna(record.sequence, references_);
    return true;
}

std::vector<Gene> read_germline_fasta(const std::string& path) {
    if (path.empty()) return {};
    SequenceReader reader(path, true);
    SequenceRecord record;
    std::vector<Gene> genes;
    std::unordered_set<std::string> seen;
    while (reader.next(record)) {
        auto name = germline_name(record);
        if (!seen.insert(name).second) {
            throw std::runtime_error("duplicate germline identifier: " + name + " in " + path);
        }
        genes.push_back(Gene{
            std::move(name), record.sequence, infer_locus(record.id + " " + record.description),
            -1, -1, {}});
        if (genes.back().locus.empty()) genes.back().locus = infer_locus(genes.back().name);
    }
    if (genes.empty()) throw std::runtime_error("no germline sequences found in: " + path);
    return genes;
}

}  // namespace swiftig
