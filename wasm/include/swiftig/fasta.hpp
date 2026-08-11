#pragma once

#include <fstream>
#include <istream>
#include <memory>
#include <string>

#include "swiftig/types.hpp"

namespace swiftig {

class SequenceReader {
public:
    explicit SequenceReader(const std::string& path, bool references = false);
    SequenceReader(const SequenceReader&) = delete;
    SequenceReader& operator=(const SequenceReader&) = delete;
    bool next(SequenceRecord& record);

private:
    enum class Format { Unknown, Fasta, Fastq };
    std::unique_ptr<std::ifstream> owned_;
    std::istream* input_ = nullptr;
    Format format_ = Format::Unknown;
    bool references_ = false;
    std::string pending_;
    std::size_t record_number_ = 0;

    bool next_fasta(SequenceRecord& record);
    bool next_fastq(SequenceRecord& record);
};

std::vector<Gene> read_germline_fasta(const std::string& path);

}  // namespace swiftig
