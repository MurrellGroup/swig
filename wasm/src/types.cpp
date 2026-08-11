#include "swiftig/types.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <string_view>

namespace swiftig {
namespace {

char complement(char base) {
    switch (base) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'R': return 'Y';
        case 'Y': return 'R';
        case 'S': return 'S';
        case 'W': return 'W';
        case 'K': return 'M';
        case 'M': return 'K';
        case 'B': return 'V';
        case 'D': return 'H';
        case 'H': return 'D';
        case 'V': return 'B';
        default: return 'N';
    }
}

int base_index(char base) {
    switch (base) {
        case 'T': return 0;
        case 'C': return 1;
        case 'A': return 2;
        case 'G': return 3;
        default: return -1;
    }
}

}  // namespace

std::string normalize_dna(const std::string& sequence, bool strip_gaps) {
    std::string normalized;
    normalized.reserve(sequence.size());
    for (unsigned char raw : sequence) {
        if (std::isspace(raw)) continue;
        char base = static_cast<char>(std::toupper(raw));
        if (base == 'U') base = 'T';
        if (base == '-' || base == '.') {
            if (!strip_gaps) normalized.push_back('N');
            continue;
        }
        switch (base) {
            case 'A': case 'C': case 'G': case 'T':
            case 'R': case 'Y': case 'S': case 'W': case 'K': case 'M':
            case 'B': case 'D': case 'H': case 'V': case 'N':
                normalized.push_back(base);
                break;
            default:
                normalized.push_back('N');
                break;
        }
    }
    return normalized;
}

std::string reverse_complement(const std::string& sequence) {
    std::string result(sequence.size(), 'N');
    for (std::size_t i = 0; i < sequence.size(); ++i) {
        const auto raw = static_cast<unsigned char>(sequence[sequence.size() - i - 1]);
        result[i] = complement(static_cast<char>(std::toupper(raw)));
    }
    return result;
}

std::string infer_locus(const std::string& gene_name) {
    std::string upper = gene_name;
    std::transform(upper.begin(), upper.end(), upper.begin(), [](unsigned char c) {
        return static_cast<char>(std::toupper(c));
    });
    static constexpr std::array<std::string_view, 7> loci{
        "IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"};
    for (const auto locus : loci) {
        if (upper.find(locus) != std::string::npos) return std::string(locus);
    }
    return {};
}

std::string translate_dna(const std::string& sequence, std::size_t frame, bool keep_gaps) {
    static constexpr std::string_view table =
        "FFLLSSSSYY**CC*W"
        "LLLLPPPPHHQQRRRR"
        "IIIMTTTTNNKKSSRR"
        "VVVVAAAADDEEGGGG";
    std::string protein;
    if (frame >= sequence.size()) return protein;
    protein.reserve((sequence.size() - frame) / 3);
    for (std::size_t i = frame; i + 2 < sequence.size(); i += 3) {
        if (keep_gaps && (sequence[i] == '-' || sequence[i + 1] == '-' || sequence[i + 2] == '-')) {
            protein.push_back('-');
            continue;
        }
        const int a = base_index(static_cast<char>(std::toupper(static_cast<unsigned char>(sequence[i]))));
        const int b = base_index(static_cast<char>(std::toupper(static_cast<unsigned char>(sequence[i + 1]))));
        const int c = base_index(static_cast<char>(std::toupper(static_cast<unsigned char>(sequence[i + 2]))));
        if (a < 0 || b < 0 || c < 0) {
            protein.push_back('X');
        } else {
            protein.push_back(table[static_cast<std::size_t>(a * 16 + b * 4 + c)]);
        }
    }
    return protein;
}

}  // namespace swiftig
