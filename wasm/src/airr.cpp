#include "swiftig/airr.hpp"

#include <algorithm>
#include <charconv>
#include <cstdio>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

namespace swiftig {
namespace {

const std::vector<std::string>& columns() {
    static const std::vector<std::string> names{
        "sequence_id", "sequence", "quality", "sequence_aa", "rev_comp", "productive",
        "vj_in_frame", "stop_codon", "complete_vdj", "locus",
        "v_call", "d_call", "j_call", "c_call",
        "sequence_alignment", "sequence_alignment_aa", "germline_alignment",
        "germline_alignment_aa", "junction", "junction_aa", "np1", "np2",
        "cdr1", "cdr1_aa", "cdr2", "cdr2_aa", "cdr3", "cdr3_aa",
        "fwr1", "fwr1_aa", "fwr2", "fwr2_aa", "fwr3", "fwr3_aa", "fwr4", "fwr4_aa",
        "v_score", "v_identity", "v_support", "v_cigar",
        "d_score", "d_identity", "d_support", "d_cigar",
        "j_score", "j_identity", "j_support", "j_cigar",
        "c_score", "c_identity", "c_support", "c_cigar",
        "v_sequence_start", "v_sequence_end", "v_germline_start", "v_germline_end",
        "d_sequence_start", "d_sequence_end", "d_germline_start", "d_germline_end",
        "j_sequence_start", "j_sequence_end", "j_germline_start", "j_germline_end",
        "c_sequence_start", "c_sequence_end", "c_germline_start", "c_germline_end",
        "cdr1_start", "cdr1_end", "cdr2_start", "cdr2_end", "cdr3_start", "cdr3_end",
        "fwr1_start", "fwr1_end", "fwr2_start", "fwr2_end",
        "fwr3_start", "fwr3_end", "fwr4_start", "fwr4_end",
        "v_sequence_alignment", "d_sequence_alignment", "j_sequence_alignment", "c_sequence_alignment",
        "v_germline_alignment", "d_germline_alignment", "j_germline_alignment", "c_germline_alignment",
        "junction_length", "junction_aa_length", "np1_length", "np2_length",
        "v_frameshift", "j_frameshift", "d_frame",
        "sequence_frame", "region_definition", "v_annotation_source", "j_annotation_source",
        "v_alternatives", "d_alternatives", "j_alternatives", "c_alternatives"};
    return names;
}

std::string boolean(bool value) { return value ? "T" : "F"; }

std::string optional_boolean(const std::optional<bool>& value) {
    return value ? boolean(*value) : std::string{};
}

std::string number(double value) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(6) << value;
    auto result = stream.str();
    while (!result.empty() && result.back() == '0') result.pop_back();
    if (!result.empty() && result.back() == '.') result.pop_back();
    return result;
}

std::string optional_position(const std::optional<std::size_t>& value) {
    return value ? std::to_string(*value) : std::string{};
}

std::string hit_name(const std::optional<SegmentHit>& hit) {
    return hit ? (hit->call.empty() ? hit->gene->name : hit->call) : std::string{};
}

std::string hit_score(const std::optional<SegmentHit>& hit) {
    return hit ? std::to_string(hit->alignment.score) : std::string{};
}

std::string hit_identity(const std::optional<SegmentHit>& hit) {
    return hit ? number(hit->alignment.identity()) : std::string{};
}

std::string hit_support(const std::optional<SegmentHit>& hit) {
    if (!hit || !hit->support) return {};
    char buffer[48];
    const int written = std::snprintf(buffer, sizeof(buffer), "%.6e", *hit->support);
    return written > 0 && static_cast<std::size_t>(written) < sizeof(buffer)
        ? std::string(buffer, static_cast<std::size_t>(written)) : std::string{};
}

std::string hit_cigar(const std::optional<SegmentHit>& hit) {
    return hit ? hit->alignment.cigar : std::string{};
}

std::string alternative_evidence(const std::vector<SegmentHit>& hits) {
    std::ostringstream stream;
    for (std::size_t i = 0; i < hits.size(); ++i) {
        if (i) stream << ';';
        const auto& hit = hits[i];
        stream << hit.gene->name << '|' << hit.alignment.score << '|'
               << number(hit.alignment.identity()) << '|'
               << hit.alignment.query_start + 1 << '|' << hit.alignment.query_end << '|'
               << hit.alignment.reference_start + 1 << '|' << hit.alignment.reference_end;
    }
    return stream.str();
}

std::string query_start(const std::optional<SegmentHit>& hit) {
    return hit ? std::to_string(hit->alignment.query_start + 1) : std::string{};
}

std::string query_end(const std::optional<SegmentHit>& hit) {
    return hit ? std::to_string(hit->alignment.query_end) : std::string{};
}

std::string reference_start(const std::optional<SegmentHit>& hit) {
    return hit ? std::to_string(hit->alignment.reference_start + 1) : std::string{};
}

std::string reference_end(const std::optional<SegmentHit>& hit) {
    return hit ? std::to_string(hit->alignment.reference_end) : std::string{};
}

std::string query_alignment(const std::optional<SegmentHit>& hit) {
    return hit ? hit->alignment.aligned_query : std::string{};
}

std::string germline_alignment(const std::optional<SegmentHit>& hit) {
    return hit ? hit->alignment.aligned_reference : std::string{};
}

std::string np1_length(const Annotation& annotation) {
    if (!annotation.v || !annotation.j) return {};
    if (annotation.d) {
        const auto a = annotation.v->alignment.query_end;
        const auto b = annotation.d->alignment.query_start;
        return std::to_string(b > a ? b - a : 0);
    }
    const auto a = annotation.v->alignment.query_end;
    const auto b = annotation.j->alignment.query_start;
    return std::to_string(b > a ? b - a : 0);
}

std::string np2_length(const Annotation& annotation) {
    if (!annotation.d || !annotation.j) return {};
    const auto a = annotation.d->alignment.query_end;
    const auto b = annotation.j->alignment.query_start;
    return std::to_string(b > a ? b - a : 0);
}

void write_row(std::ostream& output, const std::vector<std::string>& values) {
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i) output.put('\t');
        output << values[i];
    }
    output.put('\n');
}

class FastRow {
public:
    explicit FastRow(std::string& output) : output_(output) {}

    void field(std::string_view value) {
        separator();
        output_.append(value.data(), value.size());
    }

    void field(bool value) { field(value ? std::string_view{"T"} : std::string_view{"F"}); }

    template <typename Integer>
    void integer(Integer value) {
        separator();
        append_integer(value);
    }

    void optional_boolean(const std::optional<bool>& value) {
        if (value) field(*value);
        else field(std::string_view{});
    }

    template <typename Integer>
    void optional_integer(const std::optional<Integer>& value) {
        if (value) integer(*value);
        else field(std::string_view{});
    }

    void identity(double value) {
        separator();
        append_identity(value);
    }

    void alternatives(const std::vector<SegmentHit>& hits) {
        separator();
        for (std::size_t i = 0; i < hits.size(); ++i) {
            if (i) output_.push_back(';');
            const auto& hit = hits[i];
            output_.append(hit.gene->name);
            output_.push_back('|');
            append_integer(hit.alignment.score);
            output_.push_back('|');
            append_identity(hit.alignment.identity());
            output_.push_back('|');
            append_integer(hit.alignment.query_start + 1);
            output_.push_back('|');
            append_integer(hit.alignment.query_end);
            output_.push_back('|');
            append_integer(hit.alignment.reference_start + 1);
            output_.push_back('|');
            append_integer(hit.alignment.reference_end);
        }
    }

    void support(double value) {
        separator();
        char buffer[48];
        const int written = std::snprintf(buffer, sizeof(buffer), "%.6e", value);
        if (written > 0 && static_cast<std::size_t>(written) < sizeof(buffer)) {
            output_.append(buffer, buffer + written);
        }
    }

    void finish() { output_.push_back('\n'); }

private:
    template <typename Integer>
    void append_integer(Integer value) {
        char buffer[32];
        const auto converted = std::to_chars(std::begin(buffer), std::end(buffer), value);
        if (converted.ec == std::errc{}) output_.append(buffer, converted.ptr);
    }

    void append_identity(double value) {
        char buffer[64];
        const int written = std::snprintf(buffer, sizeof(buffer), "%.6f", value);
        if (written <= 0 || static_cast<std::size_t>(written) >= sizeof(buffer)) return;
        char* end = buffer + written;
        while (end > buffer && end[-1] == '0') --end;
        if (end > buffer && end[-1] == '.') --end;
        output_.append(buffer, end);
    }
    void separator() {
        if (!first_) output_.push_back('\t');
        first_ = false;
    }

    std::string& output_;
    bool first_ = true;
};

void fast_hit_name(FastRow& row, const std::optional<SegmentHit>& hit) {
    if (hit) row.field(hit->call.empty() ? hit->gene->name : hit->call);
    else row.field(std::string_view{});
}

void fast_hit_score(FastRow& row, const std::optional<SegmentHit>& hit) {
    if (hit) row.integer(hit->alignment.score);
    else row.field(std::string_view{});
}

void fast_hit_identity(FastRow& row, const std::optional<SegmentHit>& hit) {
    if (hit) row.identity(hit->alignment.identity());
    else row.field(std::string_view{});
}

void fast_hit_support(FastRow& row, const std::optional<SegmentHit>& hit) {
    if (hit && hit->support) row.support(*hit->support);
    else row.field(std::string_view{});
}

void fast_hit_cigar(FastRow& row, const std::optional<SegmentHit>& hit) {
    if (hit) row.field(hit->alignment.cigar);
    else row.field(std::string_view{});
}

void fast_query_start(FastRow& row, const std::optional<SegmentHit>& hit) {
    if (hit) row.integer(hit->alignment.query_start + 1);
    else row.field(std::string_view{});
}

void fast_query_end(FastRow& row, const std::optional<SegmentHit>& hit) {
    if (hit) row.integer(hit->alignment.query_end);
    else row.field(std::string_view{});
}

void fast_reference_start(FastRow& row, const std::optional<SegmentHit>& hit) {
    if (hit) row.integer(hit->alignment.reference_start + 1);
    else row.field(std::string_view{});
}

void fast_reference_end(FastRow& row, const std::optional<SegmentHit>& hit) {
    if (hit) row.integer(hit->alignment.reference_end);
    else row.field(std::string_view{});
}

void fast_query_alignment(FastRow& row, const std::optional<SegmentHit>& hit) {
    if (hit) row.field(hit->alignment.aligned_query);
    else row.field(std::string_view{});
}

void fast_germline_alignment(FastRow& row, const std::optional<SegmentHit>& hit) {
    if (hit) row.field(hit->alignment.aligned_reference);
    else row.field(std::string_view{});
}

void fast_alternatives(FastRow& row, const std::vector<SegmentHit>& hits) {
    row.alternatives(hits);
}

}  // namespace

void write_airr_header(std::ostream& output) { write_row(output, columns()); }

void write_airr_record(std::ostream& output, const Annotation& a) {
    std::vector<std::string> values{
        a.sequence_id, a.sequence, a.quality, a.sequence_aa, boolean(a.rev_comp), optional_boolean(a.productive),
        optional_boolean(a.vj_in_frame), optional_boolean(a.stop_codon), optional_boolean(a.complete_vdj), a.locus,
        hit_name(a.v), hit_name(a.d), hit_name(a.j), hit_name(a.c),
        a.sequence_alignment, a.sequence_alignment_aa, a.germline_alignment,
        a.germline_alignment_aa, a.junction, a.junction_aa, a.np1, a.np2,
        a.cdr1.sequence, a.cdr1.sequence_aa, a.cdr2.sequence, a.cdr2.sequence_aa,
        a.cdr3, a.cdr3_aa,
        a.fwr1.sequence, a.fwr1.sequence_aa, a.fwr2.sequence, a.fwr2.sequence_aa,
        a.fwr3.sequence, a.fwr3.sequence_aa, a.fwr4.sequence, a.fwr4.sequence_aa,
        hit_score(a.v), hit_identity(a.v), hit_support(a.v), hit_cigar(a.v),
        hit_score(a.d), hit_identity(a.d), hit_support(a.d), hit_cigar(a.d),
        hit_score(a.j), hit_identity(a.j), hit_support(a.j), hit_cigar(a.j),
        hit_score(a.c), hit_identity(a.c), hit_support(a.c), hit_cigar(a.c),
        query_start(a.v), query_end(a.v), reference_start(a.v), reference_end(a.v),
        query_start(a.d), query_end(a.d), reference_start(a.d), reference_end(a.d),
        query_start(a.j), query_end(a.j), reference_start(a.j), reference_end(a.j),
        query_start(a.c), query_end(a.c), reference_start(a.c), reference_end(a.c),
        optional_position(a.cdr1.start), optional_position(a.cdr1.end),
        optional_position(a.cdr2.start), optional_position(a.cdr2.end),
        optional_position(a.cdr3_start), optional_position(a.cdr3_end),
        optional_position(a.fwr1.start), optional_position(a.fwr1.end),
        optional_position(a.fwr2.start), optional_position(a.fwr2.end),
        optional_position(a.fwr3.start), optional_position(a.fwr3.end),
        optional_position(a.fwr4.start), optional_position(a.fwr4.end),
        query_alignment(a.v), query_alignment(a.d), query_alignment(a.j), query_alignment(a.c),
        germline_alignment(a.v), germline_alignment(a.d), germline_alignment(a.j), germline_alignment(a.c),
        a.junction.empty() ? std::string{} : std::to_string(a.junction.size()),
        a.junction_aa.empty() ? std::string{} : std::to_string(a.junction_aa.size()),
        np1_length(a), np2_length(a),
        optional_boolean(a.v_frameshift), optional_boolean(a.j_frameshift),
        a.d_frame ? std::to_string(*a.d_frame) : std::string{},
        a.sequence_frame ? std::to_string(*a.sequence_frame) : std::string{},
        a.region_definition, a.v_annotation_source, a.j_annotation_source,
        alternative_evidence(a.v_alternatives), alternative_evidence(a.d_alternatives),
        alternative_evidence(a.j_alternatives), alternative_evidence(a.c_alternatives)};
    write_row(output, values);
}

void append_airr_header(std::string& output) {
    FastRow row(output);
    for (const auto& column : columns()) row.field(column);
    row.finish();
}

void append_airr_record(std::string& output, const Annotation& a) {
    FastRow row(output);
    row.field(a.sequence_id); row.field(a.sequence); row.field(a.quality); row.field(a.sequence_aa);
    row.field(a.rev_comp); row.optional_boolean(a.productive); row.optional_boolean(a.vj_in_frame);
    row.optional_boolean(a.stop_codon); row.optional_boolean(a.complete_vdj); row.field(a.locus);
    fast_hit_name(row, a.v); fast_hit_name(row, a.d); fast_hit_name(row, a.j); fast_hit_name(row, a.c);
    row.field(a.sequence_alignment); row.field(a.sequence_alignment_aa); row.field(a.germline_alignment);
    row.field(a.germline_alignment_aa); row.field(a.junction); row.field(a.junction_aa);
    row.field(a.np1); row.field(a.np2);
    row.field(a.cdr1.sequence); row.field(a.cdr1.sequence_aa);
    row.field(a.cdr2.sequence); row.field(a.cdr2.sequence_aa); row.field(a.cdr3); row.field(a.cdr3_aa);
    row.field(a.fwr1.sequence); row.field(a.fwr1.sequence_aa);
    row.field(a.fwr2.sequence); row.field(a.fwr2.sequence_aa);
    row.field(a.fwr3.sequence); row.field(a.fwr3.sequence_aa);
    row.field(a.fwr4.sequence); row.field(a.fwr4.sequence_aa);
    fast_hit_score(row, a.v); fast_hit_identity(row, a.v); fast_hit_support(row, a.v); fast_hit_cigar(row, a.v);
    fast_hit_score(row, a.d); fast_hit_identity(row, a.d); fast_hit_support(row, a.d); fast_hit_cigar(row, a.d);
    fast_hit_score(row, a.j); fast_hit_identity(row, a.j); fast_hit_support(row, a.j); fast_hit_cigar(row, a.j);
    fast_hit_score(row, a.c); fast_hit_identity(row, a.c); fast_hit_support(row, a.c); fast_hit_cigar(row, a.c);
    fast_query_start(row, a.v); fast_query_end(row, a.v); fast_reference_start(row, a.v); fast_reference_end(row, a.v);
    fast_query_start(row, a.d); fast_query_end(row, a.d); fast_reference_start(row, a.d); fast_reference_end(row, a.d);
    fast_query_start(row, a.j); fast_query_end(row, a.j); fast_reference_start(row, a.j); fast_reference_end(row, a.j);
    fast_query_start(row, a.c); fast_query_end(row, a.c); fast_reference_start(row, a.c); fast_reference_end(row, a.c);
    row.optional_integer(a.cdr1.start); row.optional_integer(a.cdr1.end);
    row.optional_integer(a.cdr2.start); row.optional_integer(a.cdr2.end);
    row.optional_integer(a.cdr3_start); row.optional_integer(a.cdr3_end);
    row.optional_integer(a.fwr1.start); row.optional_integer(a.fwr1.end);
    row.optional_integer(a.fwr2.start); row.optional_integer(a.fwr2.end);
    row.optional_integer(a.fwr3.start); row.optional_integer(a.fwr3.end);
    row.optional_integer(a.fwr4.start); row.optional_integer(a.fwr4.end);
    fast_query_alignment(row, a.v); fast_query_alignment(row, a.d);
    fast_query_alignment(row, a.j); fast_query_alignment(row, a.c);
    fast_germline_alignment(row, a.v); fast_germline_alignment(row, a.d);
    fast_germline_alignment(row, a.j); fast_germline_alignment(row, a.c);
    if (a.junction.empty()) row.field(std::string_view{}); else row.integer(a.junction.size());
    if (a.junction_aa.empty()) row.field(std::string_view{}); else row.integer(a.junction_aa.size());
    if (!a.v || !a.j) row.field(std::string_view{});
    else if (a.d) {
        const auto left = a.v->alignment.query_end;
        const auto right = a.d->alignment.query_start;
        row.integer(right > left ? right - left : 0);
    } else {
        const auto left = a.v->alignment.query_end;
        const auto right = a.j->alignment.query_start;
        row.integer(right > left ? right - left : 0);
    }
    if (!a.d || !a.j) row.field(std::string_view{});
    else {
        const auto left = a.d->alignment.query_end;
        const auto right = a.j->alignment.query_start;
        row.integer(right > left ? right - left : 0);
    }
    row.optional_boolean(a.v_frameshift); row.optional_boolean(a.j_frameshift);
    row.optional_integer(a.d_frame); row.optional_integer(a.sequence_frame);
    row.field(a.region_definition); row.field(a.v_annotation_source); row.field(a.j_annotation_source);
    fast_alternatives(row, a.v_alternatives); fast_alternatives(row, a.d_alternatives);
    fast_alternatives(row, a.j_alternatives); fast_alternatives(row, a.c_alternatives);
    row.finish();
}

}  // namespace swiftig
