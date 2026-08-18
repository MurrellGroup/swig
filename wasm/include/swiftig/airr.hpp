#pragma once

#include <iosfwd>
#include <string>

#include "swiftig/types.hpp"

namespace swiftig {

void write_airr_header(std::ostream& output);
void write_airr_record(std::ostream& output, const Annotation& annotation);

// Allocation-light byte-identical writer used by web/CLI batch annotation.
// The ostream functions remain the reference implementation.
void append_airr_header(std::string& output);
void append_airr_record(std::string& output, const Annotation& annotation);

}  // namespace swiftig
