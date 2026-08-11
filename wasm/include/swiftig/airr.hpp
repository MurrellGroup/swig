#pragma once

#include <iosfwd>

#include "swiftig/types.hpp"

namespace swiftig {

void write_airr_header(std::ostream& output);
void write_airr_record(std::ostream& output, const Annotation& annotation);

}  // namespace swiftig
