#pragma once

#include <istream>
#include <string>
#include <vector>

#include "predtmr/sequence.h"

namespace predtmr {

class Parser {
public:
    // Reads the whole stream and detects the format from the first non-empty line.
    // The stream is consumed by this call.
    static SequenceFormat detect_format(std::istream& input);

    static Sequence parse(std::istream& input,
                          SequenceFormat format = SequenceFormat::unknown);
    static Sequence parse_string(const std::string& text,
                                 SequenceFormat format = SequenceFormat::unknown);

    static std::vector<Sequence> parse_all(std::istream& input,
                                           SequenceFormat format = SequenceFormat::unknown);
    static std::vector<Sequence> parse_all_string(const std::string& text,
                                                  SequenceFormat format = SequenceFormat::unknown);
};

}  // namespace predtmr
