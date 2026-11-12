#include "predtmr/sequence.h"

namespace predtmr {

const char* format_name(SequenceFormat format) noexcept {
    switch (format) {
        case SequenceFormat::plain:
            return "plain";
        case SequenceFormat::fasta:
            return "fasta";
        case SequenceFormat::swissprot:
            return "swissprot";
        case SequenceFormat::pdb:
            return "pdb";
        case SequenceFormat::unknown:
        default:
            return "unknown";
    }
}

}  // namespace predtmr
