#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "predtmr/parser.h"
#include "predtmr/predictor.h"
#include "predtmr/sequence.h"

namespace {

struct Options {
    predtmr::SequenceFormat format{predtmr::SequenceFormat::unknown};
    bool per_residue{false};
    std::string input_path;
};

void print_usage(std::ostream& output) {
    output << "Usage: predtmr [options] [file]\n"
           << "\n"
           << "If no file is given (or if file is '-'), the sequence is read from stdin.\n"
           << "\n"
           << "Options:\n"
           << "  --format <auto|plain|fasta|swissprot|pdb>  Input format (default: auto)\n"
           << "  --output <text|L>                          Output format (default: text)\n"
           << "  --per-residue                              Alias for --output L\n"
           << "  -h, --help                                 Show this help\n";
}

bool parse_format(const std::string& value, predtmr::SequenceFormat& format) {
    if (value == "auto") {
        format = predtmr::SequenceFormat::unknown;
    } else if (value == "plain") {
        format = predtmr::SequenceFormat::plain;
    } else if (value == "fasta") {
        format = predtmr::SequenceFormat::fasta;
    } else if (value == "swissprot" || value == "swiss") {
        format = predtmr::SequenceFormat::swissprot;
    } else if (value == "pdb") {
        format = predtmr::SequenceFormat::pdb;
    } else {
        return false;
    }
    return true;
}

std::vector<int> observed_flags(const predtmr::Sequence& sequence) {
    std::vector<int> flags(sequence.data().size(), 0);
    for (const predtmr::TransmembraneSegment& segment : sequence.observed_transmembrane()) {
        for (int position = segment.start - 1; position <= segment.end - 1; ++position) {
            if (position >= 0 && position < static_cast<int>(flags.size())) {
                flags[static_cast<std::size_t>(position)] = 1;
            }
        }
    }
    return flags;
}

void print_text(const predtmr::Sequence& sequence, const predtmr::Prediction& prediction) {
    std::cout << "Sequence: " << sequence.name() << "\n";
    std::cout << "Length: " << sequence.data().size() << "\n";
    std::cout << "Predicted transmembrane segments:\n";

    if (prediction.segments.empty()) {
        std::cout << "  none\n";
    } else {
        std::cout << "  #  start  end  score\n";
        int index = 0;
        for (const predtmr::Segment& segment : prediction.segments) {
            ++index;
            std::cout << std::setw(3) << index << "  " << std::setw(5) << segment.start
                      << "  " << std::setw(4) << segment.end << "  " << std::setprecision(4)
                      << segment.score << "\n";
        }
    }

    if (prediction.rating.valid) {
        std::cout << "Rating: Q=" << std::setprecision(4) << prediction.rating.q
                  << "% C=" << std::setprecision(3) << prediction.rating.c
                  << " SM=" << std::setprecision(3) << prediction.rating.sm << "\n";
    } else {
        std::cout << "Rating: na\n";
    }
}

void print_per_residue(const predtmr::Sequence& sequence,
                       const predtmr::Prediction& prediction,
                       const std::vector<int>& observed) {
    if (!sequence.name().empty()) {
        std::cout << sequence.name() << "\n";
    }
    for (std::size_t i = 0; i < sequence.data().size(); ++i) {
        std::cout << std::setw(4) << static_cast<int>(i + 1) << "  " << sequence.data()[i]
                  << "  " << prediction.tm[i] << "  " << observed[i] << "\n";
    }
}

}  // namespace

int main(int argc, char* argv[]) {
    Options options;

    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        if (argument == "-h" || argument == "--help") {
            print_usage(std::cout);
            return 0;
        }
        if (argument == "--per-residue") {
            options.per_residue = true;
            continue;
        }
        if (argument == "--format") {
            if (i + 1 >= argc) {
                std::cerr << "predtmr: --format requires a value\n";
                return 1;
            }
            if (!parse_format(argv[++i], options.format)) {
                std::cerr << "predtmr: unknown format '" << argv[i] << "'\n";
                return 1;
            }
            continue;
        }
        if (argument == "--output") {
            if (i + 1 >= argc) {
                std::cerr << "predtmr: --output requires a value\n";
                return 1;
            }
            const std::string value = argv[++i];
            if (value == "text") {
                options.per_residue = false;
            } else if (value == "L" || value == "l") {
                options.per_residue = true;
            } else {
                std::cerr << "predtmr: unknown output format '" << value << "'\n";
                return 1;
            }
            continue;
        }
        if (!argument.empty() && argument[0] == '-') {
            std::cerr << "predtmr: unknown option '" << argument << "'\n";
            return 1;
        }
        if (!options.input_path.empty()) {
            std::cerr << "predtmr: only one input file is supported\n";
            return 1;
        }
        options.input_path = argument;
    }

    std::ifstream file;
    std::istream* input = &std::cin;
    if (!options.input_path.empty() && options.input_path != "-") {
        file.open(options.input_path);
        if (!file.is_open()) {
            std::cerr << "predtmr: cannot open '" << options.input_path << "'\n";
            return 1;
        }
        input = &file;
    }

    const std::vector<predtmr::Sequence> sequences =
        predtmr::Parser::parse_all(*input, options.format);
    if (sequences.empty()) {
        std::cerr << "predtmr: no sequence found in input\n";
        return 1;
    }

    bool first = true;
    for (const predtmr::Sequence& sequence : sequences) {
        const std::vector<int> observed = observed_flags(sequence);
        const predtmr::Prediction prediction = predtmr::predict(sequence.data(), observed);

        if (options.per_residue) {
            if (!first) {
                std::cout << "\n";
            }
            print_per_residue(sequence, prediction, observed);
        } else {
            if (!first) {
                std::cout << "\n";
            }
            print_text(sequence, prediction);
        }
        first = false;
    }

    return 0;
}
