#include "predtmr/parser.h"

#include <cctype>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace predtmr {
namespace {

constexpr char kAllowedResidues[] = "ACDEFGHIKLMNPQRSTVWYX";

std::string read_all(std::istream& input) {
    std::ostringstream buffer;
    buffer << input.rdbuf();
    return buffer.str();
}

std::vector<std::string> split_lines(const std::string& text) {
    std::vector<std::string> lines;
    std::string current;
    for (char character : text) {
        if (character == '\n') {
            lines.push_back(current);
            current.clear();
        } else if (character != '\r') {
            current.push_back(character);
        }
    }
    if (!current.empty()) {
        lines.push_back(current);
    }
    return lines;
}

bool starts_with(const std::string& value, const char* prefix, std::size_t length) {
    return value.size() >= length && value.compare(0, length, prefix, length) == 0;
}

// Replicates the historical cleaning step: alphabetic residues outside the
// allowed alphabet become 'X'; non-alphabetic characters are dropped.
std::string clean_sequence(const std::string& raw) {
    std::string cleaned;
    cleaned.reserve(raw.size());
    for (char character : raw) {
        if (std::isalpha(static_cast<unsigned char>(character)) != 0) {
            if (std::strchr(kAllowedResidues, character) != nullptr) {
                cleaned.push_back(character);
            } else {
                cleaned.push_back('X');
            }
        }
    }
    return cleaned;
}

SequenceFormat detect_format_string(const std::string& text) {
    for (const std::string& line : split_lines(text)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '>') {
            return SequenceFormat::fasta;
        }
        if (starts_with(line, "ID   ", 5)) {
            return SequenceFormat::swissprot;
        }
        if (starts_with(line, "HEADER    ", 10)) {
            return SequenceFormat::pdb;
        }
        return SequenceFormat::plain;
    }
    return SequenceFormat::unknown;
}

std::string first_token(const std::string& value) {
    const std::size_t start = value.find_first_not_of(" \t");
    if (start == std::string::npos) {
        return {};
    }
    const std::size_t end = value.find_first_of(" \t", start);
    return value.substr(start, end == std::string::npos ? std::string::npos : end - start);
}

std::vector<Sequence> parse_fasta(const std::vector<std::string>& lines) {
    std::vector<Sequence> sequences;
    Sequence current;
    std::string raw;
    bool in_record = false;

    auto finalize = [&]() {
        current.set_data(clean_sequence(raw));
        sequences.push_back(current);
        raw.clear();
    };

    for (const std::string& line : lines) {
        if (!line.empty() && line[0] == '>') {
            if (in_record) {
                finalize();
            }
            current = Sequence();
            current.set_format(SequenceFormat::fasta);
            current.set_name(first_token(line.substr(1)));
            in_record = true;
        } else if (in_record) {
            raw += line;
        }
    }
    if (in_record) {
        finalize();
    }
    return sequences;
}

void parse_ft_transmembrane(const std::string& feature, Sequence& sequence) {
    // The historical parser used fixed offsets 14 and 21 on the text that
    // follows the five-character "FT   " prefix.
    int from = 0;
    int to = 0;
    int extent = 0;

    int index = 14;
    while (index >= 0 && index < static_cast<int>(feature.size()) &&
           std::isdigit(static_cast<unsigned char>(feature[static_cast<std::size_t>(index)])) != 0) {
        --index;
    }
    if (index >= 0 && index < static_cast<int>(feature.size()) &&
        feature[static_cast<std::size_t>(index)] == '<') {
        extent |= 1;  // EXT_LEFT
    }
    if (index == 14) {
        from = 0;
    } else {
        from = std::atoi(feature.c_str() + index + 1);
    }

    index = 21;
    while (index >= 0 && index < static_cast<int>(feature.size()) &&
           std::isdigit(static_cast<unsigned char>(feature[static_cast<std::size_t>(index)])) != 0) {
        --index;
    }
    if (index >= 0 && index < static_cast<int>(feature.size()) &&
        feature[static_cast<std::size_t>(index)] == '>') {
        extent |= 2;  // EXT_RIGHT
    }
    if (index == 21) {
        to = 0;
    } else {
        to = std::atoi(feature.c_str() + index + 1);
    }

    if ((from != 0) && (to != 0)) {
        TransmembraneSegment segment;
        segment.start = from;
        segment.end = to;
        segment.extent = extent;
        segment.direction = 0;
        sequence.add_observed_transmembrane(segment);
    }
}

std::vector<Sequence> parse_swissprot(const std::vector<std::string>& lines) {
    std::vector<Sequence> sequences;
    Sequence current;
    std::string raw;
    bool in_record = false;
    bool in_sequence = false;

    auto finalize = [&]() {
        current.set_data(clean_sequence(raw));
        sequences.push_back(current);
        raw.clear();
    };

    for (const std::string& line : lines) {
        if (starts_with(line, "ID   ", 5)) {
            if (in_record) {
                finalize();
            }
            current = Sequence();
            current.set_format(SequenceFormat::swissprot);
            current.set_name(first_token(line.substr(2)));
            in_record = true;
            in_sequence = false;
        } else if (in_record) {
            if (starts_with(line, "FT", 2)) {
                const std::string feature = line.size() > 5 ? line.substr(5) : std::string();
                if (starts_with(feature, "TRANSMEM", 8)) {
                    parse_ft_transmembrane(feature, current);
                }
            } else if (starts_with(line, "SQ", 2)) {
                in_sequence = true;
            } else if (starts_with(line, "//", 2)) {
                finalize();
                in_record = false;
                in_sequence = false;
            } else if (in_sequence && line.size() > 5 && starts_with(line, "  ", 2)) {
                raw += line.substr(5);
            }
        }
    }
    if (in_record) {
        finalize();
    }
    return sequences;
}

char three_letter_to_one(const std::string& residue) {
    if (residue == "ALA") return 'A';
    if (residue == "ARG") return 'R';
    if (residue == "ASN") return 'N';
    if (residue == "ASP") return 'D';
    if (residue == "CYS") return 'C';
    if (residue == "GLN") return 'Q';
    if (residue == "GLU") return 'E';
    if (residue == "GLY") return 'G';
    if (residue == "HIS") return 'H';
    if (residue == "ILE") return 'I';
    if (residue == "LEU") return 'L';
    if (residue == "LYS") return 'K';
    if (residue == "MET") return 'M';
    if (residue == "PHE") return 'F';
    if (residue == "PRO") return 'P';
    if (residue == "SER") return 'S';
    if (residue == "THR") return 'T';
    if (residue == "TRP") return 'W';
    if (residue == "TYR") return 'Y';
    if (residue == "VAL") return 'V';
    if (residue == "MSE") return 'M';
    return 'X';
}

std::string pdb_sequence_from_seqres(const std::vector<std::string>& lines) {
    std::string sequence;
    for (const std::string& line : lines) {
        if (!starts_with(line, "SEQRES", 6)) {
            continue;
        }
        std::size_t position = 19;
        while (position < line.size()) {
            while (position < line.size() &&
                   std::isspace(static_cast<unsigned char>(line[position])) != 0) {
                ++position;
            }
            if (position + 3 > line.size()) {
                break;
            }
            const std::string residue = line.substr(position, 3);
            sequence.push_back(three_letter_to_one(residue));
            position += 3;
        }
    }
    return sequence;
}

std::string pdb_sequence_from_atoms(const std::vector<std::string>& lines) {
    std::string sequence;
    std::string last_key;
    for (const std::string& line : lines) {
        if (!starts_with(line, "ATOM", 4) || line.size() < 26) {
            continue;
        }
        const std::string atom_name = line.substr(12, 4);
        if (atom_name.find("CA") == std::string::npos) {
            continue;
        }
        const std::string residue = line.substr(17, 3);
        const std::string key = line.substr(21, 5);
        if (key == last_key) {
            continue;
        }
        last_key = key;
        sequence.push_back(three_letter_to_one(residue));
    }
    return sequence;
}

std::vector<Sequence> parse_pdb(const std::vector<std::string>& lines) {
    Sequence sequence;
    sequence.set_format(SequenceFormat::pdb);

    for (const std::string& line : lines) {
        if (starts_with(line, "HEADER", 6)) {
            sequence.set_name(first_token(line.size() > 10 ? line.substr(10) : std::string()));
            break;
        }
    }

    std::string raw = pdb_sequence_from_seqres(lines);
    if (raw.empty()) {
        raw = pdb_sequence_from_atoms(lines);
    }
    sequence.set_data(clean_sequence(raw));

    std::vector<Sequence> sequences;
    sequences.push_back(sequence);
    return sequences;
}

std::vector<Sequence> parse_plain(const std::vector<std::string>& lines) {
    std::string raw;
    for (const std::string& line : lines) {
        raw += line;
    }

    Sequence sequence;
    sequence.set_format(SequenceFormat::plain);
    sequence.set_data(clean_sequence(raw));

    std::vector<Sequence> sequences;
    sequences.push_back(sequence);
    return sequences;
}

}  // namespace

SequenceFormat Parser::detect_format(std::istream& input) {
    return detect_format_string(read_all(input));
}

Sequence Parser::parse(std::istream& input, SequenceFormat format) {
    return parse_string(read_all(input), format);
}

Sequence Parser::parse_string(const std::string& text, SequenceFormat format) {
    const std::vector<Sequence> sequences = parse_all_string(text, format);
    if (sequences.empty()) {
        return Sequence();
    }
    return sequences.front();
}

std::vector<Sequence> Parser::parse_all(std::istream& input, SequenceFormat format) {
    return parse_all_string(read_all(input), format);
}

std::vector<Sequence> Parser::parse_all_string(const std::string& text, SequenceFormat format) {
    const std::vector<std::string> lines = split_lines(text);
    SequenceFormat resolved = format;
    if (resolved == SequenceFormat::unknown) {
        resolved = detect_format_string(text);
    }

    switch (resolved) {
        case SequenceFormat::fasta:
            return parse_fasta(lines);
        case SequenceFormat::swissprot:
            return parse_swissprot(lines);
        case SequenceFormat::pdb:
            return parse_pdb(lines);
        case SequenceFormat::plain:
            return parse_plain(lines);
        case SequenceFormat::unknown:
        default:
            return {};
    }
}

}  // namespace predtmr
