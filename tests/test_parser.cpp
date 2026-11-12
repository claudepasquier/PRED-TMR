#include "doctest.h"

#include <string>
#include <vector>

#include "predtmr/parser.h"
#include "predtmr/sequence.h"

TEST_CASE("detect_format recognises FASTA") {
    const std::string text = ">sp|P12345|TEST Protein\nMKTAYIAKQR\n";
    CHECK(predtmr::Parser::parse_string(text).format() == predtmr::SequenceFormat::fasta);
}

TEST_CASE("detect_format recognises SwissProt") {
    const std::string text =
        "ID   TEST_PROTEIN     STANDARD;      PRT;   21 AA.\n"
        "AC   P12345;\n"
        "//\n";
    CHECK(predtmr::Parser::parse_string(text).format() == predtmr::SequenceFormat::swissprot);
}

TEST_CASE("detect_format recognises PDB") {
    const std::string text = "HEADER    TEST PROTEIN\nEND\n";
    CHECK(predtmr::Parser::parse_string(text).format() == predtmr::SequenceFormat::pdb);
}

TEST_CASE("detect_format falls back to plain") {
    const std::string text = "MKTAYIAKQR\n";
    CHECK(predtmr::Parser::parse_string(text).format() == predtmr::SequenceFormat::plain);
}

TEST_CASE("FASTA parsing extracts the name and concatenates sequence lines") {
    const std::string text =
        ">sp|P12345|TEST Protein\n"
        "MKTAYIAKQR\n"
        "QISFVKSHFS\n";
    const predtmr::Sequence sequence = predtmr::Parser::parse_string(text);
    CHECK_EQ(sequence.name(), std::string("sp|P12345|TEST"));
    CHECK_EQ(sequence.data(), std::string("MKTAYIAKQRQISFVKSHFS"));
    CHECK_EQ(sequence.data().size(), static_cast<std::size_t>(20));
}

TEST_CASE("FASTA parsing handles multiple records") {
    const std::string text =
        ">first\n"
        "MKTAYIAKQR\n"
        ">second\n"
        "LLLLLLLLLLLLLLLLLLLLL\n";
    const std::vector<predtmr::Sequence> sequences =
        predtmr::Parser::parse_all_string(text);
    REQUIRE_EQ(sequences.size(), static_cast<std::size_t>(2));
    CHECK_EQ(sequences[0].name(), std::string("first"));
    CHECK_EQ(sequences[1].name(), std::string("second"));
}

TEST_CASE("SwissProt parsing extracts name, sequence and observed transmembrane segment") {
    const std::string text =
        "ID   TEST_PROTEIN     STANDARD;      PRT;   21 AA.\n"
        "AC   P12345;\n"
        "FT   TRANSMEM      5     15\n"
        "SQ   SEQUENCE   21 AA;  2000 MW;  ABCD;\n"
        "     MKTAYIAKQRQISFVKSHFSR\n"
        "//\n";
    const predtmr::Sequence sequence = predtmr::Parser::parse_string(text);
    CHECK_EQ(sequence.name(), std::string("TEST_PROTEIN"));
    CHECK_EQ(sequence.data(), std::string("MKTAYIAKQRQISFVKSHFSR"));
    REQUIRE_EQ(sequence.observed_transmembrane().size(), static_cast<std::size_t>(1));
    CHECK_EQ(sequence.observed_transmembrane()[0].start, 5);
    CHECK_EQ(sequence.observed_transmembrane()[0].end, 15);
}

TEST_CASE("plain parsing concatenates all lines") {
    const std::string text = "MKTAY\nIAKQR\n";
    const predtmr::Sequence sequence = predtmr::Parser::parse_string(text);
    CHECK_EQ(sequence.format(), predtmr::SequenceFormat::plain);
    CHECK_EQ(sequence.data(), std::string("MKTAYIAKQR"));
}

TEST_CASE("PDB parsing reads SEQRES records") {
    const std::string text =
        "HEADER    TEST PROTEIN\n"
        "SEQRES   1 A   21  MET LYS THR ALA TYR ILE ALA LYS GLN ARG GLN ILE SER PHE VAL LYS SER HIS PHE SER ARG\n"
        "END\n";
    const predtmr::Sequence sequence = predtmr::Parser::parse_string(text);
    CHECK_EQ(sequence.format(), predtmr::SequenceFormat::pdb);
    CHECK_EQ(sequence.name(), std::string("TEST"));
    CHECK_EQ(sequence.data(), std::string("MKTAYIAKQRQISFVKSHFSR"));
}
