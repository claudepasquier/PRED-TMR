// WebAssembly bindings for the PRED-TMR kernel (phase B).
//
// Exposed JavaScript API (embind, MODULARIZE + EXPORT_NAME=createPredTMR):
//
//   const module = await createPredTMR();
//   const result = module.predict("MNS...");          // raw residues
//   // -> { segments: [{start,end,score}], perResidue: [0|1], length: N }
//
//   const parsed = module.parseSequence(text, "auto"); // FASTA/SwissProt/PDB/plain
//   // -> { name, sequence, length, format }
//
// The bindings are a thin, stateless adapter over the validated C++ kernel
// (predtmr::predict / predtmr::Parser). No biological logic lives here.

#include <emscripten/bind.h>
#include <emscripten/val.h>

#include <string>

#include "predtmr/parser.h"
#include "predtmr/predictor.h"
#include "predtmr/sequence.h"

namespace {

using emscripten::val;

// Maps a JavaScript-friendly format string onto the C++ enum.
predtmr::SequenceFormat to_format(const std::string& value) {
    if (value == "plain") {
        return predtmr::SequenceFormat::plain;
    }
    if (value == "fasta") {
        return predtmr::SequenceFormat::fasta;
    }
    if (value == "swissprot" || value == "swiss") {
        return predtmr::SequenceFormat::swissprot;
    }
    if (value == "pdb") {
        return predtmr::SequenceFormat::pdb;
    }
    return predtmr::SequenceFormat::unknown;
}

// Parses an arbitrary text payload (FASTA, SwissProt, PDB or raw residues)
// using the C++ parser, so the browser and the CLI share identical semantics.
val parse_sequence(const std::string& text, const std::string& format) {
    const predtmr::Sequence sequence =
        predtmr::Parser::parse_string(text, to_format(format));

    val result = val::object();
    result.set("name", sequence.name());
    result.set("sequence", sequence.data());
    result.set("length", sequence.molecule_length());
    result.set("format", std::string(predtmr::format_name(sequence.format())));
    return result;
}

// Runs the validated kernel on a raw residue string and returns the segments
// plus the per-residue 0/1 transmembrane flag vector.
val predict_sequence(const std::string& sequence) {
    const predtmr::Prediction prediction = predtmr::predict(sequence);

    val segments = val::array();
    for (const predtmr::Segment& segment : prediction.segments) {
        val entry = val::object();
        entry.set("start", segment.start);
        entry.set("end", segment.end);
        entry.set("score", segment.score);
        segments.call<void>("push", entry);
    }

    val per_residue = val::array();
    for (const int flag : prediction.tm) {
        per_residue.call<void>("push", flag);
    }

    val result = val::object();
    result.set("segments", segments);
    result.set("perResidue", per_residue);
    result.set("length", static_cast<int>(sequence.size()));
    return result;
}

}  // namespace

EMSCRIPTEN_BINDINGS(predtmr_phase_b) {
    emscripten::function("predict", &predict_sequence);
    emscripten::function("parseSequence", &parse_sequence);
}
