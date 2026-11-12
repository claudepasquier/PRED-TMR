#pragma once

#include <string>
#include <vector>

namespace predtmr {

// Per-residue indicator arrays produced by the historical ComputeIndicators().
// The vectors are zero-initialised and padded so that the predictor can safely
// look ahead up to 40 residues past the end of the sequence, exactly like the
// original 10000-element global arrays.
struct Indicators {
    std::vector<float> start;
    std::vector<float> end;
    std::vector<float> hydro;
};

Indicators compute_indicators(const std::string& sequence);

}  // namespace predtmr
