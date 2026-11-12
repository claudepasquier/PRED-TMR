#include "predtmr/indicators.h"

#include <cstddef>
#include <cstring>

namespace predtmr {
namespace {

// Exactly the historical AA alphabet (21 residues + terminating NUL).
constexpr char kAA[] = "ACDEFGHIKLMNPQRSTVWYX";

// Values computed from SwissProt 35 (69113 sequences).
constexpr double kFreqOfAA[20] = {
    0.0757, 0.0168, 0.0529, 0.0634, 0.0408, 0.0683, 0.0223, 0.0578, 0.0595, 0.0939,
    0.0235, 0.0450, 0.0491, 0.0400, 0.0515, 0.0715, 0.0570, 0.0655, 0.0124, 0.0318};

// Values computed from all TM sequences in SP35 (excluding non-ended TM).
constexpr double kFreqOfAAInTM[20] = {
    0.1047, 0.0202, 0.0081, 0.0083, 0.0912, 0.0791, 0.0088, 0.1204, 0.0064, 0.1732,
    0.0353, 0.0175, 0.0293, 0.0109, 0.0064, 0.0576, 0.0501, 0.1150, 0.0222, 0.0342};

// Copied verbatim from PRED-TMR.cpp (float precision preserved).
constexpr float kTable1[20][10] = {
    {1.06902F, 0.989605F, 1.04213F, 0.922376F, 0.793477F, 1.30621F, 1.43313F, 1.36232F, 1.49283F, 1.52134F},
    {0.617683F, 0.657246F, 0.698084F, 0.787419F, 0.740837F, 0.63555F, 0.791247F, 0.807838F, 0.873563F, 1.05287F},
    {0.749396F, 0.664802F, 0.628002F, 0.795278F, 1.1542F, 0.101321F, 0.138839F, 0.215786F, 0.168232F, 0.154133F},
    {0.69983F, 0.639244F, 0.610494F, 0.644839F, 0.92809F, 0.0984046F, 0.146449F, 0.185811F, 0.125418F, 0.140082F},
    {1.21061F, 1.35475F, 1.24021F, 0.993785F, 0.396258F, 2.59287F, 2.33389F, 2.16701F, 2.20977F, 2.05875F},
    {1.0399F, 1.03483F, 0.976442F, 1.04985F, 0.865299F, 0.937955F, 1.10824F, 1.1364F, 1.22088F, 1.13827F},
    {0.992616F, 0.996372F, 0.993689F, 1.09724F, 1.42722F, 0.300467F, 0.379877F, 0.628835F, 0.556401F, 0.442653F},
    {0.996274F, 1.05736F, 1.0528F, 0.827607F, 0.301093F, 2.71576F, 2.18674F, 2.11153F, 2.10811F, 2.19085F},
    {0.999763F, 0.970133F, 0.941309F, 1.11062F, 1.59539F, 0.086673F, 0.123156F, 0.173749F, 0.116908F, 0.118923F},
    {1.01348F, 1.06762F, 1.09744F, 0.891112F, 0.534479F, 2.02294F, 1.90231F, 1.7696F, 1.66832F, 1.9093F},
    {1.42454F, 1.08239F, 1.091F, 1.12543F, 0.729484F, 1.6408F, 1.58593F, 1.81188F, 1.53805F, 1.48856F},
    {0.937446F, 0.986294F, 0.996174F, 1.11994F, 1.44267F, 0.276898F, 0.374594F, 0.426187F, 0.352914F, 0.417679F},
    {0.925405F, 0.902998F, 0.914081F, 0.909022F, 0.936969F, 0.514622F, 0.514863F, 0.601597F, 0.600634F, 0.563049F},
    {0.901526F, 0.902744F, 0.882041F, 0.890566F, 1.213F, 0.178418F, 0.257275F, 0.362316F, 0.268844F, 0.319081F},
    {1.31799F, 1.28724F, 1.29549F, 1.69623F, 2.38123F, 0.153249F, 0.125999F, 0.170248F, 0.157748F, 0.163248F},
    {0.987555F, 0.992227F, 0.967736F, 1.00898F, 1.19654F, 0.409588F, 0.700746F, 0.716698F, 0.780183F, 0.720726F},
    {1.00485F, 1.10271F, 1.0757F, 1.14897F, 1.22712F, 0.65175F, 0.80384F, 0.862949F, 0.944639F, 1.0053F},
    {0.931682F, 0.93327F, 0.959481F, 0.788117F, 0.419177F, 1.99204F, 1.75872F, 1.66202F, 1.75157F, 1.78374F},
    {1.66315F, 1.84874F, 1.77081F, 1.83336F, 1.04485F, 2.87411F, 2.51318F, 2.10098F, 2.09176F, 1.45397F},
    {1.14929F, 1.17764F, 1.46939F, 1.3367F, 1.14338F, 1.37371F, 1.42096F, 1.37568F, 1.45718F, 1.02841F}};

constexpr float kTable2[20][10] = {
    {1.52134F, 1.49283F, 1.36232F, 1.43313F, 1.30621F, 0.793477F, 0.922376F, 1.04213F, 0.989605F, 1.06902F},
    {1.05287F, 0.873563F, 0.807838F, 0.791247F, 0.63555F, 0.740837F, 0.787419F, 0.698084F, 0.657246F, 0.617683F},
    {0.154133F, 0.168232F, 0.215786F, 0.138839F, 0.101321F, 1.1542F, 0.795278F, 0.628002F, 0.664802F, 0.749396F},
    {0.140082F, 0.125418F, 0.185811F, 0.146449F, 0.0984046F, 0.92809F, 0.644839F, 0.610494F, 0.639244F, 0.69983F},
    {2.05875F, 2.20977F, 2.16701F, 2.33389F, 2.59287F, 0.396258F, 0.993785F, 1.24021F, 1.35475F, 1.21061F},
    {1.13827F, 1.22088F, 1.1364F, 1.10824F, 0.937955F, 0.865299F, 1.04985F, 0.976442F, 1.03483F, 1.0399F},
    {0.442653F, 0.556401F, 0.628835F, 0.379877F, 0.300467F, 1.42722F, 1.09724F, 0.993689F, 0.996372F, 0.992616F},
    {2.19085F, 2.10811F, 2.11153F, 2.18674F, 2.71576F, 0.301093F, 0.827607F, 1.0528F, 1.05736F, 0.996274F},
    {0.118923F, 0.116908F, 0.173749F, 0.123156F, 0.086673F, 1.59539F, 1.11062F, 0.941309F, 0.970133F, 0.999763F},
    {1.9093F, 1.66832F, 1.7696F, 1.90231F, 2.02294F, 0.534479F, 0.891112F, 1.09744F, 1.06762F, 1.01348F},
    {1.48856F, 1.53805F, 1.81188F, 1.58593F, 1.6408F, 0.729484F, 1.12543F, 1.091F, 1.08239F, 1.42454F},
    {0.417679F, 0.352914F, 0.426187F, 0.374594F, 0.276898F, 1.44267F, 1.11994F, 0.996174F, 0.986294F, 0.937446F},
    {0.563049F, 0.600634F, 0.601597F, 0.514863F, 0.514622F, 0.936969F, 0.909022F, 0.914081F, 0.902998F, 0.925405F},
    {0.319081F, 0.268844F, 0.362316F, 0.257275F, 0.178418F, 1.213F, 0.890566F, 0.882041F, 0.902744F, 0.901526F},
    {0.163248F, 0.157748F, 0.170248F, 0.125999F, 0.153249F, 2.38123F, 1.69623F, 1.29549F, 1.28724F, 1.31799F},
    {0.720726F, 0.780183F, 0.716698F, 0.700746F, 0.409588F, 1.19654F, 1.00898F, 0.967736F, 0.992227F, 0.987555F},
    {1.0053F, 0.944639F, 0.862949F, 0.80384F, 0.65175F, 1.22712F, 1.14897F, 1.0757F, 1.10271F, 1.00485F},
    {1.78374F, 1.75157F, 1.66202F, 1.75872F, 1.99204F, 0.419177F, 0.788117F, 0.959481F, 0.93327F, 0.931682F},
    {1.45397F, 2.09176F, 2.10098F, 2.51318F, 2.87411F, 1.04485F, 1.83336F, 1.77081F, 1.84874F, 1.66315F},
    {1.02841F, 1.45718F, 1.37568F, 1.42096F, 1.37371F, 1.14338F, 1.3367F, 1.46939F, 1.17764F, 1.14929F}};

int residue_index(char residue) noexcept {
    const char* found = std::strchr(kAA, residue);
    if (found == nullptr) {
        // The historical code would dereference a null pointer here. Unknown
        // residues are mapped to the X slot (index 20) for safety.
        return 20;
    }
    return static_cast<int>(found - kAA);
}

float hydrophobicity_ratio(char residue) noexcept {
    const int index = residue_index(residue);
    if (index >= 20) {
        // The historical code read past the end of the 20-entry arrays for X.
        // A neutral ratio keeps the port memory-safe; sequences without X are
        // bit-for-bit identical to the original.
        return 1.0F;
    }
    return static_cast<float>(kFreqOfAAInTM[index] / kFreqOfAA[index]);
}

// Factor from one of the 20-row position tables (kTable1/kTable2). Residues
// without a table row (X and anything outside the 20-letter alphabet) get a
// neutral 1.0F factor instead of an out-of-bounds read, which previously made
// the result depend on the memory layout of the target (native vs WASM).
float table_factor(const float (*table)[10], char residue, int column) noexcept {
    const int index = residue_index(residue);
    if (index >= 20) {
        return 1.0F;
    }
    return table[index][column];
}

}  // namespace

Indicators compute_indicators(const std::string& sequence) {
    const std::size_t length = sequence.size();
    // The original used 10000-element global arrays. MarkTM() reads up to 40
    // positions past the sequence end, so pad accordingly and zero-initialise.
    const std::size_t padded = length + 40;

    Indicators indicators;
    indicators.start.assign(padded, 0.0F);
    indicators.end.assign(padded, 0.0F);
    indicators.hydro.assign(padded, 0.0F);

    // ProcessFile() skips ComputeIndicators() for sequences shorter than 11.
    if (length < 11) {
        return indicators;
    }

    for (int i = 0; i < 5; ++i) {
        indicators.hydro[static_cast<std::size_t>(i)] =
            hydrophobicity_ratio(sequence[static_cast<std::size_t>(i)]);
    }

    for (std::size_t i = 0; i + 10 < length; ++i) {
        float indice1 = 1.0F;
        for (int j = 0; j < 10; ++j) {
            indice1 *= table_factor(kTable1, sequence[i + static_cast<std::size_t>(j)], j);
        }

        float indice2 = 1.0F;
        for (int j = 0; j < 10; ++j) {
            indice2 *= table_factor(kTable2, sequence[i + static_cast<std::size_t>(j)], j);
        }

        const float indice3 = hydrophobicity_ratio(sequence[i + 5]);

        indicators.start[i + 5] = indice1;
        indicators.end[i + 5] = indice2;
        indicators.hydro[i + 5] = indice3;
    }

    for (std::size_t i = length - 5; i < length; ++i) {
        indicators.hydro[i] = hydrophobicity_ratio(sequence[i]);
    }

    return indicators;
}

}  // namespace predtmr
