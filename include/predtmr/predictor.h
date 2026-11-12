#pragma once

#include <string>
#include <vector>

#include "predtmr/indicators.h"

namespace predtmr {

struct Segment {
    int start{0};   // 1-based inclusive
    int end{0};     // 1-based inclusive
    float score{0.0F};
};

struct Rating {
    float q{0.0F};
    float c{0.0F};
    float sm{0.0F};
    bool valid{false};
};

struct Prediction {
    std::vector<Segment> segments;
    std::vector<int> tm;  // per-residue 0/1 flags, length == sequence length
    Rating rating;
};

// Pure port of MarkTM(). observed_tm is an optional per-residue 0/1 vector used
// only to compute the rating; it does not influence the prediction.
Prediction predict(const std::string& sequence,
                   const std::vector<int>& observed_tm = {});

// Pure port of rating(). Returns valid == false when no observed TM exists.
Rating rate(const std::vector<int>& predicted_tm,
            const std::vector<int>& observed_tm);

}  // namespace predtmr
