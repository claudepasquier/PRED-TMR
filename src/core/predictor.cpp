#include "predtmr/predictor.h"

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <vector>

namespace predtmr {
namespace {

constexpr float kHydroCutoff = 7.5F;

}  // namespace

Prediction predict(const std::string& sequence, const std::vector<int>& observed_tm) {
    Prediction prediction;
    const int nb_aa = static_cast<int>(sequence.size());
    prediction.tm.assign(static_cast<std::size_t>(nb_aa), 0);

    // ProcessFile() bypasses ComputeIndicators()/MarkTM() for sequences < 11.
    if (nb_aa < 11) {
        return prediction;
    }

    const Indicators indicators = compute_indicators(sequence);
    const std::size_t padded = indicators.start.size();

    // The historical code used int for bestStr (truncating the float score)
    // and float for bestMiddle. Keep those exact types.
    std::vector<int> best_str(padded, 0);
    std::vector<int> end_indx(padded, 0);
    std::vector<float> best_middle(padded, 0.0F);

    for (int i = 0; i < nb_aa; ++i) {
        float start_str = 0.0F;
        float middle_str = 0.0F;
        float total_tm = 0.0F;
        int nb_aa_in_tm = 0;

        // test if it is a good start
        start_str = std::sqrt(indicators.start[static_cast<std::size_t>(i)] /
                              indicators.end[static_cast<std::size_t>(i)]);
        if (start_str > indicators.start[static_cast<std::size_t>(i)]) {
            start_str = indicators.start[static_cast<std::size_t>(i)];
        }
        start_str *= indicators.start[static_cast<std::size_t>(i)];
        start_str = std::sqrt(start_str);

        for (int j = i; j < i + 40; ++j) {
            // The original used log() (double) and fexp() (float exp).
            middle_str = static_cast<float>(
                static_cast<double>(middle_str) +
                std::log(static_cast<double>(indicators.hydro[static_cast<std::size_t>(j)])));
            total_tm += indicators.hydro[static_cast<std::size_t>(j)];
            ++nb_aa_in_tm;

            float end_str = std::sqrt(indicators.end[static_cast<std::size_t>(j + 1)] /
                                      indicators.start[static_cast<std::size_t>(j + 1)]);
            if (end_str > indicators.end[static_cast<std::size_t>(j + 1)]) {
                end_str = indicators.end[static_cast<std::size_t>(j + 1)];
            }
            end_str *= indicators.end[static_cast<std::size_t>(j + 1)];
            end_str = std::sqrt(end_str);

            const int length = j - i + 1;
            const float norm_middle =
                std::exp(middle_str * 10.0F / static_cast<float>(length));
            const float length_penalty =
                std::exp(static_cast<float>(std::abs(j - i + 1 - 21)));
            const float tot_str = start_str + norm_middle + end_str - length_penalty;

            if (start_str < 1.0F) continue;
            if (end_str < 1.0F) continue;
            if (norm_middle < 1.0F) continue;
            if (tot_str > static_cast<float>(best_str[static_cast<std::size_t>(i)])) {
                best_str[static_cast<std::size_t>(i)] = static_cast<int>(tot_str);
                end_indx[static_cast<std::size_t>(i)] = j;
                best_middle[static_cast<std::size_t>(i)] = norm_middle;
            }
        }
    }

    // Select the TM parts greedily (non-overlapping).
    std::vector<int> tm(padded, 0);
    while (true) {
        int imax = -1;
        for (int i = 0; i < nb_aa; ++i) {
            if (tm[static_cast<std::size_t>(i)] > 0) continue;
            if (best_middle[static_cast<std::size_t>(i)] < kHydroCutoff) continue;
            if ((imax == -1) ||
                (best_str[static_cast<std::size_t>(i)] >
                 best_str[static_cast<std::size_t>(imax)])) {
                if (tm[static_cast<std::size_t>(end_indx[static_cast<std::size_t>(i)])] == 0) {
                    imax = i;
                }
            }
        }
        if (imax == -1) break;
        for (int i = imax; i <= end_indx[static_cast<std::size_t>(imax)]; ++i) {
            tm[static_cast<std::size_t>(i)] = 1;
        }
    }

    for (int i = 0; i < nb_aa; ++i) {
        prediction.tm[static_cast<std::size_t>(i)] = tm[static_cast<std::size_t>(i)];
    }

    // Build the 1-based segment list.
    int in_tm = 0;
    int segment_start = 0;
    for (int i = 0; i < nb_aa; ++i) {
        if ((in_tm == 0) && (tm[static_cast<std::size_t>(i)] == 1)) {
            in_tm = 1;
            segment_start = i;
        } else if ((in_tm == 1) &&
                   ((i + 1 >= nb_aa) || (tm[static_cast<std::size_t>(i + 1)] == 0))) {
            Segment segment;
            segment.start = segment_start + 1;
            segment.end = i + 1;
            segment.score = best_middle[static_cast<std::size_t>(segment_start)];
            prediction.segments.push_back(segment);
            in_tm = 0;
        }
    }

    prediction.rating = rate(prediction.tm, observed_tm);
    return prediction;
}

Rating rate(const std::vector<int>& predicted_tm, const std::vector<int>& observed_tm) {
    Rating rating;

    float good_pred = 0.0F;      // w
    float good_non_pred = 0.0F;  // x
    float missed_pred = 0.0F;    // y
    float wrong_pred = 0.0F;     // z

    int in_pred_tm = 0;
    int in_obs_tm = 0;

    int nb_obs_tm = 0;
    int nb_exp_m = 0;
    int nb_pred_tm = 0;
    int nb_pred_m = 0;

    int exp_match = 0;
    int pred_match = 0;
    int total = 0;

    const int nb_aa = static_cast<int>(predicted_tm.size());
    for (int i = 0; i < nb_aa; ++i) {
        const int tm = predicted_tm[static_cast<std::size_t>(i)];
        const int obs = (i < static_cast<int>(observed_tm.size()))
                            ? observed_tm[static_cast<std::size_t>(i)]
                            : 0;

        if ((tm == 1) && (obs == 1)) good_pred += 1.0F;
        if ((tm == 0) && (obs == 0)) good_non_pred += 1.0F;
        if ((tm == 0) && (obs == 1)) missed_pred += 1.0F;
        if ((tm == 1) && (obs == 0)) wrong_pred += 1.0F;

        if ((in_obs_tm == 0) && (obs == 1)) {
            ++nb_obs_tm;
            in_obs_tm = 1;
            exp_match = 0;
        }
        if ((in_pred_tm == 0) && (tm == 1)) {
            ++nb_pred_tm;
            in_pred_tm = 1;
            pred_match = 0;
        }
        if ((in_obs_tm == 1) && (exp_match == 0) && (tm == 1)) {
            exp_match = 1;
            ++nb_exp_m;
        }
        if ((in_pred_tm == 1) && (pred_match == 0) && (obs == 1)) {
            pred_match = 1;
            ++nb_pred_m;
        }

        if ((in_obs_tm == 1) && (obs == 0)) in_obs_tm = 0;
        if ((in_pred_tm == 1) && (tm == 0)) in_pred_tm = 0;

        ++total;
    }

    if (nb_obs_tm == 0) {
        rating.valid = false;
        return rating;
    }

    const float pa = 100.0F * good_pred / (good_pred + missed_pred);
    const float pna = 100.0F * good_non_pred / (good_non_pred + wrong_pred);
    rating.q = (pa + pna) / 2.0F;

    const float s = (good_pred + missed_pred) / static_cast<float>(total);
    const float p = (good_pred + wrong_pred) / static_cast<float>(total);
    if ((p == 0.0F) || (s == 0.0F)) {
        rating.c = 0.0F;
    } else {
        rating.c = ((good_pred / static_cast<float>(total)) - (p * s)) /
                   std::sqrt(p * s * (1.0F - p) * (1.0F - s));
    }

    const float m = (nb_obs_tm == 0)
                        ? 0.0F
                        : static_cast<float>(nb_exp_m) / static_cast<float>(nb_obs_tm);
    const float c = (nb_pred_tm == 0)
                        ? 0.0F
                        : static_cast<float>(nb_pred_m) / static_cast<float>(nb_pred_tm);
    rating.sm = std::sqrt(m * c);
    rating.valid = true;
    return rating;
}

}  // namespace predtmr
