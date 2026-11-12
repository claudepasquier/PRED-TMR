#include "doctest.h"

#include <string>
#include <vector>

#include "predtmr/indicators.h"
#include "predtmr/predictor.h"

TEST_CASE("compute_indicators returns zeroed padded arrays for short sequences") {
    const predtmr::Indicators indicators = predtmr::compute_indicators("MKTAY");
    REQUIRE_EQ(indicators.start.size(), static_cast<std::size_t>(45));
    REQUIRE_EQ(indicators.end.size(), static_cast<std::size_t>(45));
    REQUIRE_EQ(indicators.hydro.size(), static_cast<std::size_t>(45));
    for (float value : indicators.start) CHECK_EQ(value, 0.0F);
    for (float value : indicators.end) CHECK_EQ(value, 0.0F);
    for (float value : indicators.hydro) CHECK_EQ(value, 0.0F);
}

TEST_CASE("compute_indicators fills the expected windows for an 11-residue sequence") {
    const predtmr::Indicators indicators = predtmr::compute_indicators("MKTAYIAKQRQ");
    REQUIRE_EQ(indicators.start.size(), static_cast<std::size_t>(51));
    CHECK(indicators.start[5] > 0.0F);
    CHECK(indicators.end[5] > 0.0F);
    CHECK(indicators.hydro[0] > 0.0F);
    CHECK(indicators.hydro[5] > 0.0F);
    CHECK(indicators.hydro[10] > 0.0F);
}

TEST_CASE("predict returns no segment for sequences shorter than 11") {
    const predtmr::Prediction prediction = predtmr::predict("MKTAY");
    CHECK(prediction.segments.empty());
    REQUIRE_EQ(prediction.tm.size(), static_cast<std::size_t>(5));
    for (int value : prediction.tm) CHECK_EQ(value, 0);
    CHECK(!prediction.rating.valid);
}

TEST_CASE("predict marks a strongly hydrophobic stretch") {
    const std::string sequence(60, 'L');
    const predtmr::Prediction prediction = predtmr::predict(sequence);
    REQUIRE(!prediction.segments.empty());
    for (const predtmr::Segment& segment : prediction.segments) {
        CHECK(segment.start >= 1);
        CHECK(segment.end <= static_cast<int>(sequence.size()));
        CHECK(segment.start <= segment.end);
    }
    bool has_tm = false;
    for (int value : prediction.tm) {
        if (value == 1) has_tm = true;
    }
    CHECK(has_tm);
}

TEST_CASE("rating reproduces the historical Matthews computation") {
    const std::vector<int> predicted = {1, 1, 0, 0};
    const std::vector<int> observed = {1, 0, 0, 0};
    const predtmr::Rating rating = predtmr::rate(predicted, observed);
    REQUIRE(rating.valid);
    CHECK_NEAR(rating.q, 83.3333F, 1e-3);
    CHECK_NEAR(rating.c, 0.57735F, 1e-3);
    CHECK_NEAR(rating.sm, 1.0F, 1e-3);
}

TEST_CASE("rating is invalid when no observed transmembrane segment exists") {
    const std::vector<int> predicted = {1, 1, 0, 0};
    const std::vector<int> observed = {0, 0, 0, 0};
    const predtmr::Rating rating = predtmr::rate(predicted, observed);
    CHECK(!rating.valid);
}
