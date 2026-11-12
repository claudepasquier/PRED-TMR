#pragma once

// Minimal doctest-compatible single-header harness.
// It provides the subset of macros used by the PRED-TMR phase A tests and
// keeps the build fully offline (no FetchContent / network access).

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace doctest {

struct TestCase {
    const char* name;
    void (*function)();
};

inline std::vector<TestCase>& registry() {
    static std::vector<TestCase> tests;
    return tests;
}

inline int& failure_count() {
    static int failures = 0;
    return failures;
}

struct Registrar {
    Registrar(const char* name, void (*function)()) {
        registry().push_back(TestCase{name, function});
    }
};

inline void report_failure(const char* file, int line, const std::string& message) {
    ++failure_count();
    std::cerr << file << ":" << line << ": FAILED: " << message << "\n";
}

}  // namespace doctest

#define DOCTEST_CAT_IMPL(lhs, rhs) lhs##rhs
#define DOCTEST_CAT(lhs, rhs) DOCTEST_CAT_IMPL(lhs, rhs)

#define DOCTEST_TEST_CASE(name)                                                     \
    static void DOCTEST_CAT(doctest_test_function_, __LINE__)();                    \
    [[maybe_unused]] static doctest::Registrar                                      \
        DOCTEST_CAT(doctest_registrar_, __LINE__)(                                  \
            name, &DOCTEST_CAT(doctest_test_function_, __LINE__));                  \
    static void DOCTEST_CAT(doctest_test_function_, __LINE__)()

#define TEST_CASE(name) DOCTEST_TEST_CASE(name)

#define CHECK(condition)                                                            \
    do {                                                                            \
        if (!(condition)) {                                                         \
            doctest::report_failure(__FILE__, __LINE__, #condition);                \
        }                                                                           \
    } while (false)

#define REQUIRE(condition)                                                          \
    do {                                                                            \
        if (!(condition)) {                                                         \
            doctest::report_failure(__FILE__, __LINE__, #condition);                \
            return;                                                                 \
        }                                                                           \
    } while (false)

#define CHECK_EQ(lhs, rhs)                                                          \
    do {                                                                            \
        const auto doctest_lhs = (lhs);                                             \
        const auto doctest_rhs = (rhs);                                             \
        if (!(doctest_lhs == doctest_rhs)) {                                        \
            doctest::report_failure(__FILE__, __LINE__,                             \
                                    std::string(#lhs " == " #rhs));                 \
        }                                                                           \
    } while (false)

#define REQUIRE_EQ(lhs, rhs)                                                        \
    do {                                                                            \
        const auto doctest_lhs = (lhs);                                             \
        const auto doctest_rhs = (rhs);                                             \
        if (!(doctest_lhs == doctest_rhs)) {                                        \
            doctest::report_failure(__FILE__, __LINE__,                             \
                                    std::string(#lhs " == " #rhs));                 \
            return;                                                                 \
        }                                                                           \
    } while (false)

#define CHECK_NEAR(lhs, rhs, epsilon)                                               \
    do {                                                                            \
        const double doctest_lhs = static_cast<double>(lhs);                        \
        const double doctest_rhs = static_cast<double>(rhs);                        \
        if (std::fabs(doctest_lhs - doctest_rhs) > static_cast<double>(epsilon)) {  \
            doctest::report_failure(__FILE__, __LINE__,                             \
                                    std::string(#lhs " ~= " #rhs));                 \
        }                                                                           \
    } while (false)

#define REQUIRE_NEAR(lhs, rhs, epsilon)                                             \
    do {                                                                            \
        const double doctest_lhs = static_cast<double>(lhs);                        \
        const double doctest_rhs = static_cast<double>(rhs);                        \
        if (std::fabs(doctest_lhs - doctest_rhs) > static_cast<double>(epsilon)) {  \
            doctest::report_failure(__FILE__, __LINE__,                             \
                                    std::string(#lhs " ~= " #rhs));                 \
            return;                                                                 \
        }                                                                           \
    } while (false)
