#include "doctest.h"

int main() {
    for (const doctest::TestCase& test : doctest::registry()) {
        test.function();
    }

    if (doctest::failure_count() == 0) {
        std::cout << "All tests passed\n";
        return 0;
    }

    std::cerr << doctest::failure_count() << " test(s) failed\n";
    return 1;
}
