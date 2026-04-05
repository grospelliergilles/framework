#define BOOST_TEST_MODULE TestSkylineLU
#include <boost/test/unit_test.hpp>

#include <arcane/alina/IO.h>
#include <arcane/alina/Adapters.h>
#include <arcane/alina/profiler.h>
#include "sample_problem.h"

using namespace Arcane;

namespace Arcane::Alina {
    profiler<> prof;
}

BOOST_AUTO_TEST_SUITE( test_io )

BOOST_AUTO_TEST_CASE(io_mm)
{
    std::vector<ptrdiff_t> ptr, ptr2;
    std::vector<ptrdiff_t> col, col2;
    std::vector<double>    val, val2;
    std::vector<double>    rhs, rhs2;

    size_t n = sample_problem(16, val, col, ptr, rhs);

    auto A = std::tie(n, ptr, col, val);

    Alina::IO::mm_write("test_io_crs.mm", A);
    Alina::IO::mm_write("test_io_vec.mm", rhs.data(), n, 1);

    size_t rows, cols;
    std::tie(rows, cols) = Alina::IO::mm_reader("test_io_crs.mm")(ptr2, col2, val2);

    BOOST_REQUIRE_EQUAL(n, rows);
    BOOST_REQUIRE_EQUAL(n, cols);
    BOOST_REQUIRE_EQUAL(ptr.back(), ptr2.back());
    for(size_t i = 0; i < n; ++i) {
        BOOST_REQUIRE_EQUAL(ptr[i], ptr2[i]);
        for(ptrdiff_t j = ptr[i], e = ptr[i+1]; j < e; ++j) {
            BOOST_CHECK_EQUAL(col[j], col2[j]);
            BOOST_CHECK_SMALL(val[j] - val2[j], 1e-12);
        }
    }

    std::tie(rows, cols) = Alina::IO::mm_reader("test_io_vec.mm")(rhs2);
    BOOST_REQUIRE_EQUAL(n, rows);
    BOOST_REQUIRE_EQUAL(1, cols);
    for(size_t i = 0; i < n; ++i) {
        BOOST_CHECK_SMALL(rhs[i] - rhs2[i], 1e-12);
    }
}

BOOST_AUTO_TEST_SUITE_END()

