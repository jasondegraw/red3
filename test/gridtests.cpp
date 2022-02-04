// Copyright (C) 2015-2022 Jason W. DeGraw
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
#ifdef __cpp_modules
#ifdef _MSC_VER
#define BOOST_UT_DISABLE_MODULE
#include <boost/ut.hpp>
#else
import boost.ut; // Doesn't appear to work yet with MSVC/CMake
#endif
#else
#include <boost/ut.hpp>
#endif

#include "red3/grid.hpp"

boost::ut::suite grids = [] {
  using namespace boost::ut;

  "uniform grid 1d"_test = [] {
    red3::Uniform1D uniform_x(0.25, 4);
    std::optional<red3::Grid1D> opt_x = red3::Grid1D::generate(uniform_x);
    expect((!!opt_x) >> fatal);
    red3::Grid1D x = opt_x.value();
    expect((5_i == x.size()) >> fatal);
    expect(0.25_d == x.delta(0));
    // Check grid
    expect(0.0_d == x[0]);
    expect(0.25_d == x[1]);
    expect(0.5_d == x[2]);
    expect(0.75_d == x[3]);
    expect(1.0_d == x[4]);

    red3::Uniform1D uniform_y(4, 1.0);
    std::optional<red3::Grid1D> opt_y = red3::Grid1D::generate(uniform_y);
    expect((!!opt_y) >> fatal);
    red3::Grid1D y = opt_y.value();
    expect((5_i == y.size()) >> fatal);
    expect(0.25_d == y.delta(0));
    // Check grid
    expect(0.0_d == y[0]);
    expect(0.25_d == y[1]);
    expect(0.5_d == y[2]);
    expect(0.75_d == y[3]);
    expect(1.0_d == y[4]);
  };
};


