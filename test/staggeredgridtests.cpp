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

#include "staggeredgrid.hpp"
#include <iostream>
#include <string>

boost::ut::suite staggered = [] {
  using namespace boost::ut;

  "basic 2d staggered grid"_test = [] {
    red3::Uniform1D xg(4);
    red3::Uniform1D yg(5);
    red3::StaggeredGrid grid(xg, yg);
    expect(eq(4, grid.ni) >> fatal);
    expect(eq(5, grid.nu) >> fatal);
    expect(eq(5, grid.nj) >> fatal);
    expect(eq(6, grid.nv) >> fatal);
    expect(eq(1, grid.nk) >> fatal);
    expect(eq(1, grid.nw) >> fatal);
    // Check grid
    expect(eq(5, grid.x.size()) >> fatal);
    std::vector<double> x = { 0.0, 0.25, 0.5, 0.75, 1.0 };
    for(auto i = 0; i < grid.nu; i++) {
      expect(x[i] ==  grid.x[i]) << x[i] << " == " << grid.x[i] << ", i = " << i;
    }
    expect(eq(4, grid.xm.size()) >> fatal);
    std::vector<double> xm = { 0.125, 0.375, 0.625, 0.875 };
    for(auto i = 0; i < grid.ni; i++) {
      expect(xm[i] == grid.xm[i]) << xm[i] << " == " << grid.xm[i] << ", i = " << i;
    }
    expect(eq(6, grid.y.size()) >> fatal);
    std::vector<double> y = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };
    for(auto i = 0; i < grid.nv; i++) {
      expect(close(y[i], grid.y[i], 1.0e-15)) << y[i] << " == " << grid.y[i] << ", i = " << i;
    }
    expect(eq(5, grid.ym.size()) >> fatal);
    std::vector<double> ym = { 0.1, 0.3, 0.5, 0.7, 0.9 };
    for(auto i = 0; i < grid.ni; i++) {
      expect(close(ym[i], grid.ym[i], 1.0e-15)) << ym[i] << " == " << grid.ym[i] << ", i = " << i;
    }
    // Check velocity
    grid.set_u([](double x, double y){return x;});
    grid.set_v([](double x, double y){return -y;});
    for(auto j = 0; j < grid.nj; j++) {
      for(auto i = 0; i < grid.nu; i++) {
        expect(eq(x[i], grid.u(i, j, 0))) << "Index: (" << i << "," << j << ")";
      }
    }
  for(auto j = 0; j < grid.nv; j++) {
    for(auto i = 0; i < grid.ni; i++) {
      expect(close(-y[j], grid.v(i, j, 0), 1.0e-15)) << "Index: (" << i << "," << j << ")";
    }
  }
  // Divergence, etc.
  //red3::Array<red3::StaggeredGrid> g(&grid);
  auto g = grid.p_array();
  grid.dudx(g);
  for(auto j = 0; j < grid.nj; j++) {
    for(auto i = 0; i < grid.ni; i++){
      expect(1.0_d == g(i, j, 0));
    }
  }
  grid.dvdy(g);
  for(auto j = 0; j < grid.nj; j++) {
    for(auto i = 0; i < grid.ni; i++){
      expect(-1.0_d == g(i, j, 0));
    }
  }
  grid.divg(g);
  for(auto j = 0; j < grid.nj; j++) {
    for(auto i = 0; i < grid.ni; i++){
      expect(0.0_d == g(i, j, 0));
    }
  }
  // BC Setters
  grid.set_north_u([](double x){return -x;});
  for(auto i = 0; i < grid.nu; i++) {
    expect(eq(-x[i], grid.north_u(i, 0))); // << ("Index: " + std::to_string(i));
  }

  //auto us = grid.uArray();
  };
};

