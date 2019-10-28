// Copyright (C) 2015-2019 Jason W. DeGraw
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

#include "utilities.hpp"
#include "staggeredgrid.hpp"
#include "catch.hpp"
#include <iostream>
#include <string>
#include <cmath>

TEST_CASE("Basic I Tests", "[tridai]")
{
  // 0 - 1 - 2 - 3 - 4
  //
  // a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = 0
  //
  // i = 1
  // a[1]*x[0] + b[0]*x[0] + c[0]*x[1] = 0
  // b[0]*x[0] + c[0]*x[1] = -a[1]*x[0]
  //
  // i = 3
  // a[3]*x[2] + b[3]*x[3] + c[3]*x[4] = 0
  // a[3]*x[2] + b[3]*x[3] = -c[3]*x[4]
  
  red3::StaggeredGrid grid(red3::Grid1D(5), red3::Grid1D(2), red3::Grid1D(2));
  auto g = grid.pArray();
  std::vector<double> a{ {-0.5, -0.5, -0.5, -0.5, -0.5} };
  std::vector<double> b{ { 1.0, 1.0, 1.0, 1.0, 1.0 } };
  std::vector<double> c{ {-0.5, -0.5, -0.5, -0.5, -0.5} };

  // The BCs
  for (int j = 0; j < 2; ++j) {
    for (int k = 0; k < 2; ++k) {
      // Set values
      g(0, j, k) = 0.0;
      g(4, j, k) = 4.0;

      // Modify equations
      g(1, j, k) = -a[1] * g(0, j, k);
      g(3, j, k) = -c[3] * g(4, j, k);
    }
  }

  // Finish modification of equations
  a[1] = 0.0;
  c[3] = 0.0;

  red3::tridai(a, b, c, g, (red3::index_t)1, (red3::index_t)3, (red3::index_t)0, (red3::index_t)1, (red3::index_t)0, (red3::index_t)1,
    grid.ni, grid.nj, grid.nk);

  //for (red3::index_t i = 0; i < 5; i++) {
  //  std::cout << g(i, 0, 0) << ' ' << g(i, 0, 1) << ' ' << g(i, 1, 0) << ' ' << g(i, 1, 1) << std::endl;
  //}

  //for (auto v : b) {
  //  std::cout << v << std::endl;
  //}

  for (red3::index_t j = 0; j < 2; ++j) {
    for (red3::index_t k = 0; k < 2; ++k) {
      REQUIRE(g(1, j, k) == Approx(1.0));
      REQUIRE(g(2, j, k) == Approx(2.0));
      REQUIRE(g(3, j, k) == Approx(3.0));
    }
  }

}

TEST_CASE("Basic Reverse I Tests", "[tridai]")
{
  // 0 - 1 - 2 - 3 - 4
  //
  // a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = 0
  //
  // i = 1
  // a[1]*x[0] + b[0]*x[0] + c[0]*x[1] = 0
  // b[0]*x[0] + c[0]*x[1] = -a[1]*x[0]
  //
  // i = 3
  // a[3]*x[2] + b[3]*x[3] + c[3]*x[4] = 0
  // a[3]*x[2] + b[3]*x[3] = -c[3]*x[4]

  red3::StaggeredGrid grid(red3::Grid1D(5), red3::Grid1D(2), red3::Grid1D(2));
  auto g = grid.pArray();
  std::vector<double> a{ { -0.5, -0.5, -0.5, -0.5, -0.5 } };
  std::vector<double> b{ { 1.0, 1.0, 1.0, 1.0, 1.0 } };
  std::vector<double> c{ { -0.5, -0.5, -0.5, -0.5, -0.5 } };

  // The BCs
  for (red3::index_t j = 0; j < 2; ++j) {
    for (red3::index_t k = 0; k < 2; ++k) {
      // Set values
      g(0, j, k) = 4.0;
      g(4, j, k) = 0.0;

      // Modify equations
      g(1, j, k) = -a[1] * g(0, j, k);
      g(3, j, k) = -c[3] * g(4, j, k);
    }
  }

  // Finish modification of equations
  a[1] = 0.0;
  c[3] = 0.0;

  red3::tridai(a, b, c, g, (red3::index_t)1, (red3::index_t)3, (red3::index_t)0, (red3::index_t)1, (red3::index_t)0, (red3::index_t)1,
    grid.ni, grid.nj, grid.nk);

  for (red3::index_t j = 0; j < 2; ++j) {
    for (red3::index_t k = 0; k < 2; ++k) {
      REQUIRE(g(1, j, k) == Approx(3.0));
      REQUIRE(g(2, j, k) == Approx(2.0));
      REQUIRE(g(3, j, k) == Approx(1.0));
    }
  }

}

TEST_CASE("Basic J Tests", "[tridaj]")
{
  // 0 - 1 - 2 - 3 - 4
  //
  // a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = 0
  //
  // i = 1
  // a[1]*x[0] + b[0]*x[0] + c[0]*x[1] = 0
  // b[0]*x[0] + c[0]*x[1] = -a[1]*x[0]
  //
  // i = 3
  // a[3]*x[2] + b[3]*x[3] + c[3]*x[4] = 0
  // a[3]*x[2] + b[3]*x[3] = -c[3]*x[4]

  red3::StaggeredGrid grid(red3::Grid1D(2), red3::Grid1D(5), red3::Grid1D(2));
  auto g = grid.pArray();
  std::vector<double> a{ { -0.5, -0.5, -0.5, -0.5, -0.5 } };
  std::vector<double> b{ { 1.0, 1.0, 1.0, 1.0, 1.0 } };
  std::vector<double> c{ { -0.5, -0.5, -0.5, -0.5, -0.5 } };

  // The BCs
  for (red3::index_t i = 0; i < 2; ++i) {
    for (red3::index_t k = 0; k < 2; ++k) {
      // Set values
      g(i, 0, k) = 0.0;
      g(i, 4, k) = 4.0;

      // Modify equations
      g(i, 1, k) = -a[1] * g(i, 0, k);
      g(i, 3, k) = -c[3] * g(i, 4, k);
    }
  }

  // Finish modification of equations
  a[1] = 0.0;
  c[3] = 0.0;

  red3::tridaj(a, b, c, g, (red3::index_t)0, (red3::index_t)1, (red3::index_t)1, (red3::index_t)3, (red3::index_t)0, (red3::index_t)1,
    grid.ni, grid.nj, grid.nk);

  //for (int j = 0; j < 5; ++j) {
  //  std::cout << g(0, j, 0) << ' ' << g(0, j, 1) << ' ' << g(1, j, 0) << ' ' << g(1, j, 1) << std::endl;
  //}

  //for (auto v : b) {
  //  std::cout << v << std::endl;
  //}

  for (int i = 0; i < 2; ++i) {
    for (int k = 0; k < 2; ++k) {
      REQUIRE(g(i, 1, k) == Approx(1.0));
      REQUIRE(g(i, 2, k) == Approx(2.0));
      REQUIRE(g(i, 3, k) == Approx(3.0));
    }
  }

}

TEST_CASE("Basic Reverse J Tests", "[tridaj]")
{
  // 0 - 1 - 2 - 3 - 4
  //
  // a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = 0
  //
  // i = 1
  // a[1]*x[0] + b[0]*x[0] + c[0]*x[1] = 0
  // b[0]*x[0] + c[0]*x[1] = -a[1]*x[0]
  //
  // i = 3
  // a[3]*x[2] + b[3]*x[3] + c[3]*x[4] = 0
  // a[3]*x[2] + b[3]*x[3] = -c[3]*x[4]

  red3::StaggeredGrid grid(red3::Grid1D(2), red3::Grid1D(5), red3::Grid1D(2));
  auto g = grid.pArray();
  std::vector<double> a{ { -0.5, -0.5, -0.5, -0.5, -0.5 } };
  std::vector<double> b{ { 1.0, 1.0, 1.0, 1.0, 1.0 } };
  std::vector<double> c{ { -0.5, -0.5, -0.5, -0.5, -0.5 } };

  // The BCs
  for (int i = 0; i < 2; ++i) {
    for (int k = 0; k < 2; ++k) {
      // Set values
      g(i, 0, k) = 4.0;
      g(i, 4, k) = 0.0;

      // Modify equations
      g(i, 1, k) = -a[1] * g(i, 0, k);
      g(i, 3, k) = -c[3] * g(i, 4, k);
    }
  }

  // Finish modification of equations
  a[1] = 0.0;
  c[3] = 0.0;

  red3::tridaj(a, b, c, g, (red3::index_t)0, (red3::index_t)1, (red3::index_t)1, (red3::index_t)3, (red3::index_t)0, (red3::index_t)1,
    grid.ni, grid.nj, grid.nk);

  //for (int j = 0; j < 5; ++j) {
  //  std::cout << g(0, j, 0) << ' ' << g(0, j, 1) << ' ' << g(1, j, 0) << ' ' << g(1, j, 1) << std::endl;
  //}

  //for (auto v : b) {
  //  std::cout << v << std::endl;
  //}

  for (int i = 0; i < 2; ++i) {
    for (int k = 0; k < 2; ++k) {
      REQUIRE(g(i, 1, k) == Approx(3.0));
      REQUIRE(g(i, 2, k) == Approx(2.0));
      REQUIRE(g(i, 3, k) == Approx(1.0));
    }
  }
}

TEST_CASE("Basic Periodic K Tests", "[tripak]")
{
  // 0 - 1 - 2 - 3 - 4 - 5 - 6 - 7 - 8
  //
  // d2f/dx2 = sin(pi*x/2)
  //
  // f = -(pi/2)^(-2) sin(pi*x/2)
  // df/dx = -pi/2 cos(pi*x/2)
  // d2f/dx2 = sin(pi*x/2)
  //
  //
  // a[i]*f[i-1] + b[i]*f[i] + c[i]*f[i+1] = 0
  //
  // i = 1
  // a[1]*x[0] + b[0]*x[0] + c[0]*x[1] = 0
  // b[0]*x[0] + c[0]*x[1] = -a[1]*x[0]
  //
  // i = 7
  // a[7]*x[6] + b[7]*x[7] + c[7]*x[8] = 0
  // a[7]*x[6] + b[7]*x[7] = -c[7]*x[8]

  red3::StaggeredGrid grid(red3::Grid1D(2), red3::Grid1D(2), red3::Grid1D(9));
  auto g = grid.pArray();
  std::vector<double> a{ { -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5 } };
  std::vector<double> b{ { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 } };
  std::vector<double> c{ { -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5 } };
  double pi = acos(-1.0);

  // The BCs
  for (int i = 1; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int k = 1; k < 9; ++k) {
        // Set values
        g(i, j, k) = sin(pi*((double)k)*0.5);
      }
    }
  }

  //for (int i = 0; i < 2; ++i) {
  //  for (int j = 0; j < 2; ++j) {
  ///    // Modify equations
  ///    g(i, j, 1) = -a[1] * g(i, j, 0);
  //    g(i, j, 7) = -c[7] * g(i, j, 8);
  //  }
  //}

  // Finish modification of equations
  //a[1] = 0.0;
  //c[7] = 0.0;

  //red3::utilities::tripak(a, b, c, g, (int)1, (int)3, (int)0, (int)1, (int)0, (int)1,
  //  grid.ni, grid.nj, grid.nk);

  for (int j = 0; j < 2; ++j) {
    for (int k = 0; k < 2; ++k) {
      //REQUIRE(g(1, j, k) == Approx(3.0));
      //REQUIRE(g(2, j, k) == Approx(2.0));
      //REQUIRE(g(3, j, k) == Approx(1.0));
    }
  }

}




