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
#ifndef RED3_UTILITIES_HPP
#define RED3_UTILITIES_HPP

#include "red3/defs.hpp"
#include <string>

namespace red3 {

double power2(double v);

template <int N> double power(double v)
{
  static_assert(false);
  //return v*power<N - 1>(v);
}

template <typename V, typename M, typename I> void tridai(V &a, V &b, V &c, M &r, I i1, I i2, I j1, I j2, I k1, I k2, I ni, I nj, I nk)
{
  //
  //  *******************************************************************
  //
  //     This is the new i-index array tridiagonal solver, which
  //     overwrites the LHS every time, and assumes that the arrays are
  //     aligned properly.
  //
  //     Rewritten by Jason W.DeGraw(06-20-2005)
  //
  //  *******************************************************************
  //
  //....Forward elimination
  //
  for (I i = i1; i < i2; ++i) {
    I ip1 = i + 1;
    c[i] = c[i] / b[i];
    b[ip1] = b[ip1] - a[ip1] * c[i];
    for (I k = k1; k <= k2; ++k) {
      for (I j = j1; j <= j2; ++j) {
        r(i, j, k) = r(i, j, k) / b[i];
        r(ip1, j, k) = r(ip1, j, k) - a[ip1] * r(i, j, k);
      }
    }
  }
  //
  //....Back substitution
  //
  for (I k = k1; k <= k2; ++k) {
    for (I j = j1; j <= j2; ++j) {
      r(i2, j, k) = r(i2, j, k) / b[i2];
    }
  }
  for (I i = i2 - 1; i >= i1; --i) {
    I ip1 = i + 1;
    for (I k = k1; k <= k2; ++k) {
      for (I j = j1; j <= j2; ++j) {
        r(i, j, k) = r(i, j, k) - c[i] * r(ip1, j, k);
      }
    }
  }
}

template <typename V, typename M, typename I> void tridaj(V &a, V &b, V &c, M &r, I i1, I i2, I j1, I j2, I k1, I k2, I ni, I nj, I nk)
{
  //
  //  *******************************************************************
  //
  //     called from invert
  //
  //     This is the new j - index array tridiagonal solver, which
  //     overwrites the LHS every time, and assumes that the arrays are
  //     aligned properly.
  //
  //     Rewritten by Jason W.DeGraw(06 - 20 - 2005)
  //
  //  *******************************************************************
  //
  //....Forward elimination
  //
  for (I j = j1; j < j2; ++j) {
    I jp1 = j + 1;
    c[j] = c[j] / b[j];
    b[jp1] = b[jp1] - a[jp1] * c[j];
    for (I k = k1; k <= k2; ++k) {
      for (I i = i1; i <= i2; ++i) {
        r(i, j, k) = r(i, j, k) / b[j];
        r(i, jp1, k) = r(i, jp1, k) - a[jp1]*r(i, j, k);
      }
    }
  }
  //
  //....Back substitution
  //
  for (I k = k1; k <= k2; ++k) {
    for (I i = i1; i <= i2; ++i) {
      r(i, j2, k) = r(i, j2, k) / b[j2];
    }
  }
  for (I j = j2 - 1; j >= j1; --j) {
    I jp1 = j + 1;
    for (I k = k1; k <= k2; ++k) {
      for (I i = i1; i <= i2; ++i) {
        r(i, j, k) = r(i, j, k) - c[j] * r(i, jp1, k);
      }
    }
  }
}

template <typename V, typename M, typename I> void tripak(V &a, V &b, V &c, M &f, V &q, I i1, I i2, I j1, I j2, I k1, I k2, I ni, I nj, I nk)
{
  //
  //  *******************************************************************
  //
  //     called from invert
  //
  //     New version using Sherman-Morrison Algorithm(explicitly).
  //
  //     Solve(A + u*vT)x = f
  //
  //     where
  //
  //     uT = [g 0 0 ... 0 c(j2)] and vT = [1 0 0 ... 0 a(j1)/g]
  //
  //     g is arbitrary and taken here to be -b(j1). A is nearly
  //     the tridiagonal system one would expect.
  //
  //     Then solve Ay = f and Aq = u, and compute x as
  //
  //     x = y - (v dot y) / (1 - (v dot q)) q
  //
  //     No v storage is needed, and f and q are solved in place.
  //
  //  *******************************************************************
  //
  //....Save the cyclic parts needed later
  //
  double g = -b[k1];
  //
  //....Modify the LHS
  //
  b[k2] = b[k2] + c[k2] * a[k1] / b[k1];
  b[k1] = 2.0 * b[k1];
  //
  //....Solve the two systems
  //
  //
  //....Forward elimination - do some(useless in this case) silly stuff
  //    to save ops
  //
  q[k1] = g;
  q[k2] = c[k2];
  I k = k1;
  I kp1 = k + 1;
  c[k] = c[k] / b[k];
  q[k] = -0.5;
  b[kp1] = b[kp1] - a[kp1] * c[k];
  for (I j = j1; j <= j2; ++j) {
    for (I i = i1; i <= i2; ++i) {
      f(i, j, k) = f(i, j, k) / b[k];
      f(i, j, kp1) = f(i, j, kp1) - a[kp1] * f(i, j, k);
    }
  }
  q[kp1] = -a[kp1] * q[k];
  for (k = k1 + 1; k <= k2 - 2; ++k) {
    kp1 = k + 1;
    c[k] = c[k] / b[k];
    q[k] = q[k] / b[k];
    b[kp1] = b[kp1] - a[kp1] * c[k];
    for (I j = j1; j <= j2; ++j) {
      for (I i = i1; i <= i2; ++i) {
        f(i, j, k) = f(i, j, k) / b[k];
        f(i, j, kp1) = f(i, j, kp1) - a[kp1] * f(i, j, k);
      }
    }
    q[kp1] = -a[kp1] * q[k];
  }
  k = k2 - 1;
  kp1 = k + 1;
  c[k] = c[k] / b[k];
  q[k] = q[k] / b[k];
  b[kp1] = b[kp1] - a[kp1] * c[k];
  for (I j = j1; j <= j2; ++j) {
    for (I i = i1; i <= i2; ++i) {
      f(i, j, k) = f(i, j, k) / b[k];
      f(i, j, kp1) = f(i, j, kp1) - a[kp1] * f(i, j, k);
    }
  }
  q[kp1] = q[kp1] - a[kp1] * q[k];
  //
  //....Backward substitution
  //
  for (I j = j1; j <= j2; ++j) {
    for (I i = i1; i <= i2; ++i) {
      f(i, j, k2) = f(i, j, k2) / b[k2];
    }
  }
  q[k2] = q[k2] / b[k2];
  for (k = k2 - 1; k >= k1; --k) {
    kp1 = k + 1;
    for (I j = j1; j <= j2; ++j) {
      for (I i = i1; i <= i2; ++i) {
        f(i, j, k) = f(i, j, k) - c[k] * f(i, j, kp1);
      }
    }
    q[k] = q[k] - c[k] * q[kp1];
  }
  //
  //....Combine solution vectors
  //
  double fact = 1.0 + q[k1] + a[k1] * q[k2] / g;
  for (I j = j1; j <= j2; ++j) {
    for (I i = i1; i <= i2; ++i) {
      double mult = (f(i, j, k1) + a[k1] * f(i, j, k2) / g) / fact;
      for (k = k1; k <= k2; ++k) {
        f(i, j, k) = f(i, j, k) - mult * q(k);
      }
    }
  }
}

void RED3_API fatal(const std::string& mesg);
RED3_API void* callocate(size_t num, size_t size, const std::string& name);

}

#endif
