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
#include"grid.hpp"
#include <algorithm>
#include "fmt/printf.h"

namespace red3 {

double vkruh(double ds, double i, double I)
{
  return 1.0 + tanh(ds * (i / I - 1.0)) / tanh(ds);
}

/*****************************************************************************/
/*                                                                           */
/*  hfindfactor - Given a spacing, determine the required stretching factor  */
/*                                                                           */
/*  The one sided stretching function between 0 and 1 is                     */
/*                                                                           */
/*          ds = 1 + tanh(factor * (i / I - 1)) / tanh(factor)               */
/*                                                                           */
/*  We apply Newton's Method to find the stretching factor given ds, i, and  */
/*  I. The one stretching function can be reorganized as                     */
/*                                                                           */
/*           (ds - 1)tanh(factor) = tanh(factor * (i / I - 1))               */
/*              C2 * tanh(factor) = tanh(C1 * factor)                        */
/*                                                                           */
/*  This only has solutions if |C2| > |C1|.  To see this, note that the      */
/*  limit of the LHS is C2, and the limit of the RHS is -1.  Both C1 and C2  */
/*  are negative.  The slopes at zero are C2 and C1. respectively, so the    */
/*  only possible way that the two curves can intersect is if |C2| > |C1|.   */
/*  Approximating tanh(C1*delta) as the line C1*delta for small delta, and   */
/*  approximating  C2*tanh(delta) as C2 for large delta, we obtain an first  */
/*  guess for the answer as                                                  */
/*                                                                           */
/*                             C2 = C1*factor                                */
/*                                                                           */
/*  Newton's method can now be used to determine the solution.               */
/*                                                                           */
/*****************************************************************************/
std::optional<double> OneSidedVinokur::solve(double delta, double L, size_t i, size_t I, double tolerance, size_t max_iterations)
//int vkruhff(double ds, double I, double* d, double tol, int itermax,
//  FILE* debug, char* prefix)
{
  double ds{ delta / L }; // Rescale

  //int i = 0;
  double f, fp, factor, s1, s2;
  double C1{ 1.0 / I - 1.0 };
  double C2 = ds - 1.0;

  messages.push_back("Vinokur h Stretching Factor Solution ---------------+");
  messages.push_back(fmt::sprintf("  ds = % .8e                              |", ds));
  messages.push_back(fmt::sprintf("   I = % .4e                                  |", I));
  messages.push_back(fmt::sprintf(" tol = % .3e, itermax = %5d                  |", tolerance, max_iterations));
  messages.push_back("----------------------------------------------------+");

  // Solve only for positive ds
  if (ds <= 0.0) {
    messages.push_back("No solution for negative ds.                        |");
    messages.push_back("----------------------------------------------------+");
    return {};
  }

  // Only have a nontrivial solution for |C2| > |C1|
  if (C2 > 0.0 || C2 > C1) {
    messages.push_back("No solution for given inputs.                       |");
    messages.push_back("----------------------------------------------------+");
    return {};
  }
  /*
    As an initial guess, take the intersection of the lines

      y = 1.0   (limit of tanh(d*C1)
      y = C2*d  (approximation of C2*tanh(d))
      y = C2    (limit of C2*tanh(d))
      y = C1*d  (approximation of tanh(C1*d))
  */
  size_t iter = 1;
  factor = C1 / C2;
  f = tanh(factor * C1) - C2 * tanh(factor);
  messages.push_back(" iter            d                      f           |");
  messages.push_back("----- ---------------------- ---------------------- |");
  messages.push_back(fmt::sprintf("%5d % .15e % .15e |\n", iter, delta, f));

  /* printf("%4d % .15e % .15e\n",i,delta,f); */
  while (fabs(f) > tolerance && i <= max_iterations) {
    s1 = 1.0 / cosh(delta * C1);
    s2 = 1.0 / cosh(delta);
    fp = C1 * s1 * s1 - C2 * s2 * s2;
    delta -= f / fp;
    f = tanh(delta * C1) - C2 * tanh(delta);
    i++;
    messages.push_back(fmt::sprintf("%5d % .15e % .15e |", iter, delta, f));
  }
  messages.push_back("----------------------------------------------------+");
  if (fabs(f) > tolerance) return {};
  return factor;
}


std::optional<Grid1D> Grid1D::generate(Generator1D& generator) 
{
  if (generator.failed()) {
    return {};
  }
  return Grid1D(generator.grid(), generator.midgrid(), generator.delta(), generator.uniform());
}

Grid1D::Grid1D(index_t n, double L) : uniform(true)
{
  Uniform1D generator(n, L);
  m_x = generator.grid();
  m_xm = generator.midgrid();
  m_dx = generator.delta();
}

Grid1D::Grid1D(int n) : Grid1D((index_t)std::max(1,n))
{}

Grid1D Grid1D::one()
{
  return Grid1D({ -0.5, 0.5 }, { 0.0 }, { 1.0 }, true);
}

}