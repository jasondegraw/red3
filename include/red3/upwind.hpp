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
#ifndef RED3_UPWIND_HPP
#define RED3_UPWIND_HPP
#include "staggeredgrid.hpp"
#include "array.hpp"
#include <functional>
#include "red3api.hpp"
#include <Eigen/Sparse>

namespace red3 {
namespace upwind {

enum class Differencing { FirstOrder, Hybrid, PowerLaw, Central };

class RED3_API StaggeredIncompressibleSteadyFlow : public StaggeredGrid
{
public:
  enum class BoundaryCondition {  };

  StaggeredIncompressibleSteadyFlow(double rho, double mu, Grid1D& x, Grid1D& y, Grid1D& z = Grid1D::one(),
    Differencing diff = Differencing::FirstOrder, bool xperi = false) : StaggeredGrid(x, y, z, xperi), differencing(diff),
    ae(this), aw(this), an(this), as(this), af(this), ab(this), b(this), rho(rho), mu(mu)
  {
    switch(differencing) {
    case Differencing::FirstOrder:
      A = std::bind(&StaggeredIncompressibleSteadyFlow::firstOrderUpwind, this, std::placeholders::_1);
      break;
    case Differencing::Hybrid:
      A = std::bind(&StaggeredIncompressibleSteadyFlow::hybrid, this, std::placeholders::_1);
      break;
    case Differencing::Central:
      A = std::bind(&StaggeredIncompressibleSteadyFlow::centralDifference, this, std::placeholders::_1);
      break;
    case Differencing::PowerLaw:
      A = std::bind(&StaggeredIncompressibleSteadyFlow::powerLaw, this, std::placeholders::_1);
      break;
    }

    m_Mu = Eigen::SparseMatrix<double>(nu, nu);
    m_Mv = Eigen::SparseMatrix<double>(nv, nv);
    if(nk > 1) {
      m_Mw = Eigen::SparseMatrix<double>(nw, nw);
    }

  }

  double firstOrderUpwind(double) const
  {
    return 1.0;
  }

  double centralDifference(double absP) const
  {
    return 1.0 - 0.5*absP;
  }

  double hybrid(double absP) const
  {
    return std::max(0.0, 1.0 - 0.5*absP);
  }

  double powerLaw(double absP) const
  {
    return std::max(0.0, std::pow(1.0 - 0.5*absP, 5));
  }

  void setupU();
  void setupV();
  void setupW();

  const Differencing differencing;

  ArrayU ae, aw, an, as, af, ab, b;
  const double rho;
  const double mu;

  std::function<double(double)> A;

private:
  //std::function<double(double)> m_A;
  Eigen::SparseMatrix<double> m_Mu;
  Eigen::SparseMatrix<double> m_Mv;
  Eigen::SparseMatrix<double> m_Mw;

};

}
}

#endif