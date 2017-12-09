// Copyright (C) 2015-2016 Jason W. DeGraw
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
#ifndef UPWIND_HPP
#define UPWIND_HPP
#include "staggeredgrid.hpp"
#include <functional>
#include "red3api.hpp"
#include <Eigen/Sparse>

namespace red3 {
namespace upwind {

class RED3_API IncompressibleStaggeredSolver : public StaggeredGrid
{
public:
  enum class BoundaryCondition {  };
  enum class Differencing { FirstOrder, Hybrid, PowerLaw, Central };
  IncompressibleStaggeredSolver(unsigned ni, unsigned nj, unsigned nk = 1, Differencing diff = Differencing::FirstOrder,
    bool xperi = false) : StaggeredGrid(ni, nj, nk, xperi), differencing(diff)
  {
    switch(differencing) {
    case Differencing::FirstOrder:
      m_A = std::bind(&IncompressibleStaggeredSolver::firstOrderUpwind, this, std::placeholders::_1);
      break;
    case Differencing::Hybrid:
      m_A = std::bind(&IncompressibleStaggeredSolver::hybrid, this, std::placeholders::_1);
      break;
    case Differencing::Central:
      m_A = std::bind(&IncompressibleStaggeredSolver::centralDifference, this, std::placeholders::_1);
      break;
    case Differencing::PowerLaw:
      m_A = std::bind(&IncompressibleStaggeredSolver::powerLaw, this, std::placeholders::_1);
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

  const Differencing differencing;

private:
  std::function<double(double)> m_A;
  Eigen::SparseMatrix<double> m_Mu;
  Eigen::SparseMatrix<double> m_Mv;
  Eigen::SparseMatrix<double> m_Mw;

};

}
}

#endif