// Copyright (C) 2015-2018 Jason W. DeGraw
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
#ifndef APPROXFRAC_HPP
#define APPROXFRAC_HPP
#include "staggeredgrid.hpp"
#include "array.hpp"
#include <functional>
#include "red3api.hpp"
#include <Eigen/Sparse>

namespace red3 {
namespace approxfrac {

class RED3_API IsothermalFlowSolver : public BasicStaggeredGrid<IsothermalFlowSolver>
{
public:
  enum class BoundaryCondition {  };
  enum class Differencing {  };
  IsothermalFlowSolver(double reynum, double dt, unsigned ni, unsigned nj, unsigned nk = 1, bool xperi = false)
    : BasicStaggeredGrid<IsothermalFlowSolver>(ni, nj, nk, xperi),
    g(this), //aw(this), an(this), as(this), af(this), ab(this), b(this), 
    reynum(reynum), dt(dt)
  {
   

    m_Mu = Eigen::SparseMatrix<double>(nu, nu);
    m_Mv = Eigen::SparseMatrix<double>(nv, nv);
    if(nk > 1) {
      m_Mw = Eigen::SparseMatrix<double>(nw, nw);
    }

  }


  void setupU();
  void setupV();
  void setupW();

  Array<IsothermalFlowSolver> g;
  const double reynum;
  const double dt;

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