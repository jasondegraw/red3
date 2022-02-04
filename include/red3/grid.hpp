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
#ifndef RED3_GRID_HPP
#define RED3_GRID_HPP

#include "red3/defs.hpp"
#include <vector>
#include <memory>
#include <string>
#include <optional>
#include <algorithm>

namespace red3 {

class Generator1D
{
public:

  virtual bool successful() const
  {
    return true;
  }

  virtual bool failed() const
  {
    return false;
  }

  virtual index_t iterations() const
  {
    return 0;
  }

  virtual bool uniform() const
  {
    return true;
  }

  virtual double delta(double i) const
  {
    return 1.0;
  }

  virtual std::vector<double> grid()
  {
    return { -0.5, 0.5 };
  }

  virtual std::vector<double> midgrid()
  {
    return { 0.0 };
  }

  virtual std::vector<double> delta()
  {
    return { 1.0 };
  }

  //double operator[](index_t i) const
  //{
  //  return m_x[i];
  //}

  std::vector<std::string> messages;

  //const index_t n;

};

class Uniform1D : public Generator1D
{
public:
  Uniform1D(double delta, index_t n) : m_n(std::max(n,(index_t)1)), m_delta(delta <=0 ? 1.0 : delta), m_L(m_n*m_delta)
  {}

  Uniform1D(size_t n, double L=1.0) : m_n(std::max(n,(index_t)1)), m_delta(std::abs(L)/(double)m_n), m_L(std::abs(L))
  {}

  virtual double delta(double) const
  {
    return m_delta;
  }

  std::vector<double> grid()
  {
    std::vector<double> result(m_n + 1);
    result[0] = 0.0;
    for (index_t i = 1; i < m_n; ++i) {
      result[i] = result[i - 1] + m_delta;
    }
    result[m_n] = m_L;
    return result;
  }

  std::vector<double> midgrid()
  {
    std::vector<double> result(m_n);
    result[0] = 0.5 * m_delta;
    for (index_t i = 1; i < m_n; ++i) {
      result[i] = result[i - 1] + m_delta;
    }
    return result;
  }

  std::vector<double> delta()
  {
    return std::vector<double>(m_n, m_delta);
  }

private:
  index_t m_n;
  double m_delta;
  double m_L;
};

class OneSidedVinokur : public Generator1D
{
public:
  OneSidedVinokur(double delta, double L, size_t i, size_t I)// : Generator1D(I), m_delta(delta)
  {
  }

  virtual bool successful() const
  {
    return true;
  }

  virtual bool failed() const
  {
    return false;
  }

  virtual double delta(double)
  {
    return m_delta;
  }

  std::vector<double> grid()
  {
    std::vector<double> result;
    return result;
  }

  std::vector<double> midgrid()
  {
    std::vector<double> result;
    return result;
  }

  std::vector<double> delta()
  {
    return std::vector<double>();
  }

private:
  std::optional<double> solve(double delta, double L, size_t i, size_t I, double tolerance, size_t max_iterations);

  double m_delta;

};

class RED3_API Grid1D
{
public:
  static std::optional<Grid1D> generate(Generator1D& generator);
 
  Grid1D(index_t n, double L = 1.0);
  Grid1D(int n);
 
  static Grid1D one();

  std::vector<double> grid() const
  {
    return m_x;
  }

  std::vector<double> midgrid() const
  {
    return m_xm;
  }

  std::vector<double> delta() const
  {
    return m_dx;
  }

  double& operator[](index_t i)
  {
    return m_x[i];
  }

  double& grid(index_t i)
  {
    return m_x[i];
  }

  double& midgrid(index_t i)
  {
    return m_x[i];
  }

  double& m(index_t i)
  {
    return m_x[i];
  }

  double& delta(index_t i)
  {
    return m_dx[i];
  }

  double& delta0()
  {
    return m_dx[0];
  }

  double& deltaN()
  {
    return m_dx.back();
  }

  double operator[](index_t i) const
  {
    return m_x[i];
  }

  double grid(index_t i) const
  {
    return m_x[i];
  }

  double midgrid(index_t i) const
  {
    return m_x[i];
  }

  double m(index_t i) const
  {
    return m_xm[i];
  }

  double delta(index_t i) const
  {
    return m_dx[i];
  }

  double delta0() const
  {
    return m_dx[0];
  }

  double deltaN() const
  {
    return m_dx.back();
  }

  index_t size() const
  {
    return m_x.size();
  }

  //const index_t N;
  const bool uniform;

private:
  Grid1D(std::vector<double> x, std::vector<double> xm, std::vector<double> dx, bool uniform) : uniform(uniform),
    m_x(x), m_xm(xm), m_dx(dx)
  {}

  std::vector<double> m_x;
  std::vector<double> m_xm;
  std::vector<double> m_dx;
};

}

#endif
