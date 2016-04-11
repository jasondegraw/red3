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
#ifndef ARRAY_HPP
#define ARRAY_HPP
#include <memory>
//#include "staggeredgrid.hpp"
#include "red3api.hpp"
#include "arrayops.hpp"

namespace red3 {
/*
class RED3_API Array
{
public:
  //Array() = delete;
  //Array(StaggeredGrid *parent);
  Array(const Array &other);
  virtual ~Array(){}

  Array& operator=(const Array &other);
  bool operator==(const Array &other) const;
  bool operator!=(const Array &other) const;

  double &operator[](unsigned i)
  {
    return (m_impl.get())[i];
  }

  virtual double &operator()(unsigned i, unsigned j, unsigned k)
  {
    return (m_impl.get())[INDEX(i, j, k, m_parent->nx, m_parent->ny, m_parent->nz)];
  }

protected:
  std::shared_ptr<double> m_impl;
  //StaggeredGrid *m_parent;
  void *m_parent;

};
*/

template <class T> class Array
{
public:
  Array() : m_parent(nullptr)
  {}
  Array(T *parent) : m_parent(parent)
  {
    m_impl = std::shared_ptr<double>(new double[m_parent->nx*m_parent->ny*m_parent->nz], std::default_delete<double[]>());
  }
  Array(const Array &other) : m_impl(other.m_impl), m_parent(other.m_parent)
  {}
  virtual ~Array(){}

  Array& operator=(const Array &other)
  {
    m_impl = other.m_impl;
    return *this;
  }
  bool operator==(const Array &other) const
  {
    return m_impl == other.m_impl;
  }
  bool operator!=(const Array &other) const
  {
    return m_impl != other.m_impl;
  }

  double &operator[](unsigned i)
  {
    return (m_impl.get())[i];
  }

  inline double &operator()(unsigned i, unsigned j, unsigned k)
  {
    return (m_impl.get())[INDEX(i, j, k, m_parent->nx, m_parent->ny, m_parent->nz)];
  }

protected:
  std::shared_ptr<double> m_impl;
  T *m_parent;

};

}

#endif
