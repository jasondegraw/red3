// Copyright (C) 2015-2017 Jason W. DeGraw
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

template <typename I, typename J, typename K> class ArrayX
{
public:

  bool operator==(const ArrayX<I, J, K> &other) const
  {
    return m_impl == other.m_impl;
  }

  bool operator!=(const ArrayX<I, J, K> &other) const
  {
    return m_impl != other.m_impl;
  }

  inline double &operator[](unsigned i)
  {
    return (m_impl.get())[i];
  }

  inline double &operator()(unsigned i, unsigned j, unsigned k)
  {
    return (m_impl.get())[INDEX(i, j, k, *I, *J, *K)];
  }

  int size()
  {
    return (*I)*(*J)*(*K);
  }

  void copy(ArrayX<I, J, K> &other)
  {
    for (int ijk = 0; i < (*I)*(*J)*(*K); ++ijk) {
      (m_impl.get())[ijk] = (other.m_impl.get())[ijk];
    }
  }

protected:
  ArrayX()
  {}

  void allocate()
  {
    m_impl = std::shared_ptr<double>(new double[(*I)*(*J)*(*K)], std::default_delete<double[]>());
  }

  std::shared_ptr<double> m_impl;

  friend class StaggeredGrid;

};

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
    m_impl = std::shared_ptr<double>(new double[m_parent->ni*m_parent->nj*m_parent->nk], std::default_delete<double[]>());
  }
  Array(const Array &other) : m_impl(other.m_impl), m_parent(other.m_parent)
  {}

  Array& operator=(const Array &other)
  {
    m_impl = other.m_impl;
    m_parent = other.m_parent;
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

  inline double &operator[](unsigned i)
  {
    return (m_impl.get())[i];
  }

  inline double &operator()(unsigned i)
  {
    return (m_impl.get())[i];
  }

  inline double &operator()(unsigned i, unsigned j, unsigned k)
  {
    return (m_impl.get())[INDEX(i, j, k, m_parent->ni, m_parent->nj, m_parent->nk)];
  }

  int size()
  {
    return m_parent->ni*m_parent->nj*m_parent->nk;
  }

  T* parent()
  {
    return m_parent;
  }

protected:
  std::shared_ptr<double> m_impl;
  T *m_parent;

};

template <class T> class ArrayU
{
public:
  ArrayU() : m_parent(nullptr)
  {}

  ArrayU(T *parent) : m_parent(parent)
  {
    m_impl = std::shared_ptr<double>(new double[m_parent->nu*m_parent->nj*m_parent->nk], std::default_delete<double[]>());
  }

  ArrayU(const ArrayU &other) : m_impl(other.m_impl), m_parent(other.m_parent)
  {}

  ArrayU& operator=(const ArrayU &other)
  {
    m_impl = other.m_impl;
    m_parent = other.m_parent;
    return *this;
  }

  bool operator==(const ArrayU &other) const
  {
    return m_impl == other.m_impl;
  }

  bool operator!=(const ArrayU &other) const
  {
    return m_impl != other.m_impl;
  }

  inline double &operator[](unsigned i)
  {
    return (m_impl.get())[i];
  }

  inline double &operator()(unsigned i)
  {
    return (m_impl.get())[i];
  }

  inline double &operator()(unsigned i, unsigned j, unsigned k)
  {
    return (m_impl.get())[UINDEX(i, j, k, m_parent->nu, m_parent->nj, m_parent->nk)];
  }

  inline double copy(ArrayU &other)
  {
    int nijk = m_parent->nu * m_parent->nj * m_parent->nk;

  }

  int size()
  {
    return m_parent->nu*m_parent->nj*m_parent->nk;
  }

  T* parent()
  {
    return m_parent;
  }

protected:
  std::shared_ptr<double> m_impl;
  T *m_parent;

};

template <class T> class ArrayUVW
{
public:
  ArrayUVW() : m_parent(nullptr), m_impl(nullptr)
  {}
  ArrayUVW(T *parent) : m_parent(parent)
  {
    m_impl = std::shared_ptr<double>(new double[m_parent->nu*m_parent->nv*m_parent->nw], std::default_delete<double[]>());
  }
  ArrayUVW(const ArrayUVW &other) : m_impl(other.m_impl), m_parent(other.m_parent)
  {}

  ArrayUVW& operator=(const ArrayUVW &other)
  {
    m_impl = other.m_impl;
    m_parent = other.m_parent
    return *this;
  }

  bool operator==(const ArrayUVW &other) const
  {
    return m_impl == other.m_impl;
  }

  bool operator!=(const ArrayUVW &other) const
  {
    return m_impl != other.m_impl;
  }

  inline double &operator[](unsigned i)
  {
    return (m_impl.get())[i];
  }

  inline double &operator()(unsigned i)
  {
    return (m_impl.get())[i];
  }

  inline double &operator()(unsigned i, unsigned j, unsigned k)
  {
    return (m_impl.get())[INDEX(i, j, k, m_parent->nu, m_parent->nv, m_parent->nw)];
  }

  int size()
  {
    return m_parent->nu*m_parent->nv*m_parent->nw;
  }

  T* parent()
  {
    return m_parent;
  }

protected:
  std::shared_ptr<double> m_impl;
  T *m_parent;

};

template <class T> class ArrayV
{
public:
  ArrayV() : m_parent(nullptr)
  {}
  ArrayV(T *parent) : m_parent(parent)
  {
    m_impl = std::shared_ptr<double>(new double[m_parent->ni*m_parent->nv*m_parent->nk], std::default_delete<double[]>());
  }
  ArrayV(const ArrayV &other) : m_impl(other.m_impl), m_parent(other.m_parent)
  {}

  ArrayV& operator=(const ArrayV &other)
  {
    m_impl = other.m_impl;
    m_parent = other.m_parent;
    return *this;
  }
  bool operator==(const ArrayV &other) const
  {
    return m_impl == other.m_impl;
  }
  bool operator!=(const ArrayV &other) const
  {
    return m_impl != other.m_impl;
  }

  inline double &operator[](unsigned i)
  {
    return (m_impl.get())[i];
  }

  inline double &operator()(unsigned i)
  {
    return (m_impl.get())[i];
  }

  inline double &operator()(unsigned i, unsigned j, unsigned k)
  {
    return (m_impl.get())[VINDEX(i, j, k, m_parent->ni, m_parent->nv, m_parent->nk)];
  }

protected:
  std::shared_ptr<double> m_impl;
  T *m_parent;

};

template <class T> class ArrayW
{
public:
  ArrayW() : m_parent(nullptr)
  {}
  ArrayW(T *parent) : m_parent(parent)
  {
    m_impl = std::shared_ptr<double>(new double[m_parent->ni*m_parent->nj*m_parent->nw], std::default_delete<double[]>());
  }
  ArrayW(const ArrayW &other) : m_impl(other.m_impl), m_parent(other.m_parent)
  {}

  ArrayW& operator=(const ArrayW &other)
  {
    m_impl = other.m_impl;
    m_parent = other.m_parent;
    return *this;
  }
  bool operator==(const ArrayW &other) const
  {
    return m_impl == other.m_impl;
  }
  bool operator!=(const ArrayW &other) const
  {
    return m_impl != other.m_impl;
  }

  inline double &operator[](unsigned i)
  {
    return (m_impl.get())[i];
  }

  inline double &operator()(unsigned i)
  {
    return (m_impl.get())[i];
  }

  inline double &operator()(unsigned i, unsigned j, unsigned k)
  {
    return (m_impl.get())[WINDEX(i, j, k, m_parent->ni, m_parent->nj, m_parent->nw)];
  }

protected:
  std::shared_ptr<double> m_impl;
  T *m_parent;

};

}

#endif
