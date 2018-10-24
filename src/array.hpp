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

template <class T> class ChildArray
{
public:
  ChildArray() : m_parent(nullptr)
  {}
  ChildArray(T *parent) : m_parent(parent)
  {
    m_impl = std::shared_ptr<double>(new double[m_parent->ni*m_parent->nj*m_parent->nk], std::default_delete<double[]>());
  }
  ChildArray(const ChildArray &other) : m_impl(other.m_impl), m_parent(other.m_parent)
  {}

  ChildArray& operator=(const ChildArray &other)
  {
    m_impl = other.m_impl;
    m_parent = other.m_parent;
    return *this;
  }
  bool operator==(const ChildArray &other) const
  {
    return m_impl == other.m_impl;
  }
  bool operator!=(const ChildArray &other) const
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

template <class T> class ChildArrayU
{
public:
  ChildArrayU() : m_parent(nullptr)
  {}

  ChildArrayU(T *parent) : m_parent(parent)
  {
    m_impl = std::shared_ptr<double>(new double[m_parent->nu*m_parent->nj*m_parent->nk], std::default_delete<double[]>());
  }

  ChildArrayU(const ChildArrayU &other) : m_impl(other.m_impl), m_parent(other.m_parent)
  {}

  ChildArrayU& operator=(const ChildArrayU &other)
  {
    m_impl = other.m_impl;
    m_parent = other.m_parent;
    return *this;
  }

  bool operator==(const ChildArrayU &other) const
  {
    return m_impl == other.m_impl;
  }

  bool operator!=(const ChildArrayU &other) const
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

  inline double copy(ChildArrayU &other)
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

template <class T> class ChildArrayUVW
{
public:
  ChildArrayUVW() : m_parent(nullptr), m_impl(nullptr)
  {}
  ChildArrayUVW(T *parent) : m_parent(parent)
  {
    m_impl = std::shared_ptr<double>(new double[m_parent->nu*m_parent->nv*m_parent->nw], std::default_delete<double[]>());
  }
  ChildArrayUVW(const ChildArrayUVW &other) : m_impl(other.m_impl), m_parent(other.m_parent)
  {}

  ChildArrayUVW& operator=(const ChildArrayUVW &other)
  {
    m_impl = other.m_impl;
    m_parent = other.m_parent
    return *this;
  }

  bool operator==(const ChildArrayUVW &other) const
  {
    return m_impl == other.m_impl;
  }

  bool operator!=(const ChildArrayUVW &other) const
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

template <class T> class ChildArrayV
{
public:
  ChildArrayV() : m_parent(nullptr)
  {}
  ChildArrayV(T *parent) : m_parent(parent)
  {
    m_impl = std::shared_ptr<double>(new double[m_parent->ni*m_parent->nv*m_parent->nk], std::default_delete<double[]>());
  }
  ChildArrayV(const ChildArrayV &other) : m_impl(other.m_impl), m_parent(other.m_parent)
  {}

  ChildArrayV& operator=(const ChildArrayV &other)
  {
    m_impl = other.m_impl;
    m_parent = other.m_parent;
    return *this;
  }
  bool operator==(const ChildArrayV &other) const
  {
    return m_impl == other.m_impl;
  }
  bool operator!=(const ChildArrayV &other) const
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

template <class T> class ChildArrayW
{
public:
  ChildArrayW() : m_parent(nullptr)
  {}
  ChildArrayW(T *parent) : m_parent(parent)
  {
    m_impl = std::shared_ptr<double>(new double[m_parent->ni*m_parent->nj*m_parent->nw], std::default_delete<double[]>());
  }
  ChildArrayW(const ChildArrayW &other) : m_impl(other.m_impl), m_parent(other.m_parent)
  {}

  ChildArrayW& operator=(const ChildArrayW &other)
  {
    m_impl = other.m_impl;
    m_parent = other.m_parent;
    return *this;
  }
  bool operator==(const ChildArrayW &other) const
  {
    return m_impl == other.m_impl;
  }
  bool operator!=(const ChildArrayW &other) const
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
