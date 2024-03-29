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
#ifndef RED3_ARRAY_HPP
#define RED3_ARRAY_HPP
#include <memory>
#include "red3/defs.hpp"
#include "red3/arrayops.hpp"

namespace red3 {

template <typename T, typename U> class BaseArray
{
public:
  BaseArray() : m_parent(nullptr)
  {}

  BaseArray(T* parent, index_t N) : m_parent(parent)
  {
    m_impl = std::shared_ptr<double>(new double[N], std::default_delete<double[]>());
    memset(m_impl.get(), 0, N * sizeof(double));
  }

  BaseArray(T* parent, std::shared_ptr<double> impl) : m_impl(impl), m_parent(parent)
  {}

  U& operator=(const U &other)
  {
    m_impl = other.m_impl;
    m_parent = other.m_parent;
    //return *this;
    // approxfac is triggering the need for the cast on ChildArrayW, which is strange. Why not on u or v?
    return (U&)(*this);
  }

  bool operator==(const U &other) const
  {
    return m_impl == other.m_impl;
  }

  bool operator!=(const U &other) const
  {
    return m_impl != other.m_impl;
  }

  double &operator[](index_t i)
  {
    return (m_impl.get())[i];
  }

  double &operator()(index_t i)
  {
    return (m_impl.get())[i];
  }

  virtual double &operator()(index_t i, index_t j, index_t k) = 0;

  void clear()
  {
    memset(m_impl.get(), 0, size()*sizeof(double));
  }

  virtual index_t size() const = 0;

  void swap(U &other)
  {
    if (m_parent == other.m_parent) {
      m_impl.swap(other.m_impl);
    }
  }

  T* parent() const
  {
    return m_parent;
  }

protected:
  std::shared_ptr<double> m_impl;
  T *m_parent;
};


template <class T> class ChildArray : public BaseArray<T, ChildArray<T>>
{
public:
  ChildArray() : BaseArray<T, ChildArray<T>>()
  {}

  ChildArray(T *parent) : BaseArray<T, ChildArray<T>>(parent, parent->ni*parent->nj*parent->nk)
  {}

  ChildArray(const ChildArray &other) : BaseArray<T, ChildArray<T>>(other.m_parent, other.m_impl)
  {}

  double &operator()(index_t i, index_t j, index_t k)
  {
    return (this->m_impl.get())[INDEX(i, j, k, this->m_parent->ni, this->m_parent->nj, this->m_parent->nk)];
  }

  index_t size() const
  {
    return this->m_parent->ni*this->m_parent->nj*this->m_parent->nk;
  }

};

template <class T> class ChildArrayU : public BaseArray<T, ChildArrayU<T>>
{
public:
  ChildArrayU() : BaseArray<T, ChildArrayU<T>>()
  {}

  ChildArrayU(T *parent) : BaseArray<T, ChildArrayU<T>>(parent, parent->nu*parent->nj*parent->nk)
  {}

  ChildArrayU(const ChildArrayU &other) : BaseArray<T, ChildArrayU<T>>(other.m_parent, other.m_impl)
  {}

  //ChildArrayU& operator=(const ChildArrayU &other)
  //{
  //  this->m_impl = other.m_impl;
  //  this->m_parent = other.m_parent;
  //  return *this;
  //}

  inline double &operator()(index_t i, index_t j, index_t k)
  {
    return (this->m_impl.get())[UINDEX(i, j, k, this->m_parent->nu, this->m_parent->nj, this->m_parent->nk)];
  }

  //inline double copy(ChildArrayU &other)
  //{
  //  int nijk = m_parent->nu * m_parent->nj * m_parent->nk;
  //}

  index_t size() const
  {
    return this->m_parent->nu * this->m_parent->nj * this->m_parent->nk;
  }

};

template <class T> class ChildArrayUVW : public BaseArray<T, ChildArrayUVW<T>>
{
public:
  ChildArrayUVW() : BaseArray<T, ChildArrayUVW<T>>()
  {}

  ChildArrayUVW(T *parent) : BaseArray<T, ChildArrayUVW<T>>(parent, parent->nu*parent->nv*parent->nw)
  {}

  ChildArrayUVW(const ChildArrayUVW &other) : BaseArray<T, ChildArrayUVW<T>>(other.m_parent, other.m_impl)
  {}

  double &operator()(index_t i, index_t j, index_t k)
  {
    return (this->m_impl.get())[INDEX(i, j, k, this->m_parent->nu, this->m_parent->nv, this->m_parent->nw)];
  }

  index_t size() const
  {
    return this->m_parent->nu * this->m_parent->nv * this->m_parent->nw;
  }

};

template <class T> class ChildArrayV : public BaseArray<T, ChildArrayV<T>>
{
public:
  ChildArrayV()
  {}

  ChildArrayV(T *parent) : BaseArray<T, ChildArrayV<T>>(parent, parent->ni*parent->nv*parent->nk)
  {}

  ChildArrayV(const ChildArrayV &other) : BaseArray<T, ChildArrayV<T>>(other.m_parent, other.m_impl)
  {}

  double &operator()(index_t i, index_t j, index_t k)
  {
    return (this->m_impl.get())[VINDEX(i, j, k, this->m_parent->ni, this->m_parent->nv, this->m_parent->nk)];
  }

  index_t size() const
  {
    return this->m_parent->ni * this->m_parent->nv * this->m_parent->nk;
  }

};

template <class T> class ChildArrayW : public BaseArray<T, ChildArrayW<T>>
{
public:
  ChildArrayW()
  {}

  ChildArrayW(T *parent) : BaseArray<T, ChildArrayW<T>>(parent, parent->ni*parent->nj*parent->nw)
  {}

  ChildArrayW(const ChildArrayW &other) : BaseArray<T, ChildArrayW<T>>(other.m_parent, other.m_impl)
  {}

  double &operator()(index_t i, index_t j, index_t k)
  {
    return (this->m_impl.get())[WINDEX(i, j, k, this->m_parent->ni, this->m_parent->nj, this->m_parent->nw)];
  }

  index_t size() const
  {
    return this->m_parent->ni * this->m_parent->nj * this->m_parent->nw;
  }

};

}

#endif
