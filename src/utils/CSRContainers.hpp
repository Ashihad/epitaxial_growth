#pragma once

#include <cstdlib>
#include <memory>

struct CSRMatrix {
  CSRMatrix(const std::size_t nmax, const std::size_t row_count)
      : value{new double[nmax]},
        col_index{new std::size_t[nmax]},
        row_index{new std::size_t[row_count + 1]},
        n_rows{row_count} {}
  std::unique_ptr<double[]> value;
  std::unique_ptr<std::size_t[]> col_index;
  std::unique_ptr<std::size_t[]> row_index;
  std::size_t n_rows;
};

template <typename T>
struct FastVector {
  FastVector(const std::size_t vec_size) : size{vec_size}, data{new T[size]} {}

  std::size_t size;
  std::unique_ptr<T[]> data;

  T* get() { return data.get(); };
  T& operator[](const std::size_t index) { return data[index]; }
  const T& operator[](const std::size_t index) const { return data[index]; }
};