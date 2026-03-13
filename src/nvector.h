#ifndef O4_NVECTOR_H
#define O4_NVECTOR_H
#include <array>
#include <iostream>
#include <vector>

// Simple multi-dimensionsal array class with contiguous
template <typename T, std::size_t rank> class nvector {

private:
  size_t count;
  std::array<size_t, rank> sizes;

  void init() {
    sizes[rank - 1] = 1;
    for (size_t i = rank - 1; i > 0; i--) {
      sizes[i - 1] = N[i] * sizes[i];
    }
    count = N[0];
    for (size_t i = 1; i < rank; i++) {
      count *= N[i];
    }
    v.resize(count);
  }

public:
  std::array<size_t, rank> N;
  std::vector<T> v;

  size_t size() { return count; }
  T &operator[](std::size_t i) { return v[i]; }

  nvector() {
    for (size_t i = 0; i < rank; i++) {
      N[i] = 0;
      init();
    }
  }

  nvector(const std::array<size_t, rank> &dims) : N(dims) { init(); }

  // Specialiizations
  nvector(const size_t &i) {
    N[0] = i;
    init();
  }
  nvector(const size_t &i0, const size_t &i1) {
    N[0] = i0;
    N[1] = i1;
    init();
  }
  nvector(const size_t &i0, const size_t &i1, const size_t &i2) {
    N[0] = i0;
    N[1] = i1;
    N[2] = i2;
    init();
  }
  nvector(const size_t &i0, const size_t &i1, const size_t &i2,
          const size_t i3) {
    N[0] = i0;
    N[1] = i1;
    N[2] = i2;
    N[3] = i3;
    init();
  }
  nvector(const size_t &i0, const size_t &i1, const size_t &i2, const size_t i3,
          const size_t i4) {
    N[0] = i0;
    N[1] = i1;
    N[2] = i2;
    N[3] = i3;
    N[4] = i4;
    init();
  }

  void resize(const std::array<size_t, rank> &dims) {
    for (size_t i = 0; i < rank; i++) {
      N[i] = dims[i];
    }
    init();
  }
  // Specializations
  void resize(const size_t &i) {
    N[0] = i;
    init();
  }
  void resize(const size_t &i0, const size_t &i1) {
    N[0] = i0;
    N[1] = i1;
    init();
  }
  void resize(const size_t &i0, const size_t &i1, const size_t &i2) {
    N[0] = i0;
    N[1] = i1;
    N[2] = i2;
    init();
  }
  void resize(const size_t &i0, const size_t &i1, const size_t &i2,
              const size_t i3) {
    N[0] = i0;
    N[1] = i1;
    N[2] = i2;
    N[3] = i3;
    init();
  }
  void resize(const size_t &i0, const size_t &i1, const size_t &i2,
              const size_t i3, const size_t i4) {
    N[0] = i0;
    N[1] = i1;
    N[2] = i2;
    N[3] = i3;
    N[4] = i4;
    init();
  }
  // Element access
  size_t offset(const std::array<int, rank> &idx) {
    size_t offset = idx[rank - 1];
    for (size_t i = rank - 1; i > 0; i--) {
      offset += sizes[i - 1] * idx[i - 1];
    }
    return offset;
  }

  // Element access
  T &operator()(const std::array<size_t, rank> &idx) {
    size_t offset = idx[rank - 1];
    for (size_t i = rank - 1; i > 0; i--) {
      offset += sizes[i - 1] * idx[i - 1];
    }
    return v[offset];
  }

  // Specializations
  T &operator()(const size_t &i1) { return v[sizes[0] * i1]; }
  T &operator()(const size_t &i1, const size_t &i2) {
    return v[sizes[0] * i1 + sizes[1] * i2];
  }
  T &operator()(const size_t &i1, const size_t &i2, const size_t &i3) {
    return v[sizes[0] * i1 + sizes[1] * i2 + sizes[2] * i3];
  }
  T &operator()(const size_t &i1, const size_t &i2, const size_t &i3,
                const size_t &i4) {
    return v[sizes[0] * i1 + sizes[1] * i2 + sizes[2] * i3 + sizes[3] * i4];
  }
  T &operator()(const size_t &i1, const size_t &i2, const size_t &i3,
                const size_t &i4, const size_t &i5) {
    return v[sizes[0] * i1 + sizes[1] * i2 + sizes[2] * i3 + sizes[3] * i4 +
             sizes[4] * i5];
  }

  void print() {
    std::cout << "Size: " << size() << std::endl;
    for (size_t i = 0; i < rank; i++) {
      std::cout << "Rank: " << i << " Dimension: " << N[i]
                << " Sizes: " << sizes[i] << std::endl;
    }
  }
};

#endif
