#pragma once

#include "bithelper.hpp"
#include "npy.hpp"
#include "utils.hpp"
#include <chrono>
#include <immintrin.h>
#include <iostream>
#include <string>
#include <vector>

// static inline
// __m256i blendvps_si256(__m256i& a, __m256i& b, __m256i& mask) {
// 	__m256 res = _mm256_blendv_ps(_mm256_castsi256_ps(a),
// _mm256_castsi256_ps(b), _mm256_castsi256_ps(mask)); 	return
// _mm256_castps_si256(res);
// }

template <class fT> class QD8 {
public:
  int nx, ny, nxy;
  fT dx, dy, dxy;
  std::vector<fT> *_topography = nullptr;
  std::vector<int> Sreceivers;
  std::vector<std::uint8_t> Sdistance2receivers;
  std::vector<std::uint8_t> Sdonors;
  std::vector<int> Neighbours;
  std::vector<std::uint8_t> classic_nSdonors;
  std::vector<int> classic_Sdonors;

  const std::uint8_t TopLeftMask = 0b10000000;
  const std::uint8_t TopMask = 0b01000000;
  const std::uint8_t TopRightMask = 0b00100000;
  const std::uint8_t LeftMask = 0b00010000;
  const std::uint8_t RightMask = 0b00001000;
  const std::uint8_t BottomLeftMask = 0b00000100;
  const std::uint8_t BottomMask = 0b00000010;
  const std::uint8_t BottomRightMask = 0b00000001;
  const std::vector<std::uint8_t> NeiMasks = {
      0b10000000, 0b01000000, 0b00100000, 0b00010000,
      0b00001000, 0b00000100, 0b00000010, 0b00000001};
  const std::vector<int> NeiMasksint = {
      static_cast<int>(0b10000000), static_cast<int>(0b01000000),
      static_cast<int>(0b00100000), static_cast<int>(0b00010000),
      static_cast<int>(0b00001000), static_cast<int>(0b00000100),
      static_cast<int>(0b00000010), static_cast<int>(0b00000001)};
  std::vector<int> adder;

  QD8(){};
  QD8(int nx, int ny, fT dx, fT dy) {
    this->nx = nx;
    this->ny = ny;
    this->nxy = nx * ny;
    this->dx = dx;
    this->dy = dy;
    this->dxy = std::sqrt(std::pow(dx, 2) + std::pow(dy, 2));
    this->Sreceivers = std::vector<int>(this->nxy, 0);
    this->Sdonors = std::vector<std::uint8_t>(this->nxy, 0);
    this->Neighbours = std::vector<int>(this->nxy, 0);
    this->classic_nSdonors = std::vector<std::uint8_t>(this->nxy, 0);
    this->classic_Sdonors = std::vector<int>(this->nxy * 8, 0);
    this->adder = {-this->nx - 1, -this->nx, -this->nx + 1, -1, +1,
                   +this->nx - 1, +this->nx, +this->nx + 1};
  }

  void connect_topography(std::vector<fT> &topo) { this->_topography = &topo; }

  fT &topography(int i) { return (*this->_topography)[i]; }

  void internal_neighbours(int node, std::array<int, 8> &neighbours) {
    neighbours[0] = node - this->nx - 1;
    neighbours[1] = node - this->nx;
    neighbours[2] = node - this->nx + 1;
    neighbours[3] = node - 1;
    neighbours[4] = node + 1;
    neighbours[5] = node + this->nx - 1;
    neighbours[6] = node + this->nx;
    neighbours[7] = node + this->nx + 1;
  }

  int tl_neighbour(int &i) { return i - this->nx - 1; }
  int t_neighbour(int &i) { return i - this->nx; }
  int tr_neighbour(int &i) { return i - this->nx + 1; }
  int l_neighbour(int &i) { return i - 1; }
  int r_neighbour(int &i) { return i + 1; }
  int bl_neighbour(int &i) { return i + this->nx - 1; }
  int b_neighbour(int &i) { return i + this->nx; }
  int br_neighbour(int &i) { return i + this->nx + 1; }

  std::array<fT, 8> dxneighbour() {
    return {this->dxy, this->dy,  this->dxy, this->dx,
            this->dx,  this->dxy, this->dy,  this->dxy};
  }
  __m256 dxneighbour_avx() {
    return _mm256_set_ps(this->dxy, this->dy, this->dxy, this->dx, this->dx,
                         this->dxy, this->dy, this->dxy);
  }

  __m256 fulldx() { return _mm256_set1_ps(this->dx); }
  __m256 fulldy() { return _mm256_set1_ps(this->dy); }
  __m256 fulldxy() { return _mm256_set1_ps(this->dxy); }

  void _save_SSi(__m256i &A, __m256i &B, __m256 &C, __m256 &D) {
    __m256 comparisonResult = _mm256_cmp_ps(D, C, _CMP_GT_OQ);
    C = _mm256_blendv_ps(C, D, comparisonResult);
    __m256i comparisonMask = _mm256_castps_si256(comparisonResult);
    A = blendvps_si256(A, B, comparisonMask);
  }

  void _save_SSi(__m256i &A, __m256i &B, __m256 &C, __m256 &D, __m256i &E,
                 __m256i &F) {
    __m256 comparisonResult = _mm256_cmp_ps(D, C, _CMP_GT_OQ);
    C = _mm256_blendv_ps(C, D, comparisonResult);
    __m256i comparisonMask = _mm256_castps_si256(comparisonResult);
    A = blendvps_si256(A, B, comparisonMask);
    E = blendvps_si256(E, F, comparisonMask);
  }

  __m256 _load_topo(int node) {
    return _mm256_set_ps(this->topography(node), this->topography(node + 1),
                         this->topography(node + 2), this->topography(node + 3),
                         this->topography(node + 4), this->topography(node + 5),
                         this->topography(node + 6),
                         this->topography(node + 7));
  }

  __m256i _load_topoi(int node) {
    return _mm256_set_epi32((node), (node + 1), (node + 2), (node + 3),
                            (node + 4), (node + 5), (node + 6), (node + 7));
  }

  __m256i _load_Sdonoi(int node) {
    return _mm256_set_epi32(static_cast<int>(this->Sdonors[node + 0]),
                            static_cast<int>(this->Sdonors[node + 1]),
                            static_cast<int>(this->Sdonors[node + 2]),
                            static_cast<int>(this->Sdonors[node + 3]),
                            static_cast<int>(this->Sdonors[node + 4]),
                            static_cast<int>(this->Sdonors[node + 5]),
                            static_cast<int>(this->Sdonors[node + 6]),
                            static_cast<int>(this->Sdonors[node + 7]));
  }
  __m256i _load_Sreci(int node) {
    return _mm256_set_epi32(
        (this->Sreceivers[node + 0]), (this->Sreceivers[node + 1]),
        (this->Sreceivers[node + 2]), (this->Sreceivers[node + 3]),
        (this->Sreceivers[node + 4]), (this->Sreceivers[node + 5]),
        (this->Sreceivers[node + 6]), (this->Sreceivers[node + 7]));
  }

  __m256 _load_tl_topo(int node) {
    return _mm256_set_ps(this->topography(-this->nx - 1 + node),
                         this->topography(-this->nx - 1 + node + 1),
                         this->topography(-this->nx - 1 + node + 2),
                         this->topography(-this->nx - 1 + node + 3),
                         this->topography(-this->nx - 1 + node + 4),
                         this->topography(-this->nx - 1 + node + 5),
                         this->topography(-this->nx - 1 + node + 6),
                         this->topography(-this->nx - 1 + node + 7));
  }

  __m256i _load_tl_topoi(int node) {
    return _mm256_set_epi32(
        (-this->nx - 1 + node), (-this->nx - 1 + node + 1),
        (-this->nx - 1 + node + 2), (-this->nx - 1 + node + 3),
        (-this->nx - 1 + node + 4), (-this->nx - 1 + node + 5),
        (-this->nx - 1 + node + 6), (-this->nx - 1 + node + 7));
  }

  __m256 _load_t_topo(int node) {
    return _mm256_set_ps(this->topography(-this->nx + node),
                         this->topography(-this->nx + node + 1),
                         this->topography(-this->nx + node + 2),
                         this->topography(-this->nx + node + 3),
                         this->topography(-this->nx + node + 4),
                         this->topography(-this->nx + node + 5),
                         this->topography(-this->nx + node + 6),
                         this->topography(-this->nx + node + 7));
  }

  __m256i _load_t_topoi(int node) {
    return _mm256_set_epi32((-this->nx + node), (-this->nx + node + 1),
                            (-this->nx + node + 2), (-this->nx + node + 3),
                            (-this->nx + node + 4), (-this->nx + node + 5),
                            (-this->nx + node + 6), (-this->nx + node + 7));
  }

  __m256 _load_tr_topo(int node) {
    return _mm256_set_ps(this->topography(-this->nx + 1 + node),
                         this->topography(-this->nx + 1 + node + 1),
                         this->topography(-this->nx + 1 + node + 2),
                         this->topography(-this->nx + 1 + node + 3),
                         this->topography(-this->nx + 1 + node + 4),
                         this->topography(-this->nx + 1 + node + 5),
                         this->topography(-this->nx + 1 + node + 6),
                         this->topography(-this->nx + 1 + node + 7));
  }

  __m256i _load_tr_topoi(int node) {
    return _mm256_set_epi32(
        (-this->nx + 1 + node), (-this->nx + 1 + node + 1),
        (-this->nx + 1 + node + 2), (-this->nx + 1 + node + 3),
        (-this->nx + 1 + node + 4), (-this->nx + 1 + node + 5),
        (-this->nx + 1 + node + 6), (-this->nx + 1 + node + 7));
  }

  __m256 _load_l_topo(int node) {
    return _mm256_set_ps(
        this->topography(-1 + node), this->topography(-1 + node + 1),
        this->topography(-1 + node + 2), this->topography(-1 + node + 3),
        this->topography(-1 + node + 4), this->topography(-1 + node + 5),
        this->topography(-1 + node + 6), this->topography(-1 + node + 7));
  }

  __m256i _load_l_topoi(int node) {
    return _mm256_set_epi32((-1 + node), (-1 + node + 1), (-1 + node + 2),
                            (-1 + node + 3), (-1 + node + 4), (-1 + node + 5),
                            (-1 + node + 6), (-1 + node + 7));
  }

  __m256 _load_r_topo(int node) {
    return _mm256_set_ps(
        this->topography(1 + node), this->topography(1 + node + 1),
        this->topography(1 + node + 2), this->topography(1 + node + 3),
        this->topography(1 + node + 4), this->topography(1 + node + 5),
        this->topography(1 + node + 6), this->topography(1 + node + 7));
  }

  __m256i _load_r_topoi(int node) {
    return _mm256_set_epi32((1 + node), (1 + node + 1), (1 + node + 2),
                            (1 + node + 3), (1 + node + 4), (1 + node + 5),
                            (1 + node + 6), (1 + node + 7));
  }

  __m256 _load_bl_topo(int node) {
    return _mm256_set_ps(this->topography(this->nx - 1 + node),
                         this->topography(this->nx - 1 + node + 1),
                         this->topography(this->nx - 1 + node + 2),
                         this->topography(this->nx - 1 + node + 3),
                         this->topography(this->nx - 1 + node + 4),
                         this->topography(this->nx - 1 + node + 5),
                         this->topography(this->nx - 1 + node + 6),
                         this->topography(this->nx - 1 + node + 7));
  }

  __m256i _load_bl_topoi(int node) {
    return _mm256_set_epi32(
        (this->nx - 1 + node), (this->nx - 1 + node + 1),
        (this->nx - 1 + node + 2), (this->nx - 1 + node + 3),
        (this->nx - 1 + node + 4), (this->nx - 1 + node + 5),
        (this->nx - 1 + node + 6), (this->nx - 1 + node + 7));
  }

  __m256 _load_b_topo(int node) {
    return _mm256_set_ps(this->topography(this->nx + node),
                         this->topography(this->nx + node + 1),
                         this->topography(this->nx + node + 2),
                         this->topography(this->nx + node + 3),
                         this->topography(this->nx + node + 4),
                         this->topography(this->nx + node + 5),
                         this->topography(this->nx + node + 6),
                         this->topography(this->nx + node + 7));
  }

  __m256i _load_b_topoi(int node) {
    return _mm256_set_epi32((this->nx + node), (this->nx + node + 1),
                            (this->nx + node + 2), (this->nx + node + 3),
                            (this->nx + node + 4), (this->nx + node + 5),
                            (this->nx + node + 6), (this->nx + node + 7));
  }

  __m256 _load_br_topo(int node) {
    return _mm256_set_ps(this->topography(this->nx + 1 + node),
                         this->topography(this->nx + 1 + node + 1),
                         this->topography(this->nx + 1 + node + 2),
                         this->topography(this->nx + 1 + node + 3),
                         this->topography(this->nx + 1 + node + 4),
                         this->topography(this->nx + 1 + node + 5),
                         this->topography(this->nx + 1 + node + 6),
                         this->topography(this->nx + 1 + node + 7));
  }

  __m256i _load_br_topoi(int node) {
    return _mm256_set_epi32(
        (this->nx + 1 + node), (this->nx + 1 + node + 1),
        (this->nx + 1 + node + 2), (this->nx + 1 + node + 3),
        (this->nx + 1 + node + 4), (this->nx + 1 + node + 5),
        (this->nx + 1 + node + 6), (this->nx + 1 + node + 7));
  }

  void compute_sreceivers() {
    if (this->_topography == nullptr)
      throw std::runtime_error("no topo linked sns");

    std::array<int, 8> neighbours;
    auto DXs = this->dxneighbour();
    // int tbmask = 0;

    for (int r = 1; r < this->ny - 1; ++r)
      for (int c = 1; c < this->nx - 1; ++c) {
        int i = r * this->nx + c;
        // std::cout << i << "/" << this->nxy << std::endl;
        this->internal_neighbours(i, neighbours);
        fT SS = 0;
        int ri = i;
        for (int j = 0; j < 8; ++j) {
          int tn = neighbours[j];
          fT tSS = (this->topography(i) - this->topography(tn)) / DXs[j];
          if (tSS > SS) {
            SS = tSS;
            ri = tn;
          }
        }
        this->Sreceivers[i] = ri;
        // this->Sdonors[ri] = this->Sdonors[ri] | this->NeiMasks[tbmask];
      }
  }

  void compute_sreceivers_sdonors() {
    if (this->_topography == nullptr)
      throw std::runtime_error("no topo linked sns");

    std::array<int, 8> neighbours;
    auto DXs = this->dxneighbour();

    for (int r = 1; r < this->ny - 1; ++r)
      for (int c = 1; c < this->nx - 1; ++c) {
        int i = r * this->nx + c;
        int tbmask = 0;

        // std::cout << i << "/" << this->nxy << std::endl;
        this->internal_neighbours(i, neighbours);
        fT SS = 0;
        int ri = i;
        for (int j = 0; j < 8; ++j) {
          int tn = neighbours[j];
          fT tSS = (this->topography(i) - this->topography(tn)) / DXs[j];
          if (tSS > SS) {
            SS = tSS;
            ri = tn;
            ;
            tbmask = j;
          }
        }
        this->Sreceivers[i] = ri;

        if (ri != i) {
          std::uint8_t temp = (this->Sdonors[ri] | this->NeiMasks[tbmask]);
          this->Sdonors[ri] = temp;
        }
      }
  }

  void compute_sreceivers_sdonors_i() {
    if (this->_topography == nullptr)
      throw std::runtime_error("no topo linked sns");

    std::array<int, 8> neighbours;
    auto DXs = this->dxneighbour();

    for (int r = 1; r < this->ny - 1; ++r)
      for (int c = 1; c < this->nx - 1; ++c) {
        int i = r * this->nx + c;
        int tbmask = 0;

        // std::cout << i << "/" << this->nxy << std::endl;
        this->internal_neighbours(i, neighbours);
        fT SS = 0;
        int ri = 8;
        for (int j = 0; j < 8; ++j) {
          int tn = neighbours[j];
          fT tSS = (this->topography(i) - this->topography(tn)) / DXs[j];
          if (tSS > SS) {
            SS = tSS;
            ri = j;
            ;
            tbmask = j;
          }
        }
        this->Sreceivers[i] = ri;
      }

    for (int i = 0; i < this->nxy; ++i) {
      int ri = this->Sreceivers[i];

      if (ri != 8) {
        int node = i + this->adder[ri];

        this->Neighbours[ri] = this->Neighbours[node] | this->NeiMasksint[ri];
      }
    }
  }

  void compute_sreceivers_sdonors_classic() {
    if (this->_topography == nullptr)
      throw std::runtime_error("no topo linked sns");

    std::array<int, 8> neighbours;
    auto DXs = this->dxneighbour();
    int tbmask = 0;

    std::array<int, 8> babun;
    for (int i = 0; i < 8; ++i) {
      babun[i] = this->NeiMasksint[i];
    }
    for (int i = this->nx + 1; i < this->nxy - this->nx - 1; ++i) {
      // int i = r * this->nx + c;
      // std::cout << i << "/" << this->nxy << std::endl;
      this->internal_neighbours(i, neighbours);
      fT SS = 0;
      int ri = i;
      int tmd = 0;

      for (int j = 0; j < 8; ++j) {
        int tn = neighbours[j];
        fT tSS = (this->topography(i) - this->topography(tn)) / DXs[j];
        if (tSS > SS) {
          SS = tSS;
          ri = tn;
          tmd = tn;
          ;
        }
      }

      this->Sreceivers[i] = ri;
      // this->Sdonors[ri] = (this->Sdonors[ri] | this->NeiMasks[tbmask]);
      if (i != Sreceivers[i]) {
        this->Sdonors[this->Sreceivers[i]] = tmd;
        // this->classic_Sdonors[this->Sreceivers[i]*8 +
        // this->classic_nSdonors[this->Sreceivers[i] ] ] = i ;
        // this->classic_Sdonors[this->Sreceivers[i]*8 +
        // this->classic_nSdonors[this->Sreceivers[i] ] ] = babun[tmd] ;
        // // this->classic_Sdonors[this->Sreceivers[i]*8 +
        // this->classic_nSdonors[this->Sreceivers[i] ] ] = i |
        // this->NeiMasks[tbmask];
        // // this->classic_Sdonors[this->Sreceivers[i]*8 +
        // this->classic_nSdonors[this->Sreceivers[i] ] ] = i | ri;
        // // this->classic_Sdonors[this->Sreceivers[i]*8 +
        // this->classic_nSdonors[this->Sreceivers[i] ] ] = i | tmd;
        // // this->classic_Sdonors[this->Sreceivers[i]*8 +
        // this->classic_nSdonors[this->Sreceivers[i] ] ] = i + ri;
        ++this->classic_nSdonors[this->Sreceivers[i]];
      }
    }

    // int j = 0;
    // for(int i=0; i<this->nxy ; ++i)
    // {
    // 	if(i != Sreceivers[i])
    // 	{
    // 		this->classic_Sdonors[this->Sreceivers[i]*8 +
    // this->classic_nSdonors[i]] = (this->classic_Sdonors[i] |
    // this->NeiMasksint[j]) ;
    // 		++this->classic_nSdonors[i];
    // 		++j;
    // 		if(j >= 8)j=0;
    // 	}
    // }
  }

  void compute_sreceivers_sdonors_classicv2() {
    if (this->_topography == nullptr)
      throw std::runtime_error("no topo linked sns");

    std::array<int, 8> neighbours;
    auto DXs = this->dxneighbour();
    int tbmask = 0;

    std::array<int, 8> babun;
    for (int i = 0; i < 8; ++i) {
      babun[i] = this->NeiMasksint[i];
    }
    for (int i = this->nx + 1; i < this->nxy - this->nx - 1; ++i) {
      // int i = r * this->nx + c;
      // std::cout << i << "/" << this->nxy << std::endl;
      this->internal_neighbours(i, neighbours);
      fT SS = 0;
      int ri = i;
      int tmd = 0;

      for (int j = 0; j < 8; ++j) {
        int tn = neighbours[j];
        fT tSS = (this->topography(i) - this->topography(tn)) / DXs[j];
        if (tSS > SS) {
          SS = tSS;
          ri = tn;
          tmd = j;
        }
      }

      this->Sreceivers[i] = ri;
      // this->Sdonors[ri] = (this->Sdonors[ri] | this->NeiMasks[tbmask]);
      if (i != Sreceivers[i]) {
        // this->Sdonors[this->Sreceivers[i]] =
        // this->Sdonors[this->Sreceivers[i]] | tmd;
        this->classic_Sdonors[this->Sreceivers[i] * 8 +
                              this->classic_nSdonors[this->Sreceivers[i]]] = i;
        // this->classic_Sdonors[this->Sreceivers[i]*8 +
        // this->classic_nSdonors[this->Sreceivers[i] ] ] = babun[tmd] ;
        // // this->classic_Sdonors[this->Sreceivers[i]*8 +
        // this->classic_nSdonors[this->Sreceivers[i] ] ] = i |
        // this->NeiMasks[tbmask];
        // // this->classic_Sdonors[this->Sreceivers[i]*8 +
        // this->classic_nSdonors[this->Sreceivers[i] ] ] = i | ri;
        // // this->classic_Sdonors[this->Sreceivers[i]*8 +
        // this->classic_nSdonors[this->Sreceivers[i] ] ] = i | tmd;
        // // this->classic_Sdonors[this->Sreceivers[i]*8 +
        // this->classic_nSdonors[this->Sreceivers[i] ] ] = i + ri;
        ++this->classic_nSdonors[this->Sreceivers[i]];
      }
    }

    // int j = 0;
    // for(int i=0; i<this->nxy ; ++i)
    // {
    // 	if(i != Sreceivers[i])
    // 	{
    // 		this->classic_Sdonors[this->Sreceivers[i]*8 +
    // this->classic_nSdonors[i]] = (this->classic_Sdonors[i] |
    // this->NeiMasksint[j]) ;
    // 		++this->classic_nSdonors[i];
    // 		++j;
    // 		if(j >= 8)j=0;
    // 	}
    // }
  }

  void compute_sreceivers_avx_t1() {
    if (this->_topography == nullptr)
      throw std::runtime_error("no topo linked sns");

    std::array<int, 8> neighbours;
    auto DXs = this->dxneighbour_avx();

    for (int r = 1; r < this->ny - 1; ++r)
      for (int c = 1; c < this->nx - 1; ++c) {
        int i = r * this->nx + c;
        this->internal_neighbours(i, neighbours);
        fT SS = 0;
        int ri = i;
        __m256 topos_i = _mm256_set1_ps(this->topography(i));
        __m256 topos = _mm256_set_ps(
            this->topography(neighbours[0]), this->topography(neighbours[1]),
            this->topography(neighbours[2]), this->topography(neighbours[3]),
            this->topography(neighbours[4]), this->topography(neighbours[5]),
            this->topography(neighbours[6]), this->topography(neighbours[7]));

        __m256 inter = _mm256_sub_ps(topos_i, topos);
        __m256 SSs = _mm256_div_ps(inter, DXs);
        float res[8];
        _mm256_storeu_ps(res, SSs);
        for (int j = 0; j < 8; ++j) {
          if (res[j] > SS) {
            SS = res[j];
            ri = neighbours[j];
          }
        }
        this->Sreceivers[i] = ri;
      }
  }

  void compute_sreceivers_avx_t2() {
    if (this->_topography == nullptr)
      throw std::runtime_error("no topo linked sns");

    std::array<int, 8> neighbours;

    __m256 DXs = fulldx();
    __m256 DYs = fulldy();
    __m256 DXYs = fulldxy();

    for (int i = this->nx + 1; i < this->nxy - this->nx - 1; i += 8) {

      __m256 topos = this->_load_topo(i);
      __m256i recs = this->_load_topoi(i);
      __m256 SSs = _mm256_set1_ps(0.f);

      __m256 topon = this->_load_tl_topo(i);
      __m256i recn = this->_load_tl_topoi(i);
      __m256 inter = _mm256_sub_ps(topos, topon);
      __m256 tSS = _mm256_div_ps(inter, DXYs);
      this->_save_SSi(recs, recn, SSs, tSS);

      topon = this->_load_t_topo(i);
      recn = this->_load_t_topoi(i);
      inter = _mm256_sub_ps(topos, topon);
      tSS = _mm256_div_ps(inter, DYs);
      this->_save_SSi(recs, recn, SSs, tSS);

      topon = this->_load_tr_topo(i);
      recn = this->_load_tr_topoi(i);
      inter = _mm256_sub_ps(topos, topon);
      tSS = _mm256_div_ps(inter, DXYs);
      this->_save_SSi(recs, recn, SSs, tSS);

      topon = this->_load_l_topo(i);
      recn = this->_load_l_topoi(i);
      inter = _mm256_sub_ps(topos, topon);
      tSS = _mm256_div_ps(inter, DXs);
      this->_save_SSi(recs, recn, SSs, tSS);

      topon = this->_load_r_topo(i);
      recn = this->_load_r_topoi(i);
      inter = _mm256_sub_ps(topos, topon);
      tSS = _mm256_div_ps(inter, DXs);
      this->_save_SSi(recs, recn, SSs, tSS);

      topon = this->_load_bl_topo(i);
      recn = this->_load_bl_topoi(i);
      inter = _mm256_sub_ps(topos, topon);
      tSS = _mm256_div_ps(inter, DXYs);
      this->_save_SSi(recs, recn, SSs, tSS);

      topon = this->_load_b_topo(i);
      recn = this->_load_b_topoi(i);
      inter = _mm256_sub_ps(topos, topon);
      tSS = _mm256_div_ps(inter, DYs);
      this->_save_SSi(recs, recn, SSs, tSS);

      topon = this->_load_br_topo(i);
      recn = this->_load_br_topoi(i);
      inter = _mm256_sub_ps(topos, topon);
      tSS = _mm256_div_ps(inter, DXYs);
      this->_save_SSi(recs, recn, SSs, tSS);

      int res[8];
      // _mm256_storeu_epi32(res,recn);
      _mm256_storeu_si256((__m256i *)res, recs);

      int jj = 0;
      for (int j = 7; j >= 0; --j) {
        // std::cout << res[j] << "|";
        this->Sreceivers[i + jj] = res[j];
        ++jj;
      }
      // this->Sreceivers[i] = neighbours[j];
    }
  }

  void compute_sreceivers_avx_t3() {
    if (this->_topography == nullptr)
      throw std::runtime_error("no topo linked sns");

    std::array<int, 8> neighbours;

    __m256 DXs = fulldx();
    __m256 DYs = fulldy();
    __m256 DXYs = fulldxy();

#pragma omp parallel for num_threads(8)
    for (int i = this->nx + 1; i < this->nxy - this->nx - 1; i += 8) {

      __m256 topos = this->_load_topo(i);
      __m256i recs = this->_load_topoi(i);
      __m256 SSs = _mm256_set1_ps(0.f);

      __m256 topon = this->_load_tl_topo(i);
      __m256i recn = this->_load_tl_topoi(i);
      __m256 inter = _mm256_sub_ps(topos, topon);
      __m256 tSS = _mm256_div_ps(inter, DXYs);
      this->_save_SSi(recs, recn, SSs, tSS);

      topon = this->_load_t_topo(i);
      recn = this->_load_t_topoi(i);
      inter = _mm256_sub_ps(topos, topon);
      tSS = _mm256_div_ps(inter, DYs);
      this->_save_SSi(recs, recn, SSs, tSS);

      topon = this->_load_tr_topo(i);
      recn = this->_load_tr_topoi(i);
      inter = _mm256_sub_ps(topos, topon);
      tSS = _mm256_div_ps(inter, DXYs);
      this->_save_SSi(recs, recn, SSs, tSS);

      topon = this->_load_l_topo(i);
      recn = this->_load_l_topoi(i);
      inter = _mm256_sub_ps(topos, topon);
      tSS = _mm256_div_ps(inter, DXs);
      this->_save_SSi(recs, recn, SSs, tSS);

      topon = this->_load_r_topo(i);
      recn = this->_load_r_topoi(i);
      inter = _mm256_sub_ps(topos, topon);
      tSS = _mm256_div_ps(inter, DXs);
      this->_save_SSi(recs, recn, SSs, tSS);

      topon = this->_load_bl_topo(i);
      recn = this->_load_bl_topoi(i);
      inter = _mm256_sub_ps(topos, topon);
      tSS = _mm256_div_ps(inter, DXYs);
      this->_save_SSi(recs, recn, SSs, tSS);

      topon = this->_load_b_topo(i);
      recn = this->_load_b_topoi(i);
      inter = _mm256_sub_ps(topos, topon);
      tSS = _mm256_div_ps(inter, DYs);
      this->_save_SSi(recs, recn, SSs, tSS);

      topon = this->_load_br_topo(i);
      recn = this->_load_br_topoi(i);
      inter = _mm256_sub_ps(topos, topon);
      tSS = _mm256_div_ps(inter, DXYs);
      this->_save_SSi(recs, recn, SSs, tSS);

      int res[8];
      // _mm256_storeu_epi32(res,recn);
      _mm256_storeu_si256((__m256i *)res, recs);

      int jj = 0;
      for (int j = 7; j >= 0; --j) {
        // std::cout << res[j] << "|";
        this->Sreceivers[i + jj] = res[j];
        ++jj;
      }
      // this->Sreceivers[i] = neighbours[j];
    }
  }

  void compute_sreceivers_avx_t4() {
    // if(this->_topography == nullptr) throw std::runtime_error("no topo linked
    // sns");

    // std::array<int,8> neighbours;

    // __m256 DXs = fulldx();
    // __m256 DYs = fulldy();
    // __m256 DXYs = fulldxy();

    // for (int i = this->nx + 1; i < this->nxy - this->nx - 1; i+=8)
    // {

    // 	__m256 topos = this->_load_topo(i);
    // 	__m256i recs = this->_load_topoi(i);
    // 	__m256 SSs = _mm256_set1_ps(0.f);

    // 	__m256 topon = this->_load_tl_topo(i);
    // 	__m256i recn = this->_load_tl_topoi(i);
    // 	__m256 inter = _mm256_sub_ps(topos,topon);
    // 	__m256 tSS = _mm256_div_ps(inter,DXYs);
    // 	__m256i
    // 	this->_save_SSi( recs,  recn,  SSs,  tSS);

    // 	topon = this->_load_t_topo(i);
    // 	recn = this->_load_t_topoi(i);
    // 	inter = _mm256_sub_ps(topos,topon);
    // 	tSS = _mm256_div_ps(inter,DYs);
    // 	this->_save_SSi( recs,  recn,  SSs,  tSS);

    // 	topon = this->_load_tr_topo(i);
    // 	recn = this->_load_tr_topoi(i);
    // 	inter = _mm256_sub_ps(topos,topon);
    // 	tSS = _mm256_div_ps(inter,DXYs);
    // 	this->_save_SSi( recs,  recn,  SSs,  tSS);

    // 	topon = this->_load_l_topo(i);
    // 	recn = this->_load_l_topoi(i);
    // 	inter = _mm256_sub_ps(topos,topon);
    // 	tSS = _mm256_div_ps(inter,DXs);
    // 	this->_save_SSi( recs,  recn,  SSs,  tSS);

    // 	topon = this->_load_r_topo(i);
    // 	recn = this->_load_r_topoi(i);
    // 	inter = _mm256_sub_ps(topos,topon);
    // 	tSS = _mm256_div_ps(inter,DXs);
    // 	this->_save_SSi( recs,  recn,  SSs,  tSS);

    // 	topon = this->_load_bl_topo(i);
    // 	recn = this->_load_bl_topoi(i);
    // 	inter = _mm256_sub_ps(topos,topon);
    // 	tSS = _mm256_div_ps(inter,DXYs);
    // 	this->_save_SSi( recs,  recn,  SSs,  tSS);

    // 	topon = this->_load_b_topo(i);
    // 	recn = this->_load_b_topoi(i);
    // 	inter = _mm256_sub_ps(topos,topon);
    // 	tSS = _mm256_div_ps(inter,DYs);
    // 	this->_save_SSi( recs,  recn,  SSs,  tSS);

    // 	topon = this->_load_br_topo(i);
    // 	recn = this->_load_br_topoi(i);
    // 	inter = _mm256_sub_ps(topos,topon);
    // 	tSS = _mm256_div_ps(inter,DXYs);
    // 	this->_save_SSi( recs,  recn,  SSs,  tSS);

    // 	int res[8];
    // 	// _mm256_storeu_epi32(res,recn);
    // 	_mm256_storeu_si256((__m256i*)res, recs);

    // 	int jj = 0;
    // 	for(int j = 7; j>=0; --j)
    // 	{
    // 		// std::cout << res[j] << "|";
    // 		this->Sreceivers[i+jj] = res[j];
    // 		++jj;
    // 	}
    // 	// this->Sreceivers[i] = neighbours[j];
    // }
  }

  // Utilities

  void save_topo(std::string filename) {
    const std::vector<std::uint64_t> shape{std::uint64_t(this->ny),
                                           std::uint64_t(this->nx)};
    const bool fortran_order{false};
    const std::string path{filename};
    npy::SaveArrayAsNumpy(path, fortran_order, shape.size(), shape.data(),
                          (*this->_topography));
  }

  void save_Sreceivers(std::string filename) {
    const std::vector<std::uint64_t> shape{std::uint64_t(this->ny),
                                           std::uint64_t(this->nx)};
    const bool fortran_order{false};
    const std::string path{filename};
    npy::SaveArrayAsNumpy(path, fortran_order, shape.size(), shape.data(),
                          this->Sreceivers);
  }
};
