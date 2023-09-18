#pragma once

#include "utils.hpp"

constexpr std::uint8_t TopLeftMask8 = 0b10000000;
constexpr std::uint8_t TopMask8 = 0b01000000;
constexpr std::uint8_t TopRightMask8 = 0b00100000;
constexpr std::uint8_t LeftMask8 = 0b00010000;
constexpr std::uint8_t RightMask8 = 0b00001000;
constexpr std::uint8_t BottomLeftMask8 = 0b00000100;
constexpr std::uint8_t BottomMask8 = 0b00000010;
constexpr std::uint8_t BottomRightMask8 = 0b00000001;
constexpr std::uint8_t NopeMask8 = 0b00000000;
constexpr std::uint8_t AllMask8 = 0b11111111;

template<class i_t, class f_t>
class lookup8
{
public:
	i_t nx = 0;

	std::array<i_t, 9> adder;
	std::array<f_t, 9> dxer;

	std::array<i_t, 256> NeighbourerD8;
	std::array<f_t, 256> NeighbourerD8dx;
	std::array<i_t, 256> NeighbourerNN;
	std::array<std::array<i_t>, 256> Neighbourer;
	std::array<std::array<f_t>, 256> Neighbourerdx;

	lookup8() { ; };
	lookup8(i_t nx, i_t ny, f_t dx, f_t dy)
	{

		f_t dxy = std::sqrt(std::pow(dx, 2) + std::pow(dy, 2));
		this->nx = nx;
		this->adder = { -this->nx - 1, -this->nx, -this->nx + 1, -1, 1,
										this->nx - 1,	 this->nx,	this->nx + 1,	 0. };
		this->dxer = { dxy, dy, dxy, dx, dx, dxy, dy, dxy, 0. };

		this->_compute_lookup_tables();
	};

	void _local_lookup(uint8_t indices,
										 std::array<int, 8>& arr,
										 std::array<f_t, 8>& arrdx,
										 uint8_t& nn)
	{

		// Retrieve the indices specified by the bits set in the `indices` value
		for (uint8_t i = 0; i < 8; ++i) {
			if (indices & (1 << i)) {
				// Index `i` is set, process the corresponding value
				arr[nn] = this->adder[7 - i];
				arrdx[nn] = this->dxer[7 - i];
				++nn;
			}
		}
	}

	void _compute_lookup_tables()
	{
		for (int i = 0; i < LOOKUPSIZE; ++i) {
			uint8_t ti = static_cast<uint8_t>(i);
			this->NeighbourerNN[ti] = 0;
			this->NeighbourerD8[ti] = 0;
			for (auto& v : this->Neighbourer[ti])
				v = 0;
			this->_local_lookup(ti, this->Neighbourer[ti], this->NeighbourerNN[ti]);
			if (this->NeighbourerNN[ti] == 1) {
				this->NeighbourerD8[ti] = this->Neighbourer[ti][0];
				this->NeighbourerD8dx[ti] = this->Neighbourerdx[ti][0];
			}
		}
	}
};
