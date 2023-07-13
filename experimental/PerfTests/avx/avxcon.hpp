#pragma once

#include "utils.hpp"
#include "npy.hpp" 
#include "bithelper.hpp" 

#include <bitset>
#include <chrono>
#include <cstdint>
#include <immintrin.h>
#include <iostream>
#include <string>
#include <vector>

constexpr int LOOKUPSIZE = 256;


// These are the neighbour masks
// 8 zeros:  00000000
// 16 zeros: 0000000000000000
// 24 zeros: 000000000000000000000000

// No neighbours at all, also false
constexpr uint32_t NoNeigh            = 0;

// Masks for the Donors
constexpr uint32_t DonTopLeftMask     = 0b10000000000000000000000000000000; // woud be better to have 1 << 32
constexpr uint32_t DonTopMask         = 0b01000000000000000000000000000000;
constexpr uint32_t DonTopRightMask    = 0b00100000000000000000000000000000;
constexpr uint32_t DonLeftMask        = 0b00010000000000000000000000000000;
constexpr uint32_t DonRightMask       = 0b00001000000000000000000000000000;
constexpr uint32_t DonBottomLeftMask  = 0b00000100000000000000000000000000;
constexpr uint32_t DonBottomMask      = 0b00000010000000000000000000000000;
constexpr uint32_t DonBottomRightMask = 0b00000001000000000000000000000000;

// Masks for the receivers
constexpr uint32_t RecTopLeftMask     = 0b00000000100000000000000000000000;
constexpr uint32_t RecTopMask         = 0b00000000010000000000000000000000;
constexpr uint32_t RecTopRightMask    = 0b00000000001000000000000000000000;
constexpr uint32_t RecLeftMask        = 0b00000000000100000000000000000000;
constexpr uint32_t RecRightMask       = 0b00000000000010000000000000000000;
constexpr uint32_t RecBottomLeftMask  = 0b00000000000001000000000000000000;
constexpr uint32_t RecBottomMask      = 0b00000000000000100000000000000000;
constexpr uint32_t RecBottomRightMask = 0b00000000000000010000000000000000;

// Masks for the SD donors
constexpr uint32_t SDonTopLeftMask     = 0b00000000000000001000000000000000;
constexpr uint32_t SDonTopMask         = 0b00000000000000000100000000000000;
constexpr uint32_t SDonTopRightMask    = 0b00000000000000000010000000000000;
constexpr uint32_t SDonLeftMask        = 0b00000000000000000001000000000000;
constexpr uint32_t SDonRightMask       = 0b00000000000000000000100000000000;
constexpr uint32_t SDonBottomLeftMask  = 0b00000000000000000000010000000000;
constexpr uint32_t SDonBottomMask      = 0b00000000000000000000001000000000;
constexpr uint32_t SDonBottomRightMask = 0b00000000000000000000000100000000;

// Masks for the SD receivers
constexpr uint32_t SRecTopLeftMask     = 0b00000000000000000000000010000000;// 1
constexpr uint32_t SRecTopMask         = 0b00000000000000000000000001000000;// 3
constexpr uint32_t SRecTopRightMask    = 0b00000000000000000000000000100000;// 5
constexpr uint32_t SRecLeftMask        = 0b00000000000000000000000000010000;// 7
constexpr uint32_t SRecRightMask       = 0b00000000000000000000000000001000;// 9
constexpr uint32_t SRecBottomLeftMask  = 0b00000000000000000000000000000100;// 11
constexpr uint32_t SRecBottomMask      = 0b00000000000000000000000000000010;// 13
constexpr uint32_t SRecBottomRightMask = 0b00000000000000000000000000000001;// 15

constexpr std::array<uint32_t, 8> DonMasks{DonTopLeftMask,DonTopMask,DonTopRightMask,DonLeftMask,DonRightMask,DonBottomLeftMask,DonBottomMask,DonBottomRightMask};
constexpr std::array<uint32_t, 8> RecMasks{RecTopLeftMask,RecTopMask,RecTopRightMask,RecLeftMask,RecRightMask,RecBottomLeftMask,RecBottomMask,RecBottomRightMask};
constexpr std::array<uint32_t, 8> SDonMasks{SDonTopLeftMask,SDonTopMask,SDonTopRightMask,SDonLeftMask,SDonRightMask,SDonBottomLeftMask,SDonBottomMask,SDonBottomRightMask};
constexpr std::array<uint32_t, 8> SRecMasks{SRecTopLeftMask,SRecTopMask,SRecTopRightMask,SRecLeftMask,SRecRightMask,SRecBottomLeftMask,SRecBottomMask,SRecBottomRightMask};
constexpr std::array<uint32_t, 8> MirrorIndicesMasks{7,6,5,4,3,2,1,0};

// static inline
// __m256i blendvps_si256(__m256i& a, __m256i& b, __m256i& mask) {
// 	__m256 res = _mm256_blendv_ps(_mm256_castsi256_ps(a), _mm256_castsi256_ps(b), _mm256_castsi256_ps(mask));
// 	return _mm256_castps_si256(res);
// }


template<class fT>
class conavx
{
public:

	// dimensions in number of nodes
	int nx,ny,nxy;
	// Dimensions in spatial directions
	fT dx,dy,dxy, cellarea;
	// Topography to be glued to the connector
	std::vector<fT>* _topography = nullptr;
	// Contains all the neighbours in successive bit sets:
	// - the first 8 are the donors
	// - the next 8 the receivers
	// - the next 8 the Sdonors
	// - and the final 8 the Sreceivers 	
	std::vector<int> Neighbours;

	std::vector<uint8_t> boundaries;


	



	// Adder is the regular "normal" neighbourer containing the value to add to node i to get it's respective neighbours
	std::array<int,9> adder;
	std::array<float,9> dxer;
	std::array<std::array<int,8>, LOOKUPSIZE> lookup_table;
	std::array<uint8_t, LOOKUPSIZE> lookup_nn;
	std::array<int, LOOKUPSIZE> lookup_table_SS;

	std::vector<int> EdgeNodes;
	std::vector<uint8_t> EdgeNodes_Neighbours;
	// std::vector<std::array<int,8> > IntaNodes;
	std::vector<int> IntaNodes;

	conavx() = default;
	conavx(int nx, int ny, fT dx, fT dy)
	{
		this->nx = nx; this->ny = ny; this->nxy = nx*ny;
		this->dx = dx; this->dy = dy; this->dxy = std::sqrt(std::pow(dx,2) + std::pow(dy,2));
		this->Neighbours = std::vector<int>(this->nxy,0);
		this->adder = {- this->nx - 1, - this->nx, - this->nx + 1, - 1, + 1, + this->nx - 1, + this->nx, + this->nx + 1, 0};
		this->_compute_lookup_table();
		this->boundaries = std::vector<uint8_t>(this->nxy,1);
		this->_4edges();
		this->dxer = {this->dxy,this->dy,this->dxy, this->dx,this->dx,this->dxy,this->dy,this->dxy};
		this->cellarea = this->dx * this->dy;

	}

	void connect_topography(std::vector<fT>&topo)
	{
		this->_topography = &topo;
	}

	fT& topography(int i){return (*this->_topography)[i];}


	// return the receiver at node i
	int Sreceivers(int i)
	{
		return i + this->lookup_table_SS[ui8fromui32(this->Neighbours[i],0)];
	}


	int Donors(int node, std::array<int,8>& arr)
	{
		return this->_dneighbours(node, arr, 3);
	}

	int SDonors(int node, std::array<int,8>& arr)
	{
		return this->_dneighbours(node, arr, 1);
	}

	int Receivers(int node, std::array<int,8>& arr)
	{
		return this->_dneighbours(node, arr, 2);
	}

	fT get_cell_area(int node){return this->cellarea;}


	int nDonors(const int node) const // no effect on the data class members
	{
		uint8_t idx = ui8fromui32(this->Neighbours[node],3);
		return static_cast<int>(this->lookup_nn[idx]);
	}

	int nSDonors(const int node) const
	{
		uint8_t idx = ui8fromui32(this->Neighbours[node],1);
		return static_cast<int>(this->lookup_nn[idx]);
	}

	int nReceivers(const int node) const
	{
		uint8_t idx = ui8fromui32(this->Neighbours[node],2);
		return static_cast<int>(this->lookup_nn[idx]);
	}







	int _dneighbours(int node, std::array<int,8>& arr, int bytenum)
	{
		uint8_t idx = ui8fromui32(this->Neighbours[node],bytenum);
		int nn = static_cast<int>(this->lookup_nn[idx]);
		arr = this->lookup_table[idx];
		
		for(size_t i=0;i<static_cast<size_t>(nn);++i)
			arr[i]+=node;

		return nn;
	}

	int _neighbours_from_ui8(int node, uint8_t idx, std::array<int,8>& arr)
	{
		int nn = static_cast<int>(this->lookup_nn[idx]);
		arr = this->lookup_table[idx];
		for(int i=0; i < nn; ++i)
			arr[i]+=node;

		return nn;
	}

	void _neighbours_from_ui8_nocheck(int node, uint8_t idx, std::array<int,8>& arr)
	{
		for (int i = 0; i < 8; ++i) 
		{
			if (idx & (1 << i)) 
			{
				// Index `i` is set, process the corresponding value
				arr[7 - i] = node + this->adder[7-i];
			}
			else
				arr[7-i] = -1;
		}
	}



	int Receivers_archive(int node, std::array<int,8>& arr)
	{
		int nn = 0;
		uint8_t indices = ui8fromui32(this->Neighbours[node],2);
		// if(indices>0)
		// 	std::cout << std::bitset<8>(indices) << std::endl;
		// Retrieve the indices specified by the bits set in the `indices` value
		for (int i = 0; i < 8; ++i) 
		{
			if (indices & (1 << i)) 
			{
				// Index `i` is set, process the corresponding value
				// std::cout << 7 - i << " ";
				arr[nn] = node + this->adder[7 - i];
				++nn;
			}
		}
			// if(nn>0)
			// 	std::cout << std::endl;;

		return nn;

	}

	int Receivers_avx(int node, std::array<int,8>& arr)
	{
		int nn = 0;
		int& ty = this->Neighbours[node] ;
		// Retrieve the indices specified by the bits set in the `indices` value
		for (int i = 16; i < 24; ++i) 
		{
			if ( ty & (1 << i)) 
			{
				// Index `i` is set, process the corresponding value
				// std::cout << 7 - i << " ";
				arr[nn] = node + this->adder[23 - i];
				++nn;
			}
		}

		return nn;

	}

	int Receivers_avx_2(int node, std::array<int,8>& arr)
	{
		uint8_t idx = ui8fromui32(this->Neighbours[node],2);
		int nn = static_cast<int>(this->lookup_nn[idx]);
		arr = this->lookup_table[idx];
		for(size_t i=0;i<static_cast<size_t>(nn);++i)
			arr[i]+=node;

		return nn;
	}


	int Receivers_classics(int node, std::array<int,8>& arr)
	{
		int nn = 0;
		this->internal_neighbours(node, arr);
		// Retrieve the indices specified by the bits set in the `indices` value
		for (int i = 0; i < 8; ++i) 
		{
			int tn = arr[i];
			if(tn<0 || tn >= this->nxy) continue;

			if(this->topography(tn) < this->topography(node))
			{
				arr[nn] = tn;
				++nn;
			}
		}
			// if(nn>0)
			// 	std::cout << std::endl;;

		return nn;

	}



















	void internal_neighbours(int node, std::array<int,8>& neighbours)
	{
		neighbours[0] = node - this->nx - 1;
		neighbours[1] = node - this->nx;
		neighbours[2] = node - this->nx + 1;
		neighbours[3] = node - 1;
		neighbours[4] = node + 1;
		neighbours[5] = node + this->nx - 1;
		neighbours[6] = node + this->nx;
		neighbours[7] = node + this->nx + 1;
	}

	// void internal_neighbours_with_check(int node, std::array<int,8>& neighbours)
	// {
	// 	int nn = 0;
	// 	for(int i=0;i<8;++i)
	// 	{
	// 		int tn = node + adder[i];
	// 			neighbours[i] = tn;
	// 	}
	// }

	int tl_neighbour(int& i){return i - this->nx - 1;}
	int t_neighbour(int& i){return i - this->nx;}
	int tr_neighbour(int& i){return i - this->nx + 1;}
	int l_neighbour(int& i){return i - 1;}
	int r_neighbour(int& i){return i + 1;}
	int bl_neighbour(int& i){return i + this->nx - 1;}
	int b_neighbour(int& i){return i + this->nx;}
	int br_neighbour(int& i){return i + this->nx + 1;}

	std::array<fT,8> dxneighbour(){return {this->dxy,this->dy,this->dxy,this->dx,this->dx,this->dxy,this->dy,this->dxy};}
	__m256 dxneighbour_avx(){return _mm256_set_ps(this->dxy,this->dy,this->dxy,this->dx,this->dx,this->dxy,this->dy,this->dxy);}
	
	__m256 fulldx(){return _mm256_set1_ps(this->dx);}
	__m256 fulldy(){return _mm256_set1_ps(this->dy);}
	__m256 fulldxy(){return _mm256_set1_ps(this->dxy);}

	void _save_SSi(__m256i& A, __m256i& B, __m256& C, __m256& D)
	{
		__m256 comparisonResult = _mm256_cmp_ps(D, C, _CMP_GT_OQ);
		C = _mm256_blendv_ps(C, D, comparisonResult);
		__m256i comparisonMask = _mm256_castps_si256(comparisonResult);
		A = blendvps_si256(A, B, comparisonMask);
	}

	void _save_SSi(__m256i& A, __m256i& B, __m256& C, __m256& D, __m256i& E, __m256i& F)
	{
		__m256 comparisonResult = _mm256_cmp_ps(D, C, _CMP_GT_OQ);
		C = _mm256_blendv_ps(C, D, comparisonResult);
		__m256i comparisonMask = _mm256_castps_si256(comparisonResult);
		A = blendvps_si256(A, B, comparisonMask);
		E = blendvps_si256(E, F, comparisonMask);
	}

	
	__m256 _load_topo(int node)
	{
		return _mm256_set_ps(
			this->topography(node),this->topography(node + 1), 
			this->topography(node + 2), this->topography(node + 3),
			this->topography(node + 4),this->topography(node + 5),
			this->topography(node + 6),this->topography(node + 7)
		);
	}

	__m256i _load_topoi(int node)
	{
		return _mm256_set_epi32(
			(node),(node + 1), 
			(node + 2), (node + 3),
			(node + 4),(node + 5),
			(node + 6),(node + 7)
		);
	}

	__m256 _load_topo(std::array<int,8>& arrnodes)
	{
		return _mm256_set_ps(
			this->topography(arrnodes[0]),this->topography(arrnodes[1]), 
			this->topography(arrnodes[2]), this->topography(arrnodes[3]),
			this->topography(arrnodes[4]),this->topography(arrnodes[5]),
			this->topography(arrnodes[6]),this->topography(arrnodes[7])
		);
	}

	__m256i _load_topoi(std::array<int,8>& arrnodes)
	{
		return _mm256_set_epi32(
			(arrnodes[0]),(arrnodes[1]), 
			(arrnodes[2]), (arrnodes[3]),
			(arrnodes[4]),(arrnodes[5]),
			(arrnodes[6]),(arrnodes[7])
		);
	}

	__m256 _load_topo(std::array<int,8>& arrnodes, int tadder)
	{
		return _mm256_set_ps(
			this->topography(arrnodes[0]+tadder),this->topography(arrnodes[1]+tadder), 
			this->topography(arrnodes[2]+tadder), this->topography(arrnodes[3]+tadder),
			this->topography(arrnodes[4]+tadder),this->topography(arrnodes[5]+tadder),
			this->topography(arrnodes[6]+tadder),this->topography(arrnodes[7]+tadder)
		);
	}

	__m256i _load_topoi(std::array<int,8>& arrnodes, int tadder)
	{
		return _mm256_set_epi32(
			(arrnodes[0]+tadder),(arrnodes[1]+tadder), 
			(arrnodes[2]+tadder), (arrnodes[3]+tadder),
			(arrnodes[4]+tadder),(arrnodes[5]+tadder),
			(arrnodes[6]+tadder),(arrnodes[7]+tadder)
		);
	}

	__m256 _load_topo(int i,std::vector<int>& vecnodes)
	{
		return _mm256_set_ps(
			this->topography(vecnodes[i+0]),this->topography(vecnodes[i+1]), 
			this->topography(vecnodes[i+2]), this->topography(vecnodes[i+3]),
			this->topography(vecnodes[i+4]),this->topography(vecnodes[i+5]),
			this->topography(vecnodes[i+6]),this->topography(vecnodes[i+7])
		);
	}

	__m256i _load_topoi(int i,std::vector<int>& vecnodes)
	{
		return _mm256_set_epi32(
			(vecnodes[i+0]),(vecnodes[i+1]), 
			(vecnodes[i+2]), (vecnodes[i+3]),
			(vecnodes[i+4]),(vecnodes[i+5]),
			(vecnodes[i+6]),(vecnodes[i+7])
		);
	}

	__m256 _load_topo(int i,std::vector<int>& vecnodes, int tadder)
	{
		return _mm256_set_ps(
			this->topography(vecnodes[i+0]+tadder),this->topography(vecnodes[i+1]+tadder), 
			this->topography(vecnodes[i+2]+tadder), this->topography(vecnodes[i+3]+tadder),
			this->topography(vecnodes[i+4]+tadder),this->topography(vecnodes[i+5]+tadder),
			this->topography(vecnodes[i+6]+tadder),this->topography(vecnodes[i+7]+tadder)
		);
	}

	__m256i _load_topoi(int i,std::vector<int>& vecnodes, int tadder)
	{
		return _mm256_set_epi32(
			(vecnodes[i+0]+tadder),(vecnodes[i+1]+tadder), 
			(vecnodes[i+2]+tadder), (vecnodes[i+3]+tadder),
			(vecnodes[i+4]+tadder),(vecnodes[i+5]+tadder),
			(vecnodes[i+6]+tadder),(vecnodes[i+7]+tadder)
		);
	}

	__m256i _load_neighbours(int node)
	{
		return _mm256_set_epi32(
			this->Neighbours[(node)],this->Neighbours[(node + 1)], 
			this->Neighbours[(node + 2)], this->Neighbours[(node + 3)],
			this->Neighbours[(node + 4)],this->Neighbours[(node + 5)],
			this->Neighbours[(node + 6)],this->Neighbours[(node + 7)]
		);
	}
	__m256i _load_neighbours(std::array<int,8>& arrnodes)
	{
		return _mm256_set_epi32(
			this->Neighbours[arrnodes[0]],this->Neighbours[arrnodes[1]], 
			this->Neighbours[arrnodes[2]], this->Neighbours[arrnodes[3]],
			this->Neighbours[arrnodes[4]],this->Neighbours[arrnodes[5]],
			this->Neighbours[arrnodes[6]],this->Neighbours[arrnodes[7]]
		);
	}

	__m256i _load_neighbours(std::array<int,8>& arrnodes, int tadder)
	{
		return _mm256_set_epi32(
			this->Neighbours[arrnodes[0] + tadder],this->Neighbours[arrnodes[1] + tadder], 
			this->Neighbours[arrnodes[2] + tadder], this->Neighbours[arrnodes[3] + tadder],
			this->Neighbours[arrnodes[4] + tadder],this->Neighbours[arrnodes[5] + tadder],
			this->Neighbours[arrnodes[6] + tadder],this->Neighbours[arrnodes[7] + tadder]
		);
	}

	__m256i _load_neighbours(int i, std::vector<int>& vecnodes)
	{
		return _mm256_set_epi32(
			this->Neighbours[vecnodes[0+i]],this->Neighbours[vecnodes[1+i]], 
			this->Neighbours[vecnodes[2+i]], this->Neighbours[vecnodes[3+i]],
			this->Neighbours[vecnodes[4+i]],this->Neighbours[vecnodes[5+i]],
			this->Neighbours[vecnodes[6+i]],this->Neighbours[vecnodes[7+i]]
		);
	}

	__m256i _load_neighbours(int i, std::vector<int>& vecnodes, int tadder)
	{
		return _mm256_set_epi32(
			this->Neighbours[vecnodes[0+i] + tadder],this->Neighbours[vecnodes[1+i] + tadder], 
			this->Neighbours[vecnodes[2+i] + tadder], this->Neighbours[vecnodes[3+i] + tadder],
			this->Neighbours[vecnodes[4+i] + tadder],this->Neighbours[vecnodes[5+i] + tadder],
			this->Neighbours[vecnodes[6+i] + tadder],this->Neighbours[vecnodes[7+i] + tadder]
		);
	}


	__m256 _load_tl_topo(int node)
	{
		return _mm256_set_ps(
			this->topography(- this->nx - 1 + node),this->topography(- this->nx - 1 + node + 1), 
			this->topography(- this->nx - 1 + node + 2), this->topography(- this->nx - 1 + node + 3),
			this->topography(- this->nx - 1 + node + 4),this->topography(- this->nx - 1 + node + 5),
			this->topography(- this->nx - 1 + node + 6),this->topography(- this->nx - 1 + node + 7)
		);
	} 

	__m256i _load_tl_topoi(int node)
	{
		return _mm256_set_epi32(
			(- this->nx - 1 + node),(- this->nx - 1 + node + 1), 
			(- this->nx - 1 + node + 2), (- this->nx - 1 + node + 3),
			(- this->nx - 1 + node + 4),(- this->nx - 1 + node + 5),
			(- this->nx - 1 + node + 6),(- this->nx - 1 + node + 7)
		);
	} 

	__m256 _load_t_topo(int node)
	{
		return _mm256_set_ps(
			this->topography(- this->nx + node),this->topography(- this->nx + node + 1), 
			this->topography(- this->nx + node + 2), this->topography(- this->nx + node + 3),
			this->topography(- this->nx + node + 4),this->topography(- this->nx + node + 5),
			this->topography(- this->nx + node + 6),this->topography(- this->nx + node + 7)
		);
	} 

	__m256i _load_t_topoi(int node)
	{
		return _mm256_set_epi32(
			(- this->nx + node),(- this->nx + node + 1), 
			(- this->nx + node + 2), (- this->nx + node + 3),
			(- this->nx + node + 4),(- this->nx + node + 5),
			(- this->nx + node + 6),(- this->nx + node + 7)
		);
	}  

	__m256 _load_tr_topo(int node)
	{
		return _mm256_set_ps(
			this->topography(- this->nx + 1 + node),this->topography(- this->nx + 1 + node + 1), 
			this->topography(- this->nx + 1 + node + 2), this->topography(- this->nx + 1 + node + 3),
			this->topography(- this->nx + 1 + node + 4),this->topography(- this->nx + 1 + node + 5),
			this->topography(- this->nx + 1 + node + 6),this->topography(- this->nx + 1 + node + 7)
		);
	}    

	__m256i _load_tr_topoi(int node)
	{
		return _mm256_set_epi32(
			(- this->nx + 1 + node),(- this->nx + 1 + node + 1), 
			(- this->nx + 1 + node + 2), (- this->nx + 1 + node + 3),
			(- this->nx + 1 + node + 4),(- this->nx + 1 + node + 5),
			(- this->nx + 1 + node + 6),(- this->nx + 1 + node + 7)
		);
	}    

	__m256 _load_l_topo(int node)
	{
		return _mm256_set_ps(
			this->topography(- 1 + node),this->topography(- 1 + node + 1), 
			this->topography(- 1 + node + 2), this->topography(- 1 + node + 3),
			this->topography(- 1 + node + 4),this->topography(- 1 + node + 5),
			this->topography(- 1 + node + 6),this->topography(- 1 + node + 7)
		);
	} 

	__m256i _load_l_topoi(int node)
	{
		return _mm256_set_epi32(
			(- 1 + node),(- 1 + node + 1), 
			(- 1 + node + 2), (- 1 + node + 3),
			(- 1 + node + 4),(- 1 + node + 5),
			(- 1 + node + 6),(- 1 + node + 7)
		);
	}  

	__m256 _load_r_topo(int node)
	{
		return _mm256_set_ps(
			this->topography( 1 + node),this->topography( 1 + node + 1), 
			this->topography( 1 + node + 2), this->topography( 1 + node + 3),
			this->topography( 1 + node + 4),this->topography( 1 + node + 5),
			this->topography( 1 + node + 6),this->topography( 1 + node + 7)
		);
	}

	__m256i _load_r_topoi(int node)
	{
		return _mm256_set_epi32(
			( 1 + node),( 1 + node + 1), 
			( 1 + node + 2), ( 1 + node + 3),
			( 1 + node + 4),( 1 + node + 5),
			( 1 + node + 6),( 1 + node + 7)
		);
	}  

	__m256 _load_bl_topo(int node)
	{
		return _mm256_set_ps(
			this->topography(this->nx - 1 + node),this->topography(this->nx - 1 + node + 1), 
			this->topography(this->nx - 1 + node + 2), this->topography(this->nx - 1 + node + 3),
			this->topography(this->nx - 1 + node + 4),this->topography(this->nx - 1 + node + 5),
			this->topography(this->nx - 1 + node + 6),this->topography(this->nx - 1 + node + 7)
		);
	}

	__m256i _load_bl_topoi(int node)
	{
		return _mm256_set_epi32(
			(this->nx - 1 + node),(this->nx - 1 + node + 1), 
			(this->nx - 1 + node + 2), (this->nx - 1 + node + 3),
			(this->nx - 1 + node + 4),(this->nx - 1 + node + 5),
			(this->nx - 1 + node + 6),(this->nx - 1 + node + 7)
		);
	}

	__m256 _load_b_topo(int node)
	{
		return _mm256_set_ps(
			this->topography(this->nx + node),this->topography(this->nx + node + 1), 
			this->topography(this->nx + node + 2), this->topography(this->nx + node + 3),
			this->topography(this->nx + node + 4),this->topography(this->nx + node + 5),
			this->topography(this->nx + node + 6),this->topography(this->nx + node + 7)
		);
	}

	__m256i _load_b_topoi(int node)
	{
		return _mm256_set_epi32(
			(this->nx + node),(this->nx + node + 1), 
			(this->nx + node + 2), (this->nx + node + 3),
			(this->nx + node + 4),(this->nx + node + 5),
			(this->nx + node + 6),(this->nx + node + 7)
		);
	}

	__m256 _load_br_topo(int node)
	{
		return _mm256_set_ps(
			this->topography(this->nx + 1 + node),this->topography(this->nx + 1 + node + 1), 
			this->topography(this->nx + 1 + node + 2), this->topography(this->nx + 1 + node + 3),
			this->topography(this->nx + 1 + node + 4),this->topography(this->nx + 1 + node + 5),
			this->topography(this->nx + 1 + node + 6),this->topography(this->nx + 1 + node + 7)
		);
	}

	__m256i _load_br_topoi(int node)
	{
		return _mm256_set_epi32(
			(this->nx + 1 + node),(this->nx + 1 + node + 1), 
			(this->nx + 1 + node + 2), (this->nx + 1 + node + 3),
			(this->nx + 1 + node + 4),(this->nx + 1 + node + 5),
			(this->nx + 1 + node + 6),(this->nx + 1 + node + 7)
		);
	}


	// Does not work
	void compute_archive_2()
	{
		if(this->_topography == nullptr) throw std::runtime_error("no topo linked sns");

		std::array<int,8> neighbours;

		__m256 DXs = fulldx();
		__m256 DYs = fulldy();
		__m256 DXYs = fulldxy();
		__m256 zeroVec = _mm256_setzero_ps();
		__m256i zeroVeci =  _mm256_set1_epi32(0);
		__m256i oneVeci =  _mm256_set1_epi32(-1);


		for (int i = this->nx + 1; i < this->nxy - this->nx - 1; i+=8)
		{
			// storing initial topo and initialising neighbours to no neighbours
			__m256 topos = this->_load_topo(i);
			__m256i neigh = _mm256_set1_epi32(this->NoNeigh);
			// Tracking steepest descent (init to 0.)
			__m256 SSs = _mm256_set1_ps(0.f);
			__m256i SSsi = _mm256_set1_epi32(0);


			// Loading topography from the top left corner
			__m256 topon = this->_load_tl_topo(i);
			// setting up the bitmask to donor-toward-top-left
			__m256i bitma = _mm256_set1_epi32(DonTopLeftMask);
			//## Calculating Slope
			// subtracting current topo with neighbours'
			__m256 inter = _mm256_sub_ps(topos, topon);
			// Dividing diff by diagonals
			__m256 tSS = _mm256_div_ps(inter, DXYs);
			// printM256(tSS);
			//## mask comparing to 0 for normal recs and dons
			// mask where slope < 0 (donors)
			__m256 comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_LT_OQ);
			// recasting mask to int
			__m256i comparisonMask = _mm256_castps_si256(comparisonResult);
			// printM256i(comparisonMask);
			// applying neighbour to temp vector
			__m256i addMask = _mm256_or_si256(neigh, bitma);
			// registering donor where slope < 0
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// printM256i_bitset(neigh);
			// comparisonMask = _mm256_not_si256(comparisonMask); // need to avoid flats becominglinks
			// Switching to 
			bitma = _mm256_srli_epi32(bitma,8);
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_GT_OQ);
			// printM256(comparisonResult);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			// printM256i(comparisonMask);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// printM256i_bitset(neigh);

			bitma = _mm256_srli_epi32(bitma,16);
			comparisonResult = _mm256_cmp_ps(tSS, SSs, _CMP_GT_OQ);
			SSs = _mm256_blendv_ps(SSs, tSS, comparisonResult);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			SSsi = blendvps_si256(SSsi, zeroVeci,comparisonMask);
			SSsi = blendvps_si256(SSsi, bitma,comparisonMask);


			// Loading top-left infos
			topon = this->_load_t_topo(i);
			bitma = _mm256_set1_epi32(DonTopMask);
			//## Calculating Slope
			inter = _mm256_sub_ps(topos, topon);
			tSS = _mm256_div_ps(inter, DYs);
			//## mask comparing to 0 for normal recs and dons
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_LT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonMask = _mm256_not_si256(comparisonMask); // need to avoid flats becominglinks
			bitma = _mm256_srli_epi32(bitma,8);
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_GT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_srli_epi32(bitma,16);
			comparisonResult = _mm256_cmp_ps(tSS, SSs, _CMP_GT_OQ);
			SSs = _mm256_blendv_ps(SSs, tSS, comparisonResult);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			SSsi = blendvps_si256(SSsi, zeroVeci,comparisonMask);
			SSsi = blendvps_si256(SSsi, bitma,comparisonMask);


			// Loading top-left infos
			topon = this->_load_tr_topo(i);
			bitma = _mm256_set1_epi32(DonTopRightMask);
			//## Calculating Slope
			inter = _mm256_sub_ps(topos, topon);
			tSS = _mm256_div_ps(inter, DXYs);
			//## mask comparing to 0 for normal recs and dons
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_LT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonMask = _mm256_not_si256(comparisonMask); // need to avoid flats becominglinks
			bitma = _mm256_srli_epi32(bitma,8);
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_GT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_srli_epi32(bitma,16);
			comparisonResult = _mm256_cmp_ps(tSS, SSs, _CMP_GT_OQ);
			SSs = _mm256_blendv_ps(SSs, tSS, comparisonResult);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			SSsi = blendvps_si256(SSsi, zeroVeci,comparisonMask);
			SSsi = blendvps_si256(SSsi, bitma,comparisonMask);


			// Loading top-left infos
			topon = this->_load_l_topo(i);
			bitma = _mm256_set1_epi32(DonLeftMask);
			//## Calculating Slope
			inter = _mm256_sub_ps(topos, topon);
			tSS = _mm256_div_ps(inter, DXs);
			//## mask comparing to 0 for normal recs and dons
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_LT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonMask = _mm256_not_si256(comparisonMask); // need to avoid flats becominglinks
			bitma = _mm256_srli_epi32(bitma,8);
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_GT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_srli_epi32(bitma,16);
			comparisonResult = _mm256_cmp_ps(tSS, SSs, _CMP_GT_OQ);
			SSs = _mm256_blendv_ps(SSs, tSS, comparisonResult);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			SSsi = blendvps_si256(SSsi, zeroVeci,comparisonMask);
			SSsi = blendvps_si256(SSsi, bitma,comparisonMask);


			// Loading top-left infos
			topon = this->_load_r_topo(i);
			bitma = _mm256_set1_epi32(DonRightMask);
			//## Calculating Slope
			inter = _mm256_sub_ps(topos, topon);
			tSS = _mm256_div_ps(inter, DXs);
			//## mask comparing to 0 for normal recs and dons
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_LT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonMask = _mm256_not_si256(comparisonMask); // need to avoid flats becominglinks
			bitma = _mm256_srli_epi32(bitma,8);
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_GT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_srli_epi32(bitma,16);
			comparisonResult = _mm256_cmp_ps(tSS, SSs, _CMP_GT_OQ);
			SSs = _mm256_blendv_ps(SSs, tSS, comparisonResult);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			SSsi = blendvps_si256(SSsi, zeroVeci,comparisonMask);
			SSsi = blendvps_si256(SSsi, bitma,comparisonMask);

			// Loading top-left infos
			topon = this->_load_bl_topo(i);
			bitma = _mm256_set1_epi32(DonBottomLeftMask);
			//## Calculating Slope
			inter = _mm256_sub_ps(topos, topon);
			tSS = _mm256_div_ps(inter, DXYs);
			//## mask comparing to 0 for normal recs and dons
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_LT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonMask = _mm256_not_si256(comparisonMask); // need to avoid flats becominglinks
			bitma = _mm256_srli_epi32(bitma,8);
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_GT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_srli_epi32(bitma,16);
			comparisonResult = _mm256_cmp_ps(tSS, SSs, _CMP_GT_OQ);
			SSs = _mm256_blendv_ps(SSs, tSS, comparisonResult);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			SSsi = blendvps_si256(SSsi, zeroVeci,comparisonMask);
			SSsi = blendvps_si256(SSsi, bitma,comparisonMask);


			// Loading top-left infos
			topon = this->_load_b_topo(i);
			bitma = _mm256_set1_epi32(DonBottomMask);
			//## Calculating Slope
			inter = _mm256_sub_ps(topos, topon);
			tSS = _mm256_div_ps(inter, DYs);
			//## mask comparing to 0 for normal recs and dons
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_LT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonMask = _mm256_not_si256(comparisonMask); // need to avoid flats becominglinks
			bitma = _mm256_srli_epi32(bitma,8);
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_GT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_srli_epi32(bitma,16);
			comparisonResult = _mm256_cmp_ps(tSS, SSs, _CMP_GT_OQ);
			SSs = _mm256_blendv_ps(SSs, tSS, comparisonResult);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			SSsi = blendvps_si256(SSsi, zeroVeci,comparisonMask);
			SSsi = blendvps_si256(SSsi, bitma,comparisonMask);


			// Loading top-left infos
			topon = this->_load_br_topo(i);
			bitma = _mm256_set1_epi32(DonBottomRightMask);
			//## Calculating Slope
			inter = _mm256_sub_ps(topos, topon);
			tSS = _mm256_div_ps(inter, DXYs);
			//## mask comparing to 0 for normal recs and dons
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_LT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonMask = _mm256_not_si256(comparisonMask); // need to avoid flats becominglinks
			bitma = _mm256_srli_epi32(bitma,8);
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_GT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_srli_epi32(bitma,16);
			comparisonResult = _mm256_cmp_ps(tSS, SSs, _CMP_GT_OQ);
			SSs = _mm256_blendv_ps(SSs, tSS, comparisonResult);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			SSsi = blendvps_si256(SSsi, zeroVeci,comparisonMask);
			SSsi = blendvps_si256(SSsi, bitma,comparisonMask);
			// std::cout << "~" << std::endl;
			// printM256i_bitset(neigh);


			neigh = _mm256_or_si256(neigh, SSsi);
			// printM256i_bitset(neigh);


			int res[8];
			_mm256_storeu_si256((__m256i*)res, neigh);

			int jj = 0;
			for(int j = 7; j>=0; --j)
			{
				// std::cout << res[j] << "|";
				this->Neighbours[i+jj] = res[j];
				// if(i > this->nx * 20 && this->nx * 50 > i )
				// {
				// 	std::cout << "~" << std::endl;
					// std::cout << std::bitset<32>(res[j]) << std::endl;	
				// }
				
				++jj;
			}

		}


		for (int i = this->nx + 1; i < this->nxy - this->nx - 1; i+=8)
		{
			// storing initial topo and initialising neighbours to no neighbours
			__m256i neigh = this->_load_neighbours(i);
			// printM256i_bitset_8(neigh,0);
			// Loading topography from the top left corner
			__m256i tneigh = this->_load_neighbours(i + this->adder[0]);

			// setting up the bitmask to donor-toward-top-left
			__m256i bitmatg = _mm256_set1_epi32(SRecBottomRightMask);
			__m256i bitmata = _mm256_set1_epi32(SDonTopLeftMask);
			__m256i comparisonResult = _mm256_and_si256(tneigh, bitmatg);
			__m256i comparisonMask = _mm256_cmpeq_epi32(comparisonResult, bitmatg); //_mm256_and_si256(tneigh, bitmatg);
			__m256i addMask = _mm256_or_si256(neigh, bitmata);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			 // neigh = _mm256_or_si256(neigh, comparisonResult);
			// printM256i_bitset_8(neigh,0);
			// // printM256i_bitset_8(neigh,1);
			// printM256i_bitset(neigh);
    	// __m256i comparisonMask = _mm256_set1_epi32(_mm256_movemask_epi8(comparisonResult));
			// // // applying neighbour to temp vector
			// __m256i addMask = _mm256_or_si256(tneigh, bitmata);
			// // // registering donor where slope < 0
			// // neigh = blendvps_si256(neigh, addMask,comparisonMask);


			// Loading topography from the top left corner
			tneigh = this->_load_neighbours(i + this->adder[1]);
			// setting up the bitmask to donor-toward-top-left
			bitmatg = _mm256_set1_epi32(SRecBottomMask);
			bitmata = _mm256_set1_epi32(SDonTopMask);
			comparisonResult = _mm256_and_si256(tneigh, bitmatg);
			comparisonMask = _mm256_cmpeq_epi32(comparisonResult, bitmatg); //_mm256_and_si256(tneigh, bitmatg);
			addMask = _mm256_or_si256(neigh, bitmata);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonResult = _mm256_and_si256(tneigh, bitmatg);
    	// comparisonMask = _mm256_set1_epi32(_mm256_movemask_epi8(comparisonResult));
			// // applying neighbour to temp vector
			// addMask = _mm256_or_si256(tneigh, bitmata);
			// // registering donor where slope < 0
			// neigh = blendvps_si256(neigh, addMask,comparisonMask);

			// Loading topography from the top left corner
			tneigh = this->_load_neighbours(i + this->adder[2]);
			// setting up the bitmask to donor-toward-top-left
			bitmatg = _mm256_set1_epi32(SRecBottomLeftMask);
			bitmata = _mm256_set1_epi32(SDonTopRightMask);
			comparisonResult = _mm256_and_si256(tneigh, bitmatg);
			comparisonMask = _mm256_cmpeq_epi32(comparisonResult, bitmatg); //_mm256_and_si256(tneigh, bitmatg);
			addMask = _mm256_or_si256(neigh, bitmata);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonResult = _mm256_and_si256(tneigh, bitmatg);
    	// comparisonMask = _mm256_set1_epi32(_mm256_movemask_epi8(comparisonResult));
			// // applying neighbour to temp vector
			// addMask = _mm256_or_si256(tneigh, bitmata);
			// // registering donor where slope < 0
			// neigh = blendvps_si256(neigh, addMask,comparisonMask);

			// Loading topography from the top left corner
			tneigh = this->_load_neighbours(i + this->adder[3]);
			// setting up the bitmask to donor-toward-top-left
			bitmatg = _mm256_set1_epi32(SRecRightMask);
			bitmata = _mm256_set1_epi32(SDonLeftMask);
			comparisonResult = _mm256_and_si256(tneigh, bitmatg);
			comparisonMask = _mm256_cmpeq_epi32(comparisonResult, bitmatg); //_mm256_and_si256(tneigh, bitmatg);
			addMask = _mm256_or_si256(neigh, bitmata);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);

			// comparisonResult = _mm256_and_si256(tneigh, bitmatg);
    	// comparisonMask = _mm256_set1_epi32(_mm256_movemask_epi8(comparisonResult));
			// // applying neighbour to temp vector
			// addMask = _mm256_or_si256(tneigh, bitmata);
			// // registering donor where slope < 0
			// neigh = blendvps_si256(neigh, addMask,comparisonMask);

			// Loading topography from the top left corner
			tneigh = this->_load_neighbours(i + this->adder[4]);
			// setting up the bitmask to donor-toward-top-left
			bitmatg = _mm256_set1_epi32(SRecLeftMask);
			bitmata = _mm256_set1_epi32(SDonRightMask);
			comparisonResult = _mm256_and_si256(tneigh, bitmatg);
			comparisonMask = _mm256_cmpeq_epi32(comparisonResult, bitmatg); //_mm256_and_si256(tneigh, bitmatg);
			addMask = _mm256_or_si256(neigh, bitmata);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonResult = _mm256_and_si256(tneigh, bitmatg);
    	// comparisonMask = _mm256_set1_epi32(_mm256_movemask_epi8(comparisonResult));
			// // applying neighbour to temp vector
			// addMask = _mm256_or_si256(tneigh, bitmata);
			// // registering donor where slope < 0
			// neigh = blendvps_si256(neigh, addMask,comparisonMask);

			// Loading topography from the top left corner
			tneigh = this->_load_neighbours(i + this->adder[5]);
			// setting up the bitmask to donor-toward-top-left
			bitmatg = _mm256_set1_epi32(SRecTopRightMask);
			bitmata = _mm256_set1_epi32(SDonBottomLeftMask);
			comparisonResult = _mm256_and_si256(tneigh, bitmatg);
			comparisonMask = _mm256_cmpeq_epi32(comparisonResult, bitmatg); //_mm256_and_si256(tneigh, bitmatg);
			addMask = _mm256_or_si256(neigh, bitmata);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonResult = _mm256_and_si256(tneigh, bitmatg);
    	// comparisonMask = _mm256_set1_epi32(_mm256_movemask_epi8(comparisonResult));
			// // applying neighbour to temp vector
			// addMask = _mm256_or_si256(tneigh, bitmata);
			// // registering donor where slope < 0
			// neigh = blendvps_si256(neigh, addMask,comparisonMask);

			// Loading topography from the top left corner
			tneigh = this->_load_neighbours(i + this->adder[6]);
			// setting up the bitmask to donor-toward-top-left
			bitmatg = _mm256_set1_epi32(SRecTopMask);
			bitmata = _mm256_set1_epi32(SDonBottomMask);
			comparisonResult = _mm256_and_si256(tneigh, bitmatg);
			comparisonMask = _mm256_cmpeq_epi32(comparisonResult, bitmatg); //_mm256_and_si256(tneigh, bitmatg);
			addMask = _mm256_or_si256(neigh, bitmata);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonResult = _mm256_and_si256(tneigh, bitmatg);
    	// comparisonMask = _mm256_set1_epi32(_mm256_movemask_epi8(comparisonResult));
			// // applying neighbour to temp vector
			// addMask = _mm256_or_si256(tneigh, bitmata);
			// // registering donor where slope < 0
			// neigh = blendvps_si256(neigh, addMask,comparisonMask);

			// Loading topography from the top left corner
			tneigh = this->_load_neighbours(i + this->adder[7]);
			// setting up the bitmask to donor-toward-top-left
			bitmatg = _mm256_set1_epi32(SRecTopLeftMask);
			bitmata = _mm256_set1_epi32(SDonBottomRightMask);
			comparisonResult = _mm256_and_si256(tneigh, bitmatg);
			comparisonMask = _mm256_cmpeq_epi32(comparisonResult, bitmatg); //_mm256_and_si256(tneigh, bitmatg);
			addMask = _mm256_or_si256(neigh, bitmata);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonResult = _mm256_and_si256(tneigh, bitmatg);
    	// comparisonMask = _mm256_set1_epi32(_mm256_movemask_epi8(comparisonResult));
			// // applying neighbour to temp vector
			// addMask = _mm256_or_si256(tneigh, bitmata);
			// // registering donor where slope < 0
			// neigh = blendvps_si256(neigh, addMask,comparisonMask);

			// printM256i_bitset(neigh);
			 // printM256i_bitset_8(neigh,1);
			// printM256i_bitset_8(neigh,1);




			int res[8];
			_mm256_storeu_si256((__m256i*)res, neigh);
			// printM256i_bitset(neigh);

			int jj = 0;
			for(int j = 7; j>=0; --j)
			{
				// std::cout << res[j] << "|";
				this->Neighbours[i+jj] = res[j];

				// if(i > this->nx * 20 && this->nx * 50 > i )
				// {
				// 	std::cout << "~" << std::endl;
				// if(std::bitset<8>(ui8fromui32(res[j],1))[0])
				// 	std::cout << std::bitset<8>(ui8fromui32(res[j],1)) << std::endl;;
					// exit(0);	
				// }
				++jj;
			}
		}


		this->_compute_edge_nodes();

	}





	void compute()
	{
		if(this->_topography == nullptr) throw std::runtime_error("no topo linked sns");

		std::array<int,8> neighbours;
		int dbgcpt = 0;

		__m256 DXs = fulldx();
		__m256 DYs = fulldy();

		__m256 DXYs = fulldxy();
		std::array<__m256,8> mdxer{DXYs,DYs,DXYs,DXs,DXs,DXYs,DYs,DXYs};
		
		__m256 zeroVec = _mm256_setzero_ps();
		__m256i zeroVeci =  _mm256_set1_epi32(0);
		__m256i oneVeci =  _mm256_set1_epi32(-1);

		// auto DonMasks = DonMasks();
		// auto SDonMasks = SDonMasks();
		// auto SRecMasks = SRecMasks();
		// auto MirrorIndicesMasks = MirrorIndicesMasks();


		__m256i neigh;
		__m256i tneigh;
		__m256i bitmatg;
		__m256i bitmata;
		__m256i comparisonResulti;
		__m256i comparisonMask;
		__m256i addMask;
		__m256 topos;
		__m256 SSs;
		__m256i SSsi;
		__m256 topon;
		__m256i bitma;
		__m256 inter;
		__m256 tSS;
		__m256 comparisonResult;


		for (size_t i = 0; i < this->IntaNodes.size(); i+=8)
		{
			// std::cout << i+7 << "/" << this->IntaNodes.size() << std::endl;
			// storing initial topo and initialising neighbours to no neighbours
			topos = this->_load_topo(i, this->IntaNodes);
			neigh = _mm256_set1_epi32(NoNeigh);
			// Tracking steepest descent (init to 0.)
			SSs = _mm256_set1_ps(0.f);
			SSsi = _mm256_set1_epi32(0);

			#pragma GCC unroll(8)
			for(size_t j = 0; j<8; ++j)
			{
				// Loading topography from the top left corner
				topon = this->_load_topo(i,this->IntaNodes,this->adder[j]);
				// setting up the bitmask to donor-toward-top-left
				bitma = _mm256_set1_epi32(DonMasks[j]);
				//## Calculating Slope
				// subtracting current topo with neighbours'
				inter = _mm256_sub_ps(topos, topon);
				// Dividing diff by diagonals
				tSS = _mm256_div_ps(inter, mdxer[j]);
				// printM256(tSS);
				//## mask comparing to 0 for normal recs and dons
				// mask where slope < 0 (donors)
				comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_LT_OQ);
				// recasting mask to int
				comparisonMask = _mm256_castps_si256(comparisonResult);
				// printM256i(comparisonMask);
				// applying neighbour to temp vector
				addMask = _mm256_or_si256(neigh, bitma);
				// registering donor where slope < 0
				neigh = blendvps_si256(neigh, addMask,comparisonMask);
				// printM256i_bitset(neigh);
				// comparisonMask = _mm256_not_si256(comparisonMask); // need to avoid flats becominglinks
				// Switching to 
				bitma = _mm256_srli_epi32(bitma,8);
				comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_GT_OQ);
				// printM256(comparisonResult);
				comparisonMask = _mm256_castps_si256(comparisonResult);
				// printM256i(comparisonMask);
				addMask = _mm256_or_si256(neigh, bitma);
				neigh = blendvps_si256(neigh, addMask,comparisonMask);
				// printM256i_bitset(neigh);

				bitma = _mm256_srli_epi32(bitma,16);
				comparisonResult = _mm256_cmp_ps(tSS, SSs, _CMP_GT_OQ);
				SSs = _mm256_blendv_ps(SSs, tSS, comparisonResult);
				comparisonMask = _mm256_castps_si256(comparisonResult);
				// SSsi = blendvps_si256(SSsi, zeroVeci,comparisonMask);
				SSsi = blendvps_si256(SSsi, bitma,comparisonMask);
			}

			neigh = _mm256_or_si256(neigh, SSsi);
			// printM256i_bitset_8(neigh,0);

			int res[8];
			_mm256_storeu_si256((__m256i*)res, neigh);

			int jj = 0;
			#pragma GCC unroll(8)
			for(int j = 7; j>=0; --j)
			{
				this->Neighbours[this->IntaNodes[i+jj]] = res[j];
				// std::cout << std::bitset<8>(ui8fromui32(this->Neighbours[this->IntaNodes[i+jj]],0)) << "|";
				++jj;
			}

		}

		for (size_t i = 0; i < this->IntaNodes.size(); i+=8)
		{

			neigh = this->_load_neighbours(i,this->IntaNodes);
			// printM256i_bitset_8(neigh,1);
			#pragma GCC unroll(8)
			for(size_t j=0; j < 8; ++j)
			{
				// storing initial topo and initialising neighbours to no neighbours
				// printM256i_bitset_8(neigh,0);
				// Loading topography from the top left corner
				tneigh = this->_load_neighbours(i, this->IntaNodes, this->adder[j]);;

				// setting up the bitmask to donor-toward-top-left
				bitmatg = _mm256_set1_epi32(SRecMasks[MirrorIndicesMasks[j]]);
				bitmata = _mm256_set1_epi32(SDonMasks[j]);
				comparisonResulti = _mm256_and_si256(tneigh, bitmatg);
				comparisonMask = _mm256_cmpeq_epi32(comparisonResulti, bitmatg); //_mm256_and_si256(tneigh, bitmatg);
				// printM256i_bitset(comparisonMask);
				// addMask = _mm256_or_si256(neigh, bitmata);
				// neigh = blendvps_si256(neigh, addMask,comparisonMask);
				addMask = _mm256_and_si256(comparisonMask, bitmata);
				// printM256i_bitset(addMask);
				// printM256i_bitset_8(neigh,0);
				neigh = _mm256_or_si256(addMask, neigh);
				// printM256i_bitset_8(neigh,0);
			}
			// printM256i_bitset_8(neigh,1);
			// std::cout << std::endl;
			
			int res[8];
			_mm256_storeu_si256((__m256i*)res, neigh);
			// printM256i_bitset(neigh);
			// printM256i_bitset(neigh);

			int jj = 0;
			#pragma GCC unroll(8)
			for(int j = 7; j>=0; --j)
			{
				this->Neighbours[this->IntaNodes[i + jj]] = res[j];
				++jj;
				// if(ui8fromui32(res[j],1) != 0)
				// 	++dbgcpt;
			}
		}

		this->_compute_edge_nodes();

	}

	void compute_beta1()
	{
		if(this->_topography == nullptr) throw std::runtime_error("no topo linked sns");

		std::array<int,8> neighbours;
		int dbgcpt = 0;

		__m256 DXs = fulldx();
		__m256 DYs = fulldy();

		__m256 DXYs = fulldxy();
		std::array<__m256,8> mdxer{DXYs,DYs,DXYs,DXs,DXs,DXYs,DYs,DXYs};
		
		__m256 zeroVec = _mm256_setzero_ps();
		__m256i zeroVeci =  _mm256_set1_epi32(0);
		__m256i oneVeci =  _mm256_set1_epi32(-1);

		// auto DonMasks = DonMasks();
		// auto SDonMasks = SDonMasks();
		// auto SRecMasks = SRecMasks();
		// auto MirrorIndicesMasks = MirrorIndicesMasks();


		__m256i neigh;
		__m256i tneigh;
		__m256i gatherIndices;
		__m256i gatherIndicesN;
		__m256i adderavx;
		__m256i bitmatg;
		__m256i bitmata;
		__m256i comparisonResulti;
		__m256i comparisonMask;
		__m256i addMask;
		__m256 topos;
		__m256 SSs;
		__m256i SSsi;
		__m256 topon;
		__m256i bitma;
		__m256 inter;
		__m256 tSS;
		__m256 comparisonResult;

		const int* indicesPtr = this->IntaNodes.data();
		const float* dataPtr = this->_topography->data();

		for (size_t i = 0; i < this->IntaNodes.size(); i+=8)
		{
			// std::cout << i+7 << "/" << this->IntaNodes.size() << std::endl;
			// storing initial topo and initialising neighbours to no neighbours
			// topos = this->_load_topo(i, this->IntaNodes);


			gatherIndices = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(indicesPtr) + static_cast<int>(i/8));
			topos = _mm256_i32gather_ps(dataPtr, gatherIndices, sizeof(float));

			neigh = _mm256_set1_epi32(this->NoNeigh);
			// Tracking steepest descent (init to 0.)
			SSs = _mm256_set1_ps(0.f);
			SSsi = _mm256_set1_epi32(0);

			#pragma GCC unroll(8)
			for(size_t j = 0; j<8; ++j)
			{
				// Loading topography from the top left corner
				// topon = this->_load_topo(i,this->IntaNodes,this->adder[j]);
				adderavx = _mm256_set1_epi32(this->adder[j]);
				gatherIndicesN = _mm256_add_epi32(gatherIndices,adderavx);
				topon = _mm256_i32gather_ps(dataPtr, gatherIndicesN, sizeof(float));

				// setting up the bitmask to donor-toward-top-left
				bitma = _mm256_set1_epi32(DonMasks[j]);
				//## Calculating Slope
				// subtracting current topo with neighbours'
				inter = _mm256_sub_ps(topos, topon);
				// Dividing diff by diagonals
				tSS = _mm256_div_ps(inter, mdxer[j]);
				// printM256(tSS);
				//## mask comparing to 0 for normal recs and dons
				// mask where slope < 0 (donors)
				comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_LT_OQ);
				// recasting mask to int
				comparisonMask = _mm256_castps_si256(comparisonResult);
				// printM256i(comparisonMask);
				// applying neighbour to temp vector
				addMask = _mm256_or_si256(neigh, bitma);
				// registering donor where slope < 0
				neigh = blendvps_si256(neigh, addMask,comparisonMask);
				// printM256i_bitset(neigh);
				// comparisonMask = _mm256_not_si256(comparisonMask); // need to avoid flats becominglinks
				// Switching to 
				bitma = _mm256_srli_epi32(bitma,8);
				comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_GT_OQ);
				// printM256(comparisonResult);
				comparisonMask = _mm256_castps_si256(comparisonResult);
				// printM256i(comparisonMask);
				addMask = _mm256_or_si256(neigh, bitma);
				neigh = blendvps_si256(neigh, addMask,comparisonMask);
				// printM256i_bitset(neigh);

				bitma = _mm256_srli_epi32(bitma,16);
				comparisonResult = _mm256_cmp_ps(tSS, SSs, _CMP_GT_OQ);
				SSs = _mm256_blendv_ps(SSs, tSS, comparisonResult);
				comparisonMask = _mm256_castps_si256(comparisonResult);
				// SSsi = blendvps_si256(SSsi, zeroVeci,comparisonMask);
				SSsi = blendvps_si256(SSsi, bitma,comparisonMask);
			}

			neigh = _mm256_or_si256(neigh, SSsi);
			// printM256i_bitset_8(neigh,0);

			int res[8];
			_mm256_storeu_si256((__m256i*)res, neigh);

			int jj = 0;
			#pragma GCC unroll(8)
			for(int j = 7; j>=0; --j)
			{
				this->Neighbours[this->IntaNodes[i+jj]] = res[j];
				// std::cout << std::bitset<8>(ui8fromui32(this->Neighbours[this->IntaNodes[i+jj]],0)) << "|";
				++jj;
			}

		}

		for (size_t i = 0; i < this->IntaNodes.size(); i+=8)
		{

			neigh = this->_load_neighbours(i,this->IntaNodes);
			// printM256i_bitset_8(neigh,1);
			#pragma GCC unroll(8)
			for(size_t j=0; j < 8; ++j)
			{
				// storing initial topo and initialising neighbours to no neighbours
				// printM256i_bitset_8(neigh,0);
				// Loading topography from the top left corner
				tneigh = this->_load_neighbours(i, this->IntaNodes, this->adder[j]);;

				// setting up the bitmask to donor-toward-top-left
				bitmatg = _mm256_set1_epi32(SRecMasks[MirrorIndicesMasks[j]]);
				bitmata = _mm256_set1_epi32(SDonMasks[j]);
				comparisonResulti = _mm256_and_si256(tneigh, bitmatg);
				comparisonMask = _mm256_cmpeq_epi32(comparisonResulti, bitmatg); //_mm256_and_si256(tneigh, bitmatg);
				// printM256i_bitset(comparisonMask);
				// addMask = _mm256_or_si256(neigh, bitmata);
				// neigh = blendvps_si256(neigh, addMask,comparisonMask);
				addMask = _mm256_and_si256(comparisonMask, bitmata);
				// printM256i_bitset(addMask);
				// printM256i_bitset_8(neigh,0);
				neigh = _mm256_or_si256(addMask, neigh);
				// printM256i_bitset_8(neigh,0);
			}
			// printM256i_bitset_8(neigh,1);
			// std::cout << std::endl;
			
			int res[8];
			_mm256_storeu_si256((__m256i*)res, neigh);
			// printM256i_bitset(neigh);
			// printM256i_bitset(neigh);

			int jj = 0;
			#pragma GCC unroll(8)
			for(int j = 7; j>=0; --j)
			{
				this->Neighbours[this->IntaNodes[i + jj]] = res[j];
				++jj;
				// if(ui8fromui32(res[j],1) != 0)
				// 	++dbgcpt;
			}
		}

		this->_compute_edge_nodes();

	}


	void _compute_edge_nodes()
	{
		std::array<int,8> arr;
		// auto DonMasks = DonMasks() ;
		// auto RecMasks = RecMasks() ;
		// auto SDonMasks = SDonMasks() ;
		// auto SRecMasks = SRecMasks() ;
		// auto MirrorIndicesMasks = MirrorIndicesMasks() ;

		for(size_t ti=0; ti<this->EdgeNodes.size(); ++ti)
		{
			// if(this->boundaries[i])
			// this->internal_neighbours(i,arr);
			int i = this->EdgeNodes[ti];
			this->_neighbours_from_ui8_nocheck(i,this->EdgeNodes_Neighbours[ti],arr);
			
			float topo = this->topography(i);

			for(int j=0; j< 8; ++j)
			{
				int tn = arr[j];
				if(tn > -1)
				{
					int tnn= this->Neighbours[tn];
					if(tnn & DonMasks[MirrorIndicesMasks[j]])
						this->Neighbours[i] |= RecMasks[j];
					else if(tnn & RecMasks[MirrorIndicesMasks[j]])
						this->Neighbours[i] |= DonMasks[j];


					if(tnn & SDonMasks[MirrorIndicesMasks[j]])
						this->Neighbours[i] |= SRecMasks[j];
					else if(tnn & SRecMasks[MirrorIndicesMasks[j]])
						this->Neighbours[i] |= SDonMasks[j];
					
					
				}
				
			}
		}
			
			
	}

	void compute_classics()
	{
		
		std::array<int,8> arr;
		// auto DonMasks = DonMasks() ;
		// auto RecMasks = RecMasks() ;
		// auto SDonMasks = SDonMasks() ;
		// auto SRecMasks = SRecMasks() ;
		// auto MirrorIndicesMasks = MirrorIndicesMasks() ;
		for (int i =0; i < this->nxy ; ++i)
		{
			this->Neighbours[i] = 0;
		}

		for (int i =0; i < this->nxy ; ++i)
		{
			if(this->boundaries[i] == 1)
				this->internal_neighbours(i,arr);
			else
				continue;

			float topo = this->topography(i);
			float tSS = 0;
			int tmi = 0;
			int tSSi = i;
			for(int j=0; j< 8; ++j)
			{
				int tn = arr[j];
				if(tn > -1)
				{
					float ntopo = this->topography(tn);
					if(ntopo < topo)
					{
						this->Neighbours[i] |= RecMasks[j];
						this->Neighbours[tn] |= DonMasks[MirrorIndicesMasks[j]];
						float ttSS = (topo - ntopo)/this->dxer[j];
						if(ttSS > tSS)
						{
							tSS = ttSS;
							tmi = j ;
							tSSi = tn;
						}
					}
					
				}
				
			}
			
			if(tSSi != i)
			{
				this->Neighbours[i] |= SRecMasks[tmi];
				this->Neighbours[tSSi] |= SDonMasks[MirrorIndicesMasks[tmi]];
				// std::cout << std::bitset<32>(this->Neighbours[i])  << "|" << std::bitset<32>(this->Neighbours[tSSi]) << std::endl;
				// std::cout << std::bitset<8>(ui8fromui32(this->Neighbours[i],0))  << "|" << std::bitset<8>(ui8fromui32(this->Neighbours[tSSi],1)) << std::endl;
			}


		}
		this->_compute_edge_nodes();


	}



	void compute_archive_1()
	{
		if(this->_topography == nullptr) throw std::runtime_error("no topo linked sns");

		std::array<int,8> neighbours;

		__m256 DXs = fulldx();
		__m256 DYs = fulldy();
		__m256 DXYs = fulldxy();
		__m256 zeroVec = _mm256_setzero_ps();


		for (int i = this->nx + 1; i < this->nxy - this->nx - 1; i+=8)
		{
			// storing initial topo and initialising neighbours to no neighbours
			__m256 topos = this->_load_topo(i);
			__m256i neigh = _mm256_set1_epi32(this->NoNeigh);
			__m256i orSdon = _mm256_set1_epi32(0);
			// Tracking steepest descent (init to 0.)
			__m256 SSs = _mm256_set1_ps(0.f);


			// Loading topography from the top left corner
			__m256 topon = this->_load_tl_topo(i);
			// setting up the bitmask to donor-toward-top-left
			__m256i bitma = _mm256_set1_epi32(DonTopLeftMask);
			//## Calculating Slope
			// subtracting current topo with neighbours'
			__m256 inter = _mm256_sub_ps(topos, topon);
			// Dividing diff by diagonals
			__m256 tSS = _mm256_div_ps(inter, DXYs);
			//## mask comparing to 0 for normal recs and dons
			// mask where slope < 0 (donors)
			__m256 comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_LT_OQ);
			// recasting mask to int
			__m256i comparisonMask = _mm256_castps_si256(comparisonResult);
			// applying neighbour to temp vector
			__m256i addMask = _mm256_or_si256(neigh, bitma);
			// registering donor where slope < 0
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonMask = _mm256_not_si256(comparisonMask); // need to avoid flats becominglinks
			// Switching to 
			bitma = _mm256_srli_epi32(bitma,8);
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_GT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_srli_epi32(bitma,16);
			comparisonResult = _mm256_cmp_ps(tSS, SSs, _CMP_GT_OQ);
			SSs = _mm256_blendv_ps(SSs, tSS, comparisonResult);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_slli_epi32(bitma,1);
			addMask = _mm256_or_si256(orSdon, bitma);
			orSdon = blendvps_si256(orSdon, addMask,comparisonMask);


			// Loading top-left infos
			topon = this->_load_t_topo(i);
			bitma = _mm256_set1_epi32(DonTopMask);
			//## Calculating Slope
			inter = _mm256_sub_ps(topos, topon);
			tSS = _mm256_div_ps(inter, DXYs);
			//## mask comparing to 0 for normal recs and dons
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_LT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonMask = _mm256_not_si256(comparisonMask); // need to avoid flats becominglinks
			bitma = _mm256_srli_epi32(bitma,8);
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_GT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_srli_epi32(bitma,16);
			comparisonResult = _mm256_cmp_ps(tSS, SSs, _CMP_GT_OQ);
			SSs = _mm256_blendv_ps(SSs, tSS, comparisonResult);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_slli_epi32(bitma,3);
			addMask = _mm256_or_si256(orSdon, bitma);
			orSdon = blendvps_si256(orSdon, addMask,comparisonMask);

			// Loading top-left infos
			topon = this->_load_tr_topo(i);
			bitma = _mm256_set1_epi32(DonTopRightMask);
			//## Calculating Slope
			inter = _mm256_sub_ps(topos, topon);
			tSS = _mm256_div_ps(inter, DXYs);
			//## mask comparing to 0 for normal recs and dons
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_LT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonMask = _mm256_not_si256(comparisonMask); // need to avoid flats becominglinks
			bitma = _mm256_srli_epi32(bitma,8);
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_GT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_srli_epi32(bitma,16);
			comparisonResult = _mm256_cmp_ps(tSS, SSs, _CMP_GT_OQ);
			SSs = _mm256_blendv_ps(SSs, tSS, comparisonResult);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_slli_epi32(bitma,5);
			addMask = _mm256_or_si256(orSdon, bitma);
			orSdon = blendvps_si256(orSdon, addMask,comparisonMask);

			// Loading top-left infos
			topon = this->_load_l_topo(i);
			bitma = _mm256_set1_epi32(this->DonLeftMask);
			//## Calculating Slope
			inter = _mm256_sub_ps(topos, topon);
			tSS = _mm256_div_ps(inter, DXYs);
			//## mask comparing to 0 for normal recs and dons
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_LT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonMask = _mm256_not_si256(comparisonMask); // need to avoid flats becominglinks
			bitma = _mm256_srli_epi32(bitma,8);
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_GT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_srli_epi32(bitma,16);
			comparisonResult = _mm256_cmp_ps(tSS, SSs, _CMP_GT_OQ);
			SSs = _mm256_blendv_ps(SSs, tSS, comparisonResult);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_slli_epi32(bitma,7);
			addMask = _mm256_or_si256(orSdon, bitma);
			orSdon = blendvps_si256(orSdon, addMask,comparisonMask);

			// Loading top-left infos
			topon = this->_load_r_topo(i);
			bitma = _mm256_set1_epi32(this->DonRightMask);
			//## Calculating Slope
			inter = _mm256_sub_ps(topos, topon);
			tSS = _mm256_div_ps(inter, DXYs);
			//## mask comparing to 0 for normal recs and dons
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_LT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonMask = _mm256_not_si256(comparisonMask); // need to avoid flats becominglinks
			bitma = _mm256_srli_epi32(bitma,8);
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_GT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_srli_epi32(bitma,16);
			comparisonResult = _mm256_cmp_ps(tSS, SSs, _CMP_GT_OQ);
			SSs = _mm256_blendv_ps(SSs, tSS, comparisonResult);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_slli_epi32(bitma,9);
			addMask = _mm256_or_si256(orSdon, bitma);
			orSdon = blendvps_si256(orSdon, addMask,comparisonMask);

			// Loading top-left infos
			topon = this->_load_bl_topo(i);
			bitma = _mm256_set1_epi32(this->DonBottomLeftMask);
			//## Calculating Slope
			inter = _mm256_sub_ps(topos, topon);
			tSS = _mm256_div_ps(inter, DXYs);
			//## mask comparing to 0 for normal recs and dons
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_LT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonMask = _mm256_not_si256(comparisonMask); // need to avoid flats becominglinks
			bitma = _mm256_srli_epi32(bitma,8);
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_GT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_srli_epi32(bitma,16);
			comparisonResult = _mm256_cmp_ps(tSS, SSs, _CMP_GT_OQ);
			SSs = _mm256_blendv_ps(SSs, tSS, comparisonResult);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_slli_epi32(bitma,11);
			addMask = _mm256_or_si256(orSdon, bitma);
			orSdon = blendvps_si256(orSdon, addMask,comparisonMask);

			// Loading top-left infos
			topon = this->_load_b_topo(i);
			bitma = _mm256_set1_epi32(this->DonBottomMask);
			//## Calculating Slope
			inter = _mm256_sub_ps(topos, topon);
			tSS = _mm256_div_ps(inter, DXYs);
			//## mask comparing to 0 for normal recs and dons
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_LT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonMask = _mm256_not_si256(comparisonMask); // need to avoid flats becominglinks
			bitma = _mm256_srli_epi32(bitma,8);
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_GT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_srli_epi32(bitma,16);
			comparisonResult = _mm256_cmp_ps(tSS, SSs, _CMP_GT_OQ);
			SSs = _mm256_blendv_ps(SSs, tSS, comparisonResult);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_slli_epi32(bitma,13);
			addMask = _mm256_or_si256(orSdon, bitma);
			orSdon = blendvps_si256(orSdon, addMask,comparisonMask);

			// Loading top-left infos
			topon = this->_load_br_topo(i);
			bitma = _mm256_set1_epi32(this->DonBottomRightMask);
			//## Calculating Slope
			inter = _mm256_sub_ps(topos, topon);
			tSS = _mm256_div_ps(inter, DXYs);
			//## mask comparing to 0 for normal recs and dons
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_LT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonMask = _mm256_not_si256(comparisonMask); // need to avoid flats becominglinks
			bitma = _mm256_srli_epi32(bitma,8);
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_GT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_srli_epi32(bitma,16);
			comparisonResult = _mm256_cmp_ps(tSS, SSs, _CMP_GT_OQ);
			SSs = _mm256_blendv_ps(SSs, tSS, comparisonResult);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_slli_epi32(bitma,15);
			addMask = _mm256_or_si256(orSdon, bitma);
			orSdon = blendvps_si256(orSdon, addMask,comparisonMask);

			int res[8];
			_mm256_storeu_si256((__m256i*)res, neigh);

			int jj = 0;
			for(int j = 7; j>=0; --j)
			{
				// std::cout << res[j] << "|";
				this->Neighbours[i+jj] = res[j];
				++jj;
			}

			int tsdon[8];
			_mm256_storeu_si256((__m256i*)tsdon, neigh);

			jj = 0;
			for(int j = 7; j>=0; --j)
			{
				// std::cout << res[j] << "|";
				int k = this->Sreceivers(i+jj);
				// int k = i;
				if(k != i+jj)
					this->Neighbours[k] = this->Neighbours[k] | tsdon[j];
				++jj;
			}

		}

	}

	void _compute_omp_test()
	{
		if(this->_topography == nullptr) throw std::runtime_error("no topo linked sns");

		std::array<int,8> neighbours;

		__m256 DXs = fulldx();
		__m256 DYs = fulldy();
		__m256 DXYs = fulldxy();
		__m256 zeroVec = _mm256_setzero_ps();

		#pragma omp parallel for num_threads(4)
		for (int i = this->nx + 1; i < this->nxy - this->nx - 1; i+=8)
		{
			// storing initial topo and initialising neighbours to no neighbours
			__m256 topos = this->_load_topo(i);
			__m256i neigh = _mm256_set1_epi32(this->NoNeigh);
			// Tracking steepest descent (init to 0.)
			__m256 SSs = _mm256_set1_ps(0.f);


			// Loading top-left infos
			__m256 topon = this->_load_tl_topo(i);
			__m256i bitma = _mm256_set1_epi32(DonTopLeftMask);
			//## Calculating Slope
			__m256 inter = _mm256_sub_ps(topos, topon);
			__m256 tSS = _mm256_div_ps(inter, DXYs);
			//## mask comparing to 0 for normal recs and dons
			__m256 comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_LT_OQ);
			__m256i comparisonMask = _mm256_castps_si256(comparisonResult);
			__m256i addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonMask = _mm256_not_si256(comparisonMask); // need to avoid flats becominglinks
			bitma = _mm256_srli_epi32(bitma,8);
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_GT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_srli_epi32(bitma,16);
			comparisonResult = _mm256_cmp_ps(tSS, SSs, _CMP_GT_OQ);
			SSs = _mm256_blendv_ps(SSs, tSS, comparisonResult);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			blendvps_si256(neigh, addMask,comparisonMask);


			// Loading top-left infos
			topon = this->_load_t_topo(i);
			bitma = _mm256_set1_epi32(DonTopMask);
			//## Calculating Slope
			inter = _mm256_sub_ps(topos, topon);
			tSS = _mm256_div_ps(inter, DXYs);
			//## mask comparing to 0 for normal recs and dons
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_LT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonMask = _mm256_not_si256(comparisonMask); // need to avoid flats becominglinks
			bitma = _mm256_srli_epi32(bitma,8);
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_GT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_srli_epi32(bitma,16);
			comparisonResult = _mm256_cmp_ps(tSS, SSs, _CMP_GT_OQ);
			SSs = _mm256_blendv_ps(SSs, tSS, comparisonResult);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			blendvps_si256(neigh, addMask,comparisonMask);

			// Loading top-left infos
			topon = this->_load_tr_topo(i);
			bitma = _mm256_set1_epi32(DonTopRightMask);
			//## Calculating Slope
			inter = _mm256_sub_ps(topos, topon);
			tSS = _mm256_div_ps(inter, DXYs);
			//## mask comparing to 0 for normal recs and dons
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_LT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonMask = _mm256_not_si256(comparisonMask); // need to avoid flats becominglinks
			bitma = _mm256_srli_epi32(bitma,8);
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_GT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_srli_epi32(bitma,16);
			comparisonResult = _mm256_cmp_ps(tSS, SSs, _CMP_GT_OQ);
			SSs = _mm256_blendv_ps(SSs, tSS, comparisonResult);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			blendvps_si256(neigh, addMask,comparisonMask);

			// Loading top-left infos
			topon = this->_load_l_topo(i);
			bitma = _mm256_set1_epi32(this->DonLeftMask);
			//## Calculating Slope
			inter = _mm256_sub_ps(topos, topon);
			tSS = _mm256_div_ps(inter, DXYs);
			//## mask comparing to 0 for normal recs and dons
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_LT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonMask = _mm256_not_si256(comparisonMask); // need to avoid flats becominglinks
			bitma = _mm256_srli_epi32(bitma,8);
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_GT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_srli_epi32(bitma,16);
			comparisonResult = _mm256_cmp_ps(tSS, SSs, _CMP_GT_OQ);
			SSs = _mm256_blendv_ps(SSs, tSS, comparisonResult);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			blendvps_si256(neigh, addMask,comparisonMask);

			// Loading top-left infos
			topon = this->_load_r_topo(i);
			bitma = _mm256_set1_epi32(this->DonRightMask);
			//## Calculating Slope
			inter = _mm256_sub_ps(topos, topon);
			tSS = _mm256_div_ps(inter, DXYs);
			//## mask comparing to 0 for normal recs and dons
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_LT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonMask = _mm256_not_si256(comparisonMask); // need to avoid flats becominglinks
			bitma = _mm256_srli_epi32(bitma,8);
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_GT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_srli_epi32(bitma,16);
			comparisonResult = _mm256_cmp_ps(tSS, SSs, _CMP_GT_OQ);
			SSs = _mm256_blendv_ps(SSs, tSS, comparisonResult);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			blendvps_si256(neigh, addMask,comparisonMask);

			// Loading top-left infos
			topon = this->_load_bl_topo(i);
			bitma = _mm256_set1_epi32(this->DonBottomLeftMask);
			//## Calculating Slope
			inter = _mm256_sub_ps(topos, topon);
			tSS = _mm256_div_ps(inter, DXYs);
			//## mask comparing to 0 for normal recs and dons
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_LT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonMask = _mm256_not_si256(comparisonMask); // need to avoid flats becominglinks
			bitma = _mm256_srli_epi32(bitma,8);
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_GT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_srli_epi32(bitma,16);
			comparisonResult = _mm256_cmp_ps(tSS, SSs, _CMP_GT_OQ);
			SSs = _mm256_blendv_ps(SSs, tSS, comparisonResult);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			blendvps_si256(neigh, addMask,comparisonMask);

			// Loading top-left infos
			topon = this->_load_b_topo(i);
			bitma = _mm256_set1_epi32(this->DonBottomMask);
			//## Calculating Slope
			inter = _mm256_sub_ps(topos, topon);
			tSS = _mm256_div_ps(inter, DXYs);
			//## mask comparing to 0 for normal recs and dons
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_LT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonMask = _mm256_not_si256(comparisonMask); // need to avoid flats becominglinks
			bitma = _mm256_srli_epi32(bitma,8);
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_GT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_srli_epi32(bitma,16);
			comparisonResult = _mm256_cmp_ps(tSS, SSs, _CMP_GT_OQ);
			SSs = _mm256_blendv_ps(SSs, tSS, comparisonResult);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			blendvps_si256(neigh, addMask,comparisonMask);

			// Loading top-left infos
			topon = this->_load_br_topo(i);
			bitma = _mm256_set1_epi32(this->DonBottomRightMask);
			//## Calculating Slope
			inter = _mm256_sub_ps(topos, topon);
			tSS = _mm256_div_ps(inter, DXYs);
			//## mask comparing to 0 for normal recs and dons
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_LT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			neigh = blendvps_si256(neigh, addMask,comparisonMask);
			// comparisonMask = _mm256_not_si256(comparisonMask); // need to avoid flats becominglinks
			bitma = _mm256_srli_epi32(bitma,8);
			comparisonResult = _mm256_cmp_ps(tSS, zeroVec, _CMP_GT_OQ);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			blendvps_si256(neigh, addMask,comparisonMask);
			bitma = _mm256_srli_epi32(bitma,16);
			comparisonResult = _mm256_cmp_ps(tSS, SSs, _CMP_GT_OQ);
			SSs = _mm256_blendv_ps(SSs, tSS, comparisonResult);
			comparisonMask = _mm256_castps_si256(comparisonResult);
			addMask = _mm256_or_si256(neigh, bitma);
			blendvps_si256(neigh, addMask,comparisonMask);

			int res[8];
			// _mm256_storeu_epi32(res,nein);
			_mm256_storeu_si256((__m256i*)res, neigh);

			int jj = 0;
			for(int j = 7; j>=0; --j)
			{
				// std::cout << res[j] << "|";
				this->Neighbours[i+jj] = res[j];
				++jj;
			}



			// // this->_save_SSi( neigh,  nein,  SSs,  tSS);
			// topon = this->_load_t_topo(i);
			// nein = this->_load_t_topoi(i);
			// inter = _mm256_sub_ps(topos,topon);
			// tSS = _mm256_div_ps(inter,DYs);
			// // this->_save_SSi( neigh,  nein,  SSs,  tSS);

			// topon = this->_load_tr_topo(i);
			// nein = this->_load_tr_topoi(i);
			// inter = _mm256_sub_ps(topos,topon);
			// tSS = _mm256_div_ps(inter,DXYs);
			// // this->_save_SSi( neigh,  nein,  SSs,  tSS);

			// topon = this->_load_l_topo(i);
			// nein = this->_load_l_topoi(i);
			// inter = _mm256_sub_ps(topos,topon);
			// tSS = _mm256_div_ps(inter,DXs);
			// // this->_save_SSi( neigh,  nein,  SSs,  tSS);

			// topon = this->_load_r_topo(i);
			// nein = this->_load_r_topoi(i);
			// inter = _mm256_sub_ps(topos,topon);
			// tSS = _mm256_div_ps(inter,DXs);
			// // this->_save_SSi( neigh,  nein,  SSs,  tSS);

			// topon = this->_load_bl_topo(i);
			// nein = this->_load_bl_topoi(i);
			// inter = _mm256_sub_ps(topos,topon);
			// tSS = _mm256_div_ps(inter,DXYs);
			// // this->_save_SSi( neigh,  nein,  SSs,  tSS);

			// topon = this->_load_b_topo(i);
			// nein = this->_load_b_topoi(i);
			// inter = _mm256_sub_ps(topos,topon);
			// tSS = _mm256_div_ps(inter,DYs);
			// // this->_save_SSi( neigh,  nein,  SSs,  tSS);

			// topon = this->_load_br_topo(i);
			// nein = this->_load_br_topoi(i);
			// inter = _mm256_sub_ps(topos,topon);
			// tSS = _mm256_div_ps(inter,DXYs);
			// // this->_save_SSi( neigh,  nein,  SSs,  tSS);


			// int res[8];
			// // _mm256_storeu_epi32(res,nein);
			// _mm256_storeu_si256((__m256i*)res, neigh);

			// int jj = 0;
			// for(int j = 7; j>=0; --j)
			// {
			// 	// std::cout << res[j] << "|";
			// 	this->Sreceivers[i+jj] = res[j];
			// 	++jj;
			// }
			// // this->Sreceivers[i] = neighbours[j];
		}

	}


	void compute_sreceivers_avx_t3(){;}
	// {
	// 	if(this->_topography == nullptr) throw std::runtime_error("no topo linked sns");

	// 	std::array<int,8> neighbours;

	// 	__m256 DXs = fulldx();
	// 	__m256 DYs = fulldy();
	// 	__m256 DXYs = fulldxy();

	// 	#pragma omp parallel for num_threads(8)
	// 	for (int i = this->nx + 1; i < this->nxy - this->nx - 1; i+=8)
	// 	{

	// 		__m256 topos = this->_load_topo(i);
	// 		__m256i recs = this->_load_topoi(i);
	// 		__m256 SSs = _mm256_set1_ps(0.f);

	// 		__m256 topon = this->_load_tl_topo(i);
	// 		__m256i recn = this->_load_tl_topoi(i);
	// 		__m256 inter = _mm256_sub_ps(topos,topon);
	// 		__m256 tSS = _mm256_div_ps(inter,DXYs);
	// 		this->_save_SSi( recs,  recn,  SSs,  tSS);

	// 		topon = this->_load_t_topo(i);
	// 		recn = this->_load_t_topoi(i);
	// 		inter = _mm256_sub_ps(topos,topon);
	// 		tSS = _mm256_div_ps(inter,DYs);
	// 		this->_save_SSi( recs,  recn,  SSs,  tSS);

	// 		topon = this->_load_tr_topo(i);
	// 		recn = this->_load_tr_topoi(i);
	// 		inter = _mm256_sub_ps(topos,topon);
	// 		tSS = _mm256_div_ps(inter,DXYs);
	// 		this->_save_SSi( recs,  recn,  SSs,  tSS);

	// 		topon = this->_load_l_topo(i);
	// 		recn = this->_load_l_topoi(i);
	// 		inter = _mm256_sub_ps(topos,topon);
	// 		tSS = _mm256_div_ps(inter,DXs);
	// 		this->_save_SSi( recs,  recn,  SSs,  tSS);

	// 		topon = this->_load_r_topo(i);
	// 		recn = this->_load_r_topoi(i);
	// 		inter = _mm256_sub_ps(topos,topon);
	// 		tSS = _mm256_div_ps(inter,DXs);
	// 		this->_save_SSi( recs,  recn,  SSs,  tSS);

	// 		topon = this->_load_bl_topo(i);
	// 		recn = this->_load_bl_topoi(i);
	// 		inter = _mm256_sub_ps(topos,topon);
	// 		tSS = _mm256_div_ps(inter,DXYs);
	// 		this->_save_SSi( recs,  recn,  SSs,  tSS);

	// 		topon = this->_load_b_topo(i);
	// 		recn = this->_load_b_topoi(i);
	// 		inter = _mm256_sub_ps(topos,topon);
	// 		tSS = _mm256_div_ps(inter,DYs);
	// 		this->_save_SSi( recs,  recn,  SSs,  tSS);

	// 		topon = this->_load_br_topo(i);
	// 		recn = this->_load_br_topoi(i);
	// 		inter = _mm256_sub_ps(topos,topon);
	// 		tSS = _mm256_div_ps(inter,DXYs);
	// 		this->_save_SSi( recs,  recn,  SSs,  tSS);


	// 		int res[8];
	// 		// _mm256_storeu_epi32(res,recn);
	// 		_mm256_storeu_si256((__m256i*)res, recs);

	// 		int jj = 0;
	// 		for(int j = 7; j>=0; --j)
	// 		{
	// 			// std::cout << res[j] << "|";
	// 			this->Sreceivers[i+jj] = res[j];
	// 			++jj;
	// 		}
	// 		// this->Sreceivers[i] = neighbours[j];
	// 	}

	// }






	void _local_lookup(uint8_t indices, std::array<int,8>& arr, uint8_t& nn)
	{
		// Retrieve the indices specified by the bits set in the `indices` value
		// std::cout << "idx is " << std::bitset<8>(indices) << std::endl;
		// std::cout << "Retained are: ";
		for (int i = 0; i < 8; ++i) 
		{
			if ( indices & (1 << i)) 
			{
				// Index `i` is set, process the corresponding value
				// std::cout << 7 - i << " ";
				// std::cout << 7 - i << " (" << std::bitset<8>(1 << i)  << ") ";
				arr[nn] = this->adder[7 - i];
				// std::cout << std::to_string(arr[nn]) << " | ";
				++nn;
			}
		}
		// std::cout << std::endl ;

	}

	void _compute_lookup_table()
	{
		for(int i = 0; i < LOOKUPSIZE; ++i)
		{
			uint8_t ti = static_cast<uint8_t>(i);
			this->lookup_nn[ti] = 0;
			this->lookup_table_SS[ti] = 0;
			for(auto&v:this->lookup_table[ti]) v=0;

			// std::cout << std::bitset<8>(static_cast<uint8_t>(i)) << "|" << std::to_string(static_cast<uint8_t>(i)) << "||"; 
			this->_local_lookup(ti,this->lookup_table[ti],this->lookup_nn[ti]);
			if(this->lookup_nn[ti] == 1)
				this->lookup_table_SS[ti] = this->lookup_table[ti][0];

		// std::cout << " (" << std::to_string(this->lookup_table[i][4]) << ")  " << std::endl;;
		}
		// exit(0);
	}


	// Set of functions returning the 8 bits corresponding to the local neighbours in a given topology

	uint8_t _get_ui8_rec(int node) {return ui8fromui32(this->Neighbours[node],2);}
	uint8_t _get_ui8_Srec(int node) {return ui8fromui32(this->Neighbours[node],0);}
	uint8_t _get_ui8_Don(int node) {return ui8fromui32(this->Neighbours[node],3);}
	uint8_t _get_ui8_Sdon(int node) {return ui8fromui32(this->Neighbours[node],1);}

	// setting model boundaries last'/first row/col to outletting the model
	void _4edges()
	{
		// setting the first row to 0
		for(int i=0; i<this->nx; ++i)
		{
			this->boundaries[i] = 0;
			this->EdgeNodes.emplace_back(i);
			this->EdgeNodes_Neighbours.emplace_back((i==0)?0b00001011: (this->nx-1 == i)? 0b00010110 : 0b00011111);
		}
		// setting the last row to 0
		for(int i = this->nxy-this->nx; i<this->nxy; ++i)
		{
			this->boundaries[i] = 0;
			this->EdgeNodes.emplace_back(i);
			this->EdgeNodes_Neighbours.emplace_back((i==this->nxy-this->nx) ? 0b01101000 : (i==this->nxy-1) ? 0b11010000: 0b11111000);

		}
		for(int r = 1;r < this->ny-1;++r)
		{
			int i = r * this->nx + 0;
			this->boundaries[i] = 0;
			this->EdgeNodes.emplace_back(i);
			this->EdgeNodes_Neighbours.emplace_back(0b01101011);
			i += this->nx-1;
			this->boundaries[i] = 0;
			this->EdgeNodes.emplace_back(i);
			this->EdgeNodes_Neighbours.emplace_back(0b11010110);
		}
		this->_compute_internal_nodes();
	}

	void _compute_internal_nodes()
	{

		this->IntaNodes.reserve(this->nxy);
		for(int i=0; i<this->nxy; ++i)
		{
			if(this->boundaries[i] == 1)
				this->IntaNodes.emplace_back(i);
		}
	
		int ninternalnnode = static_cast<int>(this->IntaNodes.size());
		int nmode = ninternalnnode % 8;
		
		while(nmode > 0 && nmode < 8)
		{
			this->IntaNodes.emplace_back(this->IntaNodes[this->IntaNodes.size()-1]);
			++nmode;
		}

		std::cout << "DEBUGLOG::INTANODES::" << this->IntaNodes.size() << std::endl;;

		return;

	}


	void _archive_compute_internal_nodes()
	{

		// std::vector<int> tempnodint;
		// tempnodint.reserve(this->nxy);
		// for(int i=0; i<this->nxy; ++i)
		// {
		// 	if(this->boundaries[i] == 1)
		// 		tempnodint.emplace_back(i);
		// }
	
		// int ninternalnnode = static_cast<int>(tempnodint.size());
		// std::cout << "DEBUGLOG::ninternalnnode::" << ninternalnnode << std::endl;
		// int nmode = ninternalnnode % 8;
		// int n_right_size = ninternalnnode - nmode;
		// this->IntaNodes.reserve(static_cast<int>(ninternalnnode/8)+1);



		// for (size_t i = 0; i < static_cast<size_t>(n_right_size) ; i+=8)
		// {
		// 	this->IntaNodes.emplace_back(std::array<int,8>{tempnodint[i + 0],tempnodint[i + 1],tempnodint[i + 2],tempnodint[i + 3],tempnodint[i + 4],tempnodint[i + 5],tempnodint[i + 6],tempnodint[i + 7]});
		// }

		// if(nmode == 0)
		// 	return;

		// std::array<int,8> temp2l;
		// for(int i=0; i<8; ++i)
		// {
		// 	temp2l[i] = (i<nmode) ? tempnodint[i + n_right_size] : tempnodint[n_right_size - 1];
		// }


		// this->IntaNodes.emplace_back(temp2l);

		// return;

	}

	// Utilities
	
	void save_topo(std::string filename)
	{
		const std::vector<std::uint64_t> shape{std::uint64_t(this->ny), std::uint64_t(this->nx)};
		const bool fortran_order{false};
		const std::string path{filename};
		npy::SaveArrayAsNumpy(path, fortran_order, shape.size(), shape.data(), (*this->_topography));
	}

	void save_neighbours(std::string filename)
	{
		const std::vector<std::uint64_t> shape{std::uint64_t(this->ny), std::uint64_t(this->nx)};
		const bool fortran_order{false};
		const std::string path{filename};
		npy::SaveArrayAsNumpy(path, fortran_order, shape.size(), shape.data(), this->Neighbours);
	}

	void save_Sreceivers(std::string filename)
	{
		auto tvec = this->_get_Sreceiver_vector();
		const std::vector<std::uint64_t> shape{std::uint64_t(this->ny), std::uint64_t(this->nx)};
		const bool fortran_order{false};
		const std::string path{filename};
		npy::SaveArrayAsNumpy(path, fortran_order, shape.size(), shape.data(), tvec);
	}

	std::vector<int> _get_Sreceiver_vector()
	{
		std::vector<int> tSreceivers(this->nxy);
		for(int i=0;i<this->nxy;++i)
			tSreceivers[i] = this->Sreceivers(i);

		return tSreceivers;
	}


	void save_nSDonors(std::string filename)
	{
		auto tvec = this->_get_nSDonors_vector();
		const std::vector<std::uint64_t> shape{std::uint64_t(this->ny), std::uint64_t(this->nx)};
		const bool fortran_order{false};
		const std::string path{filename};
		npy::SaveArrayAsNumpy(path, fortran_order, shape.size(), shape.data(), tvec);
	}

	std::vector<int> _get_nSDonors_vector()
	{
		std::vector<int> tnSDonors(this->nxy);
		for(int i=0;i<this->nxy;++i)
			tnSDonors[i] = this->nSDonors(i);

		return tnSDonors;
	}

	template<class T>
	void save_vec(std::string filename, std::vector<T>& data)
	{
		const std::vector<std::uint64_t> shape{std::uint64_t(this->ny), std::uint64_t(this->nx)};
		const bool fortran_order{false};
		const std::string path{filename};
		npy::SaveArrayAsNumpy(path, fortran_order, shape.size(), shape.data(), data);
	}


















	void _check_Sneighbouring()
	{
		std::array<int,8> arr, arr2;
		auto DonMasks = DonMasks();
		auto RecMasks = RecMasks();
		auto SDonMasks = SDonMasks();
		auto SRecMasks = SRecMasks();
		auto MirrorIndicesMasks = MirrorIndicesMasks();

		int n_anomalies = 0;
		int n_anomalies_Srec = 0;
		int n_sdon = 0;
		int n_zerodon = 0;
		for(int i=0; i<this->nxy; ++i)
		{
			// checking Sndonors		
			int nn = this->SDonors(i,arr);
			if(nn==0)
			{
				++n_zerodon;
				this->internal_neighbours_with_check(i,arr2);
				for(int j=0;j<8;++j)
				{
					if(arr2[j] == -1)
						continue;
					if(this->topography(arr2[j]) > this->topography(i))
						std::cout << "JOIFSDDFLKJ::" << std::bitset<32>(this->Neighbours[i]) << "||" << std::bitset<32>(this->Neighbours[arr2[j]]) << std::endl;
				}
				// std::cout << std::bitset<32>(this->Neighbours[i])  << "|" << std::bitset<8>(ui8fromui32(this->Neighbours[i],1)) << std::endl;
			}

			n_sdon++;

			for(int j=0;j<nn;++j)
			{
				int tn = arr[j];
				if(this->Sreceivers(tn) != i)
					++n_anomalies;;
			}

			int tr = this->Sreceivers(i);
			if(tr != i)
			{
				int tnd = this->nSDonors(tr);
				if(tnd < 1)
				{
					std::cout << std::bitset<32>(this->Neighbours[i])  << "|" << std::bitset<32>(this->Neighbours[tr]) << " | " << (float(tr)/this->nx)  << std::endl;
					++n_anomalies_Srec;
				}
			}

		}
		std::cout << "Out of " << this->nxy << " nodes checked " << n_sdon << " and there was " << n_zerodon << " with 0 donors. \nSDonors validity, anomaly founds: " << n_anomalies << "\nSrecs validity, anomaly founds: " << n_anomalies_Srec << std::endl;;
	



	}


};






























































