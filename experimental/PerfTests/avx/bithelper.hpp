#ifndef BITHELPER_HPP
#define BITHELPER_HPP
#include <bitset>

static inline __m256i blendvps_si256(__m256i& a, __m256i& b, __m256i& mask) 
{
	__m256 res = _mm256_blendv_ps(_mm256_castsi256_ps(a), _mm256_castsi256_ps(b), _mm256_castsi256_ps(mask));
	return _mm256_castps_si256(res);
}

// This function has been checked on the 5th of July
// Warning! it starts from the right!
// byteindex 0 -> 8 last bits from the left
//...
// byteindex 3 -> 8th first bits
std::uint8_t ui8fromui32(std::uint32_t value, int byteIndex) 
{
	return (value >> (byteIndex * 8)) & 0xFF;
}



void AiisBiifCgtD(__m256i& A, __m256i& B, __m256& C, __m256& D)
{
	__m256 comparisonResult = _mm256_cmp_ps(D, C, _CMP_GT_OQ);
	C = _mm256_blendv_ps(C, D, comparisonResult);
	__m256i comparisonMask = _mm256_castps_si256(comparisonResult);
	A = blendvps_si256(A, B, comparisonMask);
}


std::uint8_t reverseBits(std::uint8_t num) 
{
	std::uint8_t result = 0;

	for (int i = 0; i < 8; ++i) 
	{
		result <<= 1;
		result |= (num & 1);
		num >>= 1;
	}

	return result;
}


__m256i bitwiseORWithMask(__m256i& vec1, __m256i& vec2, __m256i& mask) 
{
	__m256i maskedVec1 = _mm256_and_si256(vec1, mask);
	__m256i maskedVec2 = _mm256_and_si256(vec2, _mm256_xor_si256(mask, _mm256_set1_epi32(-1)));
	return _mm256_or_si256(maskedVec1, maskedVec2);
}


// from https://stackoverflow.com/a/42616203/7114716
inline __m256i _mm256_not_si256 (__m256i a)
{    
	//return  _mm256_xor_si256 (a, _mm256_set1_epi32(0xffffffff));
	return  _mm256_xor_si256 (a, _mm256_cmpeq_epi32(a,a));//I didn't check wich one is faster   
}


// This function has been checked on the 5th of July
int first_index_from_8bits(std::uint8_t indices, int default_return = -1) 
{

	// Retrieve the indices specified by the bits set in the `indices` value
	for (int i = 0; i < 8; ++i) 
	{
		if (indices & (1 << i)) 
		{
			// Index `i` is set, process the corresponding value
			return 7 - i;
		}
	}
	return default_return;
}



void printM256(__m256 vec) 
{
    alignas(32) float temp[8];
    _mm256_storeu_ps(temp, vec);

    for (int i = 0; i < 8; i++) {
        std::cout << temp[i] << " ";
    }
    std::cout << std::endl;
}

void printM256i(__m256i vec) 
{
    alignas(32) int temp[8];
    _mm256_storeu_si256((__m256i*)temp, vec);

    for (int i = 0; i < 8; i++) {
        std::cout << temp[i] << " ";
    }
    std::cout << std::endl;
}

void printM256i_bitset(__m256i vec) 
{
    alignas(32) int temp[8];
    _mm256_storeu_si256((__m256i*)temp, vec);

    for (int i = 0; i < 8; i++) {
        std::cout << std::bitset<32>(temp[i]) << " ";
    }
    std::cout << std::endl;
}

void printM256i_bitset_8(__m256i vec, int byteIndex) 
{
    alignas(32) int temp[8];
    _mm256_storeu_si256((__m256i*)temp, vec);

    for (int i = 0; i < 8; i++) {
        std::cout << std::bitset<8>(ui8fromui32(temp[i],byteIndex)) << " ";
    }
    std::cout << std::endl;
}



























































#endif