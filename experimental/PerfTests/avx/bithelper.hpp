#ifndef BITHELPER_HPP
#define BITHELPER_HPP
#include <bitset>
#include <cassert>

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
uint8_t ui8fromui32(uint32_t value, int byteIndex) 
{
	// assert(0<= byteIndex && byteIndex <=3);
	return (value >> (byteIndex * 8)) & 0xFF;
}



void AiisBiifCgtD(__m256i& A, __m256i& B, __m256& C, __m256& D)
{
	__m256 comparisonResult = _mm256_cmp_ps(D, C, _CMP_GT_OQ);
	C = _mm256_blendv_ps(C, D, comparisonResult);
	__m256i comparisonMask = _mm256_castps_si256(comparisonResult);
	A = blendvps_si256(A, B, comparisonMask);
}


uint8_t reverseBits(uint8_t num) 
{
	uint8_t result = 0;

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
int first_index_from_8bits(uint8_t indices, int default_return = -1) 
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



// Get the maximum value in the array
uint32_t getMax(const std::vector<uint32_t>& arr) {
		uint32_t maxVal = arr[0];
		for (uint32_t i = 1; i < arr.size(); ++i) {
				if (arr[i] > maxVal) {
						maxVal = arr[i];
				}
		}
		return maxVal;
}

// Perform counting sort based on a specific digit
void countingSort(std::vector<uint32_t>& arr, int exp) {
		const int n = arr.size();
		std::vector<uint32_t> output(n);
		std::vector<uint32_t> count(10, 0);

		// Count the occurrences of each digit
		for (int i = 0; i < n; ++i) {
				++count[(arr[i] / exp) % 10];
		}

		// Compute cumulative counts to determine positions
		for (int i = 1; i < 10; ++i) {
				count[i] += count[i - 1];
		}

		// Build the output array
		for (int i = n - 1; i >= 0; --i) {
				output[count[(arr[i] / exp) % 10] - 1] = arr[i];
				--count[(arr[i] / exp) % 10];
		}

		// Copy the sorted elements back to the input array
		for (int i = 0; i < n; ++i) {
				arr[i] = output[i];
		}
}

// Radix sort implementation
void radixSort(std::vector<uint32_t>& arr) 
{
	const uint32_t maxVal = getMax(arr);

	// Perform counting sort for each digit
	for (int exp = 1; maxVal / exp > 0; exp *= 10) {
			countingSort(arr, exp);
	}
}



// Partition the array and return the index of the pivot element
int partition(std::vector<uint32_t>& arr, int low, int high) {
    uint32_t pivot = arr[high]; // Choose the last element as the pivot
    int i = low - 1;

    for (int j = low; j < high; ++j) {
        if (arr[j] <= pivot) {
            ++i;
            std::swap(arr[i], arr[j]);
        }
    }
    std::swap(arr[i + 1], arr[high]);
    return i + 1;
}

// Quicksort algorithm implementation
void quicksort(std::vector<uint32_t>& arr, int low, int high) {
    if (low < high) {
        int pivotIndex = partition(arr, low, high);

        quicksort(arr, low, pivotIndex - 1);
        quicksort(arr, pivotIndex + 1, high);
    }
}













// Get the maximum value in the array within a given range
uint32_t getMax(const std::vector<uint32_t>& arr, int start, int end) {
    uint32_t maxVal = arr[start];
    for (int i = start + 1; i <= end; ++i) {
        if (arr[i] > maxVal) {
            maxVal = arr[i];
        }
    }
    return maxVal;
}

// Perform counting sort based on a specific digit for a subsection of the array
void countingSort(std::vector<uint32_t>& arr, int start, int end, int exp) {
    const int n = end - start + 1;
    std::vector<uint32_t> output(n);
    std::vector<uint32_t> count(10, 0);

    // Count the occurrences of each digit
    for (int i = start; i <= end; ++i) {
        ++count[(arr[i] / exp) % 10];
    }

    // Compute cumulative counts to determine positions
    for (int i = 1; i < 10; ++i) {
        count[i] += count[i - 1];
    }

    // Build the output array
    for (int i = end; i >= start; --i) {
        output[count[(arr[i] / exp) % 10] - 1] = arr[i];
        --count[(arr[i] / exp) % 10];
    }

    // Copy the sorted elements back to the input array
    for (int i = start; i <= end; ++i) {
        arr[i] = output[i - start];
    }
}

// Radix sort implementation for a subsection of the array
void radixSort(std::vector<uint32_t>& arr, int start, int end) {
    const uint32_t maxVal = getMax(arr, start, end);

    // Perform counting sort for each digit
    for (int exp = 1; maxVal / exp > 0; exp *= 10) {
        countingSort(arr, start, end, exp);
    }
}



// Perform counting sort based on a specific digit
void reversedRcountingSort(std::vector<uint32_t>& arr, int exp) {
		const int n = arr.size();
		std::vector<uint32_t> output(n);
		std::vector<uint32_t> count(10, 0);
		// std::cout << "to1" << std::endl;
		// Count the occurrences of each digit
		for (int i = 0; i < n; ++i) {
				++count[(arr[i] / exp) % 10];
		}
		// std::cout << "to1" << std::endl;

		// Compute cumulative counts to determine positions
		for (int i = 1; i < 10; ++i) {
				count[i] += count[i - 1];
		}
		// std::cout << "to1" << std::endl;

		// Build the output array
		for (int i = n - 1; i >= 0; --i) {
				output[count[(arr[i] / exp) % 10] - 1] = arr[i];
				--count[(arr[i] / exp) % 10];
		}
		// std::cout << "to1" << std::endl;

		// Copy the sorted elements back to the input array
		for (int i = 0; i < n; ++i) {
				arr[i] = output[n-i-1];
		}
}

// Radix sort implementation
void reversedRadixSort(std::vector<uint32_t>& arr) 
{
	const uint32_t maxVal = getMax(arr);

	// Perform counting sort for each digit
	for (int exp = 1; maxVal / exp > 0; exp *= 10) {
			reversedRcountingSort(arr, exp);
	}
}



// Perform counting sort based on a specific digit for a subsection of the array
void reversedCountingSort(std::vector<uint32_t>& arr, int start, int end, int exp) {
    const int n = end - start + 1;
    std::vector<uint32_t> output(n);
    std::vector<uint32_t> count(10, 0);
    // std::cout << 1 << std::endl;
    // Count the occurrences of each digit
    for (int i = start; i <= end; ++i) {
        ++count[(arr[i] / exp) % 10];
    }
    // std::cout << 2 << std::endl;

    // Compute cumulative counts to determine positions
    for (int i = 1; i < 10; ++i) {
        count[i] += count[i - 1];
    }
    // std::cout << 3 << std::endl;

    // Build the output array
    for (int i = end; i >= start; --i) {
        output[count[(arr[i] / exp) % 10] - 1] = arr[i];
        --count[(arr[i] / exp) % 10];
    }
    // std::cout << 4 << std::endl;

    // Copy the sorted elements back to the input array
    for (int i = end; i >= start; --i) {
        arr[i] = output[i-start];
    }
    // std::cout << 5 << std::endl;

}

// Radix sort implementation for a subsection of the array
void reversedRadixSort(std::vector<uint32_t>& arr, int start, int end) {
    const uint32_t maxVal = getMax(arr, start, end);

    // Perform counting sort for each digit
    for (int exp = 1; maxVal / exp > 0; exp *= 10) {
        reversedCountingSort(arr, start, end, exp);
    }
}



















#endif