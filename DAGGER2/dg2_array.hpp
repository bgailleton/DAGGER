#pragma once

#include <cstddef>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <vector>

// Conditional Python support
#ifdef DG2_WITH_PYTHON
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;
#endif

/**
 * ArrayRef: A non-owning reference to array data that can be shared between
 * Python and C++
 *
 * This class provides a unified interface to access array data whether it comes
 * from:
 * - NumPy arrays (through pybind11) [when DG2_WITH_PYTHON is defined]
 * - C++ std::vector
 * - Raw C++ arrays
 *
 * The key feature is zero-copy data sharing with automatic lifetime management.
 */
template<typename T>
class ArrayRef
{
private:
	T* data_ptr_;
	size_t size_;

#ifdef DG2_WITH_PYTHON
	// Holds the Python object reference if data comes from NumPy
	// This ensures the data stays alive as long as this ArrayRef exists
	py::object py_ref_;
#endif

	// For C++ owned data (optional)
	std::shared_ptr<T[]> cpp_data_;

public:
	// Default constructor
	ArrayRef()
		: data_ptr_(nullptr)
		, size_(0)
	{
	}

#ifdef DG2_WITH_PYTHON
	// Constructor from NumPy array (zero-copy)
	explicit ArrayRef(
		py::array_t<T, py::array::c_style | py::array::forcecast> arr)
		: py_ref_(arr)
	{

		py::buffer_info buf = arr.request();

		if (buf.ndim != 1) {
			throw std::runtime_error(
				"ArrayRef expects 1D arrays. For 2D data, flatten it first.");
		}

		data_ptr_ = static_cast<T*>(buf.ptr);
		size_ = buf.size;
	}
#endif

	// Constructor from C++ vector (makes a copy and owns it)
	explicit ArrayRef(const std::vector<T>& vec)
		: size_(vec.size())
	{

		// This is the safest approach
		if (vec.empty()) {
			data_ptr_ = nullptr;
			return;
		}

		cpp_data_ = std::shared_ptr<T[]>(new T[vec.size()]);
		data_ptr_ = cpp_data_.get();
		std::copy(vec.begin(), vec.end(), data_ptr_);
	}

	// Constructor from raw pointer (non-owning)
	ArrayRef(T* ptr, size_t size)
		: data_ptr_(ptr)
		, size_(size)
	{
	}

	// Constructor that takes ownership of raw data
	ArrayRef(std::unique_ptr<T[]> data, size_t size)
		: cpp_data_(std::move(data))
		, data_ptr_(cpp_data_.get())
		, size_(size)
	{
	}

	// Access operators
	T& operator[](size_t i)
	{
		if (i >= size_)
			throw std::out_of_range("Index out of bounds");
		return data_ptr_[i];
	}

	const T& operator[](size_t i) const
	{
		if (i >= size_)
			throw std::out_of_range("Index out of bounds");
		return data_ptr_[i];
	}

	T& at(size_t i) { return (*this)[i]; }
	const T& at(size_t i) const { return (*this)[i]; }

	// Raw data access
	T* data() { return data_ptr_; }
	const T* data() const { return data_ptr_; }

	// Size and capacity
	size_t size() const { return size_; }
	bool empty() const { return size_ == 0; }

	// Iterator support
	T* begin() { return data_ptr_; }
	T* end() { return data_ptr_ + size_; }
	const T* begin() const { return data_ptr_; }
	const T* end() const { return data_ptr_ + size_; }

	// Check if this array is backed by Python data
	bool is_python_backed() const
	{
#ifdef DG2_WITH_PYTHON
		return py_ref_.ptr() != nullptr;
#else
		return false;
#endif
	}

#ifdef DG2_WITH_PYTHON
	// Get as NumPy array (creates new array if not already Python-backed)
	py::array_t<T> as_numpy()
	{
		if (is_python_backed()) {
			return py::cast<py::array_t<T>>(py_ref_);
		} else {
			// Create a new NumPy array that views our data
			return py::array_t<T>(
				{ size_ },											// shape
				{ sizeof(T) },									// strides
				data_ptr_,											// data pointer
				py::capsule(this, [](void* p) { // capsule to manage lifetime
					// This ensures the ArrayRef stays alive as long as NumPy array exists
				}));
		}
	}
#endif
};

/**
 * Grid2D: Wrapper for 2D grid data stored in row-major order
 *
 * This provides 2D indexing on top of 1D array storage, maintaining
 * compatibility with NumPy's default row-major ordering.
 */
template<typename T>
class Grid2D
{
private:
	ArrayRef<T> data_;
	size_t rows_;
	size_t cols_;
	size_t size_;

public:
	// Default constructor
	Grid2D()
		: rows_(0)
		, cols_(0)
		, size_(0)
	{
	}

	// Constructor from vector
	Grid2D(const std::vector<T>& data, size_t rows, size_t cols)
		: data_(data)
		, rows_(rows)
		, cols_(cols)
		, size_(rows * cols)
	{
		std::cout << "DEBUG::3.7" << std::endl;

		if (data.size() != size_) {
			throw std::runtime_error("Data size doesn't match grid dimensions");
		}
	}

	// Constructor from raw data with dimensions
	Grid2D(T* data, size_t rows, size_t cols)
		: data_(data, rows * cols)
		, rows_(rows)
		, cols_(cols)
		, size_(rows * cols)
	{
	}

	// Constructor that creates and owns data
	Grid2D(size_t rows, size_t cols, const T& fill_value = T{})
		: rows_(rows)
		, cols_(cols)
		, size_(rows * cols)
	{
		std::vector<T> vec(size_, fill_value);
		data_ = ArrayRef<T>(vec);
	}

#ifdef DG2_WITH_PYTHON
	// Constructor from NumPy array
	Grid2D(py::array_t<T> array, size_t rows = 0, size_t cols = 0)
	{
		py::buffer_info buf = array.request();

		if (buf.ndim == 1) {
			// 1D array, need explicit dimensions
			if (rows == 0 || cols == 0) {
				throw std::runtime_error("Must specify rows and cols for 1D array");
			}
			if (rows * cols != buf.size) {
				throw std::runtime_error("Dimensions don't match array size");
			}
			rows_ = rows;
			cols_ = cols;
			size_ = buf.size;
			data_ = ArrayRef<T>(array);
		} else if (buf.ndim == 2) {
			// 2D array, extract dimensions
			rows_ = buf.shape[0];
			cols_ = buf.shape[1];
			size_ = buf.size;

			// Flatten the array for internal storage
			auto flat = array.template cast<py::array_t<T, py::array::c_style>>();
			flat = flat.reshape({ static_cast<py::ssize_t>(size_) });
			data_ = ArrayRef<T>(flat);
		} else {
			throw std::runtime_error("Grid2D only supports 1D or 2D arrays");
		}
	}
#endif

	// 2D access
	T& operator()(size_t row, size_t col)
	{
		if (row >= rows_ || col >= cols_) {
			throw std::out_of_range("Grid2D index out of bounds");
		}
		return data_[row * cols_ + col];
	}

	const T& operator()(size_t row, size_t col) const
	{
		if (row >= rows_ || col >= cols_) {
			throw std::out_of_range("Grid2D index out of bounds");
		}
		return data_[row * cols_ + col];
	}

	// 1D access (row-major order)
	T& operator[](size_t i) { return data_[i]; }
	const T& operator[](size_t i) const { return data_[i]; }

	// Properties
	size_t rows() const { return rows_; }
	size_t cols() const { return cols_; }
	size_t size() const { return size_; }

	// Raw data access
	T* data() { return data_.data(); }
	const T* data() const { return data_.data(); }

	// Iterator support
	T* begin() { return data_.begin(); }
	T* end() { return data_.end(); }
	const T* begin() const { return data_.begin(); }
	const T* end() const { return data_.end(); }

#ifdef DG2_WITH_PYTHON
	// Get as 2D NumPy array
	py::array_t<T> as_numpy_2d()
	{
		auto flat_array = data_.as_numpy();
		return flat_array.reshape(
			{ static_cast<py::ssize_t>(rows_), static_cast<py::ssize_t>(cols_) });
	}
#endif

	// Utility functions
	void fill(const T& value) { std::fill(begin(), end(), value); }

	void resize(size_t new_rows, size_t new_cols, const T& fill_value = T{})
	{
		rows_ = new_rows;
		cols_ = new_cols;
		size_ = rows_ * cols_;
		std::vector<T> new_data(size_, fill_value);
		data_ = ArrayRef<T>(new_data);
	}
};

/**
 * Grid3D: Wrapper for 3D grid data stored in depth-row-column order
 */
template<typename T>
class Grid3D
{
private:
	ArrayRef<T> data_;
	size_t depth_;
	size_t rows_;
	size_t cols_;
	size_t size_;

public:
	// Default constructor
	Grid3D()
		: depth_(0)
		, rows_(0)
		, cols_(0)
		, size_(0)
	{
	}

	// Constructor from vector
	Grid3D(const std::vector<T>& data, size_t depth, size_t rows, size_t cols)
		: data_(data)
		, depth_(depth)
		, rows_(rows)
		, cols_(cols)
		, size_(depth * rows * cols)
	{
		if (data.size() != size_) {
			throw std::runtime_error("Data size doesn't match grid dimensions");
		}
	}

	// Constructor that creates and owns data
	Grid3D(size_t depth, size_t rows, size_t cols, const T& fill_value = T{})
		: depth_(depth)
		, rows_(rows)
		, cols_(cols)
		, size_(depth * rows * cols)
	{
		std::vector<T> vec(size_, fill_value);
		data_ = ArrayRef<T>(vec);
	}

#ifdef DG2_WITH_PYTHON
	// Constructor from NumPy array
	Grid3D(py::array_t<T> array,
				 size_t depth = 0,
				 size_t rows = 0,
				 size_t cols = 0)
	{
		py::buffer_info buf = array.request();

		if (buf.ndim == 1) {
			// 1D array, need explicit dimensions
			if (depth == 0 || rows == 0 || cols == 0) {
				throw std::runtime_error(
					"Must specify depth, rows and cols for 1D array");
			}
			if (depth * rows * cols != buf.size) {
				throw std::runtime_error("Dimensions don't match array size");
			}
		} else if (buf.ndim == 3) {
			// 3D array, extract dimensions
			depth = buf.shape[0];
			rows = buf.shape[1];
			cols = buf.shape[2];
		} else {
			throw std::runtime_error("Grid3D only supports 1D or 3D arrays");
		}

		depth_ = depth;
		rows_ = rows;
		cols_ = cols;
		size_ = buf.size;

		// Flatten the array for internal storage
		auto flat = array.template cast<py::array_t<T, py::array::c_style>>();
		flat = flat.reshape({ static_cast<py::ssize_t>(size_) });
		data_ = ArrayRef<T>(flat);
	}
#endif

	// 3D access
	T& operator()(size_t d, size_t row, size_t col)
	{
		if (d >= depth_ || row >= rows_ || col >= cols_) {
			throw std::out_of_range("Grid3D index out of bounds");
		}
		return data_[d * rows_ * cols_ + row * cols_ + col];
	}

	const T& operator()(size_t d, size_t row, size_t col) const
	{
		if (d >= depth_ || row >= rows_ || col >= cols_) {
			throw std::out_of_range("Grid3D index out of bounds");
		}
		return data_[d * rows_ * cols_ + row * cols_ + col];
	}

	// 1D access
	T& operator[](size_t i) { return data_[i]; }
	const T& operator[](size_t i) const { return data_[i]; }

	// Properties
	size_t depth() const { return depth_; }
	size_t rows() const { return rows_; }
	size_t cols() const { return cols_; }
	size_t size() const { return size_; }

	// Raw data access
	T* data() { return data_.data(); }
	const T* data() const { return data_.data(); }

#ifdef DG2_WITH_PYTHON
	// Get as 3D NumPy array
	py::array_t<T> as_numpy_3d()
	{
		auto flat_array = data_.as_numpy();
		return flat_array.reshape({ static_cast<py::ssize_t>(depth_),
																static_cast<py::ssize_t>(rows_),
																static_cast<py::ssize_t>(cols_) });
	}
#endif

	// Slice operations
	Grid2D<T> get_slice(size_t d) const
	{
		if (d >= depth_)
			throw std::out_of_range("Depth index out of bounds");

		std::vector<T> slice_data(rows_ * cols_);
		for (size_t r = 0; r < rows_; ++r) {
			for (size_t c = 0; c < cols_; ++c) {
				slice_data[r * cols_ + c] = (*this)(d, r, c);
			}
		}
		return Grid2D<T>(slice_data, rows_, cols_);
	}

	void set_slice(size_t d, const Grid2D<T>& slice)
	{
		if (d >= depth_)
			throw std::out_of_range("Depth index out of bounds");
		if (slice.rows() != rows_ || slice.cols() != cols_) {
			throw std::runtime_error("Slice dimensions don't match grid");
		}

		for (size_t r = 0; r < rows_; ++r) {
			for (size_t c = 0; c < cols_; ++c) {
				(*this)(d, r, c) = slice(r, c);
			}
		}
	}
};

// ======================
// UTILITY FUNCTIONS
// ======================

template<typename T>
T
compute_sum(ArrayRef<T>& arr)
{
	T sum = 0;
	for (size_t i = 0; i < arr.size(); ++i) {
		sum += arr[i];
	}
	return sum;
}

template<typename T>
T
compute_grid_average(Grid2D<T>& grid)
{
	T sum = 0;
	for (size_t r = 0; r < grid.rows(); ++r) {
		for (size_t c = 0; c < grid.cols(); ++c) {
			sum += grid(r, c);
		}
	}
	return sum / static_cast<T>(grid.size());
}

// Function to create a test grid with some pattern
inline Grid2D<double>
create_test_grid(size_t rows, size_t cols)
{
	std::vector<double> data(rows * cols);
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis(0.0, 100.0);

	for (size_t i = 0; i < data.size(); ++i) {
		data[i] = dis(gen);
	}

	return Grid2D<double>(data, rows, cols);
}

// Function to apply a simple filter to a grid
inline void
apply_smoothing_filter(Grid2D<double>& grid)
{
	if (grid.rows() < 3 || grid.cols() < 3)
		return;

	// Create a copy for reading original values
	std::vector<double> original(grid.begin(), grid.end());

	// Apply 3x3 smoothing filter (avoid edges)
	for (size_t r = 1; r < grid.rows() - 1; ++r) {
		for (size_t c = 1; c < grid.cols() - 1; ++c) {
			double sum = 0.0;
			for (int dr = -1; dr <= 1; ++dr) {
				for (int dc = -1; dc <= 1; ++dc) {
					sum += original[(r + dr) * grid.cols() + (c + dc)];
				}
			}
			grid(r, c) = sum / 9.0;
		}
	}
}
