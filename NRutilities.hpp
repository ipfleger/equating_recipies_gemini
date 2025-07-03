#ifndef NRUTILITIES_HPP
#define NRUTILITIES_HPP

#include <string>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <algorithm> // For std::max

// --- Modern C++ Replacements for NR Utilities ---

/**
 * @brief Throws a std::runtime_error with a formatted message.
 * Replaces the original nrerror which terminated the program.
 */
[[noreturn]] inline void nr_error(const std::string& error_text) {
    throw std::runtime_error("Numerical Recipes run-time error: " + error_text);
}


// --- Type-Safe Replacements for C Macros ---

/**
 * @brief Returns the square of a value.
 * Using a template function is type-safe and avoids side effects of C macros.
 */
template<typename T>
constexpr T dsqr(T a) {
    return a * a;
}

/**
 * @brief Returns the maximum of two values.
 * Using std::max is the idiomatic C++ way.
 */
template<typename T>
constexpr T dmax(T a, T b) {
    return std::max(a, b);
}

/**
 * @brief Copies the sign of b to the absolute value of a.
 */
template<typename T>
constexpr T sign(T a, T b) {
    return (b >= 0.0) ? std::abs(a) : -std::abs(a);
}


// --- Deprecation Notice for Memory Management Functions ---

/*
 * NOTE on Memory Management:
 *
 * The original C functions from NRutilities.c (`dvector`, `dmatrix`, `free_dvector`, etc.)
 * for manual memory allocation are now obsolete. In modern C++, these are fully replaced
 * by the standard library container `std::vector`.
 *
 * Instead of:
 * double *v = dvector(0, n - 1);
 * ...
 * free_dvector(v, 0, n - 1);
 *
 * Use `std::vector`:
 * std::vector<double> v(n);
 *
 * Instead of:
 * double **m = dmatrix(0, rows - 1, 0, cols - 1);
 * ...
 * free_dmatrix(m, 0, rows - 1, 0, cols - 1);
 *
 * Use a `std::vector` of `std::vector`s:
 * std::vector<std::vector<double>> m(rows, std::vector<double>(cols));
 *
 * `std::vector` automatically handles all memory allocation and deallocation,
 * prevents memory leaks, provides bounds checking, and offers a rich interface
 * of member functions. All code has been converted to use this modern approach.
 */

#endif // NRUTILITIES_HPP