#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>

// --- Vector and Matrix Operations (0-indexed) ---

/**
 * @brief Multiplies a matrix by a vector: r = m * v
 * @param m Input matrix.
 * @param v Input vector.
 * @param r Output vector.
 */
void matrix_multiply(const std::vector<std::vector<double>>& m, const std::vector<double>& v, std::vector<double>& r);

/**
 * @brief Multiplies a vector by a matrix: r = v * m
 * @param v Input vector (treated as a row vector).
 * @param m Input matrix.
 * @param r Output vector.
 */
void vector_multiply(const std::vector<double>& v, const std::vector<std::vector<double>>& m, std::vector<double>& r);

/**
 * @brief Multiplies two matrices: r = m * n
 * @param m First input matrix.
 * @param n Second input matrix.
 * @param r Output matrix.
 */
void matrix_multiply(const std::vector<std::vector<double>>& m, const std::vector<std::vector<double>>& n, std::vector<std::vector<double>>& r);

/**
 * @brief Computes the Euclidean norm (L2 norm) of a vector.
 * @param v Input vector.
 * @return The norm of the vector.
 */
double vector_norm(const std::vector<double>& v);

/**
 * @brief Computes P_transpose * P.
 * @param P Input matrix.
 * @param PTP Output matrix (P_transpose * P).
 */
void p_transpose_multiply_p(const std::vector<std::vector<double>>& P, std::vector<std::vector<double>>& PTP);

/**
 * @brief Computes P_transpose * A * P.
 * @param P Input matrix P.
 * @param A Input square matrix A.
 * @param PtAP Output matrix.
 */
void pt_a_p(const std::vector<std::vector<double>>& P, const std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& PtAP);

/**
 * @brief Computes v_transpose * A * v.
 * @param v Input vector v.
 * @param A Input square matrix A.
 * @return The resulting scalar value.
 */
double vt_a_v(const std::vector<double>& v, const std::vector<std::vector<double>>& A);

/**
 * @brief Computes the determinant of a square matrix.
 * @param a Input square matrix.
 * @return The determinant.
 */
double determinant(const std::vector<std::vector<double>>& a);

/**
 * @brief Computes the inverse of a square matrix.
 * @param a The matrix to invert. On success, this matrix is replaced by its inverse.
 */
void matrix_inverse(std::vector<std::vector<double>>& a);

#endif // MATRIX_HPP