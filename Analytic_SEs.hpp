/*
  Analytic_SEs.hpp                                            
*/

#ifndef ANALYTIC_SES_HPP
#define ANALYTIC_SES_HPP

// For a more idiomatic C++ approach, consider using std::vector<double>& 
// instead of raw pointers for array arguments.

/**
 * @brief Calculates estimated standard errors for equipercentile equivalents.
 * * @param nsy Number of raw score categories for old form Y.
 * @param npy Number of persons who took old form Y.
 * @param crfdy Cumulative relative frequency distribution for old form Y.
 * @param nsx Number of raw score categories for new form X.
 * @param incx Increment between scores for x.
 * @param npx Number of persons who took new form X.
 * @param prdx Percentile rank distribution for new form X.
 * @param se_eeq Output array for standard errors of equipercentile equivalents.
 */
void SE_EPequate(int nsy, double npy, double* crfdy, int nsx, 
                 double incx, double npx, double* prdx, double* se_eeq);

#endif // ANALYTIC_SES_HPP