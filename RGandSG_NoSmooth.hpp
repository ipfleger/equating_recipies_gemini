#ifndef RGANDSG_NOSMOOTH_HPP
#define RGANDSG_NOSMOOTH_HPP

#include "ERutilities.hpp" // Assumed to exist from previous conversions
#include <vector>
#include <string>
#include <iostream>

// --- Forward declarations of structs defined elsewhere ---
struct USTATS;
struct BSTATS;
struct PDATA;
struct ERAW_RESULTS;

// --- Equating Functions ---

/**
 * @brief Wrapper for equating with Random Groups (RG) design, no smoothing.
 */
void Wrapper_RN(char design, char method, char smoothing,
                USTATS& x, USTATS& y, int rep,
                PDATA& inall, ERAW_RESULTS& r);

/**
 * @brief Wrapper for equating with Single Group (SG) design, no smoothing.
 */
void Wrapper_SN(char design, char method, char smoothing, BSTATS& xy,
                int rep, PDATA& inall, ERAW_RESULTS& r);

/**
 * @brief Performs linear or mean equating for RG and SG designs.
 */
void RGandSG_LinEq(double mnx, double sdx, double mny, double sdy,
                   char method, double min, double max, double inc,
                   double& a, double& b, std::vector<double>& eraw);

// --- Print Functions ---

/**
 * @brief Prints results from the Wrapper_RN function.
 */
void Print_RN(std::ostream& os, const std::string& title, const PDATA& inall, const ERAW_RESULTS& r);

/**
 * @brief Prints results from the Wrapper_SN function.
 */
void Print_SN(std::ostream& os, const std::string& title, const PDATA& inall, const ERAW_RESULTS& r);

#endif // RGANDSG_NOSMOOTH_HPP