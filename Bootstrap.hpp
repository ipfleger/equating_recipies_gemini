/*
  Bootstrap.hpp --- header for Bootstrap.cpp
*/

#ifndef BOOTSTRAP_HPP
#define BOOTSTRAP_HPP

#include <cstdio>
#include <string>

// Note: External dependency 'ERutilities.h' is not provided.
#include "ERutilities.h"

// Forward declarations of structs assumed to be in ERutilities.h
struct PDATA;
struct BOOT_ERAW_RESULTS;
struct BOOT_ESS_RESULTS;
struct USTATS;
struct BSTATS;
struct ERAW_RESULTS;
struct ESS_RESULTS;
struct BB_SMOOTH;
struct ULL_SMOOTH;
struct BLL_SMOOTH;


void Wrapper_Bootstrap(PDATA* inall, int nrep, long* idum,
                       BOOT_ERAW_RESULTS* t, BOOT_ESS_RESULTS* u);
void Boot_initialize_eraw(PDATA* inall, BOOT_ERAW_RESULTS* t);
void Boot_USTATS(USTATS* x, long* idum, int rep, USTATS* xb);
void Boot_BSTATS(BSTATS* xv, long* idum, int rep, BSTATS* s);
void Boot_accumulate_eraw(PDATA* inall, ERAW_RESULTS* b,
                          BOOT_ERAW_RESULTS* t);
void Boot_se_eraw(PDATA* inall, BOOT_ERAW_RESULTS* t);
void Print_Boot_se_eraw(FILE* fp, const std::string& tt, PDATA* inall,
                        ERAW_RESULTS* r,
                        BOOT_ERAW_RESULTS* t, int mdiff);
void Boot_initialize_ess(PDATA* inall, BOOT_ESS_RESULTS* u);
void Boot_accumulate_ess(PDATA* inall, ESS_RESULTS* s,
                         BOOT_ESS_RESULTS* u);
void Boot_se_ess(PDATA* inall, BOOT_ESS_RESULTS* u);
void Print_Boot_se_ess(FILE* fp, const std::string& tt, PDATA* inall,
                       ESS_RESULTS* s,
                       BOOT_ESS_RESULTS* u, int mdiff);
void Parametric_boot_univ_BB(BB_SMOOTH* x, long* idum, int rep,
                             BB_SMOOTH* btx);
void Parametric_boot_univ_ULL(ULL_SMOOTH* x, long* idum, int rep,
                              ULL_SMOOTH* btx);
void Parametric_boot_biv(BLL_SMOOTH* xv, long* idum, int rep,
                         BLL_SMOOTH* btxv);

#endif // BOOTSTRAP_HPP