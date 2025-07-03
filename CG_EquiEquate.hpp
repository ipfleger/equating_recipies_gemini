/*
  CG_EquiEquate.hpp
*/

#ifndef CG_EQUIEQUATE_HPP
#define CG_EQUIEQUATE_HPP

#include <cstdio>
#include <string>

// Note: External dependency 'ERutilities.h' is not provided.
#include "ERutilities.h"

// Forward declarations of structs assumed to be in ERutilities.h
struct PDATA;
struct ERAW_RESULTS;

void FEorMFE_EE(double w1, int internal, int nsv,
                int nsx, double minx, double maxx,
                int nsy, double miny, double maxy, double inc,
                double **bxvin, double **byvin, double rv1, double rv2,
                double *fxs, double *gys, double *eraw,
                double *a, double *b, double *erawBH);
void SyntheticDensities(double w1, int internal, int nsv, int nsx, double **bxv,
                        int nsy, double **byv, double rv1, double rv2,
                        double *fxs, double *gys);
void MixSmooth(int nsv, int nsx, double unimix, double **braw,
               double *fx, double *hv);
void CondBivDist(int nsx, int nsv, double **bxv, double *hv);
void BH_LinEq(double minx, double maxx, double miny, double maxy, double inc,
              double *fxs, double *gys, double *a, double *b);
void ModCondBivDist(int internal, int nsv, int nsx, double rv1, double rv2,
                    double muv1, double muv2, double **bxv);
void Chained_EE(int nsx, double *prx1, double minv, double maxv,
                double incv, int nsv, double *Gv1, double miny,
                double incy, int nsy, double *Gy2, double *Gv2,
                double *eraw);
void Print_SynDens(FILE *fp, const std::string& tt, PDATA* inall, ERAW_RESULTS* r);

#endif // CG_EQUIEQUATE_HPP