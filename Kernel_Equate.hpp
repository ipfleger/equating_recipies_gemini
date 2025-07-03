#ifndef KERNEL_EQUATE_HPP
#define KERNEL_EQUATE_HPP

#include <vector>
#include <functional>
#include <iostream>

// Forward declarations for structs defined in other modules
struct BLL_SMOOTH;
struct USTATS;
struct BSTATS;
struct ULL_SMOOTH;
struct PDATA;
struct ERAW_RESULTS;

// Type alias for a C++ standard function pointer
using PtrKernelPdf = std::function<double(int, const std::vector<double>&, const std::vector<double>&, double, double)>;

// --- Core Kernel Functions ---

double KernelContinuPdf(int ncat, const std::vector<double>& scores, const std::vector<double>& fd, double hx, double x);
double KernelContinuCdf(int ncat, const std::vector<double>& scores, const std::vector<double>& fd, double hx, double x);
double StdNormalCdf(double x);
double StdNormalPdf(double x);
double NormalPdf(const std::vector<double>& para, double x);
void CalcKernelMoments(PtrKernelPdf pdf, double a, double b, int nDistCat, const std::vector<double>& score,
                       const std::vector<double>& fd, double hx, std::vector<double>& moments);
double KernelInverseCdf(int ncat, const std::vector<double>& scores, const std::vector<double>& fd, double h, double cdf);
void KernelEquate(int nDistCatx, const std::vector<double>& scoresx, const std::vector<double>& fdx, double hx,
                  int nDistCaty, const std::vector<double>& scoresy, const std::vector<double>& fdy, double hy,
                  std::vector<double>& Equatedx);

// --- Bandwidth Optimization Functions ---

double Pen1(int nDistCat, const std::vector<double>& scores, const std::vector<double>& fd, double hx);
double KernelPdfDerivative(int nDistCat, const std::vector<double>& scores, const std::vector<double>& fd, double hx, double x);
double Pen2(int nDistCat, const std::vector<double>& scores, const std::vector<double>& fd, double hx);
double Pen(int nDistCat, const std::vector<double>& scores, const std::vector<double>& fd, double hx, double K);
double Optimalh(int nDistCat, const std::vector<double>& scores, const std::vector<double>& fd, double K);

// --- Matrix and SEE Calculation Functions ---

void ComputeCmatrixGen(int ncat, int degree, long np, const std::vector<std::vector<double>>& B, const std::vector<double>& fd, std::vector<std::vector<double>>& Cr);
void PartialFPartialr(int ncat, const std::vector<double>& scores, const std::vector<double>& fd, double hx, std::vector<double>& Fr, double x);
double FrCrSqNorm(int ncat, int degree, const std::vector<double>& Fr, const std::vector<std::vector<double>>& Cr);

// --- Vector and Matrix Operations ---

void vPMN(int ncatx, int ncaty, const std::vector<std::vector<double>>& bdist, std::vector<double>& vP, std::vector<std::vector<double>>& M, std::vector<std::vector<double>>& N);
void MatrixMultiVector(const std::vector<std::vector<double>>& m, const std::vector<double>& v, std::vector<double>& r);
void VectorMultiMatrix(const std::vector<double>& v, const std::vector<std::vector<double>>& m, std::vector<double>& r);
double VectorNormSq(const std::vector<double>& v);

// --- High-Level Equating Functions for Different Designs ---

void KernelEquateSEERG(int nDistCatx, int degreex, long npx, const std::vector<double>& scoresx, const std::vector<double>& fdx,
                       int nDistCaty, int degreey, long npy, const std::vector<double>& scoresy,
                       const std::vector<double>& fdy, std::vector<double>& Equatedx, std::vector<double>& SEE);

void KernelEquateSEESG(BLL_SMOOTH& bivar, std::vector<double>& Equatedx, std::vector<double>& SEE);
void KernelEquateSG(BLL_SMOOTH& bivar, std::vector<double>& Equatedx);

void KernelEquateSEECB(BLL_SMOOTH& bivar1, BLL_SMOOTH& bivar2, const std::vector<double>& wts, std::vector<double>& Equatedx,
                       std::vector<double>& SEE);

void KernelEquateSEENEATPS(BLL_SMOOTH& bivar1, BLL_SMOOTH& bivar2, double wts,
                           std::vector<double>& Equatedx, std::vector<double>& SEE);

void KernelEquateNEATPS(BLL_SMOOTH& bivar1, BLL_SMOOTH& bivar2, double wts,
                        std::vector<double>& Equatedx);

void KernelEquateSEENEATChn(BLL_SMOOTH& bivar1, BLL_SMOOTH& bivar2,
                            std::vector<double>& Equatedx, std::vector<double>& SEE);

void KernelEquateNEATChn(BLL_SMOOTH& bivar1, BLL_SMOOTH& bivar2, std::vector<double>& Equatedx);


// --- Wrappers and Print Functions ---

void Wrapper_RK(char design, char method, char smoothing,
                USTATS& x, USTATS& y,
                ULL_SMOOTH& ullx, ULL_SMOOTH& ully, int rep,
                PDATA& inall, ERAW_RESULTS& r);

void Print_RK(std::ostream& os, const std::string& title, const PDATA& inall, const ERAW_RESULTS& r);

void Wrapper_SK(char design, char method, char smoothing,  BSTATS& xy,
                BLL_SMOOTH& bllxy, int rep, PDATA& inall,
                ERAW_RESULTS& r);

void Print_SK(std::ostream& os, const std::string& title, const PDATA& inall, const ERAW_RESULTS& r);

void Wrapper_CK(char design, char method, char smoothing,
                double w1, int anchor, double rv1, double rv2,
                BSTATS& xv, BSTATS& yv,
                BLL_SMOOTH& bllxv, BLL_SMOOTH& bllyv,
                int rep, PDATA& inall, ERAW_RESULTS& r);

void Print_CK(std::ostream& os, const std::string& title, const PDATA& inall, const ERAW_RESULTS& r);

#endif // KERNEL_EQUATE_HPP