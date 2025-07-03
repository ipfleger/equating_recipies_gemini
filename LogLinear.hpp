#ifndef LOGLINEAR_HPP
#define LOGLINEAR_HPP

#include "ERutilities.hpp" // Assumed to exist from previous conversions
#include <vector>
#include <string>
#include <array>
#include <iostream>

// --- Forward declarations of structs defined elsewhere ---
struct USTATS;
struct BSTATS;
struct ULL_SMOOTH;
struct BLL_SMOOTH;
struct PDATA;
struct ERAW_RESULTS;

// --- Matrix and Vector Operations ---

void mtranspose(const std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& at);
void mmult_a_b(const std::vector<std::vector<double>>& a, const std::vector<std::vector<double>>& b, std::vector<std::vector<double>>& s);
void mmult_a_v(const std::vector<std::vector<double>>& a, const std::vector<double>& v, std::vector<double>& s);

// --- Core Log-Linear Functions ---

void design_matrix(int nsu, double minu, double incu,
                   int nsv, double minv, double incv,
                   int cu, int cv, int cuv, const std::vector<std::array<int, 2>>& cpm, int scale,
                   std::vector<std::vector<double>>& B_raw, std::vector<std::vector<double>>& B);

void get_nct_bfd(int anchor, int nsx, int nsv, const std::vector<std::vector<double>>& bfd,
                 std::vector<double>& nct);

void get_bfd_mct(int anchor, int nsx, int nsv, const std::vector<double>& mct,
                 std::vector<std::vector<double>>& bfd);

void get_BtSmB(const std::vector<std::vector<double>>& B, const std::vector<double>& m, double N,
               std::vector<std::vector<double>>& BtSmB);

void get_Btnm(const std::vector<std::vector<double>>& B, const std::vector<double>& n, const std::vector<double>& m,
              std::vector<double>& Btnm);

void get_Beta0(const std::vector<std::vector<double>>& B, const std::vector<double>& n, double N,
               std::vector<double>& Beta0, std::ostream* os_debug = nullptr);

double get_mct(const std::vector<std::vector<double>>& B, const std::vector<double>& Beta, const std::vector<double>* uin, double N,
               std::vector<double>& m, std::ostream* os_debug = nullptr);

int iteration(std::ostream* os, const std::vector<std::vector<double>>& B, const std::vector<std::vector<double>>& B_raw,
              const std::vector<double>& nct, double N, const std::vector<double>* uin,
              int cu, int cv, int cuv, const std::vector<std::array<int, 2>>& cpm,
              int max_nit, int ctype, int Btype, double crit,
              std::vector<double>& Beta, std::vector<double>& mct,
              std::vector<double>& n_mts, std::vector<double>& m_mts,
              std::vector<double>& n_mts_raw, std::vector<double>& m_mts_raw,
              double& lrc, int& nzero, double& ap);

void get_LLmoments(const std::vector<std::vector<double>>& B, const std::vector<std::vector<double>>& B_raw,
                   const std::vector<double>& f, double N,
                   int cu, int cv, int cuv, const std::vector<std::array<int, 2>>& cpm,
                   std::vector<double>& mts, std::vector<double>& mts_raw);

bool crit_mts(int nc, int cu, int ctype, int Btype,
              const std::vector<double>& n_mts, const std::vector<double>& m_mts, double crit);


// --- Smoothing Wrappers ---

void Wrapper_Smooth_ULL(USTATS& x, int c, int scale, int Btype, int ctype, double crit,
                        std::ostream* os, ULL_SMOOTH& s);

void Smooth_ULL(int n, int ns, double min, double inc,
                const std::vector<double>& fd, int c,
                int scale, int Btype, int ctype, double crit,
                std::ostream* os, ULL_SMOOTH& s);

void Wrapper_Smooth_BLL(BSTATS& xv, int anchor,
                        int cu, int cv, int cuv, const std::vector<std::array<int, 2>>& cpm,
                        int scale, int Btype, int ctype, double crit,
                        std::ostream* os, BLL_SMOOTH& s);

void Smooth_BLL(int n, int nsu, double minu, double incu,
                int nsv, double minv, double incv, const std::vector<double>& nct,
                int anchor, int cu, int cv, int cuv, const std::vector<std::array<int, 2>>& cpm,
                int scale, int Btype, int ctype, double crit,
                std::ostream* os, BLL_SMOOTH& s);


// --- Equating Wrappers ---

void Wrapper_RL(char design, char method, char smoothing,
                USTATS& x, USTATS& y,
                ULL_SMOOTH& ullx, ULL_SMOOTH& ully, int rep,
                PDATA& inall, ERAW_RESULTS& r);

void Wrapper_SL(char design, char method, char smoothing, BSTATS& xy,
                BLL_SMOOTH& bllxy, int rep, PDATA& inall,
                ERAW_RESULTS& r);

void Wrapper_CL(char design, char method, char smoothing,
                double w1, int anchor, double rv1, double rv2,
                BSTATS& xv, BSTATS& yv,
                BLL_SMOOTH& bllxv, BLL_SMOOTH& bllyv,
                int rep, PDATA& inall, ERAW_RESULTS& r);


// --- Print Functions ---

void Print_ULL(std::ostream& os, const std::string& title, const USTATS& x,
               const ULL_SMOOTH& s, int print_dm, int print_mts);

void Print_BLL(std::ostream& os, const std::string& title, const BSTATS& xv, const BLL_SMOOTH& s,
               int print_dm, int print_mts, int print_freq, int print_bfd);

void Print_RL(std::ostream& os, const std::string& title, const PDATA& inall, const ERAW_RESULTS& r);
void Print_SL(std::ostream& os, const std::string& title, const PDATA& inall, const ERAW_RESULTS& r);
void Print_CL(std::ostream& os, const std::string& title, const PDATA& inall, const ERAW_RESULTS& r);

#endif // LOGLINEAR_HPP