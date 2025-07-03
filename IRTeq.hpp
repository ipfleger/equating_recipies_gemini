#ifndef IRTEQ_HPP
#define IRTEQ_HPP

#include "IRTst.hpp" // Assumes IRTst.hpp from previous conversion exists
#include "ERutilities.hpp" // Assumes ERutilities.hpp from previous conversion exists

#include <vector>
#include <string>
#include <iostream>

// --- Data Structures ---

// Stores IRT fitted distributions
struct RawFitDist {
    int nRaws = 0;
    std::vector<double> rawScrs;
    std::vector<double> newFits;
    std::vector<double> oldFits;
    std::vector<double> synFits;
    std::vector<double> mtsng; // size 4: moments for new group
    std::vector<double> mtsog; // size 4: moments for old group
    std::vector<double> mtssg; // size 4: moments for synthetic group

    RawFitDist() : mtsng(4), mtsog(4), mtssg(4) {}
};

// Stores IRT true-score and observed-score equating results
struct RawTruObsEquiv {
    int nRaws = 0;
    double TCCnewMin = 0.0;
    double TCColdMin = 0.0;
    std::vector<double> thetaTru;
    std::vector<double> unroundedEqTru;
    std::vector<double> roundedEqTru;
    std::vector<double> unroundedEqObs;
    std::vector<double> roundedEqObs;
    std::vector<double> mtsTru; // size 4: unrounded true score equivalent moments
    std::vector<double> mtsObs; // size 4: unrounded observed score equivalent moments

    RawTruObsEquiv() : mtsTru(4), mtsObs(4) {}
};

// Contains all input for IRT equating
struct IRT_INPUT {
    char method = ' ';
    std::vector<int> NewFD;
    std::string ItemNewFile;
    std::string ItemOldFile;
    std::string DistNewFile;
    std::string DistOldFile;
    // Using pointers to avoid copying large vector<ItemSpec> objects
    const std::vector<ItemSpec>* NewItems = nullptr;
    const std::vector<ItemSpec>* OldItems = nullptr;
    const IRTstControl* stControl = nullptr;
    RawFitDist* NewForm = nullptr;
    RawFitDist* OldForm = nullptr;
    RawTruObsEquiv* RawEq = nullptr;
};

// Forward declaration for a struct defined in ERutilities
struct PDATA;
struct ERAW_RESULTS;
struct ESS_RESULTS;


// --- Function Declarations ---

// IRT true score equating
void trueScoreEq(const IRTstControl& Handle, const std::vector<ItemSpec>& NewItems,
                 const std::vector<ItemSpec>& OldItems, const std::vector<double>& newScores,
                 std::vector<double>& eqvOld, std::vector<double>& theta,
                 double& newMin, double& OldMin);

double trueScore(const std::vector<ItemSpec>& Items, int n, double theta);

// IRT observed score equating
void IRTmixObsEq(const IRTstControl& Handle, const std::vector<ItemSpec>& NewItems,
                 const std::vector<ItemSpec>& OldItems, double wNew,
                 RawFitDist& newForm, RawFitDist& oldForm, RawTruObsEquiv& RawEq);

void IRTmixObsDist(const std::vector<ItemSpec>& Items, int n, int MaxScrP, int nq,
                   const std::vector<double>& xqpts, const std::vector<double>& xqwts,
                   int& nscr, std::vector<double>& xscr, std::vector<double>& xmarg);

void ObsDistGivenTheta(double theta, const std::vector<ItemSpec>& Items, int n,
                       int& nscr, std::vector<double>& xscr, std::vector<double>& xnew);

void recurs(int mino, int maxo, const std::vector<double>& xold,
            int mitem, const std::vector<double>& iitem, const std::vector<double>& xitem,
            int& minn, int& maxn, std::vector<double>& xnew);

// Functions for different IRT models
double ProbCCC(const ItemSpec& Item, int CatID, double theta);
double PdCCCoverTheta(const ItemSpec& Item, int CatID, double theta);
double Pd3PLoverTheta(int CatID, double theta, double D, double a, double b, double c);
double PdLGRoverTheta(int CatNum, int CatID, double theta, double D, double a, const std::vector<double>& b);
double PdGPCoverTheta(int CatNum, int CatID, double theta, double D, double a, const std::vector<double>& b);
double PdNRMoverTheta(int CatNum, int CatID, double theta, const std::vector<double>& a, const std::vector<double>& c);

// Wrapper and associated print functions
void Wrapper_IRTeq(char design, char method, double w1,
                   const std::string& ItemNewFile, const std::string& ItemOldFile,
                   const std::string& DistNewFile, const std::string& DistOldFile,
                   int* NewFD_raw,
                   IRT_INPUT& irtall, PDATA& pinall, ERAW_RESULTS& r);

void Print_IRTeq(std::ostream& os, const std::string& title, const PDATA& pinall,
                 const ERAW_RESULTS& r, int PrintFiles);

void Print_ESS_QD(std::ostream& os, const std::string& title, const PDATA& pinall, const ESS_RESULTS& s,
                  const RawFitDist& NewForm, const RawFitDist& OldForm);

#endif // IRTEQ_HPP