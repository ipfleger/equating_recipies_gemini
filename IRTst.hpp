/*
  IRTst.hpp
  Header for IRT scale score transformation
*/

#ifndef IRTST_HPP
#define IRTST_HPP

#include <cstdio>
#include <string>

// Note: External dependency 'ERutilities.h' must be included before this file
// in the .cpp files that use these structs, because of the struct definitions.
// To make this header more self-contained, forward declarations are used.
struct PDATA;
struct ERAW_RESULTS;

enum class symmetry { old_scale, new_scale, symmetric };
enum class ModelSpec { l3=1, gr, pc, nr };
enum class LossSpec { mix_ha, mix_sl };

struct ItemSpec {
	int ItemID;
	int CatNum;
	double ScaleConst;
	double *ScoreFunc;
	double *a;
	double *b;
	double *c;
	double *d;
	ModelSpec model;
};

struct CommonItemSpec {
	int NewID;
	int OldID;
	int CatNum;
	double ScaleConst;
	double *ScoreFunc;
	double *Na;
	double *Nb;
	double *Nc;
	double *Nd;
	double *Oa;
	double *Ob;
	double *Oc;
	double *Od;
	ModelSpec model;
};

struct IRTstControl {
	int NewItemNum;
	int OldItemNum;
	int ComItemNum;
	int NewThetaNum;
	int OldThetaNum;
	double NewRawMin;
	double NewRawMax;
	double NewRawInc;
	double OldRawMin;
	double OldRawMax;
	double OldRawInc;
	double *NewThetaValues;
	double *NewThetaWeights;
	double *OldThetaValues;
	double *OldThetaWeights;
};


/* Functions for IRT scale transformation */

ItemSpec* ItemInfoRead(FILE *inf, const char *oldOrnew, IRTstControl *Handle);

CommonItemSpec* ItemPairsRead(FILE *inf, ItemSpec *NewItem, ItemSpec *OldItem, IRTstControl *Handle);

double ProbOld(CommonItemSpec *Item, int CatID, double theta, bool original, double S, double I);
double ProbNew(CommonItemSpec *Item, int CatID, double theta, bool original, double S, double I);
// ... Other ProbDeriv function prototypes
void ScaleTransform(const char *ItemOutF, const char *DistOutF, double slope,
		double intercept, ItemSpec *NewItem, IRTstControl *Handle);

void StHaebara(IRTstControl *Handle, CommonItemSpec *ComItem,
		symmetry SYM, bool FuncStd, double S0, double I0,
		double *slope, double *intercept);
double FuncHaebara(double x[]);
void GradHaebara(double x[], double grad[]);

void StMeanMean(IRTstControl *Handle, CommonItemSpec *ComItem,
		double *slope, double *intercept);

void StMeanSigma(IRTstControl *Handle, CommonItemSpec *ComItem,
		double *slope, double *intercept);

void StItemDeAlloc(ItemSpec *Item, const char *oldOrnew, IRTstControl *Handle);
void StComItemDeAlloc(CommonItemSpec *ComItem, IRTstControl *Handle);
void StContDeAlloc(IRTstControl *Handle);

void StStockingLord(IRTstControl *Handle, CommonItemSpec *ComItem,
		symmetry SYM, bool FuncStd, double S0, double I0,
		double *slope, double *intercept);
double FuncStockingLord(double x[]);
void GradStockingLord(double x[], double grad[]);

void ThetaInfoRead(FILE *inf, const char *oldOrnew, IRTstControl *Handle);

void Wrapper_IRTst(FILE *outf, const std::string& tt,
     const char* ItemNewFile, const char* ItemOldFile, const char* ItemCommonFile,
     const char* DistNewFile, const char* DistOldFile,
     int HA, symmetry HAsym, bool HAfs, double HAs, double HAi,
     int SL, symmetry SLsym, bool SLfs, double SLs, double SLi,
     const char* ST, int PrintFiles);

#endif // IRTST_HPP