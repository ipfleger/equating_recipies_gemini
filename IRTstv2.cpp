/* IRT scale score transformation

This file, which is part of Equating Recipes, is free softrware.
You can distribute it and/or modify it under the terms of the
GNU Lesser General Public License, version 3, as published by 
the Free Software Foundation.

This file is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License, version 3, for more details.

You should have received a copy of the GNU Lesser General Public 
License, version 3, in a ReadMe file distributed along with this
file.  If not, see <http://www.gnu.org/licenses/>   

Copyright 2009 
Center for Advanced Studies in Measurement and Assessment (CASMA)
University of Iowa

   Note: All pointers are 0-offset unless otherwise indicated.
  
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "IRTst.h" 

// Forward declare runerror if not in a shared header
void runerror(const char* error_text);

static struct IRTstControl *ContHandle;      /* global setting to control the method*/
static struct CommonItemSpec *ContComItem;
static enum symmetry ContSym;
static enum OnOff ContFuncStd; 

/***********************************************************************************/           

void ScaleTransform(const char *ItemOutF, const char *DistOutF, double slope,
		double intercept, struct ItemSpec *NewItem, struct IRTstControl *Handle)
/*--------------------------------------------------------------------------------------
  Functionality:
    Use the values of the slope and intercept of the new-to-old transformation
    to convert both the parameter estimates of the new form items and
    the ability points for the new group distribution.
    
    Input:
    
      ItemOutF  The name of a file in which the transformed output for items on
                the new form is saved; The format of output is, in essence, the
                same as that of input, so the output file can be read by the function
                ItemInfoRead without any syntax error.
      DistOutF  The name of a file in which the transformed ability points for the
                new group distribution are saved along with the original weights;
                As with ItemOutF, the output is saved using the format of input for
                the ability distribution, so the file DistOutF can be read by the
                function ThetaInfoRead without problems.
      slope     The value of A (slope) for a chosen linking method
      intercept The value of B (intercept) for a chosen linking method
      NewItem   A pointer to an array of the ItemSpec structure, which is for the new form
      Handle    A pointer to a variable of the IRTstControl structure
    
    Output:
      ItemOutF  The transformed output file for the new form
      DistOutF  The transformed output file for the new group's ability distribution

      * Note: The original input remains intact for both the items and distribution.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
--------------------------------------------------------------------------------------*/   
{
	int i, k;
	FILE *outf;

	/* open the output file for the transformed item parameters */
	outf=std::fopen(ItemOutF, "w");
	if(!outf) 
		runerror("can't open output file for the transformed item parameters\n");

	/* save the transformed information for the new form items */
	std::fprintf(outf, "%d\n", Handle->NewItemNum);

	for (i = 1; i <= Handle->NewItemNum; i++) {
		std::fprintf(outf, "%3d", NewItem[i].ItemID); /* write Item ID */	
		switch(NewItem[i].model)
		{
			case l3:
				std::fprintf(outf, "  L3 2 0 1"); /* write model, category number and score functions */
				std::fprintf(outf, " %g\n", NewItem[i].ScaleConst); /* write scaling constant */
				std::fprintf(outf, " %9.5f", NewItem[i].a[2]/slope);    /* write t-a (tranformed-a) */
				std::fprintf(outf, " %9.5f", NewItem[i].b[2]*slope + intercept); /* write t-b */
				std::fprintf(outf, " %9.5f", NewItem[i].c[2]);         /* write original c */        
				break;
			case gr:
				std::fprintf(outf, "  GR %d", NewItem[i].CatNum); /* write model and category number */
				for(k=1; k <= NewItem[i].CatNum; k++) {
					std::fprintf(outf, " %g", NewItem[i].ScoreFunc[k]); /* write scoring function */
				}
				std::fprintf(outf, " %g\n", NewItem[i].ScaleConst); /* write scaling constant */
				std::fprintf(outf, " %9.5f", NewItem[i].a[2]/slope); /* write t-a */
				for(k=2; k <= NewItem[i].CatNum; k++) {
					if(k%10 == 1) std::fprintf(outf, "\n"); // BUG FIX: Was fscanf
					std::fprintf(outf, " %9.5f", NewItem[i].b[k]*slope + intercept); /* write t-b */
				}
				break;
			case pc:
				std::fprintf(outf, "  PC %d", NewItem[i].CatNum); /* write model and category number */
				for(k=1; k <= NewItem[i].CatNum; k++) {
					std::fprintf(outf, " %g", NewItem[i].ScoreFunc[k]); /* write scoring function */
				}
				std::fprintf(outf, " %g\n", NewItem[i].ScaleConst); /* write scaling constant */
				std::fprintf(outf, " %9.5f", NewItem[i].a[2]/slope); /* write t-a */
				std::fprintf(outf, " %9.5f", NewItem[i].b[0]*slope + intercept); /* write t-b */
				for(k=1; k <= NewItem[i].CatNum; k++) {
					if ((k+2)%10 == 1) std::fprintf(outf, "\n");
					std::fprintf(outf, " %9.5f", NewItem[i].d[k]*slope); /* write t-d */
				}
				break;
			case nr:
				std::fprintf(outf, "  NR %d", NewItem[i].CatNum); /* write model and category number */
				for(k=1; k <= NewItem[i].CatNum; k++) {
					std::fprintf(outf, " %g", NewItem[i].ScoreFunc[k]); /* write scoring function */
				}
				std::fprintf(outf, " 1.0\n"); /* write scaling constant */
				for(k=1; k <= NewItem[i].CatNum; k++) {
					if (k > 1 && k%10 == 1) std::fprintf(outf, "\n");
					std::fprintf(outf, " %9.5f", NewItem[i].a[k]/slope); /* write t-a */
				}
				for(k=1; k <= NewItem[i].CatNum; k++) {
					if (k%10 == 1) std::fprintf(outf, "\n");
					/* write t-c */
					std::fprintf(outf, " %9.5f", NewItem[i].c[k] - (intercept/slope)*NewItem[i].a[k]);
				}
				break;
			default:
				break;
		} /* end of switch */
		std::fprintf(outf, "\n");
	} /* end of for loop */
	std::fclose(outf);

	/* open the output file for the transformed distribution for the new group */
	outf=std::fopen(DistOutF, "w");
	if(!outf) 
		runerror("can't open output file for the transformed"
                 " distribution for the new group\n");

	/* save the transformed information for the new group */
	std::fprintf(outf, "%d\n", Handle->NewThetaNum);
	
	for (i = 1; i <= Handle->NewThetaNum; i++) {
		std::fprintf(outf, "%E", Handle->NewThetaValues[i]*slope + intercept);
		std::fprintf(outf, "   %E\n", Handle->NewThetaWeights[i]);
	}
	std::fclose(outf);

	return;
}

/***********************************************************************************/
 
struct ItemSpec *ItemInfoRead(FILE *inf, const char *oldOrnew, struct IRTstControl *Handle)
/*-----------------------------------------------------------------------------------------*/
{
	char buff[3];
	int i, k, ItemNum;
	struct ItemSpec *Item;

	std::fscanf(inf, "%d",   &ItemNum); /* reading the number of items */
	if (std::strcmp(oldOrnew, "old") == 0) Handle->OldItemNum = ItemNum;
	else                              Handle->NewItemNum = ItemNum;

	Item = static_cast<struct ItemSpec *>(std::malloc((ItemNum+1)*sizeof(struct ItemSpec)));
	if (!Item) 
		runerror("memory allocation failure 1 in ItemInfoRead()\n");

	for (i = 1; i <= ItemNum; i++) {
		std::fscanf(inf, "%d", &(Item[i].ItemID)); /* read ID */
		
		std::fscanf(inf, "%s", buff);               /* read model */
		if (std::strcmp(buff, "L3")==0 || std::strcmp(buff, "l3")==0)      Item[i].model = l3;
		else if (std::strcmp(buff, "GR")==0 || std::strcmp(buff, "gr")==0) Item[i].model = gr;
		else if (std::strcmp(buff, "PC")==0 || std::strcmp(buff, "pc")==0) Item[i].model = pc;
		else if (std::strcmp(buff, "NR")==0 || std::strcmp(buff, "nr")==0) Item[i].model = nr;
        else runerror("Invalid value for type of IRT model");

		std::fscanf(inf, "%d", &(Item[i].CatNum)); /* read CatNum */
                
		Item[i].ScoreFunc = static_cast<double *>(std::malloc((Item[i].CatNum+1) * sizeof(double)));
		if (!Item[i].ScoreFunc) 
			runerror("memory allocation failure 2 in ItemInfoRead()\n");
                
		for (k = 1; k <= Item[i].CatNum; k++) { /* read scoring function */
			std::fscanf(inf, "%lf", &(Item[i].ScoreFunc[k]));
		}
		
                
		Item[i].a = static_cast<double *>(std::malloc((Item[i].CatNum+1)*sizeof(double)));
		if (!Item[i].a) 
			runerror("memory allocation failure 3 in ItemInfoRead()\n");

		Item[i].b = static_cast<double *>(std::malloc((Item[i].CatNum+1)*sizeof(double)));
		if (!Item[i].b) 
			runerror("memory allocation failure 4 in ItemInfoRead()\n");

		Item[i].c = static_cast<double *>(std::malloc((Item[i].CatNum+1)*sizeof(double)));
		if (!Item[i].c) 
			runerror("memory allocation failure 5 in ItemInfoRead()\n");

		Item[i].d = static_cast<double *>(std::malloc((Item[i].CatNum+1)*sizeof(double)));
		if (!Item[i].d) 
			runerror("memory allocation failure 6 in ItemInfoRead()\n");

                		
		switch(Item[i].model)
		{
			case l3:
				std::fscanf(inf, "%lf", &(Item[i].ScaleConst));   /* read D */
				std::fscanf(inf, "%lf", &(Item[i].a[2]));         /* read a parameter */
				std::fscanf(inf, "%lf", &(Item[i].b[2]));         /* read b parameter */
				std::fscanf(inf, "%lf", &(Item[i].c[2]));         /* read c parameter */        
				break;
			case gr:
				std::fscanf(inf, "%lf", &(Item[i].ScaleConst));   /* read D */
				std::fscanf(inf, "%lf", &(Item[i].a[2]));         /* read a parameter */
				for(k=2; k <= Item[i].CatNum; k++) {
					std::fscanf(inf, "%lf", &(Item[i].b[k])); /* read b parameters */
				}
				break;
			case pc:
				std::fscanf(inf, "%lf", &(Item[i].ScaleConst));   /* read D */
				std::fscanf(inf, "%lf", &(Item[i].a[2]));         /* read a parameter */
				std::fscanf(inf, "%lf", &(Item[i].b[0]));         /* read b parameter */
				for(k=1; k <= Item[i].CatNum; k++) {
					std::fscanf(inf, "%lf", &(Item[i].d[k])); /* read d parameters */
					Item[i].b[k] = Item[i].b[0] - Item[i].d[k];
				}
				break;
			case nr:
				std::fscanf(inf, "%lf", &(Item[i].ScaleConst));   /* read D = 1.0 */
				if (Item[i].ScaleConst != 1.0) Item[i].ScaleConst = 1.0;
				for(k=1; k <= Item[i].CatNum; k++) {
					std::fscanf(inf, "%lf", &(Item[i].a[k])); /* read a parameters */
				}
				for(k=1; k <= Item[i].CatNum; k++) {
					std::fscanf(inf, "%lf", &(Item[i].c[k])); /* read c parameters */
				}
				break;
			default:
				break;
		} /* end of switch */
	} /* end of for loop */      
	return Item;
}

/***********************************************************************************/
 
struct CommonItemSpec *ItemPairsRead(FILE *inf, struct ItemSpec *NewItem,
		struct ItemSpec *OldItem, struct IRTstControl *Handle)
/*------------------------------------------------------------------------------*/
{
	int i, k, l, it1, it2, ItemNum, id_finding_error, Nid, Oid;
	struct CommonItemSpec *ComItem;

	std::fscanf(inf, "%d",   &ItemNum); /* reading the number of common items */
	if (ItemNum <= Handle->NewItemNum && ItemNum <= Handle->OldItemNum)
		Handle->ComItemNum = ItemNum;
	else 
		runerror("number of common items must not exceed that of new or old items");

	ComItem = static_cast<struct CommonItemSpec *>(std::malloc((ItemNum+1)*sizeof(struct CommonItemSpec)));
	if (!ComItem) 
		runerror("memory allocation failure 1 in ItemPairsRead()");
        
	for (i = 1; i <= ItemNum; i++) {
		std::fscanf(inf, "%d", &it1); /* read new item ID */
		std::fscanf(inf, "%d", &it2); /* read old item ID */
		
		ComItem[i].NewID = it1;
		ComItem[i].OldID = it2;

       /* finding memory indices for it1 and it2 */
		id_finding_error = 0;
		for(l = 1; l <= Handle->NewItemNum; l++) {
			if (NewItem[l].ItemID == it1) {
				Nid = l;
				id_finding_error ++;
			}
		}
		if (id_finding_error != 1) {
			std::fprintf(stderr, "matching error in common item %d: new item's ID cannot be found or is used more than once\n", i);
			runerror("Check the input data for common items\n");
		}

		id_finding_error = 0;
		for(l = 1; l <= Handle->OldItemNum; l++) {
			if (OldItem[l].ItemID == it2) {
				Oid = l;
				id_finding_error ++;
			}
		}
		if (id_finding_error != 1) {
			std::fprintf(stderr, "matching error in common item %d: old item's ID cannot be found or is used more than once\n", i);
			runerror("Check the input data for common items\n");
		}

		/* minimum checking for the item matching */
		if (NewItem[Nid].model == OldItem[Oid].model)
			ComItem[i].model = NewItem[Nid].model;
		else {
			std::fprintf(stderr, "matching error in IRT model between new item %d and old item %d\n", it1, it2);
			runerror("Check the input data for new and old items\n");
		}

		if (NewItem[Nid].CatNum == OldItem[Oid].CatNum)
			ComItem[i].CatNum = NewItem[Nid].CatNum;
		else {
			std::fprintf(stderr, "matching error in category number between new item %d and old item %d\n", it1, it2);
			runerror("Check the input data for new and old items\n");
		}
		
		if (NewItem[Nid].ScaleConst == OldItem[Oid].ScaleConst)
			ComItem[i].ScaleConst = NewItem[Nid].ScaleConst;
		else {
			std::fprintf(stderr, "matching error in scoring function between new item %d and old item %d\n", it1, it2);
			runerror("Check the input data for new and old items\n");
		}

		/* memory allocation for scoring function */                
		ComItem[i].ScoreFunc = static_cast<double *>(std::malloc((ComItem[i].CatNum+1)*sizeof(double)));
		if (!ComItem[i].ScoreFunc) 
			runerror("memory allocation failure 2 in ItemPairsRead()");
                
		for (k = 1; k <= ComItem[i].CatNum; k++) { /* assigning scoring function */
			ComItem[i].ScoreFunc[k] = NewItem[Nid].ScoreFunc[k];
		}
		
		/* memory allocation for item parameters from new and old forms */
		ComItem[i].Na = static_cast<double *>(std::malloc((ComItem[i].CatNum+1)*sizeof(double)));
		if (!ComItem[i].Na) 
			runerror("memory allocation failure 3 in ItemPairsRead()");

		ComItem[i].Nb = static_cast<double *>(std::malloc((ComItem[i].CatNum+1)*sizeof(double)));
		if (!ComItem[i].Nb) 
			runerror("memory allocation failure 4 in ItemPairsRead()");

		ComItem[i].Nc = static_cast<double *>(std::malloc((ComItem[i].CatNum+1)*sizeof(double)));
		if (!ComItem[i].Nc) 
			runerror("memory allocation failure 5 in ItemPairsRead()");

		ComItem[i].Nd = static_cast<double *>(std::malloc((ComItem[i].CatNum+1)*sizeof(double)));
		if (!ComItem[i].Nd) 
			runerror("memory allocation failure 6 in ItemPairsRead()");
                
		ComItem[i].Oa = static_cast<double *>(std::malloc((ComItem[i].CatNum+1)*sizeof(double)));
		if (!ComItem[i].Oa) 
			runerror("memory allocation failure 7 in ItemPairsRead()");

		ComItem[i].Ob = static_cast<double *>(std::malloc((ComItem[i].CatNum+1)*sizeof(double)));
		if (!ComItem[i].Ob) 
			runerror("memory allocation failure 8 in ItemPairsRead()");

		ComItem[i].Oc = static_cast<double *>(std::malloc((ComItem[i].CatNum+1)*sizeof(double)));
		if (!ComItem[i].Oc) 
			runerror("memory allocation failure 9 in ItemPairsRead()");

		ComItem[i].Od = static_cast<double *>(std::malloc((ComItem[i].CatNum+1)*sizeof(double)));
		if (!ComItem[i].Od) 
			runerror("memory allocation failure 10 in ItemPairsRead()");
                		
		switch(ComItem[i].model)
		{
			case l3:
				ComItem[i].Na[2] = NewItem[Nid].a[2];
				ComItem[i].Nb[2] = NewItem[Nid].b[2];
				ComItem[i].Nc[2] = NewItem[Nid].c[2];
				ComItem[i].Oa[2] = OldItem[Oid].a[2];
				ComItem[i].Ob[2] = OldItem[Oid].b[2];
				ComItem[i].Oc[2] = OldItem[Oid].c[2];
				break;
			case gr:
				ComItem[i].Na[2] = NewItem[Nid].a[2];
				ComItem[i].Oa[2] = OldItem[Oid].a[2];
				for (k=2; k <= ComItem[i].CatNum; k++) {
					ComItem[i].Nb[k] = NewItem[Nid].b[k];
					ComItem[i].Ob[k] = OldItem[Oid].b[k];
				}
				break;
			case pc:
				ComItem[i].Na[2] = NewItem[Nid].a[2];
				ComItem[i].Nb[0] = NewItem[Nid].b[0];
				ComItem[i].Oa[2] = OldItem[Oid].a[2];
				ComItem[i].Ob[0] = OldItem[Oid].b[0];
				for (k=1; k <= ComItem[i].CatNum; k++) {
					ComItem[i].Nb[k] = NewItem[Nid].b[k];
					ComItem[i].Nd[k] = NewItem[Nid].d[k];
					ComItem[i].Ob[k] = OldItem[Oid].b[k];
					ComItem[i].Od[k] = OldItem[Oid].d[k];
				}
				break;
			case nr:
				for (k=1; k <= ComItem[i].CatNum; k++) {
					ComItem[i].Na[k] = NewItem[Nid].a[k];
					ComItem[i].Nc[k] = NewItem[Nid].c[k];
					ComItem[i].Oa[k] = OldItem[Oid].a[k];
					ComItem[i].Oc[k] = OldItem[Oid].c[k];
				}
				break;
			default:
				break;
		} /* end of switch */
	} /* end of for loop */      
	return ComItem;
}

/***********************************************************************************/

double ProbOld(struct CommonItemSpec *Item, int CatID, double theta, enum OnOff original,
		double S, double I)
/*------------------------------------------------------------------------------*/
{
	double prob = 0.0;

	switch(Item->model)
	{
		case l3:
			if (original == on)
				prob = Prob3PL(CatID, theta, Item->ScaleConst,
					Item->Oa[2], Item->Ob[2], Item->Oc[2], "old", 1, 0);
			else
				prob = Prob3PL(CatID, theta, Item->ScaleConst,
					Item->Na[2], Item->Nb[2], Item->Nc[2], "old", S, I);
			break;
		case gr:
			if (original == on)
				prob = ProbLGR(Item->CatNum, CatID, theta, Item->ScaleConst,
					Item->Oa[2], Item->Ob, "old", 1, 0);
			else
				prob = ProbLGR(Item->CatNum, CatID, theta, Item->ScaleConst,
					Item->Na[2], Item->Nb, "old", S, I);
			break;
		case pc:
			if (original == on)
				prob = ProbGPC(Item->CatNum, CatID, theta, Item->ScaleConst,
					Item->Oa[2], Item->Ob, "old", 1, 0);
			else
				prob = ProbGPC(Item->CatNum, CatID, theta, Item->ScaleConst,
					Item->Na[2], Item->Nb, "old", S, I);
			break;
		case nr:
			if (original == on)
				prob = ProbNRM(Item->CatNum, CatID, theta,
					Item->Oa, Item->Oc, "old", 1, 0);
			else
				prob = ProbNRM(Item->CatNum, CatID, theta,
					Item->Na, Item->Nc, "old", S, I);
			break;
		default:
			break;
	}
	return prob;
}

/***********************************************************************************/

double ProbNew(struct CommonItemSpec *Item, int CatID, double theta, enum OnOff original,
		double S, double I)
/*------------------------------------------------------------------------------*/
{
	double prob = 0.0;

	switch(Item->model)
	{
		case l3:
			if (original == on)
				prob = Prob3PL(CatID, theta, Item->ScaleConst,
					Item->Na[2], Item->Nb[2], Item->Nc[2], "new", 1, 0);
			else
				prob = Prob3PL(CatID, theta, Item->ScaleConst,
					Item->Oa[2], Item->Ob[2], Item->Oc[2], "new", S, I);
			break;
		case gr:
			if (original == on)
				prob = ProbLGR(Item->CatNum, CatID, theta, Item->ScaleConst,
					Item->Na[2], Item->Nb, "new", 1, 0);
			else
				prob = ProbLGR(Item->CatNum, CatID, theta, Item->ScaleConst,
					Item->Oa[2], Item->Ob, "new", S, I);
			break;
		case pc:
			if (original == on)
				prob = ProbGPC(Item->CatNum, CatID, theta, Item->ScaleConst,
					Item->Na[2], Item->Nb, "new", 1, 0);
			else
				prob = ProbGPC(Item->CatNum, CatID, theta, Item->ScaleConst,
                        		Item->Oa[2], Item->Ob, "new", S, I);
			break;
		case nr:
			if (original == on)
				prob = ProbNRM(Item->CatNum, CatID, theta,
					Item->Na, Item->Nc, "new", 1, 0);
			else
				prob = ProbNRM(Item->CatNum, CatID, theta,
					Item->Oa, Item->Oc, "new", S, I);
			break;
		default:
			break;
	}
	return prob;
}

/***********************************************************************************/

double PdOldOverS(struct CommonItemSpec *Item, int CatID, double theta, double S, double I)
/*------------------------------------------------------------------------------*/
{
	double pd = 0.0;

	switch(Item->model)
	{
		case l3:
			pd = Pd3PLOldOverS(CatID, theta, Item->ScaleConst,
				Item->Na[2], Item->Nb[2], Item->Nc[2], S, I);
			break;
		case gr:
			pd = PdLGROldOverS(Item->CatNum, CatID, theta, Item->ScaleConst,
                                Item->Na[2], Item->Nb, S, I);
			break;
		case pc:
			pd = PdGPCOldOverS(Item->CatNum, CatID, theta, Item->ScaleConst,
                                Item->Na[2], Item->Nb, S, I);
			break;
		case nr:
			pd = PdNRMOldOverS(Item->CatNum, CatID, theta,
				Item->Na, Item->Nc, S, I);
			break;
		default:
			break;
	}
	return pd;
}

/***********************************************************************************/

double PdOldOverI(struct CommonItemSpec *Item, int CatID, double theta, double S, double I)
/*------------------------------------------------------------------------------*/
{
	double pd = 0.0;

	switch(Item->model)
	{
		case l3:
			pd = Pd3PLOldOverI(CatID, theta, Item->ScaleConst,
				Item->Na[2], Item->Nb[2], Item->Nc[2], S, I);
			break;
		case gr:
			pd = PdLGROldOverI(Item->CatNum, CatID, theta, Item->ScaleConst,
				Item->Na[2], Item->Nb, S, I);
			break;
		case pc:
			pd = PdGPCOldOverI(Item->CatNum, CatID, theta, Item->ScaleConst,
				Item->Na[2], Item->Nb, S, I);
			break;
		case nr:
			pd = PdNRMOldOverI(Item->CatNum, CatID, theta,
				Item->Na, Item->Nc, S, I);
			break;
		default:
			break;
	}
	return pd;
}


/***********************************************************************************/

double PdNewOverS(struct CommonItemSpec *Item, int CatID, double theta, double S, double I)
/*------------------------------------------------------------------------------*/
{
	double pd = 0.0;
   
	switch(Item->model)
	{
		case l3:
			pd = Pd3PLNewOverS(CatID, theta, Item->ScaleConst,
				Item->Oa[2], Item->Ob[2], Item->Oc[2], S, I);
			break;
		case gr:
			pd = PdLGRNewOverS(Item->CatNum, CatID, theta, Item->ScaleConst,
				Item->Oa[2], Item->Ob, S, I);
			break;
		case pc:
			pd = PdGPCNewOverS(Item->CatNum, CatID, theta, Item->ScaleConst,
				Item->Oa[2], Item->Ob, S, I);
			break;
		case nr:
			pd = PdNRMNewOverS(Item->CatNum, CatID, theta,
				Item->Oa, Item->Oc, S, I);
			break;
		default:
			break;
	}
	return pd;
}

/***********************************************************************************/

double PdNewOverI(struct CommonItemSpec *Item, int CatID, double theta, double S, double I)
/*------------------------------------------------------------------------------*/
{
	double pd = 0.0;

	switch(Item->model)
	{
		case l3:
			pd = Pd3PLNewOverI(CatID, theta, Item->ScaleConst,
				Item->Oa[2], Item->Ob[2], Item->Oc[2], S, I);
			break;
		case gr:
			pd = PdLGRNewOverI(Item->CatNum, CatID, theta, Item->ScaleConst,
				Item->Oa[2], Item->Ob, S, I);
			break;
		case pc:
			pd = PdGPCNewOverI(Item->CatNum, CatID, theta, Item->ScaleConst,
				Item->Oa[2], Item->Ob, S, I);
			break;
		case nr:
			pd = PdNRMNewOverI(Item->CatNum, CatID, theta,
				Item->Oa, Item->Oc, S, I);
			break;
		default:
			break;
	}
	return pd;
}

double Prob3PL(int CatID, double theta, double D, double a, double b, double c,
		char scale[], double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Calculate the value of the characteristic curve for the 3PL model.
    Input: 
	CatID: response category, 1 or 2 (1 for incorrect; 2 for correct)
	theta: ability value
	D: scaling constant (typically 1.7)
	a, b, c: item parameters
	scale: "old" or "new" ability scale, on which item and ability parameter
		estimates are placed on.
	S: slope of the linear tranformation
	I: intercept of the linear transformation

   Output:
	Return either original or transformed probability
	For the original probability, S = 1 and I = 0 with parameters on the
	reference scale.
	For the transformed probability, S and I with parameters on the
	transformed scale.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	double as, bs, cs; /* parameters from new-to-old transformation */
	double ar, br, cr; /* parameters from old-to-new transformation */
	double devs, devr;

	if (std::strcmp(scale, "old") == 0) {
		/* new-to-old scale transformation */
		as = a/S;
		bs = S*b + I;
		cs = c;
		devs = D*as*(theta-bs);
   		return ( cs*(CatID-1) + (1.0-cs)*std::exp(-devs*(2-CatID))/(1.0+std::exp(-devs)) );
	}
	else {
		/* old-to-new scale transformation */
		ar = S*a;
		br = (b-I)/S;
		cr = c;
		devr = D*ar*(theta - br);
		return ( cr*(CatID-1) + (1.0-cr)*std::exp(-devr*(2-CatID))/(1.0+std::exp(-devr)) );
	}
}

/***********************************************************************************/

double Pd3PLOldOverS(int CatID, double theta, double D, double na, double nb, double nc,
		double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Under the 3PL model, calculate partial derivative of P* (from new-to-old
    transformation) with respect to S.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	double as, cs;
	double ps, qs; /* Probability with transformed item parameters */
	double ps_over_S;

	/* Scale Transformation */
	as = na / S;
	cs = nc;

	ps = Prob3PL(2, theta, D, na, nb, nc, "old", S, I);
	qs = 1.0 - ps;
	ps_over_S = -D*as*( (theta - I)/S )*(ps - cs)*qs / (1.0 - cs);
   
	if (CatID == 1) return -ps_over_S;
	else return ps_over_S;
}

/***********************************************************************************/

double Pd3PLOldOverI(int CatID, double theta, double D, double na, double nb, double nc,
		double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Under the 3PL model, calculate partial derivative of P* (from new-to-old
    transformation) with respect to I.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	double as, cs;
	double ps, qs; /* Probability with transformed item parameters */
	double ps_over_I;

	/* Scale Transformation */
	as = na / S;
	cs = nc;

	ps = Prob3PL(2, theta, D, na, nb, nc, "old", S, I);
	qs = 1.0 - ps;
	ps_over_I = -D*as*(ps - cs)*qs / (1.0 - cs);
   
	if (CatID == 1) return -ps_over_I;
	else return ps_over_I;
}

/***********************************************************************************/

double Pd3PLNewOverS(int CatID, double theta, double D, double oa, double ob, double oc,
		double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Under the 3PL model, calculate partial derivative of P# (from old-to-new
    transformation) with respect to S.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
	double cr;
	double pr, qr; /* Probability with transformed item parameter estimates */
	double pr_over_S;

	/* Scale Transformation */
	cr = oc;

	pr = Prob3PL(2, theta, D, oa, ob, oc, "new", S, I);
	qr = 1.0 - pr;

	pr_over_S = D*oa*theta*(pr - cr)*qr / (1.0 - cr);
   
	if (CatID == 1) return -pr_over_S;
	else return pr_over_S;
}

/***********************************************************************************/

double Pd3PLNewOverI(int CatID, double theta, double D, double oa, double ob, double oc,
		double S, double I)
/*------------------------------------------------------------------------------
  Functionality:
    Under the 3PL model, calculate partial derivative of P# (from old-to-new
    transformation) with respect to I.

  Author: Seonghoon Kim
  Date of last revision 9/25/08
------------------------------------------------------------------------------*/
{
   double cr;
   double pr, qr; /* Probability with transformed item parameter estimates */
   double pr_over_I;

	/* Scale Transformation */
	cr = oc;

	pr = Prob3PL(2, theta, D, oa, ob, oc, "new", S, I);
	qr = 1.0 - pr;

	pr_over_I = D*oa*(pr - cr)*qr / (1.0 - cr);

	if (CatID == 1) return -pr_over_I;
	else return pr_over_I;
}

/***********************************************************************************/

double CumProbLGR(int CatID, double theta, double D, double a, double b,
		char scale[], double S, double I)
/*------------------------------------------------------------------------------*/
{
	double as, bs, ar, br;
	double cprob;

	if (CatID == 1) cprob = 1.0;
	else {
		if (std::strcmp(scale, "old") == 0) {
			/* new-to-old transformation */
			as = a/S;
			bs = S*b + I;
			cprob = 1.0/(1.0 + std::exp(-D*as*(theta - bs)));
		}
		else {
			/* old-to-new transformation */
			ar = S*a;
			br = (b - I)/S;
			cprob = 1.0/(1.0 + std::exp(-D*ar*(theta - br)));
		}
	}
	return cprob;
}

/***********************************************************************************/

double ProbLGR(int CatNum, int CatID, double theta, double D, double a, double b[],
               char scale[], double S, double I)
/*------------------------------------------------------------------------------*/
{
	double pre_cp, pos_cp;
  
	if (std::strcmp(scale, "old") == 0) {
		pre_cp = CumProbLGR(CatID, theta, D, a, b[CatID], "old", S, I);
		if (CatID < CatNum)
			pos_cp = CumProbLGR(CatID+1, theta, D, a, b[CatID+1], "old", S, I);
		else
			pos_cp = 0.0;
	}
	else {
		pre_cp = CumProbLGR(CatID, theta, D, a, b[CatID], "new", S, I);
		if (CatID < CatNum) 
			pos_cp = CumProbLGR(CatID+1, theta, D, a, b[CatID+1], "new", S, I);
		else 
			pos_cp = 0.0;
	}
	return (pre_cp - pos_cp);
}

/***********************************************************************************/

double PdLGROldOverS(int CatNum, int CatID, double theta, double D, double na, double nb[],
		double S, double I)
/*------------------------------------------------------------------------------*/
{
	double as;
	double pre_cps, pos_cps;
	double pre_cps_over_S, pos_cps_over_S;
	double ps_over_S;

	as = na/S;

	if (CatID == 1) {
		pos_cps = CumProbLGR(CatID+1, theta, D, na, nb[CatID+1], "old", S, I);
		pos_cps_over_S = -D*as*((theta - I)/S)*pos_cps*(1.0 - pos_cps);
		ps_over_S = 0.0 - pos_cps_over_S;
	}
	else {
		if (CatID < CatNum) {
			pre_cps = CumProbLGR(CatID, theta, D, na, nb[CatID], "old", S, I);
			pos_cps = CumProbLGR(CatID+1, theta, D, na, nb[CatID+1], "old", S, I);
			pre_cps_over_S = -D*as*((theta - I)/S)*pre_cps*(1.0 - pre_cps);
			pos_cps_over_S = -D*as*((theta - I)/S)*pos_cps*(1.0 - pos_cps);
			ps_over_S = pre_cps_over_S - pos_cps_over_S;
		}
		else { /* CatId == CatNum */
			pre_cps = CumProbLGR(CatID, theta, D, na, nb[CatID], "old", S, I);
			pre_cps_over_S = -D*as*((theta - I)/S)*pre_cps*(1.0 - pre_cps);
			ps_over_S = pre_cps_over_S - 0.0;
		}
	}
	return ps_over_S;
}

/***********************************************************************************/

double PdLGROldOverI(int CatNum, int CatID, double theta, double D, double na, double nb[],
		double S, double I)
/*------------------------------------------------------------------------------*/
{
	double as;
	double pre_cps, pos_cps;
	double pre_cps_over_I, pos_cps_over_I;
	double ps_over_I;

	as = na / S;

	if (CatID == 1) {
		pos_cps = CumProbLGR(CatID+1, theta, D, na, nb[CatID+1], "old", S, I);
		pos_cps_over_I = -D*as*pos_cps*(1.0 - pos_cps);
		ps_over_I = 0.0 - pos_cps_over_I;
	}
	else {
		if (CatID < CatNum) {
			pre_cps = CumProbLGR(CatID, theta, D, na, nb[CatID], "old", S, I);
			pos_cps = CumProbLGR(CatID+1, theta, D, na, nb[CatID+1], "old", S, I);
			pre_cps_over_I = -D*as*pre_cps*(1.0 - pre_cps);
			pos_cps_over_I = -D*as*pos_cps*(1.0 - pos_cps);
			ps_over_I = pre_cps_over_I - pos_cps_over_I;
		}
		else { /* CatID == CatNum */
			pre_cps = CumProbLGR(CatID, theta, D, na, nb[CatID], "old", S, I);
			pre_cps_over_I = -D*as*pre_cps*(1.0 - pre_cps);
			ps_over_I = pre_cps_over_I - 0.0;
		}
	}
	return ps_over_I;
}

/***********************************************************************************/

double PdLGRNewOverS(int CatNum, int CatID, double theta, double D, double oa, double ob[],
		double S, double I)
/*------------------------------------------------------------------------------*/
{
	double pre_cpr, pos_cpr;
	double pre_cpr_over_S, pos_cpr_over_S;
	double pr_over_S;
  
	if (CatID == 1) {
		pos_cpr = CumProbLGR(CatID+1, theta, D, oa, ob[CatID+1], "new", S, I);
		pos_cpr_over_S = D*oa*theta*pos_cpr*(1.0 - pos_cpr);
		pr_over_S = 0.0 - pos_cpr_over_S;
	}
	else {
		if (CatID < CatNum) {
			pre_cpr = CumProbLGR(CatID, theta, D, oa, ob[CatID], "new", S, I);
			pos_cpr = CumProbLGR(CatID+1, theta, D, oa, ob[CatID+1], "new", S, I);
			pre_cpr_over_S = D*oa*theta*pre_cpr*(1.0 - pre_cpr);
			pos_cpr_over_S = D*oa*theta*pos_cpr*(1.0 - pos_cpr);
			pr_over_S = pre_cpr_over_S - pos_cpr_over_S;
		}
		else { /* CatID == CatNum */
			pre_cpr = CumProbLGR(CatID, theta, D, oa, ob[CatID], "new", S, I);
			pre_cpr_over_S = D*oa*theta*pre_cpr*(1.0 - pre_cpr);
			pr_over_S = pre_cpr_over_S - 0.0;
		}
	}
	return pr_over_S;
}

/***********************************************************************************/

double PdLGRNewOverI(int CatNum, int CatID, double theta, double D, double oa, double ob[],
		double S, double I)
/*------------------------------------------------------------------------------*/
{
	double pre_cpr, pos_cpr;
	double pre_cpr_over_I, pos_cpr_over_I;
	double pr_over_I;
  
	if (CatID == 1) {
		pos_cpr = CumProbLGR(CatID+1, theta, D, oa, ob[CatID+1], "new", S, I);
		pos_cpr_over_I = D*oa*pos_cpr*(1.0 - pos_cpr);
		pr_over_I = 0.0 - pos_cpr_over_I;
	}
	else {
		if (CatID < CatNum) {
			pre_cpr = CumProbLGR(CatID, theta, D, oa, ob[CatID], "new", S, I);
			pos_cpr = CumProbLGR(CatID+1, theta, D, oa, ob[CatID+1], "new", S, I);
			pre_cpr_over_I = D*oa*pre_cpr*(1.0 - pre_cpr);
			pos_cpr_over_I = D*oa*pos_cpr*(1.0 - pos_cpr);
			pr_over_I = pre_cpr_over_I - pos_cpr_over_I;
		}
		else { /* CatID == CatNum */
			pre_cpr = CumProbLGR(CatID, theta, D, oa, ob[CatID], "new", S, I);
			pre_cpr_over_I = D*oa*pre_cpr*(1.0 - pre_cpr);
			pr_over_I = pre_cpr_over_I - 0.0;
		}
	}
	return pr_over_I;
}

/***********************************************************************************/

double ProbGPC(int CatNum, int CatID, double theta, double D, double a, double b[],
		char scale[], double S, double I)
/*------------------------------------------------------------------------------*/
{
	int k, l;
	double vjs = 0.0, vjr = 0.0;
	double a_sum, b_sum, as, bs, ar, br;

	if (std::strcmp(scale, "old") == 0) {
		for (k = 1; k <= CatNum ; k++) {
			a_sum = static_cast<double>(k)*D*a;
			b_sum = 0.0;
			for (l = 2; l <= k ; l++) {  /* b[1] = 0 */
				b_sum += b[l];
			}
			b_sum *= -D*a;

			/* new-to-old transformation */
			as = a_sum/S;
			bs = b_sum - (I/S)*a_sum;
			vjs += std::exp(as*theta + bs);  
		}
		a_sum = static_cast<double>(CatID)*D*a;
		b_sum = 0.0;
		for (l = 2; l <= CatID ; l++) {
			b_sum += b[l];
		}
		b_sum *= -D*a;

		/* new-to-old transformation */
		as = a_sum/S;
		bs = b_sum - (I/S)*a_sum;
		return std::exp(as*theta + bs)/vjs;
	}
	else {
		for (k = 1; k <= CatNum ; k++) {
			a_sum = static_cast<double>(k)*D*a;
			b_sum = 0.0;
			for (l = 2; l <= k ; l++) {
				b_sum += b[l];
			}
			b_sum *= -D*a;
			
			/* old-to-new transformation */
			ar = S*a_sum;
			br = b_sum + I*a_sum;
			vjr += std::exp(ar*theta + br);  
		}
		a_sum = static_cast<double>(CatID)*D*a;
		b_sum = 0.0;
		for (l = 2; l <= CatID ; l++) {
			b_sum += b[l];
		}
		b_sum *= -D*a;

		/* old-to-new transformation */
		ar = S*a_sum;
		br = b_sum + I*a_sum;
		return std::exp(ar*theta + br)/vjr;  
	}
}

/***********************************************************************************/

double PdGPCOldOverS(int CatNum, int CatID, double theta, double D, double na, double nb[],
		double S, double I)
/*------------------------------------------------------------------------------*/
{
	int k;
	double as, ps;
	double na_sum, as_ps_sum = 0.0;
	double ps_over_S;

	for (k = 1; k <= CatNum; k++) {
		na_sum = static_cast<double>(k)*D*na;
		as = na_sum/S;      /* new-to-old transformation */
		ps = ProbGPC(CatNum, k, theta, D, na, nb, "old", S, I);
		as_ps_sum += as * ps;
	}
	   
	na_sum = static_cast<double>(CatID)*D*na; /* for the category in question */
	as = na_sum/S;       /* new-to-old transformation */
	ps = ProbGPC(CatNum, CatID, theta, D, na, nb, "old", S, I);

	ps_over_S = -ps *((theta - I)/S)*(as - as_ps_sum);
	return ps_over_S;
}

/***********************************************************************************/

double PdGPCOldOverI(int CatNum, int CatID, double theta, double D, double na, double nb[],
		double S, double I)
/*------------------------------------------------------------------------------*/
{
	int k;
	double as, ps;
	double na_sum, as_ps_sum = 0.0;
	double ps_over_I;

	for (k = 1; k <= CatNum; k++) {
		na_sum = static_cast<double>(k)*D*na;
		as = na_sum/S;      /* new-to-old transformation */
		ps = ProbGPC(CatNum, k, theta, D, na, nb, "old", S, I);
		as_ps_sum += as * ps;
	}

	na_sum = static_cast<double>(CatID)*D*na; /* for the category in question */
	as = na_sum/S;       /* new-to-old transformation */
	ps = ProbGPC(CatNum, CatID, theta, D, na, nb, "old", S, I);

	ps_over_I = -ps * (as - as_ps_sum);
	return ps_over_I;
}

/***********************************************************************************/

double PdGPCNewOverS(int CatNum, int CatID, double theta, double D, double oa, double ob[],
		double S, double I)
/*------------------------------------------------------------------------------*/
{
	int k;
	double pr;
	double oa_sum, oa_pr_sum = 0.0;
	double pr_over_S;

	for (k = 1; k <= CatNum; k++) {
		oa_sum = static_cast<double>(k)*D*oa;
		pr = ProbGPC(CatNum, k, theta, D, oa, ob, "new", S, I);
		oa_pr_sum += oa_sum * pr;
	}
	oa_sum = static_cast<double>(CatID)*D*oa; /* for the category in question */
	pr = ProbGPC(CatNum, CatID, theta, D, oa, ob, "new", S, I);
	pr_over_S = pr * theta *(oa_sum - oa_pr_sum);
	return pr_over_S;
}

/***********************************************************************************/

double PdGPCNewOverI(int CatNum, int CatID, double theta, double D, double oa, double ob[],
		double S, double I)
/*------------------------------------------------------------------------------*/
{
	int k;
	double pr;
	double oa_sum, oa_pr_sum = 0.0;
	double pr_over_I;

	for (k = 1; k <= CatNum ; k++) {
		oa_sum = static_cast<double>(k)*D*oa;
		pr = ProbGPC(CatNum, k, theta, D, oa, ob, "new", S, I);
		oa_pr_sum += oa_sum * pr;
	}
	oa_sum = static_cast<double>(CatID)*D*oa; /* for the category in question */
	pr = ProbGPC(CatNum, CatID, theta, D, oa, ob, "new", S, I);
	pr_over_I = pr * (oa_sum - oa_pr_sum);
	return pr_over_I;
}

/***********************************************************************************/

double ProbNRM(int CatNum, int CatID, double theta, double a[], double c[],
		char scale[], double S, double I)
/*------------------------------------------------------------------------------*/
{
	int k;
	double vjs = 0.0, vjr = 0.0;
	double as, cs, ar, cr;

	if(std::strcmp(scale, "old") == 0) {
		for (k = 1; k <= CatNum; k++) {
			as = a[k]/S;
			cs = c[k] - (I/S)*a[k];
			vjs += std::exp( as*theta + cs );
		}
		as = a[CatID]/S;
		cs = c[CatID] - (I/S)*a[CatID];
		return std::exp(as*theta + cs)/vjs;
	}
	else {
		for (k = 1; k <= CatNum; k++) {
			ar = S*a[k];
			cr = c[k] + I*a[k];
			vjr += std::exp( ar*theta + cr );
		}
		ar = S*a[CatID];
		cr = c[CatID] + I*a[CatID];
		return std::exp(ar*theta + cr)/vjr;
	}
}

/***********************************************************************************/

double PdNRMOldOverS(int CatNum, int CatID, double theta, double na[], double nc[],
		double S, double I)
/*------------------------------------------------------------------------------*/
{
	int k;
	double as, ps;
	double as_ps_sum = 0.0;
	double ps_over_S;

	for (k = 1; k <= CatNum; k++) {
		as = na[k]/S;
		ps = ProbNRM(CatNum, k, theta, na, nc, "old", S, I);
		as_ps_sum += as * ps;
	}   
	as = na[CatID]/S;
	ps = ProbNRM(CatNum, CatID, theta, na, nc, "old", S, I);
	ps_over_S = - ps*((theta - I)/S)*(as - as_ps_sum);
	return ps_over_S;
}

/***********************************************************************************/

double PdNRMOldOverI(int CatNum, int CatID, double theta, double na[], double nc[],
		double S, double I)
/*------------------------------------------------------------------------------*/
{
	int k;
	double as, ps;
	double as_ps_sum = 0.0;
	double ps_over_I;

	for (k = 1; k <= CatNum ; k++) {
		as = na[k]/S;
		ps = ProbNRM(CatNum, k, theta, na, nc, "old", S, I);
		as_ps_sum += as * ps;
	}   
		
	as = na[CatID]/S;
	ps = ProbNRM(CatNum, CatID, theta, na, nc, "old", S, I);
	ps_over_I = -ps * (as - as_ps_sum);
	return ps_over_I;
}

/***********************************************************************************/

double PdNRMNewOverS(int CatNum, int CatID, double theta, double oa[], double oc[],
		double S, double I)
/*------------------------------------------------------------------------------*/
{
	int k;
	double pr;
	double oa_pr_sum = 0.0;
	double pr_over_S;

	for (k = 1; k <= CatNum ; k++) {
		pr = ProbNRM(CatNum, k, theta, oa, oc, "new", S, I);
		oa_pr_sum += oa[k] * pr;
	}   
	pr = ProbNRM(CatNum, CatID, theta, oa, oc, "new", S, I);
	pr_over_S = pr*theta*(oa[CatID] - oa_pr_sum);
	return pr_over_S;
}

/***********************************************************************************/

double PdNRMNewOverI(int CatNum, int CatID, double theta, double oa[], double oc[],
		double S, double I)
/*------------------------------------------------------------------------------*/
{
	int k;
	double pr;
	double oa_pr_sum = 0.0;
	double pr_over_I;

	for (k = 1; k <= CatNum ; k++) {
		pr = ProbNRM(CatNum, k, theta, oa, oc, "new", S, I);
		oa_pr_sum += oa[k] * pr;
	}   
		
	pr = ProbNRM(CatNum, CatID, theta, oa, oc, "new", S, I);
	pr_over_I = pr*(oa[CatID] - oa_pr_sum);
	return pr_over_I;
}          

/***********************************************************************************/

void StHaebara(struct IRTstControl *Handle, struct CommonItemSpec *ComItem,
		enum symmetry SYM, enum OnOff FuncStd, double S0, double I0,
		double *slope, double *intercept)
/*------------------------------------------------------------------------------*/
{
	int n, iter;
	double x[3];
	double ftol=0.0000000001, fret;

	ContHandle = Handle;
	ContComItem = ComItem;
	ContSym = SYM;
	ContFuncStd = FuncStd;

	x[1] = S0;
	x[2] = I0;
	n = 2;
	er_dfpmin(x, n, ftol, &iter, &fret, FuncHaebara, GradHaebara);
	*slope = x[1];
	*intercept = x[2];
	return;
}

/***********************************************************************************/

double FuncHaebara(double x[])
/*------------------------------------------------------------------------------*/
{
	int i, j, k;
	int cat_sum = 0;
	double w1_sum = 0.0, w2_sum = 0.0;
	double sym_f1, sym_f2;
	double theta, th_weight;
	double p_origi, p_trans;
	double func1, func2;
	double func1_sum = 0.0, func2_sum = 0.0;
 
	/* Q1: old scale */
	for (i = 1; i <= ContHandle->OldThetaNum; i++) {
		theta = ContHandle->OldThetaValues[i];
		th_weight = ContHandle->OldThetaWeights[i];
		w1_sum += th_weight;

		func1 = 0.0;
		for (j = 1; j <= ContHandle->ComItemNum; j++) {
			for (k = 1; k <= ContComItem[j].CatNum; k++) {
				p_origi = ProbOld(&ContComItem[j], k, theta, on, x[1], x[2]);
				p_trans = ProbOld(&ContComItem[j], k, theta, off, x[1], x[2]);
				func1 += (p_origi - p_trans)*(p_origi - p_trans);
			}
		}
		func1_sum += func1 * th_weight;
	}

	/* Q2: new scale */
	for (i = 1; i <= ContHandle->NewThetaNum; i++) {
		theta = ContHandle->NewThetaValues[i];
		th_weight = ContHandle->NewThetaWeights[i];
		w2_sum += th_weight;
		
		func2 = 0.0;
		for (j = 1; j <= ContHandle->ComItemNum; j++) {
			for (k = 1; k <= ContComItem[j].CatNum; k++) {
				p_origi = ProbNew(&ContComItem[j], k, theta, on, x[1], x[2]);
				p_trans = ProbNew(&ContComItem[j], k, theta, off, x[1], x[2]);
				func2 += (p_origi - p_trans)*(p_origi - p_trans);
			}
			if (i == 1) cat_sum += ContComItem[j].CatNum;
		}
		func2_sum += func2 * th_weight;
	}

	/* symmetric or non-symmetric setting*/
	if (ContSym == symmetric) {sym_f1 = 1.0; sym_f2 = 1.0;}
	else {
		if (ContSym == old_scale) {sym_f1 = 1.0; sym_f2 = 0.0;}
		else                      {sym_f1 = 0.0; sym_f2 = 1.0;}
	}

	/* function standardization */
	if (ContFuncStd == on) {
		return (sym_f1*func1_sum/(cat_sum*w1_sum) + sym_f2*func2_sum/(cat_sum*w2_sum));
	}
	else {
		return (sym_f1*func1_sum + sym_f2*func2_sum);
	}
}

/***********************************************************************************/

void GradHaebara(double x[], double grad[])
/*------------------------------------------------------------------------------*/
{
	int i, j, k;
	int cat_sum = 0;
	double w1_sum = 0.0, w2_sum = 0.0;
	double sym_f1, sym_f2;
	double theta, th_weight;
	double p_origi, p_trans, ps, pi;
	double ps_f1, pi_f1, ps_f2, pi_f2;
	double ps_f1_sum = 0.0, pi_f1_sum = 0.0;
	double ps_f2_sum = 0.0, pi_f2_sum = 0.0;

	/* Q1: old scale */
	for (i = 1; i <= ContHandle->OldThetaNum; i++) {
		theta = ContHandle->OldThetaValues[i];
		th_weight = ContHandle->OldThetaWeights[i];
		w1_sum += th_weight;

		ps_f1 = 0.0; pi_f1 = 0.0;
		for (j = 1; j <= ContHandle->ComItemNum; j++) {
			for (k = 1; k <= ContComItem[j].CatNum; k++) {
				p_origi = ProbOld(&ContComItem[j], k, theta, on, x[1], x[2]);
				p_trans = ProbOld(&ContComItem[j], k, theta, off, x[1], x[2]);
				ps = PdOldOverS(&ContComItem[j], k, theta, x[1], x[2]);
				pi = PdOldOverI(&ContComItem[j], k, theta, x[1], x[2]);

				ps_f1 += (p_origi - p_trans)*ps;
				pi_f1 += (p_origi - p_trans)*pi;
			}
			if (i == 1) cat_sum += ContComItem[j].CatNum;
		}
		ps_f1_sum += ps_f1 * th_weight;
		pi_f1_sum += pi_f1 * th_weight;
	}

	/* Q2: new scale */

	for (i = 1; i <= ContHandle->NewThetaNum; i++) {
		theta = ContHandle->NewThetaValues[i];
		th_weight = ContHandle->NewThetaWeights[i];
		w2_sum += th_weight;

		ps_f2 = 0.0; pi_f2 = 0.0;
		for (j = 1; j <= ContHandle->ComItemNum; j++) {
			for (k = 1; k <= ContComItem[j].CatNum; k++) {
				p_origi = ProbNew(&ContComItem[j], k, theta, on, x[1], x[2]);
				p_trans = ProbNew(&ContComItem[j], k, theta, off, x[1], x[2]);
				ps = PdNewOverS(&ContComItem[j], k, theta, x[1], x[2]);
				pi = PdNewOverI(&ContComItem[j], k, theta, x[1], x[2]);

				ps_f2 += (p_origi - p_trans)*ps;
				pi_f2 += (p_origi - p_trans)*pi;
			}
		}
		ps_f2_sum += ps_f2 * th_weight;
		pi_f2_sum += pi_f2 * th_weight;
	}

	/* symmetric or non-symmetric setting*/
	if (ContSym == symmetric) {sym_f1 = 1.0; sym_f2 = 1.0;}
	else {
		if (ContSym == old_scale) {sym_f1 = 1.0; sym_f2 = 0.0;}
		else                      {sym_f1 = 0.0; sym_f2 = 1.0;}
	}

	/* function standardization */
	if (ContFuncStd == on) {
		grad[1] = -2.0*(sym_f1*ps_f1_sum/(cat_sum*w1_sum) + sym_f2*ps_f2_sum/(cat_sum*w2_sum));
		grad[2] = -2.0*(sym_f1*pi_f1_sum/(cat_sum*w1_sum) + sym_f2*pi_f2_sum/(cat_sum*w2_sum));
	}
	else {
		grad[1] = -2.0*(sym_f1*ps_f1_sum + sym_f2*ps_f2_sum);
		grad[2] = -2.0*(sym_f1*pi_f1_sum + sym_f2*pi_f2_sum);
	}
	return;
}

/***********************************************************************************/

void StMeanMean(struct IRTstControl *Handle, struct CommonItemSpec *ComItem,
		double *slope, double *intercept)
/*------------------------------------------------------------------------------*/
{
	int j, k, a_num=0, b_num=0, index_a=0, index_b=0;
	double new_mu_a=0.0, old_mu_a=0.0;
	double new_mu_b=0.0, old_mu_b=0.0;
	double *na_vec, *oa_vec, *nb_vec, *ob_vec;

	/* counting valid a- and b-parameters by model */

	for(j=1; j <= Handle->ComItemNum; j++) {
		if(ComItem[j].model==l3 || ComItem[j].model==gr ||ComItem[j].model==pc) {
			a_num++;
			b_num += (ComItem[j].CatNum-1);
		}
		else {
			a_num += ComItem[j].CatNum;
			b_num += ComItem[j].CatNum;
		}
	}

	/* memory allocation */
	na_vec = static_cast<double *>(std::malloc((a_num+1)*sizeof(double)));
	if (!na_vec) 
		runerror("memory allocation failure 1 in StMeanMean()\n");

	oa_vec = static_cast<double *>(std::malloc((a_num+1)*sizeof(double)));
	if (!oa_vec) 
		runerror("memory allocation failure 2 in StMeanMean()\n");

	nb_vec = static_cast<double *>(std::malloc((b_num+1)*sizeof(double)));
	if (!nb_vec)
		runerror("memory allocation failure 3 in StMeanMean()\n");

	ob_vec = static_cast<double *>(std::malloc((b_num+1)*sizeof(double)));
	if (!ob_vec) 
		runerror("memory allocation failure 4 in StMeanMean()\n");


	/* copying valid a- and b-parameters by model into nb_vec and ob_vec */
	for(j=1; j <= Handle->ComItemNum; j++) {
		if(ComItem[j].model==l3 || ComItem[j].model==gr ||ComItem[j].model==pc) {
			index_a++;
			na_vec[index_a] = ComItem[j].Na[2];
			oa_vec[index_a] = ComItem[j].Oa[2];
			for(k=2; k <= ComItem[j].CatNum; k++) {
				index_b++;
				nb_vec[index_b] = ComItem[j].Nb[k];
				ob_vec[index_b] = ComItem[j].Ob[k];
			}
		}
	 	else {
			for(k=1; k <= ComItem[j].CatNum; k++) {
				index_a++;
				index_b++;
				na_vec[index_a] = ComItem[j].Na[k];
				oa_vec[index_a] = ComItem[j].Oa[k];
				nb_vec[index_b] = -(ComItem[j].Nc[k]/ComItem[j].Na[k]);
				ob_vec[index_b] = -(ComItem[j].Oc[k]/ComItem[j].Oa[k]);
			}
		}
	}

	for (j=1; j <= a_num; j++) {
		new_mu_a += na_vec[j];
		old_mu_a += oa_vec[j];
	}
	for (j=1; j <= b_num; j++) {
		new_mu_b += nb_vec[j];
		old_mu_b += ob_vec[j];
	}
	if (a_num != 0 && b_num != 0) {
		new_mu_a /= a_num;
		old_mu_a /= a_num;
		new_mu_b /= b_num;
		old_mu_b /= b_num;

		*slope = new_mu_a/old_mu_a;
		*intercept = old_mu_b - (*slope)*new_mu_b;
	}
	else {
		*slope = 1.0;
		*intercept = 0.0;
	}

	std::free(na_vec);
	std::free(oa_vec);
	std::free(nb_vec);
	std::free(ob_vec);
	return;
}

/***********************************************************************************/

void StMeanSigma(struct IRTstControl *Handle, struct CommonItemSpec *ComItem,
		double *slope, double *intercept)
/*------------------------------------------------------------------------------*/
{
	int j, k, b_num=0, index=0;
	double new_mu_b=0.0, old_mu_b=0.0;
	double new_si_b=0.0, old_si_b=0.0;
	double *nb_vec, *ob_vec;

	/* counting valid b-parameters by model */

	for(j=1; j <= Handle->ComItemNum; j++) {
		if(ComItem[j].model==l3 || ComItem[j].model==pc ||ComItem[j].model==gr) {
			b_num += (ComItem[j].CatNum-1);
		}
		else b_num += ComItem[j].CatNum;
	}

	/* memory allocation */
	nb_vec = static_cast<double *>(std::malloc((b_num+1)*sizeof(double)));
	if (!nb_vec) 
		runerror("memory allocation failure 1 in StMeanSigma()\n");

	ob_vec = static_cast<double *>(std::malloc((b_num+1)*sizeof(double)));
	if (!ob_vec) 
		runerror("memory allocation failure 2 in StMeanSigma()\n");


	/* copying valid b-parameters by model into nb_vec and ob_vec */
	for(j=1; j <= Handle->ComItemNum; j++) {
		if(ComItem[j].model==l3 || ComItem[j].model==pc ||ComItem[j].model==gr) {
			for(k=2; k <= ComItem[j].CatNum; k++) {
				index++;
				nb_vec[index] = ComItem[j].Nb[k];
				ob_vec[index] = ComItem[j].Ob[k];
			}
		}
		else {
			for(k=1; k <= ComItem[j].CatNum; k++) {
				index++;
				nb_vec[index] = -(ComItem[j].Nc[k]/ComItem[j].Na[k]);
				ob_vec[index] = -(ComItem[j].Oc[k]/ComItem[j].Oa[k]);
			}
		}
	}

	for (j=1; j <= b_num; j++) {
		new_mu_b += nb_vec[j];
		old_mu_b += ob_vec[j];
		new_si_b += nb_vec[j]*nb_vec[j];
		old_si_b += ob_vec[j]*ob_vec[j];
	}
	if (b_num != 0) {
		new_mu_b /= b_num;
		old_mu_b /= b_num;
		new_si_b = std::sqrt(new_si_b/b_num - new_mu_b*new_mu_b);
		old_si_b = std::sqrt(old_si_b/b_num - old_mu_b*old_mu_b);

		*slope = old_si_b/new_si_b;
		*intercept = old_mu_b - (*slope)*new_mu_b;
	}
	else {
		*slope = 1.0;
		*intercept = 0.0;
	}

	std::free(nb_vec);
	std::free(ob_vec);
	return;
}

/***********************************************************************************/

void StStockingLord(struct IRTstControl *Handle, struct CommonItemSpec *ComItem,
		enum symmetry SYM, enum OnOff FuncStd, double S0, double I0,
		double *slope, double *intercept)
/*------------------------------------------------------------------------------*/
{
	int n, iter;
	double x[3];
	double ftol=0.0000000001, fret;

	ContHandle = Handle;
	ContComItem = ComItem;
	ContSym = SYM;
	ContFuncStd = FuncStd;

	x[1] = S0;
	x[2] = I0;
	n = 2;

	er_dfpmin(x, n, ftol, &iter, &fret, FuncStockingLord, GradStockingLord); 

	*slope = x[1];
	*intercept = x[2];
	return;
}

/***********************************************************************************/

double FuncStockingLord(double x[])
/*------------------------------------------------------------------------------*/
{
	int i, j, k;
	double w1_sum = 0.0, w2_sum = 0.0;
	double sym_f1, sym_f2;
	double theta, th_weight;
	double p_origi, p_trans, Ujk;
	double func1, func2;
	double func1_sum = 0.0, func2_sum = 0.0;

	/* F1: old scale */
	for (i = 1; i <= ContHandle->OldThetaNum; i++) {
		theta = ContHandle->OldThetaValues[i];
		th_weight = ContHandle->OldThetaWeights[i];
		w1_sum += th_weight;

		func1 = 0.0;
		for (j = 1; j <= ContHandle->ComItemNum; j++) {
			for (k = 1; k <= ContComItem[j].CatNum; k++) {
				Ujk = ContComItem[j].ScoreFunc[k];
				p_origi = ProbOld(&ContComItem[j], k, theta, on, x[1], x[2]);
				p_trans = ProbOld(&ContComItem[j], k, theta, off, x[1], x[2]);

				func1 += Ujk*(p_origi - p_trans);
			}
		}
		func1_sum += (func1 * func1) * th_weight;
	}

	/* F2: new scale */
	for (i = 1; i <= ContHandle->NewThetaNum; i++) {
		theta = ContHandle->NewThetaValues[i];
		th_weight = ContHandle->NewThetaWeights[i];
		w2_sum += th_weight;

		func2 = 0.0;
		for (j = 1; j <= ContHandle->ComItemNum; j++) {
			for (k = 1; k <= ContComItem[j].CatNum; k++) {
				Ujk = ContComItem[j].ScoreFunc[k];
				p_origi = ProbNew(&ContComItem[j], k, theta, on, x[1], x[2]);
				p_trans = ProbNew(&ContComItem[j], k, theta, off, x[1], x[2]);

				func2 += Ujk*(p_origi - p_trans);
			}
		}
		func2_sum += (func2 * func2) * th_weight;
	}

	/* symmetric or non-symmetric setting*/
	if (ContSym == symmetric) {sym_f1 = 1.0; sym_f2 = 1.0;}
	else {
		if (ContSym == old_scale) {sym_f1 = 1.0; sym_f2 = 0.0;}
		else                      {sym_f1 = 0.0; sym_f2 = 1.0;}
	}

	/* function standardization */
	if (ContFuncStd == on) {
		return (sym_f1*func1_sum/w1_sum + sym_f2*func2_sum/w2_sum);
	}
	else {
		return (sym_f1*func1_sum + sym_f2*func2_sum);
	}
}

/***********************************************************************************/

void GradStockingLord(double x[], double grad[])
/*------------------------------------------------------------------------------*/
{
	int i, j, k;
	double w1_sum = 0.0, w2_sum = 0.0;
	double sym_f1, sym_f2;

	double theta, th_weight;
	double p_origi, p_trans, ps, pi, Ujk;
	double func1, func2, ps_sum, pi_sum;
	double ps_f1_sum = 0.0, pi_f1_sum = 0.0;
	double ps_f2_sum = 0.0, pi_f2_sum = 0.0;


	/* F1: old scale */
	for (i = 1; i <= ContHandle->OldThetaNum; i++) {
		theta = ContHandle->OldThetaValues[i];
		th_weight = ContHandle->OldThetaWeights[i];
		w1_sum += th_weight;

		func1 = 0.0;
		ps_sum = 0.0; pi_sum = 0.0;
		for (j = 1; j <= ContHandle->ComItemNum; j++) {
			for (k = 1; k <= ContComItem[j].CatNum; k++) {
				Ujk = ContComItem[j].ScoreFunc[k];
				p_origi = ProbOld(&ContComItem[j], k, theta, on, x[1], x[2]);
				p_trans = ProbOld(&ContComItem[j], k, theta, off, x[1], x[2]);
				ps = PdOldOverS(&ContComItem[j], k, theta, x[1], x[2]);
				pi = PdOldOverI(&ContComItem[j], k, theta, x[1], x[2]);

				func1 += Ujk * (p_origi - p_trans);
				ps_sum += Ujk*ps;
				pi_sum += Ujk*pi;
			}
		}
		ps_f1_sum += func1 * ps_sum * th_weight;
		pi_f1_sum += func1 * pi_sum * th_weight;
	}

	/* F2: new scale */
	for (i = 1; i <= ContHandle->NewThetaNum; i++) {
		theta = ContHandle->NewThetaValues[i];
		th_weight = ContHandle->NewThetaWeights[i];
		w2_sum += th_weight;

		func2 = 0.0;
		ps_sum = 0.0; pi_sum = 0.0;

		for (j = 1; j <= ContHandle->ComItemNum; j++) {
			for (k = 1; k <= ContComItem[j].CatNum; k++) {
				Ujk = ContComItem[j].ScoreFunc[k];
				p_origi = ProbNew(&ContComItem[j], k, theta, on, x[1], x[2]);
				p_trans = ProbNew(&ContComItem[j], k, theta, off, x[1], x[2]);
				ps = PdNewOverS(&ContComItem[j], k, theta, x[1], x[2]);
				pi = PdNewOverI(&ContComItem[j], k, theta, x[1], x[2]);

				func2 += Ujk * (p_origi - p_trans);
				ps_sum += Ujk*ps;
				pi_sum += Ujk*pi;
			}
		}
		ps_f2_sum += func2 * ps_sum * th_weight;
		pi_f2_sum += func2 * pi_sum * th_weight;
	}

	/* symmetric or non-symmetric setting*/
	if (ContSym == symmetric) {sym_f1 = 1.0; sym_f2 = 1.0;}
	else {
		if (ContSym == old_scale) {sym_f1 = 1.0; sym_f2 = 0.0;}
		else                      {sym_f1 = 0.0; sym_f2 = 1.0;}
	}

	/* function standardization */
	if (ContFuncStd == on) {
		grad[1] = -2.0*(sym_f1*ps_f1_sum/w1_sum + sym_f2*ps_f2_sum/w2_sum);
		grad[2] = -2.0*(sym_f1*pi_f1_sum/w1_sum + sym_f2*pi_f2_sum/w2_sum);
	}
	else {
		grad[1] = -2.0*(sym_f1*ps_f1_sum + sym_f2*ps_f2_sum);
		grad[2] = -2.0*(sym_f1*pi_f1_sum + sym_f2*pi_f2_sum);
	}
	return;
}

/***********************************************************************************/

void ThetaInfoRead(FILE *inf, const char *oldOrnew, struct IRTstControl *Handle)
/*------------------------------------------------------------------------------*/
{
	int i, ThetaNum;

	std::fscanf(inf, "%d",   &ThetaNum); /* reading the number of ability points */
	if (std::strcmp(oldOrnew, "old") == 0) Handle->OldThetaNum = ThetaNum;
	else                              Handle->NewThetaNum = ThetaNum;

	/* memory allocation */
	if (std::strcmp(oldOrnew, "old") == 0) {
	        Handle->OldThetaValues = static_cast<double *>(std::malloc((ThetaNum+1)*sizeof(double)));
		if (!Handle->OldThetaValues) 
			runerror("memory allocation failure 1 in ThetaInfoRead()\n");

	        Handle->OldThetaWeights = static_cast<double *>(std::malloc((ThetaNum+1)*sizeof(double)));
		if (!Handle->OldThetaWeights) 
			runerror("memory allocation failure 2 in ThetaInfoRead()\n");
	}
	else {
	        Handle->NewThetaValues = static_cast<double *>(std::malloc((ThetaNum+1)*sizeof(double)));
		if (!Handle->NewThetaValues) 
			runerror("memory allocation failure 3 in ThetaInfoRead()\n");

	        Handle->NewThetaWeights = static_cast<double *>(std::malloc((ThetaNum+1)*sizeof(double)));
		if (!Handle->NewThetaWeights) 
			runerror("memory allocation failure 4 in ThetaInfoRead()\n");
	}
        
	for (i=1; i <= ThetaNum; i++) {
		if (std::strcmp(oldOrnew, "old") == 0) {
	        	std::fscanf(inf, "%lf", &(Handle->OldThetaValues[i]));
	        	std::fscanf(inf, "%lf", &(Handle->OldThetaWeights[i]));
		}
		else {
	        	std::fscanf(inf, "%lf", &(Handle->NewThetaValues[i]));
	        	std::fscanf(inf, "%lf", &(Handle->NewThetaWeights[i]));
		}
	}
	return;
}

/***********************************************************************************/

void StItemDeAlloc(struct ItemSpec *Item, const char *oldOrnew, struct IRTstControl *Handle)
/*------------------------------------------------------------------------------*/
{
	int j, ItemNum;
	if (std::strcmp(oldOrnew, "old")==0) ItemNum = Handle->OldItemNum;
	else                            ItemNum = Handle->NewItemNum;
	
	for (j=1; j <= ItemNum; j++) {
		std::free(Item[j].ScoreFunc);
		std::free(Item[j].a);
		std::free(Item[j].b);
		std::free(Item[j].c);
		std::free(Item[j].d);
	}

	std::free(Item);
	
	return;
}

/***********************************************************************************/

void StComItemDeAlloc(struct CommonItemSpec *ComItem, struct IRTstControl *Handle)
/*------------------------------------------------------------------------------*/
{
	int j;
	for (j=1; j <= Handle->ComItemNum; j++) {
		std::free(ComItem[j].ScoreFunc);
		std::free(ComItem[j].Na);
		std::free(ComItem[j].Nb);
		std::free(ComItem[j].Nc);
		std::free(ComItem[j].Nd);
		std::free(ComItem[j].Oa);
		std::free(ComItem[j].Ob);
		std::free(ComItem[j].Oc);
		std::free(ComItem[j].Od);
	}
	std::free(ComItem);
	
	return;
}

/***********************************************************************************/

void StContDeAlloc(struct IRTstControl *Handle)
/*------------------------------------------------------------------------------*/
{
	std::free(Handle->NewThetaValues);
	std::free(Handle->NewThetaWeights);
	std::free(Handle->OldThetaValues);
	std::free(Handle->OldThetaWeights);

	return;
}

/***********************************************************************************/

 void Wrapper_IRTst(FILE *outf, char tt[],
     char ItemNewFile[], char ItemOldFile[], char ItemCommonFile[], 
     char DistNewFile[], char DistOldFile[], 
     int HA, enum symmetry HAsym, enum OnOff HAfs, double HAs, double HAi,
     int SL, enum symmetry SLsym, enum OnOff SLfs, double SLs, double SLi,
     char ST[], int PrintFiles)         
/*------------------------------------------------------------------------------*/
{ 
  FILE *inf;
  struct ItemSpec *NewItems, *OldItems;
  struct CommonItemSpec *ComItems;
  struct IRTstControl StInfo; 
  double MSa, MSb, MMa, MMb, HAa = 0, HAb = 0, SLa = 0, SLb = 0, a = 0, b = 0; 
  char ItemNewFileST[103],    /* file name for transformed new item parameters */
       DistNewFileST[103]; /* file name for transformed new group ability dist */
  char line1[1000],                            /* first line; \n replace by \t */
       line2[1000];                                             /* second line */                           

  /* read the input for the new form items */

  inf = std::fopen(ItemNewFile, "r"); 
  if(!inf) runerror("can't open file for new form items \n");
  NewItems = ItemInfoRead(inf, "new", &StInfo);
  std::fclose(inf);

  /* read the input for the old form items */

  inf = std::fopen(ItemOldFile, "r");
  if(!inf) runerror("can't open file for old form items \n");
  OldItems = ItemInfoRead(inf, "old", &StInfo);
  std::fclose(inf);  

  /* read the input for the common items */

  inf = std::fopen(ItemCommonFile, "r");
  if(!inf) runerror("can't open file for common items\n");
  ComItems = ItemPairsRead(inf, NewItems, OldItems, &StInfo);
  std::fclose(inf);

  if(std::strcmp(DistNewFile,"NULL") != 0 && std::strcmp(DistOldFile,"NULL") != 0){

    /* read the input for the new group's ability distribution */

    inf = std::fopen(DistNewFile, "r");
    if(!inf) runerror("can't open file for new form ability distribution \n");
    ThetaInfoRead(inf, "new", &StInfo);
    std::fclose(inf);

    /* read the input for the old group's ability distribution */

    inf = std::fopen(DistOldFile, "r");
    if(!inf) runerror("can't open file for old form ability distribution \n");
    ThetaInfoRead(inf, "old", &StInfo);
    std::fclose(inf);
  }

  /* initial print statements */

  std::fprintf(outf,"\n\n%s\n\n",tt);
  std::fprintf(outf,"IRT Scale Transformation\n\n");

  std::fprintf(outf,"Input files:\n\n");

  std::fprintf(outf,"   %s (item parameters for new form X)\n",ItemNewFile);
  if(PrintFiles == 1) Print_file(ItemNewFile,outf);

  std::fprintf(outf,"   %s (item parameters for new form Y)\n",ItemOldFile);
  if(PrintFiles == 1) Print_file(ItemOldFile,outf);

  std::fprintf(outf,"   %s (pairs of ID's for common items)\n",ItemCommonFile);
  if(PrintFiles == 1) Print_file(ItemCommonFile,outf);

  if(std::strcmp(DistNewFile,"NULL") != 0 && std::strcmp(DistOldFile,"NULL") != 0){
    std::fprintf(outf,"   %s (ability distribution for new form X)\n",DistNewFile);
    if(PrintFiles == 1) Print_file(DistNewFile,outf);
    std::fprintf(outf,"   %s (ability distribution for old form Y)\n\n",DistOldFile);
    if(PrintFiles == 1) Print_file(DistOldFile,outf);
  }

  if(HA == 1){
    std::fprintf(outf,"Haebara symmetric option = %s\n",
             (HAsym==0 ? "old_scale" :(HAsym==1 ? "new_scale" : "symmetric")));
    std::fprintf(outf,"Haebara function standardization option = %s\n",
             HAfs ? "on" : "off"); 
    std::fprintf(outf,"Haebara starting value for slope (A) = %9.5f\n",HAs);
    std::fprintf(outf,"Haebara starting value for intercept (B) = %9.5f\n\n",HAi);
  }

  if(SL == 1){
    std::fprintf(outf,"Stocking-Lord symmetric option = %s\n",
             (SLsym==0 ? "old_scale" :(SLsym==1 ? "new_scale" : "symmetric")));
    std::fprintf(outf,"Stocking-Lord function standardization option = %s\n",
             SLfs ? "on" : "off"); 
    std::fprintf(outf,"Stocking-Lord starting value for slope (A) = %9.5f\n",SLs);
    std::fprintf(outf,"Stocking-Lord starting value for intercept (B) = %9.5f\n\n",SLi);
  }

  std::fprintf(outf,"------------------------------------------------\n");
  std::fprintf(outf,"                         A=slope     B=intercept\n");
  
  /* mean/sigma method */

  StMeanSigma(&StInfo, ComItems, &MSa, &MSb);
  std::fprintf(outf,"Mean/Sigma Method:    %10.5f      %10.5f\n", MSa, MSb);
	
  /* mean/mean method */

  StMeanMean(&StInfo, ComItems, &MMa, &MMb);
  std::fprintf(outf,"Mean/Mean Method:     %10.5f      %10.5f\n", MMa, MMb);

  if(HA != 1 && SL != 1)
    std::fprintf(outf,"------------------------------------------------\n");

  /* Haebara method */

  if(HA == 1){                                 /* Haebara method requested */
    if(std::strcmp(DistNewFile,"NULL") == 0 || std::strcmp(DistOldFile,"NULL") == 0)
       runerror("Haebara results cannot be computed because at least\n"
                "one of the ability distribution files is not present");
    StHaebara(&StInfo, ComItems, HAsym, HAfs, HAs, HAi, &HAa, &HAb);
    std::fprintf(outf,"Haebara Method:       %10.5f      %10.5f\n", HAa, HAb);
  }

  if(HA == 1 && SL != 1)
    std::fprintf(outf,"------------------------------------------------\n");

  /* Stocking-Lord method */

  if(SL == 1){                            /* Stocking-Lord method requested */
    if(std::strcmp(DistNewFile,"NULL") == 0 || std::strcmp(DistOldFile,"NULL") == 0)
      runerror("Stocking-Lord results cannot be computed because at least\n"
               "one of the ability distribution files is not present");
	StStockingLord(&StInfo, ComItems, SLsym, SLfs, SLs, SLi, &SLa, &SLb);
	std::fprintf(outf,"Stocking-Lord Method: %10.5f      %10.5f\n", SLa, SLb);
  }
  if(HA == 1 && SL == 1)
    std::fprintf(outf,"------------------------------------------------\n");

  /* scale transformation */ 
   
  if(std::strcmp(ST,"NL") != 0){

    std::strcpy(ItemNewFileST, ItemNewFile);
	std::strcat(ItemNewFileST,".ST");
    std::strcpy(DistNewFileST, DistNewFile);
	std::strcat(DistNewFileST,".ST");

    if(std::strcmp(ST,"MS") == 0)      {a = MSa; b = MSb;}
    else if(std::strcmp(ST,"MM") == 0) {a = MMa; b = MMb;}
    else if(std::strcmp(ST,"HA") == 0) {a = HAa; b = HAb;}
    else if(std::strcmp(ST,"SL") == 0) {a = SLa; b = SLb;}
    else runerror("invalid scale transformation method"); 

    ScaleTransform(ItemNewFileST, DistNewFileST, a, b, NewItems, &StInfo);

    /* print Transformed Item Paramter Estimates for New Form */

    std::fprintf(outf,"\nTransformed Item Paramter Estimates for New Form\n"
                 "Based on %s Method\n"
                 "Results stored in file named: %s\n\n",
                 (!std::strcmp(ST,"MS") ? "Mean/Sigma" :
                 (!std::strcmp(ST,"MM") ? "Mean/Mean" : 
                 (!std::strcmp(ST,"HA") ? "Haebara" : "Stocking-Lord"))),  
           ItemNewFileST);

    inf = std::fopen(ItemNewFileST,"r");
    std::fgets(line1,998,inf); /* skip the first line */
    while(std::fgets(line1,998,inf) != nullptr){
      line1[std::strlen(line1)-1] = '\t';                      /* replace \n with \t */
      std::fprintf(outf,"%s",std::strcat(line1,std::fgets(line2,998,inf)));
    }  
    std::fclose(inf);
    std::fprintf(outf,"------------------------------------------------\n");

    /* print Transformed Ability Distribution for New Group */

    std::fprintf(outf,"\nTransformed Ability distribution for New Group\n"
                 "Based on %s Method\n"
                 "Results stored in file named: %s\n\n",
                 (!std::strcmp(ST,"MS") ? "Mean/Sigma" :
                 (!std::strcmp(ST,"MM") ? "Mean/Mean" : 
                 (!std::strcmp(ST,"HA") ? "Haebara" : "Stocking-Lord"))),  
           DistNewFileST);

    inf = std::fopen(DistNewFileST,"r");
    std::fgets(line1,998,inf); /* skip the first line */
    while(std::fgets(line1,998,inf) != nullptr)
      std::fprintf(outf,"%s",std::strcat(line1,std::fgets(line2,998,inf))); 
    std::fclose(inf);
    std::fprintf(outf,"------------------------------------------------\n");
    
  }

  /* deallocate memory */

  StComItemDeAlloc(ComItems, &StInfo);
  StItemDeAlloc(NewItems, "new", &StInfo);
  StItemDeAlloc(OldItems, "old", &StInfo);
  StContDeAlloc(&StInfo);

  return;
}

/***********************************************************************************/