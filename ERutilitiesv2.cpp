/*
  ERutilities.cpp 

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

  NOTES:

       Unless otherwise noted, raw scores can differ by any constant 
       positive real number, and raw scores can be negative. (An obvious
       exception is the beta-binomial model.)  In general, the user specifies
       the lowest raw score (min), the highest raw score (max), and the 
       increment (inc) between any two adjacent raw scores. Freqencies and 
       related statistics for raw scores are stored in zero-offset vectors 
       or matrices. So, if fd[] is a frequency distribution, then fd[0] is
       associated with min, and fd[(int) ((max-min)/inc + .5)] is associated
       with max.
 
       In these utilities, mind and maxd are the minimum and maximum scores,
       respectively, in the data, whereas min and max are the user-specified
       minimum and maximum scores. It must always be true that min <= mind 
       and max >= maxd.  Often, directly or indirectly, the user will set 
       min = mind and max = maxd, especially for number correct scores.  
       This distinction allows zero frequencies at the low and high ends of
       distributions.       
       
       minp and maxp are the minimum and maximum possible raw scores.
       It must be true that minp<=min and maxp>=max.  minp and maxp are
       used exclusively with scale scales; whereas min and max are used 
       exclusively with raw scores.  Distinguishing between min and minp
       (as well as between max and maxp) can be necessary (or at least
       desirable) in cases such as the following:
         (a) a particular new form has fewer items than the old form,
             but the scale score conversion for the new form must be
             provided for all of the possible old-form scores;
         (b) min is set equal to mind and max is set equal to maxd 
       For more information about minp and maxp, see comments in 
       ReadSSConvTableForY().
       
       No function that involves scale scores is
       design/methodology/smoothing dependent. 

       The single-group design (S or SG) requires one struct BSTATS.
       The random-groups design (R or RG) requires two struct USTATS.
       The common-item non-equivalent groups design (C or CG)
         requires two struct BSTATS. 

  R. L. Brennan

  March 18, 2009                           
*/
#include <cstdio>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "ERutilities.hpp"
#include "NRutilities.hpp"

/**************************************************************************/

int loc(double x, double min, double inc)
/*
  Location of score x in zero-offset vector; 

  Input:
    x = score
    min = min score
    inc = score increment

    NOTE:  best if x, min, and inc are specified very precisely 
           (say at least eight digits)

  Returns location of score x in zero-offset vector
  
  Function calls other than C or NR utilities: None
                                               
  R. L. Brennan

  Date of last revision: 6/30/08         
*/
{
  return static_cast<int>((x - min) / inc + 0.5);
}

/**************************************************************************/

int nscores(double max, double min, double inc)
/*
  Number of scores (or categories in zero-offset vector; 

  Input:
    max = max score
    min = min score
    inc = score increment

    NOTE:  best if max, min, and inc are specified very precisely 
           (say at least eight digits)

  Returns number of scores (or categories) in zero-offset vector
  with (implicit) scores of min(inc)max;   

  Function calls other than C or NR utilities: None
                                               
  R. L. Brennan

  Date of last revision: 6/30/08         

*/
{
  return loc(max, min, inc) + 1;
}

/**************************************************************************/

double score(int loc, double min, double inc)
/*
  Score associated with location in zero-offset vector 

  Input:
    loc = location in zero-offset vector
    min = min score
    inc = score increment

    NOTE:  best if min and inc are specified very precisely 
           (say at least eight digits)

  Returns score associated with location in zero-offset vector

  Function calls other than C or NR utilities: None
                                               
  R. L. Brennan

  Date of last revision: 6/30/08  
*/
{
  return min + loc * inc;
}

/**************************************************************************/

int skipcols(FILE *fp, int k)
/*
  Skip k non-whitespace-delimited columns (i.e., charcters) 
  from file fp; return k or EOF if end-of-file or error 
     
  Function calls other than C or NR utilities: None
                                               
  R. L. Brennan

  Date of last revision: 6/30/08     
*/
{
  int i;
  for (i = 1; i <= k; i++)
    if (std::fgetc(fp) == EOF) return EOF;
  return k;
}

/***************************************************************************/

int atodouble(FILE *fp, char *str, int begin, int end, double *x)
/* Read end-begin+1 characters from file fp; 
  store in str[]; convert to double

  Input:
    fp = pointer to file
    str = string (storage already allocated)
    begin = beginning column of x in a record
    end = end columnn of x in a record

  Output
    x = double version of x as contained in str[]

  Return string length or EOF if end-of-file or error 
    
  Function calls other than C or NR utilities: None
                                               
  R. L. Brennan

  Date of last revision: 6/30/08     
*/
{
  int i;
  str[end - begin + 1] = '\0';
  for (i = begin; i <= end; i++)
    if ((str[i - begin] = std::fgetc(fp)) == EOF) return EOF;
  *x = std::atof(str);
  return end - begin + 1;
}

/***************************************************************************/

int atointeger(FILE *fp, char *str, int begin, int end, int *x)
/* Read end-begin+1 characters from file fp; 
  store in str[]; convert to integer

  Input:
    fp = pointer to file
    str = string (storage already allocated)
    begin = beginning column of x in a record
    end = end columnn of x in a record

  Output
    x = integer version of x as contained in str[]

  Return string length or EOF if end-of-file or error 
    
  Function calls other than C or NR utilities: None
                                               
  R. L. Brennan

  Date of last revision: 6/30/08         
*/
{
  int i;
  str[end - begin + 1] = '\0';
  for (i = begin; i <= end; i++)
    if ((str[i - begin] = std::fgetc(fp)) == EOF) return EOF;
  *x = std::atoi(str);
  return end - begin + 1;
}

/***************************************************************************/

void flushline(FILE *fp)
/* For file pointed to by fp, delete all characters 
   from current file position to newline or return character.

  Function calls other than C or NR utilities: None
                                               
  R. L. Brennan

  Date of last revision: 6/30/08  
*/
{
  char c;

  while (((c = std::fgetc(fp)) != '\n' && c != '\r'));

}

/***************************************************************************/

void runerror(char error_text[])
/* standard error handler
    
  Function calls other than C or NR utilities: None
                                               
  R. L. Brennan

  Date of last revision: 6/30/08   
*/
{

  std::printf("%s\n", error_text);
  std::printf("\n\nType return to exit\n");
  std::getchar();
  std::exit(EXIT_FAILURE);

}

/*******************************************************************/

void convertFtoW(char fname[], int nv, int fields[][3], char tname[])
/* Convert nv selected variables in a fixed format file (fp)
  to a file (named tname[]) with exactly nv tab-delimited 
  variables that can be read directly using other utility functions
    
  Input:
    fname[] = name of fixed-format file
    nv = number of variables to be converted 
    fields[][3] has nv rows (0...nv-1) that represent each variable
      to be converted. For each row,   
        col 1 :  = 0 --> write variable as integer
                 = 1 --> write variable as  double
        col 2: beginning column in file fp for variable
        col 3: ending column in file fp for variable
    tname[] = name of tab-delimited file with selected variables
              file is opened and closed here

    NOTE: not all variables in fname[] need to be converted 

  Output
    Function creates and closes file named tname[] that 
    contains the tab-delimited variables
               
  E.G., suppose file fp named "fixed" contains 
        variables in cols 3-4, 5-6, 8, 10-14:
        
                         11111
        columns 12345678901234
        
                  8614 4 184.9  
                   234 3  19.3   
                  32 1 2 230.4
                        
        If the only needed variables are in cols 3-4, 5-6, 10-14, 
        then could use:                      
           
        int fields[3][3] = {{0,3,4},{0,5,6},{1,10,14}};         
    
        convertFtoW("fixed",3,fields,"tab-delimited);
        
        The file white now contains three variables:
        
                 >>86>>14>>184.9   where >> stands for a single tab
                 >>2>>34>>19.3   
                 >>32>>1>>230.4
                 
  Function calls other than C or NR utilities:
    skipcols()
    atointeger(0
    atodouble()
    flushline()
                                               
  R. L. Brennan

  Date of last revision: 6/30/08  
*/
{
  int i,
    e,                                     /* ending column of a number */
    ii;                                                   /* an integer */
  double dd;                                                  /* a double */
  char str[20];                           /* should-be-long-enough length */
  FILE *fp,                                     /* pointer for input file */
    *ftp;                            /* pointer for tab-delimited file */

  fp = std::fopen(fname, "r");
  if (fp == nullptr) runerror("No file open\n");
  ftp = std::fopen(tname, "w");
  if (ftp == nullptr) runerror("No file open\n");

  for (;;) {
    e = 0;                                                 /* new record */
    for (i = 0; i <= nv - 1; i++) {
      if (skipcols(fp, fields[i][1] - e - 1) == EOF) goto done;   /* skip cols */
      if (fields[i][0] == 0) {
        if (atointeger(fp, str, fields[i][1], fields[i][2], &ii) == EOF)
          goto done;
        std::fprintf(ftp, "\t%d", ii);
      }
      else {
        if (atodouble(fp, str, fields[i][1], fields[i][2], &dd) == EOF)
          goto done;
        std::fprintf(ftp, "\t%lf", dd);
      }
      e = fields[i][2];
    }
    std::fprintf(ftp, "\n");
    flushline(fp);                             /* end of work for record */
  }

done:
  std::fclose(fp);
  std::fclose(ftp);
}

/*************************************************************************/

void cum_rel_freqs(double min, double max, double inc,
  double *rfd, double *crfd)
/*
  Compute cumulative relative frequencies from relative frequencies

  Input
    min = min score
    max = max score
    inc = increment
    rfd[] = rel freq dist

  Output
    crfd[] = cum rel freq dist (space already allocated)

  Function calls other than C or NR utilities:
    loc()
                                               
  R. L. Brennan

  Date of last revision: 6/30/08  
*/
{
  int i;

  crfd[0] = rfd[0];
  for (i = 1; i <= loc(max, min, inc); i++)
    crfd[i] = crfd[i - 1] + rfd[i];

}

/***************************************************************************/

double perc_rank(double min, double max, double inc, double *crfd, double x)
/*
   Compute percentile rank (pr) given cumulative relative frequencies crfd[]

   pr can be computed for any score x, not simply those scores that are 
     actually achieved (i.e, associated with a distribution). Code would be
     simpler for only "achievable" scores, and simpler still for only 
     integer scores.  

   Formula used is the analogue of Equation 2.14 in Kolen nd Brennan (2004) 
     for the more general case of x scores being any real numbers 
     differing by a constant amount (inc) and ranging from min to max

   In effect, the inverse of this function is the perc_point() function.

  Input
    min = min score
    max = max score
    inc = increment
    crfd[] = cum rel freq dist
    x = score for which pr is desired 

  Returns percentile rank

  Function calls other than C or NR utilities:
    loc()
    score()
                                               
  R. L. Brennan

  Date of last revision: 6/30/08  
*/
{
  int i;
  double pr, xstar;

  if (x < min - inc / 2) pr = 0.;

  else if (x < min + inc / 2)
    pr = 100 * ((x - (min - inc / 2)) / inc) * crfd[0];

  else if (x >= max + inc / 2) pr = 100.;

  else {
    for (i = 1; i <= loc(max, min, inc); i++) {
      xstar = score(i, min, inc);
      if (x < xstar + inc / 2) break;
    }
    pr = 100 * (crfd[i - 1] + ((x - (xstar - inc / 2)) / inc) * (crfd[i] - crfd[i - 1]));
  }

  return pr;

}

/***********************************************************/

double interpolate(double x, int ns, double *f)
/* Interpolated value of f at x assuming there are ns score
   categories.  Function treats first category as having 
   a score of 0; last category has score of ns-1. 
   For x < 0 or x > ns-1, extrapolated value is provided.

  Input
    x = score
    ns = number of score categories
    f = f[0]...f[ns-1]: usually (relative) frequencies

  Returned: interpolated value of f at score of x

  Function calls other than C or NR utilities: None
                                               
  R. L. Brennan

  Date of last revision: 6/30/08  
*/
{
  int xi = static_cast<int>(x),    /* truncated (integer) value for x */
    max = ns - 1;
  double value;

  if (x <= 0)                           /* extrapolate low */
    value = f[0] + x * (f[1] - f[0]);
  else if (x >= max)                   /* extrapolate high */
    value = f[max] + (x - max) * (f[max] - f[max - 1]);
  else                                     /* interpolate */
    value = f[xi] + (x - xi) * (f[xi + 1] - f[xi]);

  /* an extrapolated value can be less than 0.  Hence ... */

  return (value > 0) ? value : 1.0e-10;
}

/***************************************************************************/

int ReadRawGet_moments(char fname[], int scol,
  double *moments, double *mind, double *maxd)
/*
    Compute moments from raw scores in file;
    scores need not be integers or positive;
    frequency distribution NOT computed or output
    assumes space allocated for moments[4];
    assumes data are whitespaced delimited
  
  Input
    fname[] = name of file for input
    scol =  column for scores (whitespace delimited columns)
  
  Output
    moments[] = mean, sd, skew, kurt (space already allocated)
    mind = lowest score (double) in file
    maxd = highest score (double) in file
    
  Returns n = number of persons

  Function calls other than C or NR utilities:
    flushline()
                                               
  R. L. Brennan

  Date of last revision: 6/30/08  
*/
{
  int i,
    n = 0;                                        /* number of examinees */
  double x,                                                   /* score */
    low = 10000.,                          /* initial value of mind */
    high = -10000.,                        /* initial value of maxd */
    mean = 0., var = 0., skew = 0., kurt = 0., dev, dev2;
  FILE *fp;

  fp = std::fopen(fname, "r");
  if (fp == nullptr) runerror("No file open\n");

  /* first pass to compute mean, n, low, high */

  for (;;) {
    for (i = 1; i <= scol - 1; i++) if (std::fscanf(fp, "%*s") == EOF) break;
    if (std::fscanf(fp, "%lf", &x) == EOF) break;
    n++;
    if (x < low) low = x;
    if (x > high) high = x;
    mean += x;
    flushline(fp);
  }
  *mind = low;
  *maxd = high;
  mean /= n;
  moments[0] = mean;

  /* second pass to get sum of powers of deviation scores */

  std::rewind(fp);
  for (;;) {
    for (i = 1; i <= scol - 1; i++) if (std::fscanf(fp, "%*s") == EOF) break;
    if (std::fscanf(fp, "%lf", &x) == EOF) break;
    dev = x - mean;
    dev2 = dev * dev;
    var += dev2;
    dev *= dev2;
    skew += dev;
    dev2 = dev2 * dev2;
    kurt += dev2;
    flushline(fp);
  }
  var /= n;
  skew /= n;
  kurt /= n;

  /* compute s.d., skewness, and kurtosis */

  moments[1] = std::sqrt(var);
  var *= moments[1];
  moments[2] = skew / var;
  var *= moments[1];
  moments[3] = kurt / var;

  std::fclose(fp);

  return n;
}

/**************************************************************************/

int ReadRawGet_mind_maxd(char fname[], int scol, double *mind, double *maxd)

/*
    get minimum and maximum raw score in a file;
    assumes data are whitespace-delimited
  
  Input
    fname[] = name of input file
    scol =  score column (whitespace delimited columns) 
    
  Output 
    mind = minimum score in file 
    maxd = maximum score in file
    
  Returns n = number of examinees in file
    
  Function calls other than C or NR utilities:
    flushline()
                                               
  R. L. Brennan

  Date of last revision: 6/30/08  
*/

{
  int i,
    n = 0;                                        /* number of examinees */
  double x,                                                     /* score */
    low = 10000,                          /* minimum starting value */
    high = -10000;                        /* maximum starting value */
  FILE *fp;

  fp = std::fopen(fname, "r");
  if (fp == nullptr) runerror("No file open\n");

  for (;;) {
    for (i = 1; i <= scol - 1; i++)
      if (std::fscanf(fp, "%*s") == EOF) break;                    /* skip to x */
    if (std::fscanf(fp, "%lf", &x) == EOF) break;                      /* read x */

    n++;
    if (x < low) low = x;
    if (x > high) high = x;
    flushline(fp);
  }
  *mind = low;
  *maxd = high;

  std::fclose(fp);
  return n;
}

/**************************************************************************/

void ReadRawGet_USTATS(char fname[], int scol, double min,
  double max, double inc, char id, struct USTATS *s)
/*
  Read raw data and get or assign all elements for struct s;
  assumes data are whitespace-delimited;
  space allocated here for fd[] 
  
  Input
    fname[] = name of input file
    scol =  score column
    min = minimum raw score
    max = maximum raw score 
    inc = increment for raw scores
	id = single-character id for variable
        	
    NOTE: min can be lower than the lowest frequency in the file (s->mind),
          and max can be higher than the highest frequency (s->maxd).   
          Assuming inc = 1, this feature might be used to force s->fd[] to 
          range from s->fd[0] to s->fd[k] (for a test with k dich scored 
          items) even though s->fd[0] = s->fd[k] = 0

    NOTE: Raw scores in file, min, max, and inc should contain
          sufficient precision whenever non-integer raw scores
          are possible.  This is especially true when the fractional
          part of inc is a repeating decimal (e.g., .33333...) in which
          case it is suggested that at least eight decimal digits be used.
           
  Output
    populates values for struct USTATS s
     
    NOTE: space for s->fd[] must be allocated before calling function
        number of cells = number of scores = nscores(max,min,inc)

  E.G.: if inc = 1, min = -5, and max = 20, then
  
        scores              -5 -4 -3 -2 -1 0 1 2 3 ... 20
        position in s->fd[]  0  1  2  3  4 5 6 7 8 ... 25
        i.e., frequency for score of x is stored in x-min
              position p has freq for score of p+min 
              highest position is max-min
              s->fd[] has max-min+1 cells
    
  Function calls other than C or NR utilities:
    nscores()
    loc()
    flushline()
    MomentsFromFD()
    cum_rel_freqs()
    perc_rank()
    runerror()
                                               
  R. L. Brennan

  Date of last revision: 6/30/08  
*/
{
  int i,
    n = 0;                                         /* number of examinees */

  double x,                                                      /* score */
    low = 10000.,                         /* minimum starting value */
    high = -10000.;                       /* maximum starting value */
  FILE *fp;

  fp = std::fopen(fname, "r");
  if (fp == nullptr) runerror("No file open\n");

  std::strcpy(s->fname, fname);
  s->id = id;
  s->min = min;
  s->max = max;
  s->inc = inc;

  s->ns = nscores(max, min, inc);                            /* allocations */
  s->fd = ivector(0, s->ns - 1);
  s->dbl_fd = dvector(0, s->ns - 1);
  s->cfd = ivector(0, s->ns - 1);
  s->rfd = dvector(0, s->ns - 1);
  s->crfd = dvector(0, s->ns - 1);
  s->prd = dvector(0, s->ns - 1);

  for (i = 0; i <= s->ns - 1; i++) s->fd[i] = 0;                     /* initialize */
  for (;;) {
    for (i = 1; i <= scol - 1; i++)
      if (std::fscanf(fp, "%*s") == EOF) break;                     /* skip to x */
    if (std::fscanf(fp, "%lf", &x) == EOF) break;                       /* read x */

    if (x < min) {
      std::printf("In file %s, for person %d, score = %lf", fname, n + 1, x);
      runerror(" < min");
    }
    if (x > max) {
      std::printf("In file %s, for person %d, score = %lf", fname, n + 1, x);
      runerror(" > max");
    }

    if (x < low) low = x;
    if (x > high) high = x;
    s->fd[loc(x, min, inc)]++;                        /* increment s->fd[] */
    n++;
    flushline(fp);
  }
  s->mind = low;
  s->maxd = high;
  if ((s->n = MomentsFromFD(min, max, inc, nullptr, s->fd, s->mts)) != n)
    runerror("\nError somewhere in ReadRawGet_USTATS()");

  s->cfd[0] = s->fd[0];
  for (i = 1; i <= loc(max, min, inc); i++) s->cfd[i] = s->cfd[i - 1] + s->fd[i];
  for (i = 0; i <= loc(max, min, inc); i++) s->rfd[i] = static_cast<double>(s->fd[i]) / n;
  cum_rel_freqs(min, max, inc, s->rfd, s->crfd);

  for (i = 0; i <= loc(max, min, inc); i++) {
    s->prd[i] = perc_rank(min, max, inc, s->crfd, score(i, min, inc));
    s->dbl_fd[i] = s->fd[i];
  }

  std::fclose(fp);
}

/**************************************************************************/

void ReadFdGet_USTATS(char fname[], int scol, int fcol, double min,
  double max, double inc, char id, struct USTATS *s)
/* Read frequency distribution from file fp;
  get or assign all elements for struct s
  assumes data are whitespace-delimited;
  space allocated here for s->fd[];
  assumes frequency column follows score column;
  scores with zero frequencies need not be in file;
  scores need not be ordered low-to-high, nor high-to-low
  
  Input
    fname[] = name of input file
    scol = score column (whitespace-delimited columns) 
    fcol = frequency (whitespace-delimited columns)
    min = minumum raw score
    max = maximum raw score 
    inc = increment for raw scores
	id = single-character id for variable
	
  See notes and examples in comments for ReadRawGet_USTATS()
    
  Function calls other than C or NR utilities:
    nscores()
    loc()
    flushline()
    MomentsFromFD()
    cum_rel_freqs()
    perc_rank()
    runerror()
                                               
  R. L. Brennan

  Date of last revision: 6/30/08     
*/

{
  int i,
    f,                                                    /* frequency */
    n = 0;                                        /* number of examinees */
  double x;                                                     /* score */
  FILE *fp;

  std::strcpy(s->fname, fname);
  s->id = id;
  s->min = min;
  s->max = max;
  s->inc = inc;
  s->ns = nscores(max, min, inc);

  s->fd = ivector(0, s->ns - 1);                            /* allocations */
  s->dbl_fd = dvector(0, s->ns - 1);
  s->cfd = ivector(0, s->ns - 1);
  s->rfd = dvector(0, s->ns - 1);
  s->crfd = dvector(0, s->ns - 1);
  s->prd = dvector(0, s->ns - 1);

  fp = std::fopen(fname, "r");
  if (fp == nullptr) runerror("No file open\n");

  for (i = 0; i <= s->ns - 1; i++) s->fd[i] = 0;                     /* initialize */
  for (;;) {
    for (i = 1; i <= scol - 1; i++)
      if (std::fscanf(fp, "%*s") == EOF) break;                     /* skip to x */
    if (std::fscanf(fp, "%lf", &x) == EOF) break;                       /* read x */

    for (i = 1; i <= fcol - scol - 1; i++)
      if (std::fscanf(fp, "%*s") == EOF) break;                     /* skip to f */
    if (std::fscanf(fp, "%d", &f) == EOF) break;                        /* read f */

    if (x < min) {
      std::printf("In file %s, there is at least one ", fname);
      runerror("score < min");
    }
    if (x > max) {
      std::printf("In file %s, there is at least one ", fname);
      runerror("score > max");
    }

    s->fd[loc(x, min, inc)] = f;                    /* freq for score of x */
    n += f;
    flushline(fp);
  }
  for (i = 0; i <= s->ns - 1; i++) if (s->fd[i]) break;
  s->mind = score(i, min, inc);
  for (i = s->ns - 1; i >= 0; i--) if (s->fd[i]) break;
  s->maxd = score(i, min, inc);
  if ((s->n = MomentsFromFD(min, max, inc, nullptr, s->fd, s->mts)) != n)
    runerror("\nError somewhere in ReadFdGet_USTATS()");

  s->cfd[0] = s->fd[0];
  for (i = 1; i <= loc(max, min, inc); i++) s->cfd[i] = s->cfd[i - 1] + s->fd[i];
  for (i = 0; i <= s->ns - 1; i++) s->rfd[i] = s->fd[i] / (static_cast<double>(n));
  cum_rel_freqs(min, max, inc, s->rfd, s->crfd);

  for (i = 0; i <= loc(max, min, inc); i++) {
    s->prd[i] = perc_rank(min, max, inc, s->crfd, score(i, min, inc));
    s->dbl_fd[i] = s->fd[i];
  }

  std::fclose(fp);
}

/*************************************************************************/

void ReadRawGet_BSTATS(char fname[], int rows, int cols,
  double rmin, double rmax, double rinc,
  double cmin, double cmax, double cinc,
  char rid, char cid, struct BSTATS *s)
/* Reads raw data and get or assign all elements for bivariate struct s;
  assumes data are whitespace-delimited
   
  Input

    fname[] = name of input file
    rows = column (in file) for score for rows in s->bfd[][] 
    cols = column (in file) for score for columns in s->bfd[][]
    rmin = min value for row scores in s->bfd[][]
    rmax = max value for row scores in s->bfd[][]
    rinc = increment for row scores in s->bfd[][]
    cmin = min value for column scores in s->bfd[][]
    cmax = max value for column scores in s->bfd[][]
    cinc = increment for column scores in s->bfd[][]
	rid = single-character id for rows
	cid = single-character id for columns
      
    
    NOTE: Almost always rinc = cinc
    NOTE: It is not necessary that rows<cols. This permits rows/cols to be
          associated with either variable. 
    NOTE: Most other functions assume that the cols var is associated with
          common items, v.
    NOTE: Raw scores in file, rmin, rmax, rinc, cmin, cmax, and cinc should 
          contain sufficient precision whenever non-integer raw scores
          are possible.  This is especially true when the fractional
          part of rinc (or cinc) is a repeating decimal (e.g., .33333...) 
          in which case it is suggested that at least eight decimal digits
          be used.

  Output
    populates struct BSTATS s
  
    NOTE: space allocated here for fd1[], fd2[], bfd[][], and bp12[][]
 
  Example      Suppose rinc = cinc = 1.  Also suppose that
               rows for matrix are X with scores ranging from
               0 to highest # correct raw score = nx; 
               columns for matrix are V with scores ranging from 0 to
               highest # correct raw score = nv; column (in file)
               for matrix rows (X) is 3; column (in file) for matrix
               columns (V) is 1.  Then we might use
               
               struct BSTATS xv;
               char infname[] = "input data", ;
               
               ReadRawGet_BSTATS(infname,3,1,0,nx,1,0,nv,1,'X','V',&xv);  
               
  Function calls other than C or NR utilities:
    nscores()
    loc()
    flushline()
    MomentsFromFD()
    cum_rel_freqs()
    perc_rank()
    runerror()
                                               
  R. L. Brennan

  Date of last revision: 6/30/08      
*/
{
  int i, j, v,
    n = 0,                                           /* number of persons */
    onec = rows,
    twoc = cols,
    h1;                                                   /* hold value */

  double one,                                      /* row score --> var 1 */
    two;                                      /* col score --> var 2 */

  FILE *fp;

  fp = std::fopen(fname, "r");
  if (fp == nullptr) runerror("No file open\n");

  s->ns1 = nscores(rmax, rmin, rinc);
  s->ns2 = nscores(cmax, cmin, cinc);

  s->fd1 = ivector(0, s->ns1 - 1);                           /* allocation */
  s->fd2 = ivector(0, s->ns2 - 1);
  s->bfd = imatrix(0, s->ns1 - 1, 0, s->ns2 - 1);

  s->dbl_fd1 = dvector(0, s->ns1 - 1);
  s->dbl_fd2 = dvector(0, s->ns2 - 1);
  s->dbl_bfd = dmatrix(0, s->ns1 - 1, 0, s->ns2 - 1);

  s->cfd1 = ivector(0, s->ns1 - 1);
  s->rfd1 = dvector(0, s->ns1 - 1);
  s->crfd1 = dvector(0, s->ns1 - 1);
  s->prd1 = dvector(0, s->ns1 - 1);

  s->cfd2 = ivector(0, s->ns2 - 1);
  s->rfd2 = dvector(0, s->ns2 - 1);
  s->crfd2 = dvector(0, s->ns2 - 1);
  s->prd2 = dvector(0, s->ns2 - 1);

  s->bp12 = dmatrix(0, s->ns1 - 1, 0, s->ns2 - 1);

  n = 0;                                                    /* initialize */
  s->mind1 = s->mind2 = 10000;
  s->maxd1 = s->maxd2 = -10000;
  s->cov = 0.;
  for (i = 0; i <= s->ns1 - 1; i++) s->fd1[i] = 0;
  for (j = 0; j <= s->ns2 - 1; j++) s->fd2[j] = 0;
  for (i = 0; i <= s->ns1 - 1; i++) for (j = 0; j <= s->ns2 - 1; j++) s->bfd[i][j] = 0;

  s->min1 = rmin; s->max1 = rmax; s->inc1 = rinc; s->id1 = rid;
  s->min2 = cmin; s->max2 = cmax; s->inc2 = cinc; s->id2 = cid;
  std::strcpy(s->fname, fname);

  if (cols < rows) {                      /* interchange values if cols<rows */
    h1 = onec;
    onec = twoc; twoc = h1;
  }

  for (;;) {

    for (i = 1; i <= onec - 1; i++)
      if (std::fscanf(fp, "%*s") == EOF) break;    /* skip to 1st score in record */
    if (std::fscanf(fp, "%lf", (cols > rows) ? &one : &two) == EOF) break; /* row sc */

    for (i = 1; i <= twoc - onec - 1; i++)
      if (std::fscanf(fp, "%*s") == EOF) break;    /* skip to 2nd score in record */
    if (std::fscanf(fp, "%lf", (cols > rows) ? &two : &one) == EOF) break; /* col sc */

    if (one < rmin) {
      std::printf("In file %s, for person %d, row score = %f", fname, n + 1, one);
      runerror(" < rmin");
    }
    if (one > rmax) {
      std::printf("In file %s, for person %d, row score = %f", fname, n + 1, one);
      runerror(" > rmax");
    }
    if (two < cmin) {
      std::printf("In file %s, for person %d, col score = %f", fname, n + 1, two);
      runerror(" < cmin");
    }
    if (two > cmax) {
      std::printf("In file %s, for person %d, col score = %f", fname, n + 1, two);
      runerror(" > cmax");
    }

    if (one < s->mind1) s->mind1 = one;
    if (one > s->maxd1) s->maxd1 = one;
    if (two < s->mind2) s->mind2 = two;
    if (two > s->maxd2) s->maxd2 = two;

    s->fd1[loc(one, rmin, rinc)]++;
    s->fd2[loc(two, cmin, cinc)]++;
    s->bfd[loc(one, rmin, rinc)][loc(two, cmin, cinc)]++;
    s->cov += one * two;
    n++;
    flushline(fp);
  }

  if ((s->n = MomentsFromFD(s->min1, s->max1, s->inc1, nullptr,
    s->fd1, s->mts1)) != n)
    runerror("\nError somewhere in ReadRawGet_BSTATS()");
  if ((s->n = MomentsFromFD(s->min2, s->max2, s->inc2, nullptr,
    s->fd2, s->mts2)) != n)
    runerror("\nError somewhere in ReadRawGet_BSTATS()");
  s->cov = s->cov / s->n - (s->mts1[0]) * (s->mts2[0]);
  s->corr = s->cov / (s->mts1[1] * s->mts2[1]);

  s->cfd1[0] = s->fd1[0];
  for (i = 1; i <= loc(rmax, rmin, rinc); i++) s->cfd1[i] = s->cfd1[i - 1] + s->fd1[i];
  for (i = 0; i <= loc(rmax, rmin, rinc); i++) s->rfd1[i] = static_cast<double>(s->fd1[i]) / n;
  cum_rel_freqs(rmin, rmax, rinc, s->rfd1, s->crfd1);

  for (i = 0; i <= loc(rmax, rmin, rinc); i++) {
    s->prd1[i] = perc_rank(rmin, rmax, rinc, s->crfd1, score(i, rmin, rinc));
    s->dbl_fd1[i] = s->fd1[i];
  }

  s->cfd2[0] = s->fd2[0];
  for (i = 1; i <= loc(cmax, cmin, cinc); i++) s->cfd2[i] = s->cfd2[i - 1] + s->fd2[i];
  for (i = 0; i <= loc(cmax, cmin, cinc); i++) s->rfd2[i] = static_cast<double>(s->fd2[i]) / n;
  cum_rel_freqs(cmin, cmax, cinc, s->rfd2, s->crfd2);

  for (i = 0; i <= loc(cmax, cmin, cinc); i++) {
    s->prd2[i] = perc_rank(cmin, cmax, cinc, s->crfd2, score(i, cmin, cinc));
    s->dbl_fd2[i] = s->fd2[i];
  }

  /* bivariate proportions for frequency estimation; as well as
     double verions of bfd[][] */

  for (v = 0; v <= s->ns2 - 1; v++)               /* v indexes second variable here */                                            /* v external to x */
    for (i = 0; i <= s->ns1 - 1; i++) {
      s->bp12[i][v] = static_cast<double>(s->bfd[i][v]) / s->n;
      s->dbl_bfd[i][v] = s->bfd[i][v];
    }

  std::fclose(fp);
}

/****************************************************************************/

void ReadSSConvTableForY(char nameyct[], double minp, double maxp, double inc,
  double **yct)
/* Read raw-to-scale score conversion table for base form y
      from file fp;
    assumes data are white-space delimited with raw scores in first col
      and scale scores in second column --- almost always unrounded scale
      scores should be used;
    assumes that space already allocated for yct[][] with
        #rows = 2 and  #cols = nscores(maxp,minp,inc)+2;
    assumes there are nscores+2 records where
            first record is for minp-(inc/2), 
            second record is for min,
            third record is for min+inc,
            ...,
            3rd last record is for max-inc
            2nd last record is for max,
            last record is for maxp+(inc/2)
    
  Input
  
    nameyct[] = name of file containing conversion table
                NOTE: all records (including last) must end with a 
                single newline (i.e., return)
    minp = minimum possible raw score on base form, not simply the minimum
           score in a particular data set 
           (should be set to 0 for number correct)
    maxp = maximum possible raw score on base form, not simply the maximum
           score in a particular data set 
           (should be set to # items for number correct)
    inc = raw score increment
           
  Output
    
    yct[][] = scale scores associated with raw scores for form Y, where, in 
              accordance with the conventions in K&B (2004, pp. 56-59),         
            
              yct[0][0] = raw score of minp-(inc/2)
              yct[0][1] = raw score of minp
              yct[0][2] = raw score of minp+inc
              ...
              yct[0][nscores(maxp,minp,inc)-1] = raw score of maxp-inc
              yct[0][nscores(maxp,minp,inc)] = raw score of maxp
              yct[0][nscores(maxp,minp,inc)+1] = raw score of maxp+(inc/2)
         
              Scale scores are located in yct[1][]
            
              NOTE: highest location is nscores(maxp,minp,inc)+1
                    number of scores is nscores(maxp,minp,inc)+2                
              NOTE: Assumes storage already allocated using
                    yct = dmatrix(0,1,0,nscores(maxp,minp,inc)+1)
                    This is done in Wrapper_ESS() 
              NOTE: The value scored in yct[0][i] is the raw
                    score value in the file.  This could be more
                    or less precise than score(i,minp,inc) which
                    depends on the precision of inc. The yct[0][i]
                    values are used by Equated_ss() to get equated scale
                    scores, which means that yct[0][i] should be specified
                    with considerable precision (say five decimal places) 
                    when minp, maxp, and inc are not all whole numbers.                       
           
  NOTE: this function is called from Wrapper_ESS()
        
  Function calls other than C or NR utilities:
    nscores()
    flushline()
    runerror()
                                               
  R. L. Brennan

  Date of last revision: 6/30/08      
*/
{
  int i;
  double x,                                                  /* raw score */
    ss;                                                /* scale score */
  FILE *fp;

  fp = std::fopen(nameyct, "r");
  if (fp == nullptr) runerror("No file open\n");

  for (i = 0; i <= nscores(maxp, minp, inc) + 1; i++) {
    if (std::fscanf(fp, "%lf%lf", &x, &ss) == EOF)
      runerror("end-of-file encountered before expected in"
        " ReadRawToSSConvTable()");
    yct[0][i] = x;
    yct[1][i] = ss;
    flushline(fp);
  }

  std::fclose(fp);
}

/**************************************************************************/

int MomentsFromFD(double min, double max, double inc, double *scores,
  int *fd, double *moments)
/*
  Compute moments from a freqency distribution, fd[] with
  frequencies stored as integers; scores are either
       (i) min(inc)max raw scores created using score(i,min,inc)
       (ii) scores[]
  If scores == NULL, then (i) is used
  If scores != NULL, assumes scores are in 
    scores[0]...scores[loc(max,min,inc)]=
                scores[nscores(max,min,inc)-1]

  Assumes storage allocated for moments[4]    
 
  Input
    min = minimum raw score
    max = maximum raw score
    inc = raw score increment
    scores[] = scores
    fd[] = frequencies 
          
  Output
    moments[4] = mean, sd, skew, kurt (space already allocated)
    
  Returns
    n = number of persons 

  Function calls other than C or NR utilities:
    nscores()
    score()
                                               
  R. L. Brennan

  Date of last revision: 6/30/08     
*/

{
  int i,
    ns = nscores(max, min, inc),                    /* number of scores */
    n = 0;                                       /* number of examinees */
  double mean = 0., var = 0., skew = 0., kurt = 0., dev, dev2;
  double *dscores;

  if (scores == nullptr) {
    dscores = dvector(0, ns - 1);
    for (i = 0; i <= ns - 1; i++) dscores[i] = score(i, min, inc);
  }
  else
    dscores = scores;

  /* first pass to compute mean and n */

  for (i = 0; i <= ns - 1; i++) {
    n += fd[i];
    mean += dscores[i] * fd[i];
  }

  mean /= n;
  moments[0] = mean;

  /* second pass to get sum of powers of deviation scores */

  for (i = 0; i <= ns - 1; i++) {
    dev = dscores[i] - mean;
    dev2 = dev * dev;
    var += dev2 * fd[i];
    dev *= dev2 * fd[i];
    skew += dev;
    dev2 = dev2 * dev2 * fd[i];
    kurt += dev2;
  }
  var /= n;
  skew /= n;
  kurt /= n;

  /* compute s.d., skewness, and kurtosis */

  moments[1] = std::sqrt(var);
  var *= moments[1];
  moments[2] = skew / var;
  var *= moments[1];
  moments[3] = kurt / var;

  // Bug Fix: Free memory if it was allocated in this function
  if (scores == nullptr) {
    free_dvector(dscores, 0, ns - 1);
  }

  return n;
}

/**************************************************************************/

void MomentsFromRFD(double min, double max, double inc, double *scores,
  double *rfd, double *moments)
/*
  Compute moments from a relative freqency distribution in rfd[]; 
  possbile scores are either
       (i) min(inc)max raw scores are created using score(i,min,inc)
       (ii) scores[]
  If scores == NULL, then (i) is used
  If scores != NULL, assumes scores are in 
    scores[0]...scores[loc(max,min,inc)]=
                scores[nscores(max,min,inc)-1]

  Assumes storage allocated for moments[4]    
 
  Input
    min = minimum raw score
    max = maximum raw score
    inc = raw score increment
    scores[] = scores
    rfd[] = relative frequencies 
          
  Output
    moments[4] = mean, sd, skew, kurt (space already allocated)
    
  Returns
    n = number of persons 

  Function calls other than C or NR utilities:
    nscores()
    score()
                                               
  R. L. Brennan

  Date of last revision: 6/30/08   
*/

{
  int i,
    ns = nscores(max, min, inc);                    /* number of scores */
  double mean = 0., var = 0., skew = 0., kurt = 0., dev, dev2;
  double *dscores;

  if (scores == nullptr) {
    dscores = dvector(0, ns - 1);
    for (i = 0; i <= ns - 1; i++) dscores[i] = score(i, min, inc);
  }
  else
    dscores = scores;

  /* first pass to compute mean and n */

  for (i = 0; i <= ns - 1; i++) {
    mean += dscores[i] * rfd[i];
  }

  moments[0] = mean;

  /* second pass to get sum of powers of deviation scores */

  for (i = 0; i <= ns - 1; i++) {
    dev = dscores[i] - mean;
    dev2 = dev * dev;
    var += dev2 * rfd[i];
    dev *= dev2 * rfd[i];
    skew += dev;
    dev2 = dev2 * dev2 * rfd[i];
    kurt += dev2;
  }

  /* compute s.d., skewness, and kurtosis */

  moments[1] = std::sqrt(var);
  var *= moments[1];
  moments[2] = skew / var;
  var *= moments[1];
  moments[3] = kurt / var;

  // Bug Fix: Free memory if it was allocated in this function
  if (scores == nullptr) {
      free_dvector(dscores, 0, ns - 1);
  }
}

/*******************************************************************/

void EquiEquate(int nsy, double miny, double incy,
  double *crfdy, int nsx, double *prdx, double *eraw)
/*
  Computes equipercentile equivalents on scale of y for percentile
  ranks on scale of x. See comments in perc_point() for details.
	
	Input
		nsy 	 Number of raw score categories for old form Y
		miny 	 Minimum raw score for old form Y
		incy 	 Increment between consecutive raw scores for old form Y
		crfdy[]  Cumulative rel freq dist for old form Y
        nsx      Number of raw score categories for X
		prdx[]	 Percentile rank distribution for new form X 
 
	Output 
        eraw[]   equiperentile equivalents

  Function calls other than C or NR utilities:
    perc_point()
                                               
  R. L. Brennan

  Date of last revision: 6/30/08         
*/
{
  int i;

  for (i = 0; i <= nsx - 1; i++)
    eraw[i] = perc_point(nsy, miny, incy, crfdy, prdx[i]);
}

/*******************************************************************/

double perc_point(int ns, double min, double inc,
  double *crfd, double pr)
/*
  Computes percentile point (pp) as average of "upper" pp
  (see Kolen & Brennan, 2004, Equation 2.15, p. 45; see also
  Equation 2.18 on p. 46) and "lower" pp (see Kolen & Brennan, 
  2004, p. 45, Equation 2.16)
	
	Input
		ns 		    Number of raw score categories
		min 		Minimum raw score
		inc 		Increment between consecutive raw scores
		crfd[] 		Cumulative rel freq dist 
		pr		    Percentile rank 
 
	Returns percentile point on scale implied by ns, min, and inc 

  Before getting to return statement, all work is based on getting
  a score, say x, on the piecewise continuous 0...(ns-1) scale,
  where k is an integer.  That is, possible values of x are 
  [-.5, (ns-1)+.5]. The return statement transforms x to the scale 
  implied by ns, min, and inc, which means that the transformed values 
  of x are [min - inc/2, min +(ns-1)*inc + inc/2]. 
  (Best if min and inc are specified very precisely--- say at least
  eight decimal digits.)
  
  NOTE: This function gets the percentile point (pp) associated with 
        crfd[] given any pr.  The two obvious uses are (a) obtaining 
        a pp for a pr that does not equal one of the ns values 
        immediately obtainable from crfd[]; and (b) obtaining 
        equipercentile equivalents.  In the latter case,
        pr is for x and all other variables are for y.  (See Kolen & 
        Brennan, 2004, Equation 2.18 on p. 46)

  In effect, this function is the inverse of the perc_rank() function

  Function calls other than C or NR utilities: None
                                               
  R. L. Brennan

  Date of last revision: 6/30/08     
*/
{
  double prp = pr / 100,              /* pr expressed as proportion */
    ppU,     /* "upper" perc point on scale [-.5, (ns-1)+.5] */
    ppL;     /* "lower" perc point on scale [-.5, (ns-1)+.5] */
  int i,  /* i =  x*_U = smallest integer such that crfd[i] > prp */
    j;    /* j = x*_L = largest integer such that crfd[j] < prp */

  /* Special case: PR=0 and 0 freqs at bottom of rfd[].
     First line of code means that prp <= 1.0e-8 
     is an operational definition of prp==0.  Remaining code is
     to handle possibility of 0 freq's at bottom of rfd[].
     Note that for loop starts at 0 */

  if (prp <= 1.0e-8) {
    ppL = -.5;
    for (i = 0; i <= ns - 1; i++) if (crfd[i] > 1.0e-8) break;
    ppU = i - .5;
    return  min + inc * ((ppU + ppL) / 2);
  }

  /* Special case: PR=1 and 0 freqs at top of rfd[].
     First line of code means that 1 - 1.0e-8 is 
     an operational definition of prp==1.  Remaining code is
     to handle possibility of 0 freq's at top of rfd[].
     Not that for loop starts at ns-1 */

  if (prp >= 1. - 1.0e-8) {
    ppU = ns - .5;
    for (j = ns - 1; j >= 0; j--) if (crfd[j] < 1. - 1.0e-8) break;
    ppL = 1 + j + .5;
    return  min + inc * ((ppU + ppL) / 2);
  }

  /* Special case:  crfd[0] > prp can happen in equipercentile
     equating. If this occurs, then we have the following anomalous
     circumstances: 
       (a) x*_U = 0 by defn on p. 45, in which case
           crfd[x*_U - 1] = crfd[-1] which is undefined; and 
       (b) x*_L is non-existent by defn on p. 45.
     Next two lines of code are consistent with graphical 
     procedures in Kolen and Brennan (sect. 2.5.1). 
  */

  if (crfd[0] > prp)
    return min + inc * (prp / crfd[0] - .5);

  /* majority of work occurs next */

  /* upper pp -- get x*_U */

  for (i = 1; i <= ns - 1; i++) if (crfd[i] > prp) break;
  if (crfd[i] != crfd[i - 1])
    ppU = (prp - crfd[i - 1]) / (crfd[i] - crfd[i - 1]) + (i - .5);
  else
    ppU = i - .5;

  /* lower pp -- get x*_L */

  for (j = ns - 2; j >= 0; j--) if (crfd[j] < prp) break;
  if (crfd[j + 1] != crfd[j])
    ppL = (prp - crfd[j]) / (crfd[j + 1] - crfd[j]) + (j + .5);
  else
    ppL = j + .5;

  /* return area */

  return min + inc * ((ppU + ppL) / 2);
}

/*******************************************************************************/

void Wrapper_ESS(struct PDATA *inall, struct ERAW_RESULTS *r, double minp,
  double maxp, double incp, char *nameyct, int round,
  int lprss, int hprss, struct ESS_RESULTS *s)
/*
  Wrapper used for getting equated scale scores essu[][] and essr[][]
  Assigns or computes values for all variables in struct ESS_RESULTS s
  Can be used for any design or methodology 
  
    Input:
    
      uses the following from struct PDATA inall
        nm        = number of methods 
        names[][] = names of methods
        min       = min raw score for x 
        max       = max raw score for x
        inc       = increment for x         
        fdx[]     = fd for new form x 
        rep       = replication number for bootstrap; = 0 if no bootstrap
        
      uses the following from struct ERAW_RESULTS r
        eraw[][] = equated raw scores
        
      minp    = minimum raw score for base form y conversion yct
      maxp    = maximum raw score for base form y conversion yct
	  incp    = increment in raw scores in base form y conversion yct
      nameyct = name of file with raw-to-scale score conversion for Y
      round: if round=i --> round to i-th decimal place
      lprss = lowest possible rounded scale score
      hprss = highest possible rounded scale score
       
    Output: 
      
      s = struct ESS_RESULTS:
            essu[][] = unrounded equated scale scores
            essr[][] = rounded equated scale scores  
            mtsu[][] = moments for equated unrounded scale scores 
            mtsr[][] = moments for equated rounded scale scores 
        
      NOTE: The following are added to struct PDATA inall: minp, maxp, nameyct,
            round, lprss, hprss, and **yct
        
      NOTE: struct ESS_RESULTS must be different struct's for actual 
            equating and bootstrap use of Wrapper_ESS()
        
      NOTE:       
      yct[][] = conversion table for Y; first dimension ranges from
	            0 to (nscores(maxp,minp,inc)+1
				(see comments for ReadSSConvTableForY() for details)

  Function calls other than C or NR utilities: 
    nscores()
    ReadSSConvTableForY()
    Equated_ss()
    MomentsFromFD()
                                               
  R. L. Brennan

  Date of last revision: 6/30/08                                     
*/
{
  int i,
    ns = nscores(inall->max, inall->min, inall->inc);

  inall->minp = minp;
  inall->maxp = maxp;
  inall->incp = incp;
  std::strcpy(inall->nameyct, nameyct);
  inall->round = round;
  inall->lprss = lprss;
  inall->hprss = hprss;

  if (inall->rep == 0) {                                 /* for actual equating */
    inall->yct = dmatrix(0, 1, 0, nscores(maxp, minp, inall->inc) + 1);
    ReadSSConvTableForY(nameyct, minp, maxp, inall->inc, inall->yct);
  }

  if (inall->rep <= 1) {       /* for actual equating or 1st bootstrap rep only */
    s->essu = dmatrix(0, inall->nm - 1, 0, ns - 1);
    s->essr = dmatrix(0, inall->nm - 1, 0, ns - 1);
    s->mtsu = dmatrix(0, inall->nm - 1, 0, 3);         /* 0,3 is for the 4 moments */
    s->mtsr = dmatrix(0, inall->nm - 1, 0, 3);         /* 0,3 is for the 4 moments */
  }

  for (i = 0; i <= inall->nm - 1; i++)
    Equated_ss(inall->min, inall->max, inall->inc, minp, maxp, incp, r->eraw[i],
      inall->yct, inall->round, inall->lprss, inall->hprss, s->essu[i], s->essr[i]);

  /* compute moments:  Note that when inc==1, essu[*][min-minp+1] is the
     unrounded scale score associated with fdx[0], where fdx[0]
     is associated with scores ranging from min to max.
     Recall that essu[*][0] is the unrounded scale score
     associated with minp-inc/2 = minp-.5 when inc = 1.  In this example,
     min-minp+1 = loc(min,minp,inc) + 1
     
  */
  if (inall->fdx != nullptr) {
    for (i = 0; i <= inall->nm - 1; i++) {
      MomentsFromFD(inall->min, inall->max, inall->inc,
        s->essu[i], inall->fdx, s->mtsu[i]);
      MomentsFromFD(inall->min, inall->max, inall->inc,
        s->essr[i], inall->fdx, s->mtsr[i]);
    }
  }

}
/******************************************************************************/

void Equated_ss(double min, double max, double inc, double minp, double maxp,
  double incp, double *eraw, double **yct, int round, int lprss,
  int hprss, double *essu, double *essr)
/* Get SS conversion given raw score equating results (eraw) and raw-to-ss
  conversion Y (yct).  This function is for one method only, and it 
  does not make any assumptions about the equating method used.  For the 4 
  linear results this function can be called 4 times using eraw[0]...eraw[3] 

  NOTE: It is assumed that increment for x is the same as inc in the 
        conversion table for the base form Y
   
  Input
  
    min = minimum raw score for x
    max = maximum raw score for x
    inc = raw score increment for x and y
    minp = minimum raw score for base form (y) conversion
    maxp = maximum raw score for base form (y) conversion
    incp = increment in raw score for base form (y) conversion
    eraw[] = equated raw scores (length = nscores(max,min,inc)+1 
             in positions [0]...[nscores(max,min,inc)-1 =
                                 loc(max,min,inc)]
    yct[][] = scale scores associated with raw scores for form Y
              (see comments for ReadSSConvTableForY() for details)

    round = 0 --> no rounding
          = 1 --> round to first digit before decimal point
                  (e.g., 123.45 --> 123; 123.55 --> 124)
          = 2 --> round to second digit before decimal point
                  (e.g., 123.45 --> 120; 125.55 --> 130)
          ....
    lprss = lowest possible rounded scale score
    hprss = highest possible rounded scale score
  
  Output
  
    essu[] = unrounded scale scores 
    essr[] = rounded scale scores 

  SIMPLE EXAMPLE:  assuming inc = 1 (pp. 81, 85 in K&B)

            Y-raw   Y-ess   X-raw    eraw    
  
              19  18.3403     20  19.1647
              20  19.2844     21  20.3676 
              21  20.1839
                           
           essu[] for x=20 = 18.3403 + .1647*(19.2844-18.3403) = 18.4958
         
           see Equation 2.22 in K&B (2004, p. 58)

  COMPLICATED EXAMPLE:
 
    minp = -1.0; min = -.33333; max = 4.33333; maxp = 5.0; inc = .33333 

             yct[][]                     eraw[]
   ---------------------------   ------------------------      essu[]
   loc  raw score  scale score   loc       raw             --------------
    j  yct[0][loc] yct[1][loc]    i      score  eraw[loc]  loc  essu[loc]

     0    -1.16667     10                                    
     1    -1.00000     23                                   
     2    -0.66667     26                                           
     3    -0.33333     33          0  -0.33333     -.25      0   34.00000    
     4     0.00000     37          1   0.00000     -.10      1   35.80000
    ..     ...         ..         ..   .......     ....     ..
    16     4.00000     65         13   4.00000     4.20     13   67.40000
    17     4.33333     69         14   4.33333     4.40     14   69.40000
    18     4.66667     71                                   
    19     5.00000     77                                   
    20     5.16667     90                                  
  
    e.g., essu[] for raw score of 4,
          which has equated raw score of 4.20 located in eraw[13], is

          essu[16] = 65 + (69 - 65)*[(4.20-4.00000 )/(4.33333-4.0000)]
                   = 65 + 4*(.6) = 67.4

    Note that:
      (a) eraw[] goes from eraw[0] to
          eraw[loc(max,min,inc)] = eraw[nscores(max,min.inc)-1]
      (b) yct[*][] goes from yct[*][0] to 
          yct[*][loc(maxp,minp,inc)+2] = yct[*][nscores(maxp,minp,inc)+1]
      (c) raw score associated with eraw[i] is score(i,min,inc) 

  Function calls other than C or NR utilities: 
    nscores()
    runerror
                                               
  R. L. Brennan

  Date of last revision: 6/30/08    
*/
{

  int i, j, jlast = 0,
    hyct = nscores(maxp, minp, incp) + 1,             /* highest loc in yct[*][] */
    heraw = nscores(max, min, inc) - 1;            /* highest location in eraw[] */

  double y;

  for (i = 0; i <= heraw; i++) {
    y = eraw[i];                                        /* equated raw score */

    if (y <= yct[0][0]) essu[i] = yct[1][0];
    else if (y >= yct[0][hyct]) essu[i] = yct[1][hyct];
    else {
      /* find j such that yct[0][j] > y */
      for (j = jlast; j <= hyct; j++) if (yct[0][j] > y) break;

      essu[i] = yct[1][j - 1] + (yct[1][j] - yct[1][j - 1]) *
        (y - yct[0][j - 1]) / (yct[0][j] - yct[0][j - 1]); /* interpolate */

      jlast = j;
    }
  }

  /* rounding */

  for (i = 0; i <= heraw; i++) {
    if (round) {
      essr[i] = std::pow(10.0, static_cast<double>(round - 1)) *
        (static_cast<int>(essu[i] / std::pow(10.0, static_cast<double>(round - 1)) + .5));
      if (essr[i] < lprss) essr[i] = lprss;
      if (essr[i] > hprss) essr[i] = hprss;
    }
    else
      essr[i] = essu[i];
  }

}

/*******************************************************************************/
// NOTE: The rest of the file continues below, including all Print functions
// and the functions prefixed with "er_"
/*******************************************************************************/

void Print_ESS(FILE *fp, char tt[], struct PDATA *inall, struct ESS_RESULTS *s)
/*
  print raw-to-scale score conversion table for y = yct[][],
  as well as equated scale scores (unrounded and rounded);
  can be used for any design, any methodology, or any set of nm methods;
   
  Input
    fp = file pointer for output
    tt[] = user supplied text identifier
    variables used from struct PDATA inall:
        nm = number of methods 
        name[][]: names of methods 
        minp = minimum raw score in yct
        maxp = maximum raw score in yct
        incp = increment in raw scores for yct
		min = minimum raw score for x
        max = maximum raw score for x
        inc = increment in raw scores for x
        nameyct[100]: name of file containing yct[][] 
        yct[][]: conversion table for Y  
        round: if round = i, then round to i-th place
        lprss = lowest possible rounded scale score
        hprss = highest possible rounded scale score
 
        rep = replication number for bootstrap; 0 for actual equating
        
     variables used from struct ESS_RESULTS s:
        essu[][]: unrounded equated scale scores
        essr[][]: rounded equated scale scores
        mtsu[][]: moments for equated unrounded scale scores 
        mtsr[][]: moments for equated rounded scale scores
        
   NOTE:       
        yct[][] = raw-to-scale score conversion table for y
        yct[0][] = raw scores
        yct[1][] = scale scores
        essu[][] = methods x columns matrix of 
                   unrounded equated scale scores
        essr[][] = methods x columns matrix of 
                   rounded equated scale scores
        see Equated_ss() for more details

  Function calls other than C or NR utilities: None
                                               
  R. L. Brennan

  Date of last revision: 6/30/08    

*/
{
  int i, j, k;

  /* print conversion table for Y (base form) */

  std::fprintf(fp, "\n\n%s\n\n", tt);
  if (inall->rep > 0)  std::fprintf(fp, "Bootstrap relication %d\n\n", inall->rep);
  std::fprintf(fp, "Name of file containing yct[][]: %s\n\n", inall->nameyct);

  std::fprintf(fp, "Minimum raw score in base form conversion, Y = %12.5f\n", inall->minp);
  std::fprintf(fp, "Maximum raw score in base form conversion, Y = %12.5f\n\n", inall->maxp);
  if (inall->round > 0) {
    std::fprintf(fp, "Lowest possible rounded scale score = %6d\n", inall->lprss);
    std::fprintf(fp, "Highest possible rounded scale score = %6d\n\n", inall->hprss);
  }

  std::fprintf(fp, "\nCONVERSION TABLE FOR Y (BASE FORM)\n\n");
  std::fprintf(fp, "                         Y | ");

  std::fprintf(fp, "\n         Raw         Scale | ");
  std::fprintf(fp, "\n       Score         Score | ");
  std::fprintf(fp, "\n\n");

  for (i = 0; i <= nscores(inall->maxp, inall->minp, inall->incp) + 1; i++)
    std::fprintf(fp, "%12.5f  %12.5f |\n", inall->yct[0][i], inall->yct[1][i]);
  std::fprintf(fp, "\n");
  for (i = 1; i <= 27; i++) std::fprintf(fp, "-");

  /* print raw to scale-score (unrounded) conversion table for x;
     then print rounded version if inall->round>0 */

  for (k = 1; k <= ((inall->round > 0) ? 2 : 1); k++) {         /* beginning of loop for k */
    /* k=1 (unrounded); k=2 (rounded) */
    std::fprintf(fp, "\n\n%s\n\n", tt);
    if (inall->rep > 0)  std::fprintf(fp, "Bootstrap relication %d\n\n", inall->rep);
    std::fprintf(fp, "Name of file containing yct[][]: %s\n\n", inall->nameyct);

    /* following code set up for any number of methods */

    for (j = 1; j <= 15 + (inall->nm * 12 - 18) / 2; j++) std::fprintf(fp, " ");
    if (k == 1) std::fprintf(fp, "Equated Scale Scores (unrounded)\n");
    else std::fprintf(fp, "Equated Scale Scores (rounded)\n");

    std::fprintf(fp, "\n               ");
    for (j = 0; j <= inall->nm - 1; j++) std::fprintf(fp, "   Method %d:", j);
    std::fprintf(fp, "\nRaw Score (x)  ");
    for (j = 0; j <= inall->nm - 1; j++) std::fprintf(fp, "  %s", inall->names[j]);
    std::fprintf(fp, "\n");

    for (i = 0; i <= loc(inall->max, inall->min, inall->inc); i++) {
      std::fprintf(fp, "\n %12.5f  ", score(i, inall->min, inall->inc));
      for (j = 0; j <= inall->nm - 1; j++) {
        if (k == 1) std::fprintf(fp, "%12.5f", s->essu[j][i]);
        else std::fprintf(fp, "%12.0f", s->essr[j][i]);
      }
    }

    if (inall->fdx != nullptr) {
      std::fprintf(fp, "\n\n         Mean  ");
      for (j = 0; j <= inall->nm - 1; j++) std::fprintf(fp, "%12.5f",
        (k == 1) ? s->mtsu[j][0] : s->mtsr[j][0]);
      std::fprintf(fp, "\n         S.D.  ");
      for (j = 0; j <= inall->nm - 1; j++) std::fprintf(fp, "%12.5f",
        (k == 1) ? s->mtsu[j][1] : s->mtsr[j][1]);
      std::fprintf(fp, "\n         Skew  ");
      for (j = 0; j <= inall->nm - 1; j++) std::fprintf(fp, "%12.5f",
        (k == 1) ? s->mtsu[j][2] : s->mtsr[j][2]);
      std::fprintf(fp, "\n         Kurt  ");
      for (j = 0; j <= inall->nm - 1; j++) std::fprintf(fp, "%12.5f",
        (k == 1) ? s->mtsu[j][3] : s->mtsr[j][3]);

      std::fprintf(fp, "\n\n");
      for (j = 1; j <= 63; j++) std::fprintf(fp, "*");
    }
  }
}

/*******************************************************************************/

void Print_USTATS(FILE *fp, char tt[], struct USTATS *s)
/*
  Print results in USTATS s
  
  Input
    fp = file pointer for output
    tt[] = user supplied text identifier
    s = struct USTATS

  NOTE:  raw scores are printed using score(i,s->min,s->inc),
         which could be more or less precise than scores read 
         using ReadRawGet_USTATS()

  Function calls other than C or NR utilities:
    score()
                                               
  R. L. Brennan

  Date of last revision: 6/30/08  
*/
{
  int i;

  std::fprintf(fp, "\n\n%s\n\n", tt);
  std::fprintf(fp, "Input filename:  %s\n\n", s->fname);

  std::fprintf(fp, "Variable identifier: %c\n\n", s->id);

  std::fprintf(fp, "Number of persons = %6d\n\n", s->n);
  std::fprintf(fp, "Min score in data = %12.5f\n", s->mind);
  std::fprintf(fp, "Max score in data = %12.5f\n", s->maxd);
  std::fprintf(fp, "             Mean = %12.5f\n", s->mts[0]);
  std::fprintf(fp, "             S.D. = %12.5f\n", s->mts[1]);
  std::fprintf(fp, "             Skew = %12.5f\n", s->mts[2]);
  std::fprintf(fp, "             Kurt = %12.5f\n\n", s->mts[3]);

  std::fprintf(fp, "Min score for fd[] = %10.5f\n", s->min);
  std::fprintf(fp, "Max score for fd[] = %10.5f\n", s->max);
  std::fprintf(fp, "Increment between scores = %10.5f\n", s->inc);
  std::fprintf(fp, "Number of scores = %6d\n\n", s->ns);
  std::fprintf(fp, "       Score  Freq  CFreq    Rel Freq   Cum RFreq   Perc Rank\n\n");
  for (i = 0; i <= s->ns - 1; i++)
    std::fprintf(fp, "%12.5f%6d%7d%12.5f%12.5f%12.5f\n", score(i, s->min, s->inc),
      s->fd[i], s->cfd[i], s->rfd[i], s->crfd[i], s->prd[i]);
  std::fprintf(fp, "\n");
  for (i = 1; i <= 61; i++) std::fprintf(fp, "*");
}

/*******************************************************************************/

void Print_BSTATS(FILE *fp, char tt[], struct BSTATS *s, int pbfd)
/*
  Print results in BSTATS s
  
  Input
    fp = file pointer for output
    tt[] = user supplied text identifier
    s = struct BSTATS
    pbfd: 0 --> don't print s->bfd[][]; 
          1 --> print s->bfd[][];

  Function calls other than C or NR utilities:
    score()
                                               
  R. L. Brennan

  Date of last revision: 6/30/08  
*/
{
  int i, j;
  double sum, grandsum = 0.;

  std::fprintf(fp, "\n\n%s\n\n", tt);
  std::fprintf(fp, "Input filename:  %s   \n\n", s->fname);

  std::fprintf(fp, "Number of Persons = %6d", s->n);

  std::fprintf(fp, "\n\nVariable 1 identifier (rows of bfd[][]): %c\n\n", s->id1);
  std::fprintf(fp, "Min score in data = %12.5f\n", s->mind1);
  std::fprintf(fp, "Max score in data = %12.5f\n", s->maxd1);
  std::fprintf(fp, "             Mean = %12.5f\n", s->mts1[0]);
  std::fprintf(fp, "             S.D. = %12.5f\n", s->mts1[1]);
  std::fprintf(fp, "             Skew = %12.5f\n", s->mts1[2]);
  std::fprintf(fp, "             Kurt = %12.5f\n\n", s->mts1[3]);

  std::fprintf(fp, " Min score for fd1[] = %12.5f\n", s->min1);
  std::fprintf(fp, " Max score for fd1[] = %12.5f\n\n", s->max1);
  std::fprintf(fp, "       Score  Freq  CFreq    Rel Freq   Cum RFreq   Perc Rank\n\n");
  for (i = 0; i <= s->ns1 - 1; i++)
    std::fprintf(fp, "%12.5f%6d%7d%12.5f%12.5f%12.5f\n", score(i, s->min1, s->inc1),
      s->fd1[i], s->cfd1[i], s->rfd1[i], s->crfd1[i], s->prd1[i]);

  std::fprintf(fp, "\n\nVariable 2 identifier (columns of bfd[][]): %c\n\n", s->id2);
  std::fprintf(fp, "Min score in data = %12.5f\n", s->mind2);
  std::fprintf(fp, "Max score in data = %12.5f\n", s->maxd2);
  std::fprintf(fp, "             Mean = %12.5f\n", s->mts2[0]);
  std::fprintf(fp, "             S.D. = %12.5f\n", s->mts2[1]);
  std::fprintf(fp, "             Skew = %12.5f\n", s->mts2[2]);
  std::fprintf(fp, "             Kurt = %12.5f\n\n", s->mts2[3]);

  std::fprintf(fp, " Min score for fd2[] = %12.5f\n", s->min2);
  std::fprintf(fp, " Max score for fd2[] = %12.5f\n\n", s->max2);
  std::fprintf(fp, "       Score  Freq\n\n");
  std::fprintf(fp, "       Score  Freq  CFreq    Rel Freq   Cum RFreq   Perc Rank\n\n");
  for (i = 0; i <= s->ns2 - 1; i++)
    std::fprintf(fp, "%12.5f%6d%7d%12.5f%12.5f%12.5f\n", score(i, s->min2, s->inc2),
      s->fd2[i], s->cfd2[i], s->rfd2[i], s->crfd2[i], s->prd2[i]);

  std::fprintf(fp, "\n\n Covariance(1,2) = %12.5f", s->cov);
  std::fprintf(fp, "\nCorrelation(1,2) = %12.5f", s->corr);

  if (pbfd) {

    /* print bivariate frequency distribution [var1][var2] */

    std::fprintf(fp, "\n\nBivariate Frequency Distribution bfd[var1][var2]:\n"
      "Loc means location in bfd[var1][var2]\n\n");

    std::fprintf(fp, "        Loc |");
    for (j = 0; j <= s->ns2 - 1; j++)  std::fprintf(fp, "%8d", j);
    std::fprintf(fp, " |\n");
    std::fprintf(fp, "Loc   Score |");
    for (j = 0; j <= s->ns2 - 1; j++)  std::fprintf(fp, "%8.2f", score(j, s->min2, s->inc2));
    std::fprintf(fp, " |   fd1[]\n");
    std::fprintf(fp, "-------------");
    for (j = 0; j <= s->ns2 - 1; j++)  std::fprintf(fp, "--------");
    std::fprintf(fp, "----------\n");

    for (i = 0; i <= s->ns1 - 1; i++) {
      std::fprintf(fp, "%3d%8.2f |", i, score(i, s->min1, s->inc1));
      for (j = 0; j <= s->ns2 - 1; j++)  std::fprintf(fp, "%8d", s->bfd[i][j]);
      std::fprintf(fp, " |%8d\n", s->fd1[i]);
    }

    std::fprintf(fp, "------------");
    for (j = 0; j <= s->ns2 - 1; j++)  std::fprintf(fp, "--------");
    std::fprintf(fp, "-----------\n");
    std::fprintf(fp, "      fd2[] |");
    for (j = 0; j <= s->ns2 - 1; j++)  std::fprintf(fp, "%8d", s->fd2[j]);
    std::fprintf(fp, " |%8d", s->n);
    std::fprintf(fp, "\n\n");

    std::fprintf(fp, "*************");
    for (j = 0; j <= s->ns2 - 1; j++)  std::fprintf(fp, "********");
    std::fprintf(fp, "**********\n");

    /* print bivariate proportions used for frequency estimation */


    std::fprintf(fp, "\n\nBiv Props bp12[var1][var2]\n\n");

    std::fprintf(fp, "        Loc |");
    for (j = 0; j <= s->ns2 - 1; j++)  std::fprintf(fp, "%8d", j);
    std::fprintf(fp, " |\n");
    std::fprintf(fp, "Loc   Score |");
    for (j = 0; j <= s->ns2 - 1; j++)  std::fprintf(fp, "%8.2f", score(j, s->min2, s->inc2));
    std::fprintf(fp, " |    marg\n");
    std::fprintf(fp, "-------------");
    for (j = 0; j <= s->ns2 - 1; j++)  std::fprintf(fp, "--------");
    std::fprintf(fp, "----------\n");

    for (i = 0; i <= s->ns1 - 1; i++) {
      std::fprintf(fp, "%3d%8.2f |", i, score(i, s->min1, s->inc1));
      sum = 0.;
      for (j = 0; j <= s->ns2 - 1; j++) {
        std::fprintf(fp, "%8.5f", s->bp12[i][j]);
        sum += s->bp12[i][j];
      }
      std::fprintf(fp, " |%8.5f\n", sum);
      grandsum += sum;
    }

    std::fprintf(fp, "------------");
    for (j = 0; j <= s->ns2 - 1; j++)  std::fprintf(fp, "--------");
    std::fprintf(fp, "-----------\n");
    std::fprintf(fp, "     rfd2[] |");
    for (j = 0; j <= s->ns2 - 1; j++)  std::fprintf(fp, "%8.5f", static_cast<double>(s->fd2[j]) / s->n);
    std::fprintf(fp, " |%8.5f", grandsum);
    std::fprintf(fp, "\n\n");

    std::fprintf(fp, "*************");
    for (j = 0; j <= s->ns2 - 1; j++)  std::fprintf(fp, "********");
    std::fprintf(fp, "**********\n");

  }
}

/*******************************************************************/

void Print_vector(FILE *fp, char label[], double *vector, int nrows,
  char rowhead[], char colhead[])
/*
  Print a double vector

  Input
    fp         output file pointer
    label      information printed first
    vector     vector
    nrows      number of rows (i.e., elements) in vector
    rowhead    heading for rows (often "Raw Score")
    colhead    heading for column (i.e, vector)

  NOTE: data are printed in columns of width 18.  So, to make
        rowhead and colhead align with printed data, the 
        headings should be 18 characters long.  See example below.

        Print_matrix(outf,"Dummy output",vector,40,
                     "         Raw Score",
                     "            Column") 

  Function calls other than C or NR utilities: None
                                               
  R. L. Brennan

  Date of last revision: 6/30/08    
*/
{
  int i;

  std::fprintf(fp, "\n\n%s", label);
  std::fprintf(fp, "\n\n");
  std::fprintf(fp, "%s  ", rowhead);
  std::fprintf(fp, "%s", colhead);
  std::fprintf(fp, "\n");
  for (i = 0; i < nrows; i++)
    std::fprintf(fp, "\n%12d  %12.5f", i, vector[i]);
  std::fprintf(fp, "\n\n");
  for (i = 1; i <= 26; i++) std::fprintf(fp, "-");

}

/*******************************************************************/

void Print_matrix(FILE *fp, char label[], double **matrix, int nrows,
  int ncols, char rowhead[], char colheads[])
/*
  Print a double matrix

  Input
    fp         output file pointer
    label      information printed first
    matrix     matrix
    nrows      number of rows
    ncols      number of columns
    rowhead    heading for rows (often "Raw Score")
    colheads   heading for each of j=0,...,ncols columns

  NOTE: data are printed in columns of width 18.  So, to make
        rowhead and colheads align with printed data, the 
        headings should be 12 characters long.  See example below.

        Print_matrix(outf,"Dummy output",matrix,40,3,
          "         Raw Score",
          "         First Col        Second Col         Third Col") 

  Function calls other than C or NR utilities: None
                                               
  R. L. Brennan

  Date of last revision: 6/30/08    
*/
{
  int i, j;

  std::fprintf(fp, "\n\n%s", label);
  std::fprintf(fp, "\n\n");
  std::fprintf(fp, "%s  ", rowhead);
  std::fprintf(fp, "%s", colheads);
  std::fprintf(fp, "\n");
  for (i = 0; i < nrows; i++) {
    std::fprintf(fp, "\n%18d  ", i);
    for (j = 0; j < ncols; j++)
      std::fprintf(fp, "%18.5f", matrix[i][j]);
  }
  std::fprintf(fp, "\n\n");
  for (i = 1; i <= 18 * (ncols + 1) + 2; i++) std::fprintf(fp, "-");
}

/*******************************************************************/

void Print_file(char FileName[], FILE *fp)
/*
  Print contents of file named FileName to file pointed to by fp
  Author: R. L. Brennan

  Date of last revision: 3-17-09  
*/
{
  FILE *inf;
  int i, c;
  inf = std::fopen(FileName, "r");

  std::fprintf(fp, "\n\nContents of %s", FileName);
  std::fprintf(fp, "\n\n");
  for (i = 1; i <= 40; i++) std::fprintf(fp, "-");
  std::fprintf(fp, "\n");
  while ((c = std::fgetc(inf)) != EOF) std::fputc(c, fp);
  for (i = 1; i <= 40; i++) std::fprintf(fp, "-");
  std::fprintf(fp, "\n");

  std::fclose(inf);
}

/*******************************************************************************/
/*******************************************************************************/

void print_matrix(double **matrix, int nrow, int ncol, int offset)
/*
   Purpose:
      print the content of the matrix

   Input:       
     matrix		matrix whose content will be printed
	 nrow		number of rows to print
	 ncol		number of columns to print
	 offset		starting position of the matrix

   Output:                                                           
     matrix elements on the standard output

   Author: Jaehoon Seol
   Date: August 20, 2009
   Version: 1.0
   References : 
       
   Comments: 
 */
{
  int i, j;
  std::printf("\n");
  for (i = 0; i < nrow; i++)
  {
    for (j = 0; j < ncol; j++)
      std::printf("%f ", matrix[i + offset][j + offset]);
    std::printf("\n");
  }

}

/*******************************************************************************/

void print_vector(double *vect, int nrow, int offset)
/*
   Purpose:
      print the content of the vector

   Input:       
     vect		vect whose content will be printed
	 nrow		number of rows to print 
	 offset		starting position of the vector

   Output:                                                           
     vector elements on the standard output

   Author: Jaehoon Seol
   Date: August 20, 2009
   Version: 1.0
   References : 
       
   Comments: 
 */
{
  int i;
  std::printf("\n");
  for (i = 0; i < nrow; i++)
    std::printf("%f ", vect[i + offset]);

}

/*******************************************************************************/
void er_error2(char err_source[], char err_message[])
/*
   Purpose:
      print error source and error message

   Input:       
     err_source   error source
	 err_message  error message

   Output:                                                           
     error message on the standard output

   Author: Jaehoon Seol
   Date: August 20, 2009
   Version: 1.0
   References : 
       
   Comments: This has the same functionality as nrerror(). 
 */
{
  std::fprintf(stderr, "Equating Recipes error occured\n");
  std::fprintf(stderr, "Source: %s, Error: %s\n", err_source, err_message);
  std::fprintf(stderr, "Exiting the system. Good bye!\n");
  std::exit(EXIT_FAILURE);
}

/*******************************************************************************/
void release_matrix(enum datatype mode, void **matrix, long rid_low, long cid_low)
/*
   Purpose:
      release memory allocate for matrix

   Input:       
     mode      data type
     matrix    matrix whose memory is released
	 rid_low   starting row index
	 cid_low   starting column index

   Output:                                                           
     memory is allocated

   Precondition:                                                     
   
   Author: Jaehoon Seol   
   Date: August 20, 2009
   Version: 1.0
   References : 
       
   Comments: This has the same functionality as the free_xxxxxx() functions in 
   Numerical Recipes.
 */
{
  char	**matrixc = nullptr;
  int		**matrixi = nullptr;
  float	**matrixf = nullptr;
  double	**matrixd = nullptr;

  switch (mode)
  {
  case CHARACTER:
    matrixc = static_cast<char **>(matrix);
    std::free(matrixc[rid_low] + cid_low - 1);
    std::free(matrixc + rid_low - 1);
    break;
  case INTEGER:
    matrixi = static_cast<int **>(matrix);
    std::free(matrixi[rid_low] + cid_low - 1);
    std::free(matrixi + rid_low - 1);
    break;
  case FLOAT:
    matrixf = static_cast<float **>(matrix);
    std::free(matrixf[rid_low] + cid_low - 1);
    std::free(matrixf + rid_low - 1);
    break;
  case DOUBLE:
    matrixd = static_cast<double **>(matrix);
    std::free(matrixd[rid_low] + cid_low - 1);
    std::free(matrixd + rid_low - 1);
    break;
  default:
    er_error2("er_release_matrixor", "wrong input mode");
  }
}

/*******************************************************************************/
void **allocate_matrix(enum datatype mode, long rid_low, long rid_high,
  long cid_low, long cid_high)
/*
   Purpose:
      allocate memory for a matrix whose indices are rid_low..rid_high 
	  for row and cid_low..cid_high for column

   Input:       
     mode      data type 
	 rid_low   starting row index
	 rid_high  ending row index
	 cid_low   starting column index
	 cid_high  ending column index

   Output:                                                           
     pointer to the allocated memory. 

   Precondition:    

   Author: Jaehoon Seol
   Date: August 20, 2009
   Version: 1.0
   References : 
       
   Comment:
      If mode = DOUBLE, then the return value has to be cast by (double **).
	  It applies to other data type as well.

   Comment:  This has the same functionality as the Numerical Recipes functions:
             float **matrix(long nrl, long nrh, long ncl, long nch);
             double **dmatrix(long nrl, long nrh, long ncl, long nch);
             int **imatrix(long nrl, long nrh, long ncl, long nch);
             char **cmatrix(long nrl, long nrh, long ncl, long nch); 

 */
{
  long i,
    nrow = rid_high - rid_low + 1,				  /* number of rows      */
    ncol = cid_high - cid_low + 1;				  /* number of columns   */
  float	**matrixf = nullptr;
  double	**matrixd = nullptr;
  char	**matrixc = nullptr;
  int		**matrixi = nullptr;
  void	**rmatrix = nullptr;

  switch (mode)
  {
  case CHARACTER:
    /* allocate memory for row vector pointers */
    matrixc = static_cast<char **>(std::malloc((size_t)((nrow + 1) * sizeof(char*))));
    if (!matrixc)
      er_error2("allocate_matrix",
        "can not allocate memory for char pointers");
    matrixc += 1;
    matrixc -= rid_low;
    /* allocate memory for whole matrix    */
    matrixc[rid_low] = static_cast<char *>(std::malloc((size_t)((nrow*ncol + 1) * sizeof(char))));
    /* adjust pointers to point proper row */
    if (!matrixc[rid_low])
      er_error2("allocate_matrix",
        "can not allocate memory for char elements");
    /* set values of vector to 0.0 */
    std::memset(matrixc[rid_low], 0, (nrow*ncol + 1) * sizeof(char));
    matrixc[rid_low] += 1;
    matrixc[rid_low] -= cid_low;
    for (i = rid_low + 1; i <= rid_high; i++)
      matrixc[i] = matrixc[i - 1] + ncol;
    rmatrix = matrixc;
    break;
  case INTEGER:
    /* allocate memory for row vector pointers */
    matrixi = static_cast<int **>(std::malloc((size_t)((nrow + 1) * sizeof(int*))));
    if (!matrixi)
      er_error2("allocate_matrix",
        "can not allocate memory for int pointers");
    matrixi += 1;
    matrixi -= rid_low;
    /* allocate memory for whole matrix    */
    matrixi[rid_low] = static_cast<int *>(std::malloc((size_t)((nrow*ncol + 1) * sizeof(int))));
    /* adjust pointers to point proper row */
    if (!matrixi[rid_low])
      er_error2("allocate_matrix",
        "can not allocate memory for int elements");
    /* set values of vector to 0.0 */
    std::memset(matrixi[rid_low], 0, (nrow*ncol + 1) * sizeof(int));
    matrixi[rid_low] += 1;
    matrixi[rid_low] -= cid_low;
    for (i = rid_low + 1; i <= rid_high; i++)
      matrixi[i] = matrixi[i - 1] + ncol;
    rmatrix = matrixi;
    break;
  case FLOAT:
    /* allocate memory for row vector pointers */
    matrixf = static_cast<float **>(std::malloc((size_t)((nrow + 1) * sizeof(float*))));
    if (!matrixf)
      er_error2("allocate_matrix",
        "can not allocate memory for float pointers");
    matrixf += 1;
    matrixf -= rid_low;
    /* allocate memory for whole matrix    */
    matrixf[rid_low] = static_cast<float *>(std::malloc((size_t)((nrow*ncol + 1) * sizeof(float))));
    /* adjust pointers to point proper row */
    if (!matrixf[rid_low])
      er_error2("allocate_matrix",
        "can not allocate memory for float elements");
    /* set values of vector to 0.0 */
    std::memset(matrixf[rid_low], 0, (nrow*ncol + 1) * sizeof(float));
    matrixf[rid_low] += 1;
    matrixf[rid_low] -= cid_low;
    for (i = rid_low + 1; i <= rid_high; i++)
      matrixf[i] = matrixf[i - 1] + ncol;
    rmatrix = matrixf;
    break;
  case DOUBLE:
    /* allocate memory for row vector pointers */
    matrixd = static_cast<double **>(std::malloc((size_t)((nrow + 1) * sizeof(double*))));
    if (!matrixd)
      er_error2("allocate_matrix",
        "can not allocate memory for double pointers");
    matrixd += 1;
    matrixd -= rid_low;
    /* allocate memory for whole matrix    */
    matrixd[rid_low] = static_cast<double *>(std::malloc((size_t)((nrow*ncol + 1) * sizeof(double))));
    /* adjust pointers to point proper row */
    if (!matrixd[rid_low])
      er_error2("allocate_matrix",
        "can not allocate memory for double elements");
    /* set values of vector to 0.0 */
    std::memset(matrixd[rid_low], 0, (nrow*ncol + 1) * sizeof(double));
    matrixd[rid_low] += 1;
    matrixd[rid_low] -= cid_low;
    for (i = rid_low + 1; i <= rid_high; i++)
      matrixd[i] = matrixd[i - 1] + ncol;
    return matrixd;
  default:
    er_error2("allocate_matrix", "wrong input mode");
  }
  return rmatrix;
}

/*******************************************************************************/
void er_scale(double vect1[], double scale, double vect2[], int n, int offset)
/*
   Purpose:
      Performs the scale BLAS operation:

			vect1 = scale*vect2

   Input:       
     vect1		target vector 
	 scale		scaling factor
     vect2		source vector 
	 n			size of vect1 and vect2
	 offset 	starting address of vect1[] and vect2[]
	                 = 0 if 0-offset is used.
					 = 1 if 1-offset is used.

   Output:                                                           
     vect1

   Precondition:                                                     
   
   Author: Jaehoon Seol    
   Date: August 20, 2009
   Version: 1.0
   References : 
       
   Comments:
      If mode = DOUBLE, then the return value has to be cast by (double **).
	  It applies to other data type as well.
 */
{
  int i;												   /* loop index */

  /* assert : offset == 0 or 1 */
  if (offset != 0 && offset != 1)
  {
    std::fprintf(stderr, "Equating Recipes error occured\n");
    std::fprintf(stderr, "Source: %s, Error: %s\n", "er_scale", "wrong offset info");
    std::fprintf(stderr, "Exiting the system. Good bye!\n");
    std::exit(EXIT_FAILURE);
  }

  for (i = offset; i < n + offset; i++)
    vect1[i] = scale * vect2[i];
}

/*******************************************************************************/
double er_dot(double vect1[], double vect2[], int n, int offset)
/*
   Purpose:
      Computes dot product of two vectors vect1[] and vect2[].

   Input:    
     vect1           input vector
	 vect2			 input vector
	 n               number of elements in vect1[] and vect2[]
	 offset 	     starting address of vect1[] and vect2[]
	                 = 0 if 0-offset is used.
					 = 1 if 1-offset is used.

   Output:     
	 returns dot product of vect1[] and vect2[] 

   Precondition:     

   Author: Jaehoon Seol  
   Date: August 20, 2009
   Version: 1.0
   References :  
   Comments: 
 */
{
  int i;												   /* loop index */
  double temp = 0.0;

  /* assert : offset == 0 or 1 */
  if (offset != 0 && offset != 1)
  {
    std::fprintf(stderr, "Equating Recipes error occured\n");
    std::fprintf(stderr, "Source: %s, Error: %s\n", "er_dot", "wrong offset info");
    std::fprintf(stderr, "Exiting the system. Good bye!\n");
    std::exit(EXIT_FAILURE);

  }

  for (i = offset; i < n + offset; i++)
    temp += vect1[i] * vect2[i];
  return temp;
}

/*******************************************************************************/
void er_daxpy(double vectY[], double alpha, double vectX[], int n, int offset)
/*
   Purpose:
      Implementation of daxy operation:

	           y = alpha*x +y   -------- (1) 

   Input :     
     vectY           vector Y
     alpha           scaling factor
     vectX           vector X
	 n               number of elements in vector X and Y 
	 offset 	     starting address of vect and matx
	                 = 0 if 0-offset is used.
					 = 1 if 1-offset is used.

   Output : 
     vectY           vector Y storing (1) 

   Precondition:     

   Author: Jaehoon Seol  
   Date: August 20, 2009
   Version: 1.0
   References :  
   Comments: 
 */
{
  int i;
  /* assert : offset == 0 or 1 */
  if (offset != 0 && offset != 1)
  {
    std::fprintf(stderr, "Equating Recipes error occured\n");
    std::fprintf(stderr, "Source: %s, Error: %s\n", "er_daxpy", "wrong offset info");
    std::fprintf(stderr, "Exiting the system. Good bye!\n");
    std::exit(EXIT_FAILURE);
  }

  for (i = offset; i < n + offset; i++)
    vectY[i] += alpha * vectX[i];

}

/*******************************************************************************/
void er_r1update(double **matx, double scale, double vect[], int n, int offset)
/*
   Purpose:
      Computes scaled rank 1 update of matx using scale and vect, i.e.,

	           matx = matx + scale * vect* vect^t

   Input :     
     matx            contains initial matrix information
     scale           scaling factor
	 vect			 vector used for rank 1 update
	 n               dimension infomration for matx and vect
	 offset 	     starting address of vect and matx
	                 = 0 if 0-offset is used.
					 = 1 if 1-offset is used.

   Output : 
	 matxh           final computation of (1) is stored in matxh. 

   Precondition:     

   Author: Jaehoon Seol  
   Date: August 20, 2009
   Version: 1.0
   References :  
   Comments: 
 */
{
  int i, j;											   /* loop index */

  /* assert : offset == 0 or 1 */
  if (offset != 0 && offset != 1)
  {
    std::fprintf(stderr, "Equating Recipes error occured\n");
    std::fprintf(stderr, "Source: %s, Error: %s\n", "er_r1update", "wrong offset info");
    std::fprintf(stderr, "Exiting the system. Good bye!\n");
    std::exit(EXIT_FAILURE);
  }

  /* rank 1 update */
  for (i = offset; i < n + offset; i++)
    for (j = offset; j < n + offset; j++)
      matx[i][j] += scale * vect[i] * vect[j];
}

/*******************************************************************************/
void er_mvmult(double vect2[], double **matx, double vect1[], int n, int offset)
/*
   Purpose:
      Computes matrix vector multiplication:

	           vect2 = matx * vect1    -----------------(1)

   Input :     
     matx            contains initial matrix information 
	 vect1			 vector used for multiplication
	 n               dimension of vect1
	 offset 	     starting address of vect and matx
	                 = 0 if 0-offset is used.
					 = 1 if 1-offset is used.

   Output : 
	 vect2           final computation of (1) is stored in vect2. 

   Precondition:     

   Author: Jaehoon Seol  
   Date: August 20, 2009
   Version: 1.0
   References :  
   Comments: 
 */
{
  int i, j;											   /* loop index */

  for (i = offset; i < n + offset; i++)
  {
    vect2[i] = 0.0;
    for (j = offset; j < n + offset; j++)
      vect2[i] += matx[i][j] * vect1[j];
  }

}

/*******************************************************************************/
float er_random(long *seed)
/*
   Purpose:
      Generate random numbers

   Input :     
     seed            seed for random number generation 

   Output : 
	 seed            modified seed information
	 er_random( )	 returns a random number

   Precondition:     

   Author: Jaehoon Seol  
   Date: August 20, 2009  
   Version: 1.0
   References :  
   Comments: 
 */
{
  int		i;
  static int		id = 0;
  long	q,												/* quotient  */
    r,												/* remainder */
    t,
    seed1,
    ltemp;
  static long seed2 = 1;
  static long state[32];
  float rvalue;

  seed1 = *seed;
  if (seed1 == 0)
  {
    *seed = 1454538876;
    return 0.28538090;
  }
  if (seed1 < 0) {
    seed1 = std::abs(seed1);
    seed2 = seed1;
    for (i = 39; i >= 0; i--) {
      q = seed1 / 53668;
      r = seed1 - q * 53668;
      t = 40014 * r - q * 12211;
      seed1 = (t < 0 ? t + 2147483563 : t);
      if (i < 32) state[i] = seed1;
    }
    id = state[0] / 67108862;
  }
  q = seed1 / 53668;
  r = seed1 - q * 53668;
  t = 40014 * r - q * 12211;
  seed1 = (t < 0 ? t + 2147483563 : t);
  q = seed2 / 52774;
  r = seed2 - q * 52774;
  t = 40692 * r - q * 3791;
  seed2 = (t < 0 ? t + 2147483399 : t);
  ltemp = state[id] - seed2;
  ltemp = (ltemp < 1 ? ltemp + 2147483562 : ltemp);
  state[id] = seed1;
  id = ltemp / 67108862;
  rvalue = 4.656613057391769e-10*ltemp;
  *seed = seed1;
  return rvalue;
}

/*******************************************************************************/
void er_ludcmp(double **a, int n, int *pivot, double *det)
/*
   Purpose:
      Implementation of Gaussian elimination with partial pivoting. 
      Store L and U over A.

   Input:                                  
     a             n x n matrix to be factored.
     n             dimension of the matrix

   Output:                
     a             L U matrices overwrite a matrix
     pivot         n dimensional vector containing pivot information
	 

   Precondition:                                                     
     None 
   
   Author: Jaehoon Seol
   Date: August 20, 2009
   Version: 1.0
   References : 
      J. Stoer & R. Bulirsch, Introduction to Numerical analysis, 2nd edition 
	  David S. Watkins, Fundamentals of Matrix Computations,2002
   Comments:
 */
{
  int i, k,                                               /* row index  */
    j,                                                 /* col. index */
    maxid;                            /* index of row with max value */
  double temp,
    sign = 1.0,                    /* sign of the det(a)          */
    maxvalue;

  for (k = 1; k < n; k++)
  {
    /* find max among a[k..n,k] */
    maxvalue = std::fabs(a[k][k]);
    maxid = k;
    for (i = k + 1; i <= n; i++)
      if (std::fabs(a[i][k]) > maxvalue)
      {
        maxvalue = std::fabs(a[i][k]);
        maxid = i;
      }
    if (maxvalue == 0.0)
    {
      std::fprintf(stderr, "Equating Recipes error occured\n");
      std::fprintf(stderr, "Source: %s at 1, Error: %s \n", "er_ludcmp", "singular matrix");
      std::fprintf(stderr, "k = %d \n", k);
      std::fprintf(stderr, "Exiting the system. Good bye!\n");
      std::exit(EXIT_FAILURE);
    }
    else
    {
      pivot[k] = maxid;                 /* record row interchange  */
      /* swap row k and row maxid */
      if (k != maxid)
      {
        for (j = 1; j <= n; j++)
        {
          temp = a[k][j];
          a[k][j] = a[maxid][j];
          a[maxid][j] = temp;
        }
        sign *= -1.0;
      }
      /* calculate multipliers */
      for (i = k + 1; i <= n; i++)
        a[i][k] = a[i][k] / a[k][k];
      for (i = k + 1; i <= n; i++)
        for (j = k + 1; j <= n; j++)
          a[i][j] -= a[k][j] * a[i][k];
    }
  }
  if (a[n][n] == 0.0)
  {
    std::fprintf(stderr, "Equating Recipes error occured\n");
    std::fprintf(stderr, "Source: %s at 2, Error: %s\n", "er_ludcmp", "singluar matrix");
    std::fprintf(stderr, "Exiting the system. Good bye!\n");
    std::exit(EXIT_FAILURE);
  }
  else
    pivot[n] = n;

  /* compute determinant of a */
  temp = 1.0;
  for (i = 1; i <= n; i++) temp *= a[i][i];
  *det = sign * temp;
}

/*******************************************************************************/
void    er_lubksb(double **a, int n, int pivot[], double b[])
/*
   Purpose:
     Solves the following matrix equation:
             
               a  x = b
   
     Matrix a contains L U factorization information that satisfies

               P a = L U

   Input:                                  
      a             n x n matrix to be factored.
      n             dimension of the matrix 
      pivot         n dimensional vector containing pivot information
      b             right hand side of the equation

   Output:                
      b             solution x over-write vector b
	 
   Precondition:                                                     
     None 
 
   Author: Jaehoon Seol
   Date: August 20, 2009
   Version: 1.0
   References : 
      J. Stoer & R. Bulirsch, Introduction to Numerical analysis, 2nd edition 
	  David S. Watkins, Fundamentals of Matrix Computations,2002
   Comments:
 */
{
  int i, j, k;
  double sum,
    temp;
  /* compute P^t b */
  for (k = 1; k < n; k++)
  {
    temp = b[pivot[k]];
    b[pivot[k]] = b[k];
    b[k] = temp;
  }
  /* forward substitution */
  for (i = 1; i <= n; i++)
  {
    sum = 0.0;
    for (j = 1; j < i; j++)
      sum += a[i][j] * b[j];
    b[i] -= sum;
  }
  /* backward substitution */
  for (i = n; i >= 1; i--)
  {
    temp = a[i][i];
    if (temp == 0.0)
    {
      std::fprintf(stderr, "Equating Recipes error occured\n");
      std::fprintf(stderr, "Source: %s, Error: %s\n",
        "er_solve", "input matrix is singular");
      std::fprintf(stderr, "Exiting the system. Good bye!\n");
      std::exit(EXIT_FAILURE);
    }
    sum = 0.0;
    for (j = i + 1; j <= n; j++)
      sum += a[i][j] * b[j];
    b[i] = (b[i] - sum) / temp;
  }
}

/*******************************************************************************/
void er_matrix_inverse(int n, double **a)
/*
   Purpose:
      Find the inverse of the matrix a using the Gauss-Jordan 
      algorithm with partial pivoting. The original matrix a[][]
      is over-written by its inverse. 

   Input:                                                            
     a    matrix whose inverse we are trying to find.
     n    dimension of the matrix

   Output:                                                           
     a    On exit, matrix a is overwritten by its inverse.

   Precondition:                                                     
     matrix a is non-singular. 
 
   Author: Jaehoon Seol
   Date: August 20, 2009
   Version: 1.0
   References : 
      J. Stoer & R. Bulirsch, Introduction to Numerical analysis, 2nd edition 
   Comments:
 */
{
  int i, j, k, r,                              /* loop index           */
    itemp;                                   /* temp. variable       */
  double max,                       /* max value on the current column */
    temp;
  int    *pivot = static_cast<int *>(std::malloc(n * sizeof(int)));   /* vector for pivot  */
  double *b = static_cast<double *>(std::malloc(n * sizeof(double)));/* temp. vector      */

  for (i = 0; i < n; i++) pivot[i] = i;
  for (j = 0; j < n; j++)
  {
    /* search for pivot */
    max = std::fabs(a[j][j]);
    r = j;
    for (i = j + 1; i < n; i++)
      if (std::fabs(a[i][j]) > max)
      {
        max = std::fabs(a[i][j]);
        r = i;
      }
    /* input matrix is singular */
    if (max == 0.0)
    {
      std::fprintf(stderr, "Equating Recipes error occured\n");
      std::fprintf(stderr, "Source: %s, Error: %s\n",
        "er_gauss", "input matrix is singular");
      std::fprintf(stderr, "Exiting the system. Good bye!\n");
      std::exit(EXIT_FAILURE);
    }
    /* interchange rows if r>j */
    if (r > j)
    {
      /* swap a[j][.] with a[r][.] */
      for (k = 0; k < n; k++)
      {
        temp = a[j][k];
        a[j][k] = a[r][k];
        a[r][k] = temp;
      }
      /* update pivot information */
      itemp = pivot[j];
      pivot[j] = pivot[r];
      pivot[r] = itemp;
    }
    /* transformation        */
    temp = 1.0 / a[j][j];
    for (i = 0; i < n; i++) a[i][j] = a[i][j] * temp;
    a[j][j] = temp;
    for (k = 0; k < n; k++)
    {
      if (k == j) continue;
      for (i = 0; i < n; i++)
      {
        if (i == j) continue;
        a[i][k] = a[i][k] - a[i][j] * a[j][k];
      }
      a[j][k] = -temp * a[j][k];
    }
  }
  /* column interchange */
  for (i = 0; i < n; i++)
  {
    for (k = 0; k < n; k++) b[pivot[k]] = a[i][k];
    for (k = 0; k < n; k++) a[i][k] = b[k];
  }
  std::free(pivot);
  std::free(b);
}

/*******************************************************************************/

void er_sort(float* vector, int l, int r)
/*
   Purpose:
      sort vector from index left to index right

   Input:                                                            
     vector     vector to be sorted.
	 left       starting index
	 right      ending index

   Output:                                                           
     vector

   Author: Jaehoon Seol
   Date: August 20, 2009
   Precondition:                                                     
     None 
  
   Comments:
 */
{
  int i;
  float temp;
  for (i = l + 1; i <= r; i++)
    if (vector[i] < vector[l])
    {
      temp = vector[l];
      vector[l] = vector[i];
      vector[i] = temp;
    }
  for (i = l + 2; i <= r; i++)
  {
    int j = i;
    float v = vector[i];
    while (v < vector[j - 1])
    {
      vector[j] = vector[j - 1];
      j--;
    }
    vector[j] = v;
  }
}

/*******************************************************************************/
void er_qrdcmp(double **a, int nrow, int ncol, double *coeff, double *diag)
/*
   Purpose:
      Performs QR decomposition of the matrix a.

   Input:                                                            
     a			matrix whose whose QR decomposition is to be computed
     nrow		number of rows of a
	 ncol		number of columns of a
	 coeff		vector containing 1/2*(u*u) where u is the vector used
	            to construct Q matrix
	 diag		vector containing diagonal element of R matrix

   Output:                                                           
     a			On exit, upper triangular matrix a is R

   Precondition:                 
 
   Author: Jaehoon Seol
   Date: August 20, 2009
   Version: 1.0
   References : 
      J. Stoer & R. Bulirsch, Introduction to Numerical analysis, 2nd edition 
	  J. Wilkinson and C. Reinsch, Handbook for Automatic Computation: Volume 2: Linear Algebra 
   Comments:
 */
{
  int		i, j, k,
    ncol2;
  double	delem,
    melem,
    sigma,
    sum,
    t1;
  const double tol = 1.0e-13;

  ncol2 = ncol;
  if (nrow == ncol) ncol2 = ncol - 1;
  std::memset(diag, 0, ncol2 * sizeof(double));
  std::memset(coeff, 0, ncol2 * sizeof(double));
  for (k = 1; k <= ncol2; k++) {
    /* find max. element on the column */
    delem = a[k][k];
    melem = std::fabs(delem);
    for (i = k + 1; i <= nrow; i++)
      if (melem < std::fabs(a[i][k]))
        melem = std::fabs(a[i][k]);
    if (melem > tol) {
      for (i = k; i <= nrow; i++)
        a[i][k] /= melem;
      sum = 0.0;
      for (i = k; i <= nrow; i++)
        sum += std::pow(a[i][k], 2.0);
      if (delem > 0.0)
        sigma = std::sqrt(sum);
      else
        sigma = -std::sqrt(sum);
      a[k][k] += sigma;
      coeff[k] = sigma * a[k][k];
      diag[k] = -melem * sigma;
      /* update a matrix */
      for (j = k + 1; j <= ncol; j++) {
        t1 = 0.0;
        for (i = k; i <= nrow; i++)
          t1 += a[i][k] * a[i][j] / coeff[k];
        for (i = k; i <= nrow; i++)
          a[i][j] -= t1 * a[i][k];
      }
    }
  }
  diag[nrow] = a[nrow][ncol];
}

/*******************************************************************************/
double er_rtsafe(void(*funcd)(double, double *, double *), double x0, double x1, double error)
/*
   Purpose:  
      This function implements bisection method to find x such that

           f(x)==0.0.

   Input:             
      funcd   function pointer whose root we are trying to find.
	  x0      left end of the function domain 
	  x1      right end of the function domain

   Output:                
      x that satisfies f(x)==0.0
	 
   Precondition:                                                     
     f(x) is defined on [x0,x1] with f(x0)*f(x1) < = 0.
 
   Author: Jaehoon Seol
   Date: August 20, 2009
   Version: 1.0
   References : 
      J. Stoer & R. Bulirsch, Introduction to Numerical analysis, 2nd edition 
	  David S. Watkins, Fundamentals of Matrix Computations,2002
   Comments:
 */
{
  double diff,
    left,
    right,
    side1,
    side2,
    mid,
    temp;

  (*funcd)(x0, &left, &temp);
  (*funcd)(x1, &right, &temp);
  /* assert that there is a root between x1 and x2 */
  if (left*right > 0)               /* function values have the same sign */
  {
    std::fprintf(stderr, "Equating Recipes error occured\n");
    std::fprintf(stderr, "Source: %s, Error: %s\n", "er_find_root",
      "no root found exists on the interval");
    std::fprintf(stderr, "Exiting the system. Good bye!\n");
    std::exit(EXIT_FAILURE);
  }

  diff = std::fabs(right - left);        /* difference of the interval         */
  while (diff > error)
  {
    mid = (left + right) / 2.0;
    (*funcd)(mid, &side1, &temp);
    (*funcd)(right, &side2, &temp);
    if (side1*side2 <= 0)
      left = mid;
    else
      right = mid;
    diff = std::fabs(right - left);
  }
  return mid;
}

/*******************************************************************************/
int er_lnsrch(double xold[], int n, double gk[], double sk[], double maxstep,
  double(*ftn_ptr)(double[]), double xnew[])
/*
   Purpose:
      Find lamda so that f(xold + lamda*sk) has decreased sufficiently. To check if
	  the decrease is sufficient, we use the following condition

   Input:    
     xold[1..n]      n-dimensional point
     n               dimension of the array xold, gk, sk, xnew, f
	 fvalue_old      f(xold)
	 gk[1..n]
	 sk[1..n]    the Newton direction
	 xnew[1..n]      new point along the direction sk from xold
	                 where the funtion ftn_ptr has decreased sufficiently.
	 maxstep         maximum step length used to keep the function outside regions 
	                 where it is undefined or subject to overflow.
	 ftn_ptr         pointer to an input function f( )  
	  

   Output:     
	 fvalue_new      new function value
	 return value = 0 if normal exit
	              = 1 if xnew is too close to xold
   Precondition:     

   Author: Jaehoon Seol
   Date: August 20, 2009
   Version: 1.0
   References :  
   Comments: 
      For detailed explanation of the algorithm, refer to p384-p385, NR.
 */
{
  int i,
    rcode;				 /* return code value                        */
  /* rcode = 0   on normal exit               */
  /* = 1   xnew is too close to xold    */
  double a,				 /* coefficient  of gk(lamda)                */
    b,				 /* coefficient  of gk(lamda)                */
    lamda0,			 /* temp. lamda                              */
    lamda1,
    lamda2,
    lamda1p2,
    lamda2p2,
    dgx0,
    scale,
    fvalue,
    fvalue2,
    fvalue_old,		 /* old function value                       */
    temp1,			 /* temp. used to compute a                  */
    temp2,			 /* temp. used to compute b                  */
    temp,			 /* temp. variable                           */
    nrm2 = 0.0,		 /* L2 norm of sk                            */
    test = 0.0;
  double error = 2.0*std::pow(10.0, -6.0);		 /* tolerance level           */
  const double alpha = std::pow(10.0, -7.0);	 /* constant alpha in 9.7.7   */


  rcode = 0;								 /* normal exit               */
  nrm2 = std::sqrt(er_dot(sk, sk, n, 1));			 /* compute L2-norm of sk     */
  if (nrm2 > maxstep)						 /* normalize sk vector       */
    er_scale(sk, maxstep / nrm2, sk, n, 1);
  dgx0 = er_dot(gk, sk, n, 1);
  for (i = 1; i <= n; i++)
  {
    temp = (std::fabs(xold[i]) > 1.0 ? std::fabs(xold[i]) : 1.0);
    temp = std::fabs(sk[i]) / temp;
    if (temp > test) test = temp;
  }
  error /= test;

  lamda1 = 1.0;
  fvalue_old = (*ftn_ptr)(xold);
  fvalue = fvalue_old + alpha * lamda1*dgx0 + 1.0;
  while (fvalue > fvalue_old + alpha * lamda1*dgx0)
  {
    std::memcpy(xnew + 1, xold + 1, n * sizeof(double));
    er_daxpy(xnew, lamda1, sk, n, 1);
    fvalue = (*ftn_ptr)(xnew);
    /* xnew is too close to xold   */
    if (lamda1 < error)
    {
      std::memcpy(xnew + 1, xold + 1, n * sizeof(double));
      rcode = 1;
      break;
    }
    /* backtracking: compute lamda0 */
    if (lamda1 != 1.0)
    {
      /* compute a & b in Eq. 9.7.13 */
      temp1 = fvalue - fvalue_old - lamda1 * dgx0;
      temp2 = fvalue2 - fvalue_old - lamda2 * dgx0;
      scale = 1.0 / (lamda1 - lamda2);
      lamda1p2 = std::pow(lamda1, -2.0);
      lamda2p2 = std::pow(lamda2, -2.0);
      a = scale * (temp1*lamda1p2 - temp2 * lamda2p2);
      b = scale * (temp2*lamda1*lamda2p2 - temp1 * lamda2*lamda1p2);
      if (a == 0.0)						 /* use Eq. 9.7.11       */
        lamda0 = -dgx0 / (2.0*b);
      else								 /* use Eq. 9.7.14       */
      {
        /* check b*b-3.0*a*gk'(0) >= 0.0 */
        temp = std::pow(b, 2.0) - 3.0*a*dgx0;
        if (temp < 0.0)
        {
          std::fprintf(stderr, "Equating Recipes error occured\n");
          std::fprintf(stderr, "Source: %s, Error: %s\n", "er_lnsrch", "Roundoff problem");
          std::fprintf(stderr, "Exiting the system. Good bye!\n");
          std::exit(EXIT_FAILURE);
        }
        lamda0 = -b + std::sqrt(temp);
        lamda0 /= (3.0*a);
      }
      if (lamda0 > 0.5*lamda1) lamda0 = 0.5*lamda1;
    }
    else
      lamda0 = -dgx0 / (2.0*(fvalue - fvalue_old - dgx0));
    fvalue2 = fvalue;
    lamda2 = lamda1;
    lamda1 = (lamda0 > 0.1*lamda1 ? lamda0 : 0.1*lamda1);
  }
  return rcode;
}

/*******************************************************************************/
void er_dfpmin(double xold[], int n, double error, int *numiter, double *fvalue,
  double(*ftn_ptr)(double[]), void(*dftn_ptr)(double[], double[]))
/*
   Purpose:
      This is an implementation of BFGS that solves the minimization problem for a real 
	  function f : R^n -> R of n variable.

   Input:    
     xold[1..n]      n-dimensional starting point 
     n               dimension of the array startx
	 error           tolerance level
	 ftp_ptr         function pointer to f( )
	 dftn_ptr        function pointer to D[f( )]
	  

   Output:     
	 xold[1..n]        xold is over-written by the final minimizer. 

   Precondition:     

   Author: Jaehoon Seol  
   Date: August 20, 2009
   Version: 1.0
   References :  
   Comments: 
 */
{
  int    i,
    num = 0;					 /* counts number of iterations      */
  double scale1,
    scale2,
    maxstep,
    nrm2;
  double *gk,						  /* gk = D[f(x_(k+1))-f(x_k)]       */
    *qk,						  /* qk = g_(k+1)-g_k                */
    *hddf,					  /* hddf = hessian * qk             */
    **matxh,					  /* hessian matrix                  */
    *xnew,					  /* x_(k+1)                         */
    *sk,						  /* sk = Hk * gk, search direction  */
    *pk;					      /* pk = xnew-xold                  */

  pk = static_cast<double *>(std::malloc((n + 1) * sizeof(double)));
  qk = static_cast<double *>(std::malloc((n + 1) * sizeof(double)));
  gk = static_cast<double *>(std::malloc((n + 1) * sizeof(double)));
  sk = static_cast<double *>(std::malloc((n + 1) * sizeof(double)));
  hddf = static_cast<double *>(std::malloc((n + 1) * sizeof(double)));
  xnew = static_cast<double *>(std::malloc((n + 1) * sizeof(double)));
  matxh = static_cast<double **>(allocate_matrix(DOUBLE, 1, n, 1, n));

  (*dftn_ptr)(xold, gk);
  /* set hessian matrix by identity matrix */
  for (i = 1; i <= n; i++) matxh[i][i] = 1.0;
  nrm2 = std::sqrt(er_dot(xold, xold, n, 1));
  maxstep = (nrm2 > n ? 100.0*nrm2 : 100.0*n);
  while (1)
  {
    er_mvmult(sk, matxh, gk, n, 1);              /* compute sk = Hk*gk   */
    er_scale(sk, -1.0, sk, n, 1);
    er_lnsrch(xold, n, gk, sk, maxstep, ftn_ptr, xnew); /* determine x(k+1)*/
    std::memcpy(sk + 1, xnew + 1, n * sizeof(double));         /* update pk       */
    er_daxpy(sk, -1.0, xold, n, 1);
    *fvalue = (*ftn_ptr)(xold);
    std::memcpy(xold + 1, xnew + 1, n * sizeof(double));       /* update xold     */
    std::memcpy(qk + 1, gk + 1, n * sizeof(double));          /* update qk and gk */
    er_scale(qk, -1.0, qk, n, 1);                    /* " "          */
    (*dftn_ptr)(xold, gk);                        /* " "          */
    er_daxpy(qk, 1.0, gk, n, 1);                     /* " "          */
    /* stopping condition */
    if (std::sqrt(er_dot(sk, sk, n, 1)) < std::pow(10.0, -8.0) ||
      std::sqrt(er_dot(gk, gk, n, 1)) < error)
      break;
    er_mvmult(hddf, matxh, qk, n, 1);         /* update hessian matrix   */
    scale1 = er_dot(qk, sk, n, 1);
    scale2 = er_dot(qk, hddf, n, 1);
    if (scale1*scale1 > 3.0e-8*er_dot(qk, qk, n, 1)*er_dot(sk, sk, n, 1))
    {
      er_scale(qk, 1.0 / scale1, sk, n, 1);
      er_daxpy(qk, -1.0 / scale2, hddf, n, 1);
      er_r1update(matxh, 1.0 / scale1, sk, n, 1);
      er_r1update(matxh, -1.0 / scale2, hddf, n, 1);
      er_r1update(matxh, scale2, qk, n, 1);
    }
    num++;                            /* count number of iterations */
  }
  *numiter = num;
  /* release all allocated memory */
  std::free(pk);
  std::free(qk);
  std::free(gk);
  std::free(sk);
  std::free(hddf);
  std::free(xnew);
  release_matrix(DOUBLE, matxh, 1, 1);
}