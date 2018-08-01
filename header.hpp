/**************************************************/
/**          Alexandre GAILLARD - 2017           **/
/**   the stochastic neoclassical growth model   **/
/**                    2018                      **/
/**************************************************/


/****************/
//    INCLUDE   //
/****************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <limits>
#include <assert.h>
#include <cstring>
#include <fstream>
#include <iostream>

#define deriv(val1,val2,val3,x1,x2,x3) ((1.0 - (x3 - x2)/(x3 - x1))*((val3 - val2)/(x3-x2)) + ((x3 - x2)/(x3 - x1))*((val2 - val1)/(x2-x1)))



// INDEX //
#define maxigrid 1000 // define the grid of saving (next period wealth)
#define maxygrid 10

#define ifulldim (maxigrid*maxygrid)
#define inx(igridindex,jclassindex) (((jclassindex)*(maxigrid))+(igridindex))
#define linspace(x0,xmax,n,i) ((i)*(((xmax)-(x0))/(n))+(x0))


// EXOGENOUS GRID //
double K[maxigrid];
double YY[ifulldim];
double prod[maxygrid];
double ytrans[maxygrid][maxygrid];
double yinv[maxygrid];



/***************************/
//  CALIBRATION DEFINITION //
/***************************/

// Government steady-state //
const double beta = 0.9896;
const double tau = 2;
const double theta = 0.357;
const double alpha = 0.4;
const double delta = 0.0196;
const double rho = 0.95;
const double sigma = 0.007;


// UTILITY//
#define MUc(x) (pow((x),-tau))
#define inv_MU(u) (pow((u),(-(1/tau))))
#define U(x) (pow((x),(1.0-tau))/(1.0-tau))



//Convergence criterion
const double epsilon=	0.0000001; //Convergence criterion on labor supply (and other stuff)



//Output files
const char policyfile[]="policy.out";
const char distfile[]="dist.out";
const char valuefile[]="value.out";
