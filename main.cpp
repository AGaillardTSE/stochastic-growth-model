/*************************************************/
/**          Alexandre GAILLARD - 2017          **/
/**   the stochastic neoclassical growth model  **/
/**                    2018                     **/
/*************************************************/


// DECISION RULE: endogenous grid method [CAROLL, C. (2006)] //
// SIMULATION: use of histogram method //




#include "header.hpp"
#include "useful.cpp"
#include "tauchen.cpp"
#include "POLICY.cpp"








// MAIN //
int main(int argc, char* argv[])
{

	
// INITIALIZATION //
int	i,y,k;


// STEADY STATE VALUES //
double k_ss, y_ss, c_ss;

k_ss = pow((1/beta - 1 + delta)/alpha, 1/(alpha-1));
y_ss = pow(k_ss, alpha);
c_ss = pow(k_ss, alpha) - delta*k_ss;


// K GRID //
double coverage = 0.25;
double step = 2*coverage*k_ss/(maxigrid - 1);

K[0] = (1-coverage)*k_ss;

for(i = 1; i<maxigrid; i++) {
    K[i] = K[i-1] + step;
}


// TRANSITION + STATE Z OF ENTREPRENEURS //
const double p_e = rho;
const double std_e = sigma;
const double m_e = 3;
const double sigma_e = std_e*(pow((1-p_e*p_e),0.5));
tauchenfun(p_e, m_e, 0.0, std_e, prod, ytrans);

//inv_distri(yinv, ytrans);
//
//
//for(int y = 0; y < maxygrid; y++){
//    //prod[y] = exp(prod[y]);
//    printf("%f\t", prod[y]);
//}
//printf("\n");
//printf("\n");
//for(int y = 0; y < maxygrid; y++){
//    for(int k = 0; k < maxygrid; k++){
//        printf("%f\t", ytrans[y][k]);
//    }
//    printf("\n");
//}
//printf("\n");
//
//
//for(int y = 0; y < maxygrid; y++){
//    printf("%f\t", yinv[y]);
//}
//printf("\n");



// GRID FOR Y VALUE (CASH ON HAND TODAY) //
for(i=0;i<maxigrid;i++)
{
    for(y=0;y<maxygrid;y++)
    {
        YY[inx(i,y)] = exp(prod[y])*pow(K[i], alpha) + (1-delta)*K[i];
        
        //printf("%f", YY[inx(i,y)]);getchar();
    }
}





// TIME //
time_t rawtime,timeofstart,timeofend;
struct tm * timeinfo;
time ( &rawtime );
time (&timeofstart);
timeinfo = localtime ( &rawtime );


// MARGINAL UTILITY, VALUES FUNCTION AND POLICIES //
double *VF, *cons, *optiK;          // for decision rules

// Note for users :: please, always use pointers and save your computer's memory ;) == banish all arrays //
VF = (double *) calloc((ifulldim), sizeof(double));      // value function
cons = (double *) calloc((ifulldim), sizeof(double));
optiK = (double *) calloc((ifulldim), sizeof(double));


// START BY GUESSING VF //
double VFstart = (1/(1-beta))*U(c_ss);

for(i = 0; i < maxigrid; i++){
    for(y = 0; y < maxygrid; y++){
        VF[inx(i,y)] = VFstart + i/(double)maxigrid;      // REQUIERE TO BE INCREASING IN K (the case here)
    }
}

printf("STARTING COMPUTATION\n");

POLICY_EGM(VF,cons, optiK);




time ( &rawtime );
time (&timeofend);
timeinfo = localtime ( &rawtime );


printf("=====Program ended: %s in %f seconds\n",asctime(timeinfo),difftime(timeofend,timeofstart));

return 0;

}


