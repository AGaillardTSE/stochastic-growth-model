/*************************************************/
/**          Alexandre GAILLARD - 2017          **/
/**   the stochastic neoclassical growth model  **/
/**                    2018                     **/
/*************************************************/



// POLICY ITERATION //
void POLICY_EGM(double *eVFnew, double *cons, double *OPTI_K)
{


// INTEGER //
int i,ii,y,ynext,iter;


// INITIALIZATION //
double *VFendo, *VF, *Yendo, *eVF, *deVF, critere, weight, slope, Ycase, slope2;
VFendo = (double *) calloc((ifulldim), sizeof(double));     // Value function on the endogenous grid (is defined on a current grid -> same dimension of the other space).
Yendo = (double *) calloc((ifulldim), sizeof(double));       // endogenous grid values
eVF = (double *) calloc((ifulldim), sizeof(double));         // expected value function
deVF = (double *) calloc((ifulldim), sizeof(double));        // derivative of the expected value function
VF = (double *) calloc((ifulldim), sizeof(double));        // derivative of the expected value function

// START LOOP OVER DECISION RULES //
critere=1.0;
iter=0;


while (critere > epsilon) 
{

    // copy vector //
    bascule(eVFnew,eVF,ifulldim);
 
    for(y = 0; y < maxygrid; y++){
        for(i = 0; i < maxigrid; i++){
        
            /** 1. APPROXIMATE THE DERIVATIVE **/
            // At the corner: use the next point
            // in Between: use the weighted average between the two consecutives linear slopes.
            if(i == 0){deVF[inx(i,y)] = (eVF[inx(1,y)]-eVF[inx(0,y)])/(K[1] - K[0]);}
            if(i == (maxigrid-1)){deVF[inx(i,y)] = (eVF[inx((maxigrid-1),y)]-eVF[inx((maxigrid-2),y)])/(K[(maxigrid-1)] - K[(maxigrid-2)]);}
            
            // standard case //
            //if(i > 0 && i < (maxigrid-1)){deVF[inx(i,y)] = deriv(eVF[inx((i-1),y)],eVF[inx(i,y)],eVF[inx((i+1),y)],K[(i-1)],K[i],K[(i+1)]);}
            if(i > 0 && i < (maxigrid-1)){deVF[inx(i,y)] = (eVF[inx((i+1),y)]-eVF[inx((i-1),y)])/(K[(i+1)] - K[(i-1)]);}
            
            
            /** 2 COMPUTE THE Y ENDO AND CORRESPONDING Vendo **/
            cons[inx(i,y)] = inv_MU(deVF[inx(i,y)]);     // define the consumption level tomorrow implied by deVF and K[i]. POLICIES ARE "INVARIANT" WHEN VALUE IS CONVERGED, so current policy == next policy
            Yendo[inx(i,y)] = cons[inx(i,y)] + K[i];
            VFendo[inx(i,y)] = U(cons[inx(i,y)]) + eVF[inx(i,y)];
    
        } // end igridj
    } // end ygrid


    
    /** 3. INTERPOLATE THE VALUE FUNCTION AND COMPUTE CRITERION **/
    critere=0.0;
    
    for(y = 0; y < maxygrid; y++){

        for(i = 0; i < maxigrid; i++){

            ii = 0;
            Ycase = 0;
            if((ii == 0) && (YY[inx(i,y)] < Yendo[inx(ii,y)])){Ycase = 1;} // solution is before Yendo, so extrapolate backward
            
            while((YY[inx(i,y)] > Yendo[inx(ii,y)]) && (ii < maxigrid)){
                if(ii == (maxigrid-1)){Ycase = 2;}  // solution is after Yendo, so extrapolate forward
                ii++;
            }
            
            if(Ycase == 1){
                slope = (VFendo[inx(1,y)] - VFendo[inx(0,y)])/(Yendo[inx(1,y)] - Yendo[inx(0,y)]);
                slope2 = (K[1] - K[0])/(Yendo[inx((maxigrid-1),y)] - Yendo[inx((maxigrid-2),y)]);
                VF[inx(i,y)] = VFendo[inx(0,y)] - (Yendo[inx(0,y)] - YY[inx(i,y)])*slope;
                OPTI_K[inx(i,y)] = K[0] - (Yendo[inx(0,y)] - YY[inx(i,y)])*slope2;
            }
            
            if(Ycase == 2){
                slope = (VFendo[inx((maxigrid-1),y)] - VFendo[inx((maxigrid-2),y)])/(Yendo[inx((maxigrid-1),y)] - Yendo[inx((maxigrid-2),y)]);
                slope2 = (K[(maxigrid-1)] - K[(maxigrid-2)])/(Yendo[inx((maxigrid-1),y)] - Yendo[inx((maxigrid-2),y)]);
                VF[inx(i,y)] = VFendo[inx((maxigrid-1),y)] + (YY[inx(i,y)] - Yendo[inx((maxigrid-1),y)])*slope;
                OPTI_K[inx(i,y)] = K[(maxigrid-1)] + (YY[inx(i,y)] - Yendo[inx((maxigrid-1),y)])*slope2;
            }
            
            if(Ycase == 0){
                weight = (YY[inx(i,y)] - Yendo[inx((ii-1),y)])/(Yendo[inx(ii,y)] - Yendo[inx((ii-1),y)]);
                VF[inx(i,y)] = inter1d(weight,VFendo[inx((ii-1),y)],VFendo[inx(ii,y)]);
                OPTI_K[inx(i,y)] = inter1d(weight,K[(i-1)],K[i]);
            }
            
        } // end igrid
    } // end ygrid
 


          
    /** 4. COMPUTE EXPECTED VALUE FUNCTION **/
    for(i = 0; i < maxigrid; i++){
        for(y = 0; y < maxygrid; y++){
            
            eVFnew[inx(i,y)] = 0.0;
            
            for(ynext = 0; ynext < maxygrid; ynext++){
                eVFnew[inx(i,y)] += ytrans[y][ynext]*VF[inx(i,ynext)];
                
                //if(i == 2){printf("%f %f %f %f",eVFnew[inx(i,y)], ytrans[y][ynext], ytrans[y][ynext]*VF[inx(i,ynext)], VF[inx(i,ynext)]); getchar();}
            }
            
             eVFnew[inx(i,y)] = beta*eVFnew[inx(i,y)];
            
            
            /** 4. COMPUTE THE CRITERION **/
            critere = max(critere,fabs(eVF[inx(i,y)] - eVFnew[inx(i,y)]));
            

        }
    }


  // printf("Cnvg Criterion : %30.25f\n", critere);
   
    

 
}//end of while loop



FILE *vfile;

// VALUE FUNCTION //
vfile=fopen(valuefile, "w"); setbuf (vfile, NULL );
for(i=0;i<maxigrid;i++)
{
    for(y=0;y<maxygrid;y++)
    {
        fprintf(vfile,"%5d\t%5d\t%20.15f\t%20.15f\t%20.15f\t%20.15f\t%20.15f\n",i,y,YY[inx(i,y)],eVF[inx(i,y)],eVFnew[inx(i,y)],deVF[inx(i,y)],OPTI_K[inx(i,y)]);
    }
    fprintf(vfile,"\n");
}
fclose(vfile);

free(deVF);
free(eVF);
free(VFendo);
free(Yendo);
	
}


