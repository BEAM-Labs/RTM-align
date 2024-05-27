/* These functions implement d0 normalization. The d0 for final TM-score
 * output is implemented by parameter_set4final. For both RNA alignment
 * and protein alignment, using d0 set by parameter_set4search yields
 * slightly better results during initial alignment-superposition iteration.
 */
#include <math.h>
#include "basic_fun.h"

double normalize_tm(double tm, int len){
    double mu = 0.0;
    double std = 0.0;
    mu = 23.3130/(len+60.6832)+0.10;
    std = 2.7698 / (len+20.3846)+0.0362;
    mu = mu + 2* std;
    tm = (tm - mu) / std;
    tm = exp(tm) / (1 + exp(tm));
    return tm;
}



void parameter_set4search(const int xlen, const int ylen,
    double &D0_MIN, double &Lnorm,
    double &score_d8, double &d0, double &d0_search, double &dcu0)
{
    //parameter initilization for searching: D0_MIN, Lnorm, d0, d0_search, score_d8
    double c = 0.0;
    D0_MIN=0.5; 
    dcu0=4.25;                       //update 3.85-->4.25
 
    Lnorm=getmin(xlen, ylen);        //normaliz TMscore by this in searching
    c = 0.3222  + 0.7681 * exp(-1 * Lnorm / 25.7198)  + 0.0029 * Lnorm* exp(-1 * Lnorm / 151.8844);
    d0 = pow(Lnorm, 0.3410) * sqrt(1 - c);
    D0_MIN=d0+0.8;              //this should be moved to above
    d0=D0_MIN;                  //update: best for search    

    d0_search=d0;
    if (d0_search>8)   d0_search=8;
    if (d0_search<4.5) d0_search=4.5;

    score_d8=1.5*pow(Lnorm*1.0, 0.3)+3.5; //remove pairs with dis>d8 during search & final
}


void parameter_set4final_C3prime(const double len, double &D0_MIN,
    double &Lnorm, double &d0, double &d0_search)
{
    D0_MIN=0.3; 
    double c = 0.0;
 
    Lnorm=len;            //normaliz TMscore by this in searching
    c = 0.3222  + 0.7681 * exp(-1 * Lnorm / 25.7198)  + 0.0029 * Lnorm* exp(-1 * Lnorm / 151.8844);
    d0 = pow(Lnorm, 0.3410) * sqrt(1 - c);

    d0_search=d0;
    if (d0_search>8)   d0_search=8;
    if (d0_search<4.5) d0_search=4.5;
}

void parameter_set4final(const double len, double &D0_MIN, double &Lnorm,
    double &d0, double &d0_search, const int mol_type)
{
    if (mol_type>0) // RNA
    {
        parameter_set4final_C3prime(len, D0_MIN, Lnorm,
            d0, d0_search);
        return;
    }
    D0_MIN=0.5; 
 
    Lnorm=len;            //normaliz TMscore by this in searching
    if (Lnorm<=21) d0=0.5;          
    else d0=(1.24*pow((Lnorm*1.0-15), 1.0/3)-1.8);
    if (d0<D0_MIN) d0=D0_MIN;   
    d0_search=d0;
    if (d0_search>8)   d0_search=8;
    if (d0_search<4.5) d0_search=4.5;
}

void parameter_set4scale(const int len, const double d_s, double &Lnorm,
    double &d0, double &d0_search)
{
    d0=d_s;          
    Lnorm=len;            //normaliz TMscore by this in searching
    d0_search=d0;
    if (d0_search>8)   d0_search=8;
    if (d0_search<4.5) d0_search=4.5;  
}
