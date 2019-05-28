
/* 
 * File:   main.c
 * Author: Diogo Costa, diogo.costa@usask.ca
 *
 * Created on January 9, 2019, 12:16 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*
 * 
 */
 
const int loopno = 9; // no of neighbouring pixels to identify

char errmg1[35] = "Calculation completed successfully";
char errmg2[171] = "Problem with the photo or pixels selected. Take a new photo and make sure (1) you have good light conditions and (2) select the reference and measurement pixels correctly";
char errmg3[39] = "Concentration above the test kit limit";

double median(int n, double x[], double x_unc[]) {
    double temp, temp_u;
    int i, j, m;
    
    // remove NaN values
    m = 0; // no NaN no of fields
    for(i=0; i<n; i++)
    {
      if(isnan(x[i]))
      {
          x[i] = 9999;
          x_unc[i] = 9999;
      }else{
          m = m + 1;      
        }
    }
    
    // the following two loops sort the array x in ascending order
    for(i=0; i<n-1; i++) {
        
        for(j=i+1; j<n; j++) {          
            if(x[j] < x[i]) {
                temp = x[i];
                temp_u = x_unc[i];
                // swap elements          
                x[i] = x[j];
                x[j] = temp;
                x_unc[i] = x_unc[j];
                x_unc[j] = temp_u;
            }
        }
    }
    
    if(m%2==0) {
        // if there is an even number of elements, return mean of the two elements in the middle
        return((x[m/2] + x[m/2 - 1]) / 2.0);
    } else {
        // else return the element in the middle
        return x[m/2];
    }
}

// RGB to XYZ color space conversion algorithm
double conv_rgb2xyz(double *inR, double *inG, double *inB, double *outX, double *outY, double *outZ) 
{
    
    double var_R = (*inR / 255.0f); //R from 0 to 255
    double var_G = (*inG / 255.0f); //G from 0 to 255
    double var_B = (*inB / 255.0f); //B from 0 to 255

    if (var_R > 0.04045f) {
        var_R = pow(( (var_R + 0.055f) / 1.055f), 2.4f);
    } else {
        var_R = var_R / 12.92f;
    }

    if (var_G > 0.04045f) {
        var_G = pow(( (var_G + 0.055f) / 1.055f), 2.4f);
    } else {
        var_G = var_G / 12.92f;
    }

    if (var_B > 0.04045f) {
        var_B = pow(( (var_B + 0.055f) / 1.055f), 2.4f);
    } else {
        var_B = var_B / 12.92f;
    }

    var_R = var_R * 100;
    var_G = var_G * 100;
    var_B = var_B * 100;

    //Observer. = 2°, Illuminant = D65
    *outX = var_R * 0.4124f + var_G * 0.3576f + var_B * 0.1805f;
    *outY = var_R * 0.2126f + var_G * 0.7152f + var_B * 0.0722f;
    *outZ = var_R * 0.0193f + var_G * 0.1192f + var_B * 0.9505f;
    
    return 0;
}

// Delta E calc: auxiliary function 1
double func(double x)
{
    double y;
    double a = 0.3333333333333333333;
    if (x > 0.008856f) {
        y = pow(x,a);
    } else {
        y = 7.787f * x + 16/116;
    }
    
    return y;
}

// Delta E calc: 2000 method auxiliary function 2
double f_atan(double x1, double x2)
{
    double y;
    if (x1 == x2){
        y = 0;
    } else if (x1 >= 0  && x2 < 0) {
        y = atan(x2 / x1) + 2 * M_PI;
    } else if (x1 < 0) {
        y = atan(x2 / x1) + M_PI;
    } else {
        y = atan(x2 /x1);
    }   
    return y;
}

// Delta E calc: 2000 method auxiliary function 3
double f_deltahp(double x1, double x2, double x3, double x4)
{ 
    double y;
    
    if ((x3 + x4) == 0) {
        y = 0;
    } else if (fabs(x2 - x1) <= M_PI) {
        y = x2 - x1;
    } else if ((x2 -x1) > M_PI) {
        y = x2 - x1 - 2 * M_PI;
    } else {
        y = x2 - x1 + 2 * M_PI;
    }

    return y;
}

// Delta E calc: 2000 method auxiliary function 4
double hueDiff(double x1, double x2, double x3, double x4)
{
    double y;
    
    if ((x3 * x4) == 0) {
        y = x1 + x2;
    } else if (fabs(x1 - x2) <= M_PI) {
        y = (x1+x2)/2;
    } else if (fabs(x1-x2) > M_PI && (x1+x2) < (2 * M_PI)) {
        y = (x1 + x2 + 2 * M_PI) /2 ;
    } else {
        y = (x1+x2-2*M_PI)/2;
    }
    
    return y;
}

//Delta E: 2000 method
double delta_2000(double SamplesC1, double SamplesC2, double SamplesL1, 
                  double SamplesL2, double Samplesa1, double Samplesa2,
                  double Samplesb1, double Samplesb2, double deltaCab,
                  double deltaHab) 
{       
    double C_bar, G, Lsp, Lbp, asp, abp, bsp, bbp, Csp, Cbp, hsp, hbp;
    double deltaLp, deltaCp, deltahp, deltaHp, Lp_bar, Cp_bar, hp_bar;
    double SL, SC, T, deltaTheta, RC, RT, SH, deltaE_i;

    C_bar = (SamplesC1 + SamplesC2)/2;
    G = 0.5*(1-sqrt(pow(C_bar,7)/(pow(C_bar,7)+pow(25,7))));

    Lsp = SamplesL1; // s for standard
    Lbp = SamplesL2; // b for sample
    asp = (1 + G) * Samplesa1;
    abp = (1 + G) * Samplesa2;
    bsp = Samplesb1;
    bbp = Samplesb2;
    Csp = sqrt(pow(asp,2) + pow(bsp,2));
    Cbp = sqrt(pow(abp,2) + pow(bbp,2));
    hsp = f_atan(asp,bsp);
    hbp = f_atan(abp,bbp);

    deltaLp = Lbp - Lsp;
    deltaCp = Cbp - Csp;
    deltahp = f_deltahp(hsp, hbp, Csp, Cbp);
    deltaHp = 2 * sqrt(Cbp * Csp) * sin(deltahp / 2);

    Lp_bar = (Lsp + Lbp) / 2;
    Cp_bar = (Csp + Cbp) / 2;
    hp_bar = hueDiff(hsp, hbp, Csp, Cbp);

    SL = 1 + (0.015 * pow((Lp_bar - 50),2)) / sqrt(20 + pow((Lp_bar - 50),2));
    SC = 1 + 0.045 * Cp_bar;
    T = 1 - 0.17 * cos(hp_bar - M_PI/6) + 0.24 * cos(2 * hp_bar) + 
            0.32 * cos(3 * hp_bar + M_PI/30) - 
            0.20 * cos(4 * hp_bar - M_PI/(180/63));
    deltaTheta = (30 * exp(-pow(((hp_bar - M_PI/(180/275))/
                  (25 * M_PI/180)),2))) * M_PI/180;
    RC = 2 * sqrt((pow(Cp_bar,7)) / (pow(Cp_bar,7) + pow(25,7)));
    RT = -sin(2 * deltaTheta) * RC;
    SH = 1 + 0.015 * Cp_bar * T;

    deltaE_i = sqrt(pow((deltaLp / SL),2) + pow((deltaCp / SC),2) + 
                pow((deltaHp / SH),2) +
                RT * (deltaCp / SC) * (deltaHp / SH));

    return deltaE_i;
}

// Delta E: 1994 method
double delta_1994(double SamplesC1_i, double SamplesC2_i, double SamplesL1_i,
                  double SamplesL2_i, double deltaCab_i, double deltaHab_i)
{
    double deltaE_i, a, aa, b, c;
    double Cab, SL, SC, SH;
           
    Cab = sqrt(SamplesC1_i * SamplesC2_i);
    SL = 1;
    SC = 1 + 0.045*Cab;
    SH = 1 + 0.015*Cab;
    a = SamplesL1_i-SamplesL2_i;  
    aa = a/SL;
    b = deltaCab_i/SC;
    c = deltaHab_i/SH;
            
    deltaE_i = sqrt(pow(aa,2) + pow(b,2) + pow(c,2)); // ATTENTION: removed this term because its gives a complex num)
    
    return deltaE_i;
}

// Delta E: Lab method
double delta_lab(double SamplesL1_i, double SamplesL2_i,
                 double deltaCab_i, double deltaHab_i)
{
    double deltaE_i, a, b, c;
    
    a = 1 * (SamplesL1_i-SamplesL2_i);
    b = 0.5 * deltaCab_i;
    c = 0.75 * deltaHab_i;
    deltaE_i = sqrt(pow(a,2) + pow(b,2) + pow(c,2));  
    
    return deltaE_i;
}

// Delta E calculation
double deltae_calc(double SamplesRGB1[], double SamplesRGB2[], int deltaem)
{
    double a, b, c, Yn = 100;
    float Xn = 95.05, Zn = 108.90;
    double deltaE_i;
    double SamplesXYZ1[3], SamplesXYZ2[3];
    double SamplesL1, SamplesL2, Samplesa1, Samplesa2;
    double Samplesb1, Samplesb2, Samplesc1, Samplesc2;
    double deltaEab, deltaCab,deltaHab;
    
    /*  --------  */   
    /*  Calculations for reference used for DeltaE calc: zero concentration */
    /*  --------  */   
    conv_rgb2xyz(&SamplesRGB2[0],&SamplesRGB2[1],&SamplesRGB2[2],&SamplesXYZ2[0],&SamplesXYZ2[1],&SamplesXYZ2[2]);
    //SamplesXYZ2[0] = 100 * SamplesXYZ2[0];
    //SamplesXYZ2[1] = 100 * SamplesXYZ2[1];
    //SamplesXYZ2[2] = 100 * SamplesXYZ2[2];
    SamplesL2 = 116*(func(SamplesXYZ2[1]/Yn)) - 16;
    Samplesa2 = 500*(func(SamplesXYZ2[0]/Xn) - func(SamplesXYZ2[1]/Yn));
    Samplesb2 = 200*(func(SamplesXYZ2[1]/Yn) - func(SamplesXYZ2[2]/Zn));
    Samplesc2 = sqrt(pow(Samplesa2,2) + pow(Samplesb2,2));
     
    /*  --------  */    
    /*  Calculation of Delta E */
    /*  --------  */
    
    // Convert RGB to XYZ
    conv_rgb2xyz(&SamplesRGB1[0],&SamplesRGB1[1],&SamplesRGB1[2],&SamplesXYZ1[0],&SamplesXYZ1[1],&SamplesXYZ1[2]);
    // Intermediate Calculations
    //SamplesXYZ1[0] = 100 * SamplesXYZ1[0];
    //SamplesXYZ1[1] = 100 * SamplesXYZ1[1];
    //SamplesXYZ1[2] = 100 * SamplesXYZ1[2];
    SamplesL1 = 116*(func(SamplesXYZ1[1]/Yn)) - 16;
    Samplesa1 = 500*(func(SamplesXYZ1[0]/Xn) - func(SamplesXYZ1[1]/Yn));
    Samplesb1 = 200*(func(SamplesXYZ1[1]/Yn) - func(SamplesXYZ1[2]/Zn));
    Samplesc1 = sqrt(pow(Samplesa1,2) + pow(Samplesb1,2));

    /* Calculation of deltaEab, deltaCab, deltaHab */ 
    a = SamplesL1 - SamplesL2;
    b = Samplesa1 - Samplesa2;
    c = Samplesb1 - Samplesb2;
    deltaEab = sqrt(pow(a,2) + pow(b,2) + pow(c,2));
    deltaCab = sqrt(pow(b,2) + pow(c,2));
    deltaHab = sqrt(pow(deltaEab,2) - pow(a,2) - pow(deltaCab,2));
    
    // DeltaE_1994
    if (deltaem == 0) { 
        deltaE_i = delta_1994(Samplesc1, Samplesc2, SamplesL1,
                            SamplesL2, deltaCab, deltaHab);
    // DeltaE_2000    
    } else if (deltaem == 1) {
        deltaE_i = delta_2000(Samplesc1, Samplesc2, SamplesL1, SamplesL2, 
                              Samplesa1, Samplesa2, Samplesb1, Samplesb2,
                              deltaCab,  deltaHab);
    // DeltaE_lab
    } else { 
        deltaE_i = delta_lab(SamplesL1, SamplesL2, deltaCab, 
                            deltaHab);     
    } 
        return deltaE_i;
}

// Obtain concentration: interpolation of DeltaEs
double getconc(double colchar_de[], double res_de, int nboxes, double colchar_cref[], int *outerrorno)
{
    int i;
    double con_u = 0.f;
    double diff_de_i, debel = -10000, deabo = 10000;
    double crefbe = 0, crefab = 0;
    
    // Look for upper and lower reference values
    for (i = 0; i < nboxes; i++) 
    {
        diff_de_i = colchar_de[i] - res_de;
        if (diff_de_i > 0.f && diff_de_i < deabo){
            deabo = colchar_de[i];
            crefab = colchar_cref[i];
        }
        else if (diff_de_i < 0.f && diff_de_i > debel){
            debel = colchar_de[i];
            crefbe = colchar_cref[i];
        } else if (diff_de_i == 0.f) {
            con_u = colchar_cref[i];
            break;
        }
    }
    
    // Look for errors
    if (debel == -10000 && deabo == 10000)
    {
      *outerrorno = 1;
    }
    else if (debel != -10000 && deabo == 10000) // this condition was missing
    {
        *outerrorno = 2;
    }
    
    // Interpolation
    if (diff_de_i != 0.f){
        con_u = crefbe + (crefab - crefbe) * (res_de - debel) / (deabo - debel);
    }
    
    return con_u;
    
}

// Process the image
double ProcImag(int chem, int nboxes, int deltaem, double parcal[3][3], double colchar_rgb[][9], double res_rgb[][9], double colchar_cref[], double con_c[], double con_u[], char **errmsg)
{
    
    // Declarations
    int i, j, l;
    int boxes_tot = nboxes * 3;
    double con_c_median; // uncorrected and corrected concentration estimations
    double colchar0_rgb[3], colchar_rgb_i[3],res_rgb_i[3]; // 
    double colchar_de[nboxes]; // deltaE for colorchart
    //int colchar_dep; // pointer for colchar_de
    double res_de; // deltaE for results
    int errorno[9] = {0,0,0,0,0,0,0,0,0}; // error message number
    int err0 = 0, err2 = 0; // err1 = 0

    // Loop over the different pixel combinations 
    for (l = 0; l < loopno; l++)
    {
        // RGB for zero concentration in colorchart
        colchar0_rgb[0] = colchar_rgb[boxes_tot-3][l]; 
        colchar0_rgb[1] = colchar_rgb[boxes_tot-2][l];
        colchar0_rgb[2] = colchar_rgb[boxes_tot-1][l];

        // Determine DeltaE values: colorchart 
        for (i = 0; i < nboxes; i++)
        {
            j = i * 3;
            colchar_rgb_i[0] = colchar_rgb[j][l];
            colchar_rgb_i[1] = colchar_rgb[j+1][l];
            colchar_rgb_i[2] = colchar_rgb[j+2][l];
            colchar_de[i] = deltae_calc(colchar0_rgb, colchar_rgb_i, deltaem); 
            //printf( "%f\n", colchar_cref[i]);
        }
        
        // Determine the l measurement RBS
        res_rgb_i[0] = res_rgb[0][l];
        res_rgb_i[1] = res_rgb[1][l];
        res_rgb_i[2] = res_rgb[2][l];
        
        // Determine DeltaE values: measurement 
        res_de = deltae_calc(colchar0_rgb,res_rgb_i,deltaem);

        // Determine concentration from the DeltaEs
        con_u[l] = getconc(colchar_de,res_de, nboxes, colchar_cref,&errorno[l]);
        con_c[l] = con_u[l] / parcal[chem][deltaem];          
    }
    
    // handling the error message
    for (l = 0; l < loopno; l++){
        if (errorno[l]==0) { // for errmg1
            err0 = err0 + 1;
        //}else if(errorno[l]==1) { // for errmg2 // Please leave this comment so I can follow the code in the future
        //    err1 = err1 + 1;
        }else if(errorno[l]==2) { // for errmg3
            err2 = err2 + 1;
        };
    }
   
    
    if (err0>2){ //If 3 or more pixels were non-failure then use the result.  Otherwise throw error
        con_c_median = median(loopno,con_c,con_u);
        *errmsg = errmg1;
    }else{ // the error message can be different depending on the errorno values
        con_c_median = 9999;
        if (err2>4) { // this is when the algorithm cannot find an upper color reference from the ref card - this indicates that the concentration is higher than the upper limit of the card
            *errmsg = errmg3;
        } else { // this is the condition for err1 (lots of NaN) - take another photo)
            *errmsg = errmg2;
        };          
    };
    
    return con_c_median;
}

// This subroutine is just to generate the necessary input data for the subroutine "ProcImag" and it should not be integrated in the app. 
// You should only integrate subroutine "con_c = ProcImag(nboxes, deltaem, parcal, colchar_rgb, res_rgb_c, colchar_cref)" and all the other subroutines it calls.
int main()
{
    
    // Advanced Settings
    int nboxes = 7; // number of boxes in the colorchart

    int chem = 0; // chemical to analyze: 0) NO3, 1) PO4, 2) customized
   
    int deltaem = 0; // delta E method to use: 0) 1994, 1) 2000, 2) lab
    double parcal[3][3] = {{2.1868,2.7771,2.3146}, // NO3
                           {3.6458, 2.6888, 3.2431},  //PO4
                           {1,1,1}}; // costumized (the user should be able to cahnge this in advanced settings)
    double conref_advset[21] = {50, 20, 10, 5, 2, 1, 0,
                                10, 5, 2, 1, 0.5, 0.25, 0,
                                0, 0, 0, 0, 0, 0, 0}; // this last line can be changed in advaced settings
    char *errmsg=NULL;
    
    // Declarations and initiations
    int i, j, boxes_tot;
    double con_c_median;
    //double * colchar_rgb;
    boxes_tot = nboxes * 3;
    double colchar_rgb[boxes_tot][loopno];
    double res_rgb[3][loopno];
    double colchar_cref[nboxes];
    double con_c[loopno],con_u[loopno];

    // User input
    // NO3 - reference card
    // lab standard: 0.49 mg/l (taken from MATLAB)
    double colchar_userinput_NO3_1[21] = {215, 42, 108,
                                    224, 96, 144,
                                    215, 121, 147,
                                    214, 135, 154,
                                    211, 166, 161,
                                    205, 185, 178,
                                    209, 185, 161};
    double colchar_userinput_NO3_2[21] = {215,42,108,
                                        226,98,146,
                                        216, 122, 148,
                                        213, 134 ,  153,
                                        210, 165,   160,
                                        205, 185 ,  178,
                                        209, 185,   161};
    double colchar_userinput_NO3_3[21] = {216,    43,   107,
                                        224,    96,   144,
                                        217,   123,   150,
                                        218,   139,   158,
                                        208 ,  163,   158,
                                        205,   185,   178,
                                        207,   183,   159};
    double colchar_userinput_NO3_4[21] = {215,    42,   108,
                                        223,    95,   143,
                                        216,   123,   150,
                                        214,   137,   155,
                                        209,   164,   158,
                                        205,   185,   178,
                                        207,   183,   159};
    double colchar_userinput_NO3_5[21] = {216,    43,   109,
                                        228,   100,   148,
                                        215,   122,   149,
                                        212,   135,   153,
                                        213,   164,   160,
                                        205,   185,   178,
                                        207,   183,   159};
   double colchar_userinput_NO3_6[21] = {215,    42,   108,
                                        220,    92,   140,
                                        216,   122,   148,
                                        215,   136,   155,
                                        209,   164,   159,
                                        203,   183,   176,
                                        209,   185,   161};
   double colchar_userinput_NO3_7[21] = {203,    49,   101,
                                        220,    92,   141,
                                        214,   120,   146,
                                        212,   133,   152,
                                        211,   166,   161,
                                        203,   183,   176,
                                        212,   188,   164};
   double colchar_userinput_NO3_8[21] = {215,    42 ,  108,
                                        228,   100,   148,
                                        215,   121,   147,
                                        215,   136,   155,
                                        209,   164,   158,
                                        204,   184,   177,
                                        209,   185,   161};
   double colchar_userinput_NO3_9[21] = {215,    42,   108,
                                        230,   102,   150,
                                        217,   123,   150,
                                        213,   136,   154,
                                        208,   163,   158,
                                        206,   186,   179,
                                        209,   185,   161};

   // NO3 - measurements
    double res_rgb_NO3_1[3] = {214,181,176};
    double res_rgb_NO3_2[3] = {215,   182,   177};
    double res_rgb_NO3_3[3] = {214,   181,   176};
    double res_rgb_NO3_4[3] = {215,   182,   177};
    double res_rgb_NO3_5[3] = {214,   181,   176};
    double res_rgb_NO3_6[3] = {214,   181,   176};
    double res_rgb_NO3_7[3] = {214,  181,   176};
    double res_rgb_NO3_8[3] = {215,   182,   177};
    double res_rgb_NO3_9[3] = {214,   181,   176};
   
// // PO4 ref card
// // lab standard: 0.9 mg/l (taken from MATLAB)
//   double colchar_userinput_PO4_1[21] = {223,   231,   174,
//                                    213,   224,   155,
//                                    198 ,  215,   160,
//                                    138,   180,   140,
//                                     99,   182,   154,
//                                      3,   136,   177,
//                                     15,    91,   153};
//    double colchar_userinput_PO4_2[21] = {230,   229,   175,
//                                    216,   228,   156,
//                                    202,   219,   164,
//                                    145,   197,   149,
//                                    100,   185,   154,
//                                      1 ,  137,   177,
//                                     17,    92,   157};
//    double colchar_userinput_PO4_3[21] = {225,   230,   174,
//                                        216,   227,   159,
//                                        201,   218,   163,
//                                        138,   189,   146,
//                                         97,   184,   152,
//                                          0,   135,   175,
//                                         10,    88,   152};
//    double colchar_userinput_PO4_4[21] = {228,   231,   176,
//                                        215,   225,   162,
//                                        203,   220,   165,
//                                        137,   191,   142,
//                                        108,   189,   172,
//                                          0,   133,   174,
//                                         16 ,   91,   156};
//    double colchar_userinput_PO4_5[21] = {228,   231,   178,
//                                        211,   226,   161,
//                                        199,   216,   161,
//                                        140,   185,   144,
//                                         95,   183,   158,
//                                          0,   132,   173,
//                                         13,    89,   147};
//   double colchar_userinput_PO4_6[21] = { 228,   233,   177,
//                                        215,   227,   155,
//                                        201,   218,   163,
//                                        139,   193,   144,
//                                         99,   184,   155,
//                                          0,   132,   173,
//                                         15,    90,   145};
//   double colchar_userinput_PO4_7[21] = {224,   229 ,  173,
//                                        217,   229,   157,
//                                        199 ,  216,   161,
//                                        143,   195,   146,
//                                         97,   180 ,  152,
//                                          0,   133,   174,
//                                         15,    86,   142};
//   double colchar_userinput_PO4_8[21] = {226,   229,   174,
//                                        216,   227,   158,
//                                        202,   219,   164,
//                                        142,   194,   146,
//                                         98,   188,   163,
//                                          1,   134,   175,
//                                         16,    86,   148};
//   double colchar_userinput_PO4_9[21] = {226,   231,   175,
//                                        212,   227,   162,
//                                        200,   217,   162,
//                                        137,   193,   144,
//                                         97,   181,   147,
//                                          1,   137,   177,
//                                         17,    91,   152};
//
//    // measurements PO4
//    double res_rgb_PO4_1[3] = {210,   221,   179};
//    double res_rgb_PO4_2[3] = {217,   228,   186};
//    double res_rgb_PO4_3[3] = {210,   221,   179};
//    double res_rgb_PO4_4[3] = {210,   222,   176};
//    double res_rgb_PO4_5[3] = {215,   220,  180};
//    double res_rgb_PO4_6[3] = {213,   221,   184};
//    double res_rgb_PO4_7[3] = {212,   220,   179};
//    double res_rgb_PO4_8[3] = { 214,   221,   180};
//    double res_rgb_PO4_9[3] = { 210,   221,   179};

    
     
    // Getting the reference colors (reference card)
    for (i=0;i<boxes_tot;i++){
        colchar_rgb[i][0] = colchar_userinput_NO3_1[i];
        colchar_rgb[i][1] = colchar_userinput_NO3_2[i];
        colchar_rgb[i][2] = colchar_userinput_NO3_3[i];
        colchar_rgb[i][3] = colchar_userinput_NO3_4[i];
        colchar_rgb[i][4] = colchar_userinput_NO3_5[i];
        colchar_rgb[i][5] = colchar_userinput_NO3_6[i];
        colchar_rgb[i][6] = colchar_userinput_NO3_7[i];
        colchar_rgb[i][7] = colchar_userinput_NO3_8[i];
        colchar_rgb[i][8] = colchar_userinput_NO3_9[i];
       
//        colchar_rgb[i][0] = colchar_userinput_PO4_1[i];
//        colchar_rgb[i][1] = colchar_userinput_PO4_2[i];
//        colchar_rgb[i][2] = colchar_userinput_PO4_3[i];
//        colchar_rgb[i][3] = colchar_userinput_PO4_4[i];
//        colchar_rgb[i][4] = colchar_userinput_PO4_5[i];
//        colchar_rgb[i][5] = colchar_userinput_PO4_6[i];
//        colchar_rgb[i][6] = colchar_userinput_PO4_7[i];
//        colchar_rgb[i][7] = colchar_userinput_PO4_8[i];
//        colchar_rgb[i][8] = colchar_userinput_PO4_9[i];
    }
    
    // Getting the reference concentrations
    for (i = 0; i<nboxes;i++){
    j = chem * nboxes + i;
        colchar_cref[i] = conref_advset[j];
        printf( "%f\n", colchar_cref[i]);
    }
    // Getting the measurement colors
    for (i=0;i<3;i++){
        res_rgb[i][0] = res_rgb_NO3_1[i];
        res_rgb[i][1] = res_rgb_NO3_2[i];
        res_rgb[i][2] = res_rgb_NO3_3[i];
        res_rgb[i][3] = res_rgb_NO3_4[i];
        res_rgb[i][4] = res_rgb_NO3_5[i];
        res_rgb[i][5] = res_rgb_NO3_6[i];
        res_rgb[i][6] = res_rgb_NO3_7[i];
        res_rgb[i][7] = res_rgb_NO3_8[i];
        res_rgb[i][8] = res_rgb_NO3_9[i];
       
//        res_rgb[i][0] = res_rgb_PO4_1[i];
//        res_rgb[i][1] = res_rgb_PO4_2[i];
//        res_rgb[i][2] = res_rgb_PO4_3[i];
//        res_rgb[i][3] = res_rgb_PO4_4[i];
//        res_rgb[i][4] = res_rgb_PO4_5[i];
//        res_rgb[i][5] = res_rgb_PO4_6[i];
//        res_rgb[i][6] = res_rgb_PO4_7[i];
//        res_rgb[i][7] = res_rgb_PO4_8[i];
//        res_rgb[i][8] = res_rgb_PO4_9[i];
    }

    // Process image
    // Kiya: the python interface should retrieve the arguments of ProcImag
    // this code should replace part of the variable initiation above
    con_c_median = ProcImag(chem, nboxes, deltaem, parcal, colchar_rgb, res_rgb, colchar_cref, con_c,con_u, &errmsg);
    printf( "con_c_median = %f\n", con_c_median);
    printf("%s\n", errmsg);
    
    return 0;
 
}   
