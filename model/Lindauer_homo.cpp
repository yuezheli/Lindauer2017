// this model us an implementation in mrgsolve published by Lindauer et al., 2017
// url:
// https://ascpt.onlinelibrary.wiley.com/doi/full/10.1002/psp4.12130

[PARAM]
// mouse PK parameters
  TVV1     = 2877      // Volume of the central blood compartment                 V1    ; mL
  TVV2     = 2854     // Volume of the peripheral compartment                    V2    ; mL
  TVQ      = 384      // Distribution clearance                                  Q     ; mL/day
  TVVMAX   = 144     // Maximal rate of MM-elimination                          VMAX  ; ug/day
  TVKM     = 0.078     // KM value of MM-elimination                              KM    ; ug/mL
  TVCL     = 167     // Linear clearance                                        CL    ; mL/day
  TVPLQ    = 12.7      // Tumor plasma flow                                       PLQ   ; L/h/L
  TVL      = 0.002     // Fractional Lymphflow                                    L     ;
  TVCLup   = 0.0366    // Rate of pinocytosis per unit endosomal space            CLup  ; L/h/L
  TVKdeg   = 42.9      // Degradation of free antibody from the endosomal space   Kdeg  ; 1/h
  TVv_ref  = 0.842     // vascular reflection coefficient                         v_ref ; 
  TVv_ref_is = 0.2     // lymph reflection coefficient                            v_ref_is ; 

// FcRn binding parameters
  TVFcRni      = 49.8     // Initial FcRn concentration                         FcRni ; uM
  TVFR         = 0.715    // Fraction of FcRn recycled to vascular space        FR
  TVKon_FcRn   = 792*1E6 // Association rate constant of FcRn binding          Kon_FcRn ; 1/M/h
  TVKoff_FcRn  = 23.9     //Dissociation rate constant of FcRn binding          Koff_FcRn ; 1/h
 
// Drug - Target Interaction
  N_Tcell  = 1000           // Number of Tcells per uL blood                                        N_Tcell ; 
  V_blood  = 5100e3           // Blood volume, updated for human value                             ; uL
  N_PD1_TC = 10000          // Number of PD1 receptors per Tcell                                    N_PD1_TC ;
  Tmulti  = 4.32            // Arbitrary multiplier of target concentration in tumor versus blood 
  EMAXTP  = 94.7            // (fold) Maximal increase in target production
  EC50TP  = 1.46            // Concentration of Target-Receptor complex that induces 50% of the maximal induction of target production; nmol/L
  Kon_PD1_iv  = 2880/1E3   // Kon_PD1_iv;  1/nM/h
  Koff_PD1    = 0.144     // Koff_PD1 ; 1/h
  K_IVIV      = 1           // In vitro/in vivo proportionality factor for Kon
  KdegPD1 = 0.00246          // 1/h Degradation rate from Ferl et al referring to tumor              KdegPD1 ; 1/h

// tumor growth parameters
 TVL0   = 0.0017     // exponential Growth                                   L0 ; 1/day
 TVW0   = 41.5e3       // Tumor volume at time of inoculation                  W0 ; mm3
 TVSLtg = 2.886e-7  // Slope of drug effect (ROtumor) on tumor kill rate    Sltg ; 1/%day
 Gamma   = 2.28     // Gamma of drug effect (ROtumor) on tumor kill rate

[SET]
delta = 1 

[CMT] @annotated
A1			: central compartment (blood compartment) monoclonal antibody concentration
A2			: peripheral compartment monoclonal antibody concentration
A3			: Vascular space antibody concentration
A4			: Endosomal space mAb unbound to FcRn
A5			: Endosomal space mAb bound to FcRn
A6 			: Endosomal FcRn
A7			: Interstitial antibody concentration
A8			: Drug receptor in the tumor
A9			: Drug receptor in blood
A10			: PD1 concentration in the tumor interstitium
A11			: tumor volume

[FIX]
 // Constant values
 AV = 6.0221415E23  // Avrogadro's Number
 MW = 149000    // g/mol ; Molecular weight of antibody
 pi = 3.1416 

[OMEGA] @annotated
EW0 : 0.375  : ETA on initial tumor size
ECL : 0.173  : ETA on clearance
EL0 : 0.2    : ETA on growth rate (exponential)

[SIGMA]
1

[MAIN]
double CL      = TVCL / (1000*24) * exp(ECL) ; // mL/day -> L/h
double L0      = TVL0 / 24 *exp(EL0); // 1/day   -> 1/h; adjust for human weight
double W0      = TVW0 * exp(EW0) ;


// PK Model based on mousePK
 double  V1		= TVV1 / 1000 		; // mL->L
 double  V2		= TVV2 / 1000	 	; // mL->L 
 double  Q		= TVQ  / (1000*24) 	; // mL/day -> L/h
 double  VMAX   = TVVMAX / (MW*24) * 1E3; // ug/day -> nmol/h
 double  KM		= TVKM /   (MW)    * 1E6; // mg/L -> nmol/L
 double  PLQraw	= TVPLQ;
 double  Lraw		= PLQraw*TVL ;		
 double  CLupraw	= TVCLup	;
 double  Kdeg		= TVKdeg	;
 double  v_ref		= TVv_ref ;
 double  v_ref_is 	= TVv_ref_is ;	
 double  K12		= Q/V1 ;
 double  K21		= Q/V2 ;
 double  K			= CL/V1 ;
 double  VSS		= V1+V2 ;

// FcRn binding
 double FcRni     = TVFcRni * 1000 ; // umol/L -> nmol/L 
 double FR        = TVFR;
 double Kon_FcRn  = TVKon_FcRn/1E9 ; // 1/M/h -> 1/nM/h     
 double Koff_FcRn = TVKoff_FcRn; 

// Drug - Target Interaction
 double N_PD1_b  = N_PD1_TC * N_Tcell * V_blood; //  
 double M_PD1_b  = N_PD1_b/AV*1E9      ; // // PD1 amount in nmoles in blood
 double C_PD1_b  = M_PD1_b/V1          ; // // PD1 concentration in blood
 double Kon_PD1  = Kon_PD1_iv/K_IVIV   ; // In vivo Kon 
 double kout     = KdegPD1             ; // 1/h   rate constant of target decline - Assumption: = same as Kdeg


// tumor Growth 
double PSI     = 20        ; // recommended value by Simeoni
double SLtg    = TVSLtg / 24   ; // 1/day%  -> 1/h%

double V_is_in    = 0.55 * W0/1000000     ; // Initial intertitial volume
double M_PD1_ti   = C_PD1_b*Tmulti*V_is_in; // PD1 initial amount in tumor interstiti
double kin        = M_PD1_ti * kout ; // Rate constant of target production

// set initial tumor size with variation
A11_0 = W0;

// set initial PD-1 concentration in tumor interstitial within the script
A10_0 = M_PD1_ti;

// set initial FcRn concentration
A6_0 = 49800;

[ODE]

 double V_tot   = A11/1000000;       // uL/ mm^3-> L Total tumor volume
 double V_es    = V_tot * 0.005;         // L endosomal volume
 double V_is    = V_tot * 0.55 ;        // L interstitial volume
 double V_vs    = V_tot * 0.07 ;        // L vascular space volume
 double PLQ = PLQraw  * V_vs ;
 double L   = Lraw    * V_vs ;
 double CLup    = CLupraw * V_es ;

 double C1	= A1/V1;
 double C2	= A2/V2;

 double Cvs	= A3/V_vs;
 double Ce_ub	= A4;
 double Ce_b	= A5;	
 double FcRn	= A6;
 double Cis	= A7/V_is;
 double PD1_t = A8;
 double PD1_b = A9;
 double M_PD1_t=A10;

 double C_PD1_t    = M_PD1_t/V_is	  ; // PD1 concentration in the tumor interstitium

 double ROtumor   = 100*PD1_t/C_PD1_t ; // Receptor occupancy in tumor intertitium
 double ROblood   = 100*PD1_b/C_PD1_b ;

 double DE = SLtg * pow(ROtumor, Gamma) ; // drug effect

// Central compartment C1
dxdt_A1 = -K*(C1*V1) - C1*VMAX/(KM+C1) - PLQ*C1 + PLQ * Cvs - K12*(C1*V1) + K21*(C2*V2)  - Kon_PD1*(C1*V1)*(C_PD1_b-PD1_b) + Koff_PD1*PD1_b*V1 ; 

// Peripheral compartment C2
dxdt_A2 =  K12*(C1*V1) - K21*(C2*V2);  

//  Vascular space tumor
dxdt_A3 =  (PLQ*C1) - ((PLQ - L) * Cvs) - ((1-v_ref) * L * Cvs) - (CLup *Cvs) + (CLup*FR*Ce_b); 

// Endosomal space mAb unbound to FcRn
dxdt_A4 = (CLup/V_es * (Cvs + Cis)) - Kon_FcRn*Ce_ub*FcRn + Koff_FcRn*Ce_b - Kdeg*Ce_ub;  

// Endosomal space mAb bound to FcRn
dxdt_A5 =   Kon_FcRn*Ce_ub*FcRn - Koff_FcRn*Ce_b  - CLup/V_es*Ce_b; 

// Endosomal FcRn
dxdt_A6 =  - Kon_FcRn*Ce_ub*FcRn + Koff_FcRn*Ce_b + CLup/V_es*Ce_b;

// Interstitial
dxdt_A7 =  ((1-v_ref)*L*Cvs) - ((1-v_ref_is) * L * Cis) - (CLup*Cis) + (CLup*(1-FR)*Ce_b)  - Kon_PD1*(Cis*V_is)*(C_PD1_t-PD1_t) + Koff_PD1*PD1_t*V_is;

// Drug Receptor binding
dxdt_A8 =  Kon_PD1*(Cis)*(C_PD1_t-PD1_t)  - Koff_PD1*PD1_t - KdegPD1*PD1_t;
dxdt_A9 =  Kon_PD1*(C1)*(C_PD1_b-PD1_b) - Koff_PD1*PD1_b - KdegPD1*PD1_b;
dxdt_A10 = kin * (1+EMAXTP*PD1_t/(EC50TP+PD1_t)) - kout*M_PD1_t;

// Tumor Volume
dxdt_A11 = L0 * A11 - DE*A11; // exponential growth model

[TABLE]

capture CP              = (A1/V1)*MW/1000000 ;
// the last 2 parts are to remove proportional residual error and additive residual error
capture TV              = A11      ;
capture xC_PD1_t 		    = A10 / (A11*0.55/1000000)    ; // PD1 concentration in the tumor interstitium
capture ROb             = 100 * A9 / C_PD1_b  	      ; // receptor occupancy (RO) in blood
capture ROt             = 100 * A8 / xC_PD1_t         ; // RO in tumor
capture radiuschange    = pow( (A11/ W0), 0.333) ;