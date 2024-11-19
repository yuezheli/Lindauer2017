# define power function similar to C

function pow(x,y)
return Real(x)^y
end

# define the model
function ode_Lindauer!(du, u, p,t)
A1, A2, Avs, Ce_ub, Ce_b, FcRn, Ais, PD1_t, PD1_b, M_PD1_t, TV = u

# define some fixed constants
AV = 6.0221415E23  # Avrogadro's Number
MW = 149000    # g/mol ; Molecular weight of antibody
pi = 3.1416 

# mouse PK parameters, keep the same unit as published in the Lindauer paper
TVV1     = 1.26      # Volume of the central blood compartment                 V1    ; mL
TVV2     = 0.819     # Volume of the peripheral compartment                    V2    ; mL
TVQ      = 4.82      # Distribution clearance                                  Q     ; mL/day
TVVMAX   = 0.518     # Maximal rate of MM-elimination                          VMAX  ; ug/day
TVKM     = 0.366     # KM value of MM-elimination                              KM    ; ug/mL
TVPLQ    = 12.7      # Tumor plasma flow                                       PLQ   ; L/h/L
TVL      = 0.002     # Fractional Lymphflow                                    L     ;
TVCLup   = 0.0366    # Rate of pinocytosis per unit endosomal space            CLup  ; L/h/L
TVKdeg   = 42.9      # Degradation of free antibody from the endosomal space   Kdeg  ; 1/h
TVv_ref  = 0.842     # vascular reflection coefficient                         v_ref ; 
TVv_ref_is = 0.2     # lymph reflection coefficient                            v_ref_is ; 

# FcRn binding parameters, keep the same unit as published in the Lindauer paper
TVFcRni      = 49.8     # Initial FcRn concentration                           FcRni ; uM
TVFR         = 0.715    # Fraction of FcRn recycled to vascular space          FR
TVKon_FcRn   = 80.6*1E6 # Association rate constant of FcRn binding            Kon_FcRn ; 1/M/h
TVKoff_FcRn  = 6.55     #Dissociation rate constant of FcRn binding            Koff_FcRn ; 1/h

# Drug - Target Interaction, keep the same unit as published in the Lindauer paper
N_Tcell  = 1000           # Number of Tcells per uL blood                                        N_Tcell ; 
V_blood  = 1400           # Blood volume                                                         V_blood ; uL
N_PD1_TC = 10000          # Number of PD1 receptors per Tcell                                    N_PD1_TC ;
Tmulti  = 4.32            # Arbitrary multiplier of target concentration in tumor versus blood 
EMAXTP  = 94.7            # (fold) Maximal increase in target production
EC50TP  = 1.46            # Concentration of Target-Receptor complex that induces 50% of the maximal induction of target production; nmol/L
Kon_PD1_iv  = 340.2/1E3   # Kon_PD1_iv;  1/nM/h
Koff_PD1    = 0.10584     # Koff_PD1 ; 1/h
K_IVIV      = 1           # In vitro/in vivo proportionality factor for Kon
KdegPD1 = 0.0194          # 1/h Degradation rate from Ferl et al referring to tumor              KdegPD1 ; 1/h

# tumor growth parameters, keep the same unit as published in the Lindauer paper
TVL0   = 0.113     # exponential Growth                                   L0 ; 1/day
TVSLtg = 1.98e-05  # Slope of drug effect (ROtumor) on tumor kill rate    Sltg ; 1/%day
Gamma   = 2.28     # Gamma of drug effect (ROtumor) on tumor kill rate

# fix variable parameters 
TVCL     = 0.334     # Linear clearance                              CL    ; mL/day
TVL1   = 187         # Linear Growth                                 L1 ; mm3/day
W0   = 170           # Tumor volume at time of inoculation           W0 ; mm3
CL = TVCL / (1000*24)
L0 = TVL0/24

# PK Model based on mousePK
V1   = TVV1 / 1000     ; # mL->L
V2   = TVV2 / 1000   ; # mL->L 
Q    = TVQ  / (1000*24)  ; # mL/day -> L/h
VMAX   = TVVMAX / (MW*24) * 1E3; # ug/day -> nmol/h
KM   = TVKM /   (MW)    * 1E6; # mg/L -> nmol/L
PLQraw = TVPLQ;
Lraw   = PLQraw*TVL ;    
CLupraw  = TVCLup  ;
Kdeg   = TVKdeg  ;
v_ref    = TVv_ref ;
v_ref_is   = TVv_ref_is ;  
K12    = Q/V1 ;
K21    = Q/V2 ;
K      = CL/V1 ;
VSS    = V1+V2 ;

# FcRn binding
FcRni     = TVFcRni * 1000 ; # umol/L -> nmol/L 
FR        = TVFR;
Kon_FcRn  = TVKon_FcRn/1E9 ; # 1/M/h -> 1/nM/h     
Koff_FcRn = TVKoff_FcRn; 

# Drug - Target Interaction
N_PD1_b  = N_PD1_TC * N_Tcell * V_blood; #  
M_PD1_b  = N_PD1_b/AV*1E9      ; # # PD1 amount in nmoles in blood
C_PD1_b  = M_PD1_b/V1          ; # # PD1 concentration in blood
Kon_PD1  = Kon_PD1_iv/K_IVIV   ; # In vivo Kon 
kout     = KdegPD1             ; # 1/h   rate constant of target decline - Assumption: = same as Kdeg


# tumor Growth 
L1      = TVL1 / 24 ; # mm3/day -> mm3/h
PSI     = 20        ; # recommended value by Simeoni
SLtg    = TVSLtg / 24   ; # 1/day%  -> 1/h%

V_is_in    = 0.55 * W0/1000000     ; # Initial intertitial volume
M_PD1_ti   = C_PD1_b*Tmulti*V_is_in; # PD1 initial amount in tumor interstiti
kin        = M_PD1_ti * kout ; # Rate constant of target production

V_tot   = TV/1000000;       # uL/ mm^3-> L Total tumor volume
V_es    = V_tot * 0.005;         # L endosomal volume
V_is    = V_tot * 0.55 ;        # L interstitial volume
V_vs    = V_tot * 0.07 ;        # L vascular space volume
PLQ = PLQraw  * V_vs ;
L   = Lraw    * V_vs ;
CLup    = CLupraw * V_es ;

C1  = A1/V1;
C2  = A2/V2;

Cvs = Avs/V_vs;
Cis = Ais/V_is;

C_PD1_t    = M_PD1_t/V_is   ; # PD1 concentration in the tumor interstitium

ROtumor   = 100*PD1_t/C_PD1_t ; # Receptor occupancy in tumor intertitium; percentage
ROblood   = 100*PD1_b/C_PD1_b ;

DE = SLtg * pow(ROtumor, Gamma) ; # drug effect

du[1] = -K*(C1*V1) - C1*VMAX/(KM+C1) - PLQ*C1 + PLQ * Cvs - K12*(C1*V1) + K21*(C2*V2)  - Kon_PD1*(C1*V1)*(C_PD1_b-PD1_b) + Koff_PD1*PD1_b*V1

du[2] = K12*(C1*V1) - K21*(C2*V2)

du[3] = (PLQ*C1) - ((PLQ - L) * Cvs) - ((1-v_ref) * L * Cvs) - (CLup *Cvs) + (CLup*FR*Ce_b)

du[4] = (CLup/V_es * (Cvs + Cis)) - Kon_FcRn*Ce_ub*FcRn + Koff_FcRn*Ce_b - Kdeg*Ce_ub

du[5] = Kon_FcRn*Ce_ub*FcRn - Koff_FcRn*Ce_b  - CLup/V_es*Ce_b

du[6] =  -Kon_FcRn*Ce_ub*FcRn + Koff_FcRn*Ce_b + CLup/V_es*Ce_b

du[7] = ((1-v_ref)*L*Cvs) - ((1-v_ref_is) * L * Cis) - (CLup*Cis) + (CLup*(1-FR)*Ce_b)  - Kon_PD1*(Cis*V_is)*(C_PD1_t-PD1_t) + Koff_PD1*PD1_t*V_is

du[8] = Kon_PD1*(Cis)*(C_PD1_t-PD1_t)  - Koff_PD1*PD1_t - KdegPD1*PD1_t

du[9] = Kon_PD1*(C1)*(C_PD1_b-PD1_b) - Koff_PD1*PD1_b - KdegPD1*PD1_b

du[10] = kin * (1+EMAXTP*PD1_t/(EC50TP+PD1_t)) - kout*M_PD1_t

du[11] = L0 * TV/ pow( 1+ pow( (L0/L1*TV) , PSI), 1/PSI) - DE*TV

end
