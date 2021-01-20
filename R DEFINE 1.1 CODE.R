####################################################################################
# DEFINE 1.1
####################################################################################
#STRUCTURE OF THE CODE
#1. DEFINING THE ENDOGENOUS VARIABLES 
#Non-Monte Carlo variables
#Monte Carlo variables
#2. VALUES FOR PARAMETERS & VARIABLES 
#Parameters necessary for the scenarios
#Parameters per scenario
#Non-scenario parameters
#Parameters/endogenous variables for the steady state
#Initial values for endogenous variables
#3. MODEL EQUATIONS
#4. FILLING IN THE MONTE-CARLO VARIABLES
#5. PRINTING THE RESULTS IN R
#6. EXCEL FILES WITH VARIABLES AND PARAMETER VALUES
#7. FIGURES
#8. TABLES
####################################################################################

rm(list=ls(all=T))
T<-83  #number of years
loops<-28 #number of scenarios that are run (minimum 6)
tables_dummy<-1 #should be set equal to 0 when we do not run all the scenarios
figures_dummy<-1 #should be set equal to 0 when we do not run all the scenarios 
scenarios_print<-1 #number of scenarios that are printed in R
monte_sce<-200 #number of Monte Carlo simulations for scenarios 2 and 3
Data<- read.csv("_DEFINE_Carbontaxes.csv") #this file includes carbon taxes and carbon intensities

######################################
#1. DEFINING THE ENDOGENOUS VARIABLES 
######################################

########################
#Non-Monte Carlo variables
########################

#2.1.1 Matter, recycling and waste
MY<- vector(length=T)
M<- vector(length=T)
REC<- vector(length=T)
DEM<- vector(length=T)
SES<- vector(length=T)
W<- vector(length=T)
CEN<- vector(length=T)
O2<- vector(length=T)	
HW_CUM<- vector(length=T)
hazratio<- vector(length=T)
REV_M<- vector(length=T)
CON_M<- vector(length=T)
RES_M<- vector(length=T)
dep_M<- vector(length=T)

#2.1.2 Energy
E<- vector(length=T)
E_NF<- vector(length=T)
E_F<- vector(length=T)
ED<- vector(length=T)
REV_E<- vector(length=T)
CON_E<- vector(length=T)
RES_E<- vector(length=T)
dep_E<- vector(length=T)

#2.1.3 Emissions and climate change 
EMIS_IN<- vector(length=T)
g_EMIS_L<- vector(length=T)
EMIS_L<- vector(length=T)
EMIS<- vector(length=T)
CO2_CUM<- vector(length=T)
T_AT<- vector(length=T)

#2.1.4 Ecological efficiency and technology
omega<- vector(length=T)
g_omega<- vector(length=T)
mu<- vector(length=T)
rho<- vector(length=T)
epsilon<- vector(length=T)
theta<- vector(length=T)
seq<- vector(length=T)

#2.2.1 Output determination and climate damages
Y_star_M<- vector(length=T)
Y_star_E<- vector(length=T)
Y_star_K<- vector(length=T)
Y_star_N<- vector(length=T)
Y_star<- vector(length=T)
Y<- vector(length=T)
Y_N<- vector(length=T)
um<- vector(length=T)
ue<- vector(length=T)
u<- vector(length=T)
re<- vector(length=T)
D_T<- vector(length=T)
D_TP<- vector(length=T)
D_TF<- vector(length=T)

#2.2.2 Firms
TP_G<- vector(length=T)
TP<- vector(length=T)
RP<- vector(length=T)
DP<- vector(length=T)
r<- vector(length=T)
I_PRI_D<- vector(length=T)
I_PRI_D_S1<- vector(length=T) 
I_PRI_D_S2<- vector(length=T) 
I_PRI_D_S3<- vector(length=T)
I_PRI_D_S4<- vector(length=T)
IG_PRI_D_S1<- vector(length=T) 
IG_PRI_D_S2<- vector(length=T) 
IG_PRI_D_S3<- vector(length=T)
IG_PRI_D_S4<- vector(length=T)
beta_S1<- vector(length=T)
beta_S2<- vector(length=T)
beta_S3<- vector(length=T)
beta_S4<- vector(length=T)
beta_0_S1<- vector(length=T)
beta_0_S2<- vector(length=T)
beta_0_S3<- vector(length=T)
beta_0_S4<- vector(length=T)
g_beta_0<- vector(length=T)
tucr<- vector(length=T)
tucn<- vector(length=T)
ucn<- vector(length=T)
g_ucn<- vector(length=T)
ucr<- vector(length=T)
g_ucr<- vector(length=T)
IC_PRI_D_S1<- vector(length=T)
IC_PRI_D_S2<- vector(length=T)
IC_PRI_D_S3<- vector(length=T)
IC_PRI_D_S4<- vector(length=T)
NLG_D_S1<- vector(length=T)
NLG_D_S2<- vector(length=T)
NLG_D_S3<- vector(length=T)
NLG_D_S4<- vector(length=T)
NLC_D_S1<- vector(length=T)
NLC_D_S2<- vector(length=T)
NLC_D_S3<- vector(length=T)
NLC_D_S4<- vector(length=T)
IG_PRI_S1<- vector(length=T)
IG_PRI_S2<- vector(length=T)
IG_PRI_S3<- vector(length=T)
IG_PRI_S4<- vector(length=T)
IC_PRI_S1<- vector(length=T) 
IC_PRI_S2<- vector(length=T)
IC_PRI_S3<- vector(length=T)
IC_PRI_S4<- vector(length=T)
IG_PRI<- vector(length=T)
IC_PRI<- vector(length=T)
I_PRI<- vector(length=T)
kappa<- vector(length=T)
L<- vector(length=T)
KG_PRI_S1<- vector(length=T)
KG_PRI_S2<- vector(length=T)
KG_PRI_S3<- vector(length=T)
KG_PRI_S4<- vector(length=T)
KC_PRI_S1<- vector(length=T)
KC_PRI_S2<- vector(length=T)
KC_PRI_S3<- vector(length=T)
KC_PRI_S4<- vector(length=T)
KG_PRI<- vector(length=T)
KC_PRI<- vector(length=T)
K_PRI<- vector(length=T)
KGE_PRI_S1<- vector(length=T)
KGE_PRI_S2<- vector(length=T)
KGE_PRI_S3<- vector(length=T)
KGE_PRI_S4<- vector(length=T)
KGNE_PRI_S1<- vector(length=T)
KGNE_PRI_S2<- vector(length=T)
KGNE_PRI_S3<- vector(length=T)
KGNE_PRI_S4<- vector(length=T)
KCE_PRI_S1<- vector(length=T)
KCE_PRI_S2<- vector(length=T)
KCE_PRI_S3<- vector(length=T)
KCE_PRI_S4<- vector(length=T)
KCNE_PRI_S1<- vector(length=T)
KCNE_PRI_S2<- vector(length=T)
KCNE_PRI_S3<- vector(length=T)
KCNE_PRI_S4<- vector(length=T)
KSEQ_PRI_S1<- vector(length=T)
KSEQ_PRI_S2<- vector(length=T)
KGE<- vector(length=T)
KGNE<- vector(length=T)
KCE<- vector(length=T)
KCNE<- vector(length=T)
KSEQ<- vector(length=T)
delta<- vector(length=T)
v<- vector(length=T)
g_lambda<- vector(length=T)
sigma_0<- vector(length=T)
lambda<- vector(length=T)
wage<- vector(length=T)
N<- vector(length=T)
ur<- vector(length=T)
b_C<- vector(length=T)
b_G<- vector(length=T)
x_1<- vector(length=T)
x_2<- vector(length=T)
x_20<- vector(length=T)
g_x_20<- vector(length=T)
yield_C<- vector(length=T)
yield_G<- vector(length=T)
coupon_C<- vector(length=T)
coupon_G<- vector(length=T)
B_C<- vector(length=T)
B_G<- vector(length=T)
p_C<- vector(length=T)
p_G<- vector(length=T)
B<- vector(length=T)
DL<- vector(length=T)
def<- vector(length=T)
illiq<- vector(length=T)
dsr<- vector(length=T)

#2.2.3 Households 
Y_HG<- vector(length=T)
Y_H<- vector(length=T)
CO_PRI_N<- vector(length=T)
CO_PRI<- vector(length=T)
V_HF<- vector(length=T)
SEC_H<- vector(length=T)
B_CH<- vector(length=T)
B_GH<- vector(length=T)
D_N<- vector(length=T)
D<- vector(length=T)
lambda_30<- vector(length=T)
g_lambda_30<- vector(length=T)
lambda_40<- vector(length=T)
b_CH<- vector(length=T)
b_GH<- vector(length=T)
DC<- vector(length=T)
g_POP<- vector(length=T)
POP<- vector(length=T)
LF<- vector(length=T)
lf_1<- vector(length=T)

#2.2.4 Commercial banks 
BP<- vector(length=T)
CAP<- vector(length=T)
BP_U<- vector(length=T)
BP_D<- vector(length=T)
HPM<- vector(length=T)
SEC_B<- vector(length=T)
A<- vector(length=T)
w_G<- vector(length=T)
w_C_S1<- vector(length=T)
w_C_S2<- vector(length=T)
w_C_S3<- vector(length=T) 
w_C_S4<- vector(length=T)
CR<- vector(length=T)
CR_G<- vector(length=T)
CR_C_S1<- vector(length=T)
CR_C_S2<- vector(length=T)
CR_C_S3<- vector(length=T)
CR_C_S4<- vector(length=T)
CR_C<- vector(length=T)
LC_S1<- vector(length=T)
LC_S2<- vector(length=T)
LC_S3<- vector(length=T)
LC_S4<- vector(length=T)
LG_S1<- vector(length=T)
LG_S2<- vector(length=T)
LG_S3<- vector(length=T)
LG_S4<- vector(length=T)
LC<- vector(length=T)
LG<- vector(length=T)
lev_B<- vector(length=T)
CAR<- vector(length=T)
w_LT<- vector(length=T)
sh_NLG<- vector(length=T)
sh_NLC_S1<- vector(length=T)
sh_NLC_S2<- vector(length=T)
sh_NLC_S3<- vector(length=T)
sh_NLC_S4<- vector(length=T)
sh_LG<- vector(length=T)
sh_LC_S1<- vector(length=T)
sh_LC_S2<- vector(length=T)
sh_LC_S3<- vector(length=T)
sh_LC_S4<- vector(length=T)
int_G<- vector(length=T)
int_C_S1<- vector(length=T)
int_C_S2<- vector(length=T)
int_C_S3<- vector(length=T)
int_C_S4<- vector(length=T)
int<- vector(length=T)
spr<- vector(length=T)
spr_G<- vector(length=T)
spr_C_S1<- vector(length=T)
spr_C_S2<- vector(length=T)
spr_C_S3<- vector(length=T)
spr_C_S4<- vector(length=T)

#2.2.5 Government sector
GNS <- vector(length=T)
SEC<- vector(length=T)
IG_GOV <- vector(length=T)
IC_GOV <- vector(length=T)
I_GOV <- vector(length=T)
KG_GOV<- vector(length=T)
KC_GOV <- vector(length=T)
K_GOV <- vector(length=T)
K<- vector(length=T)
KG<- vector(length=T)
KC<- vector(length=T)
CO_GOV<- vector(length=T)
SUB<- vector(length=T)
gov_SUB<- vector(length=T)
TAX_H<- vector(length=T)
TAX_F<- vector(length=T)
TAX_C<- vector(length=T)
tau_C<- vector(length=T)
TAX<- vector(length=T)

#2.2.6 Central banks
CBP<- vector(length=T)
B_GCB<- vector(length=T)
B_CCB<- vector(length=T)
b_CCB<- vector(length=T)
b_GCB<- vector(length=T)
SEC_CB<- vector(length=T)
SEC_CBred<- vector(length=T)

#Auxiliary variables 
g_Y<- vector(length=T)
gI<- vector(length=T)
gC<- vector(length=T)
g_K<- vector(length=T)
ur_per<- vector(length=T)
g_Y_per<- vector(length=T)
theta_per<- vector(length=T)
greenI<- vector(length=T)
r_total<- vector(length=T)
Y_POP_ratio<- vector(length=T)
I_Y_ratio<- vector(length=T)
I_K_ratio<- vector(length=T)
C_Y_ratio<- vector(length=T)
C_K_ratio<- vector(length=T)
Y_K_ratio<- vector(length=T)
Y_H_Y_ratio<- vector(length=T)
Y_H_K_ratio<- vector(length=T)
A_K_ratio<- vector(length=T)
zero<- vector(length=T)
sigma_0_optimal<- vector(length=T)
h<- vector(length=T)
h_optimal<- vector(length=T)
g_h<- vector(length=T)
omega_ratio<- vector(length=T)
mu_ratio<- vector(length=T)
rho_ratio<- vector(length=T)
epsilon_ratio<- vector(length=T)
Wbill<- vector(length=T)
Interest<- vector(length=T)
Depreciation<- vector(length=T)
L_K<- vector(length=T)
I<- vector(length=T)
SEC_Y<- vector(length=T)
SEC_Y_per<- vector(length=T)
fiscal_balance_per<- vector(length=T)
fiscal_balance<- vector(length=T)
PORT_BCH<- vector(length=T)
PORT_BGH<-vector(length=T)
PORT_SECH<-vector(length=T)
PORT_D<-vector(length=T)
SEC_BN<-vector(length=T)
sh_L<-vector(length=T)
g_bC<-vector(length=T)
g_bG<-vector(length=T)
g_pG<-vector(length=T)
g_pC<-vector(length=T)
B_C_issue<-vector(length=T)
B_G_issue<-vector(length=T)
IG_cum<- vector(length=T)
IC_cum<- vector(length=T)
V_CB<- vector(length=T)
gGOVCO<- vector(length=T)
HPM_K_ratio<- vector(length=T)
D_K_ratio<- vector(length=T)
ID_K_ratio<- vector(length=T)
A_N<- vector(length=T)
g_lf<- vector(length=T)
W_POP_ratio<- vector(length=T)
E_N_ratio<- vector(length=T)
Y_E_ratio<- vector(length=T)
g_EN_ratio<- vector(length=T)
g_YE_ratio<- vector(length=T)
lambda_perworker<- vector(length=T)
random1<- vector(length=T)
random2<- vector(length=T)
BAILOUT<- vector(length=T)
BAILOUT_levB<- vector(length=T)
BAILOUT_CAR<- vector(length=T)
CAP_before<- vector(length=T)
CAR_before<- vector(length=T)
lev_B_before<- vector(length=T)
E_ratio<- vector(length=T)
CO2_ratio<- vector(length=T)
lf<- vector(length=T)
haz_flow<- vector(length=T)
NLG_DN_S1<- vector(length=T)
NLG_DN_S2<- vector(length=T)
NLG_DN_S3<- vector(length=T)
NLG_DN_S4<- vector(length=T)
NLC_DN_S1<- vector(length=T)
NLC_DN_S2<- vector(length=T)
NLC_DN_S3<- vector(length=T)
NLC_DN_S4<- vector(length=T)
V_H<- vector(length=T)
Y_HD<- vector(length=T)
KG_GOV_Y_ratio<- vector(length=T)
KC_GOV_Y_ratio<- vector(length=T)
IG_GOV_Y_ratio<- vector(length=T)
IC_GOV_Y_ratio<- vector(length=T)
KG_Y_ratio<- vector(length=T)
KC_Y_ratio<- vector(length=T)
KGE_KCE<- vector(length=T)
KGNE_KCNE<- vector(length=T)
KSEQ_KCE<- vector(length=T)
KG_K_ratio<- vector(length=T)
LG_L_ratio_pseudo<- vector(length=T)
LG_L_ratio <- vector(length=T)
g_Y_cum<- vector(length=T)
def_per<- vector(length=T)
CAR_per<- vector(length=T)
beta<- vector(length=T)
I_PRI_S1<- vector(length=T)
I_PRI_S2<- vector(length=T)
I_PRI_S3<- vector(length=T)
I_PRI_S4<- vector(length=T)
IG_PRI_D<- vector(length=T)
IC_PRI_D<- vector(length=T)
TAX_C_change<- vector(length=T)
SUB_change<- vector(length=T)
IG_GOV_change<- vector(length=T)
endog<- vector(length=T)
tau_C_baseline<- vector(length=T)
tau_C_policy_High<- vector(length=T)
tau_C_policy_Medium<- vector(length=T)
tau_C_baseline_per<- vector(length=T)
tau_C_policy_High_per<- vector(length=T)
tau_C_policy_Medium_per<- vector(length=T)
EMIS_IN_Y<- vector(length=T)
EMIS_IN_Y_SSP360<- vector(length=T)
g_Y_per_mean_max_policy<- vector(length=T)
EMIS_mean_max_policy<- vector(length=T)
T_AT_mean_max_policy<- vector(length=T)
g_Y_per_mean_min_policy<- vector(length=T)
EMIS_mean_min_policy<- vector(length=T)
T_AT_mean_min_policy<- vector(length=T)
fiscal_balance_help<- vector(length=T)
private_balance<- vector(length=T)
firms_balance<- vector(length=T)
households_balance<- vector(length=T)
banks_balance<- vector(length=T)
private_balance_per<- vector(length=T)
firms_balance_per<- vector(length=T)
households_balance_per<- vector(length=T)
banks_balance_per<- vector(length=T)
fiscal_balance_help_per<- vector(length=T)
g_POP_stag<- vector(length=T)
IG_E<- vector(length=T) 
IG_E_Y<- vector(length=T) 
IG_E_Y_cum<- vector(length=T) 
yield_C_per<- vector(length=T)
yield_G_per<- vector(length=T)

sce=1
for (j in 1:loops){
  if (j==2 | j==3) {monte<-monte_sce} else {monte<-1}
  
  ########################
  #Monte Carlo variables
  ########################
  
  assign(paste("Monte_g_Y_per",sce, sep=""), matrix(nrow=T,ncol=monte))
  assign(paste("Monte_epsilon_ratio",sce, sep=""), matrix(nrow=T,ncol=monte))
  assign(paste("Monte_mu_ratio",sce, sep=""), matrix(nrow=T,ncol=monte))
  assign(paste("Monte_IG_PRI",sce, sep=""), matrix(nrow=T,ncol=monte))
  assign(paste("Monte_IG_GOV",sce, sep=""), matrix(nrow=T,ncol=monte))
  assign(paste("Monte_IG_E_Y",sce, sep=""), matrix(nrow=T,ncol=monte))
  assign(paste("Monte_SUB",sce, sep=""), matrix(nrow=T,ncol=monte))
  assign(paste("Monte_ur_per",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_Y",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_CO_PRI",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_I_PRI",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_lambda",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_N",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_Y_POP_ratio",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_lambda_perworker",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_EMIS_IN_Y",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_EMIS_IN",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_E_N_ratio",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_I_K_ratio",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_ur_per",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_POP",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_beta_S1",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_greenI",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_kappa",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_theta_per",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_r_total",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_dsr",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_L_K",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_L",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_E",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_SEC",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_CO_GOV",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_I",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_def_per",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_lev_B",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_PORT_BCH",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_PORT_BGH",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_p_C",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_p_G",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_h",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_EMIS",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_T_AT",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_dep_E",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_dep_M",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_CO2_CUM",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_rho",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_mu",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_epsilon",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_W",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_W_POP_ratio",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_hazratio",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_CAR_per",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_CR",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_CR_C_S1",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_CR_G",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_spr",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_spr_C_S1",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_spr_G",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_SEC_Y_per",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_fiscal_balance_per",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_fiscal_balance_help_per",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_firms_balance_per",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_households_balance_per",sce, sep=""), matrix(nrow=T,ncol=monte))
  assign(paste("Monte_banks_balance_per",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_illiq",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_LF",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_lf",sce, sep=""), matrix(nrow=T,ncol=monte)) 
  assign(paste("Monte_yield_C",sce, sep=""), matrix(nrow=T,ncol=monte))
  assign(paste("Monte_yield_G",sce, sep=""), matrix(nrow=T,ncol=monte))
  assign(paste("Monte_yield_C_per",sce, sep=""), matrix(nrow=T,ncol=monte))
  assign(paste("Monte_yield_G_per",sce, sep=""), matrix(nrow=T,ncol=monte))
  
  assign(paste("g_Y_per_mean",sce, sep=""), vector(length=T))
  assign(paste("epsilon_ratio_mean",sce, sep=""), vector(length=T))
  assign(paste("mu_ratio_mean",sce, sep=""), vector(length=T))
  assign(paste("IG_PRI_mean",sce, sep=""), vector(length=T))
  assign(paste("IG_GOV_mean",sce, sep=""), vector(length=T))
  assign(paste("IG_E_Y_mean",sce, sep=""), vector(length=T))
  assign(paste("SUB_mean",sce, sep=""), vector(length=T))
  assign(paste("Y_sdplus",sce, sep=""), vector(length=T))
  assign(paste("Y_sdminus",sce, sep=""), vector(length=T))
  assign(paste("CR_C_S1_sdplus",sce, sep=""), vector(length=T))
  assign(paste("CR_C_S1_sdminus",sce, sep=""), vector(length=T))
  assign(paste("EMIS_sdplus",sce, sep=""), vector(length=T))
  assign(paste("EMIS_sdminus",sce, sep=""), vector(length=T))
  assign(paste("r_total_sdplus",sce, sep=""), vector(length=T))
  assign(paste("r_total_sdminus",sce, sep=""), vector(length=T))
  assign(paste("def_per_sdplus",sce, sep=""), vector(length=T))
  assign(paste("def_per_sdminus",sce, sep=""), vector(length=T))
  assign(paste("lev_B_sdplus",sce, sep=""), vector(length=T))
  assign(paste("lev_B_sdminus",sce, sep=""), vector(length=T))
  assign(paste("ur_per_mean",sce, sep=""), vector(length=T))
  assign(paste("Y_mean",sce, sep=""), vector(length=T))
  assign(paste("CO_PRI_mean",sce, sep=""), vector(length=T))
  assign(paste("I_PRI_mean",sce, sep=""), vector(length=T))
  assign(paste("lambda_mean",sce, sep=""), vector(length=T))
  assign(paste("N_mean",sce, sep=""), vector(length=T))
  assign(paste("Y_POP_ratio_mean",sce, sep=""), vector(length=T))
  assign(paste("lambda_perworker_mean",sce, sep=""), vector(length=T))
  assign(paste("EMIS_IN_Y_mean",sce, sep=""), vector(length=T))
  assign(paste("EMIS_IN_mean",sce, sep=""), vector(length=T))
  assign(paste("E_N_ratio_mean",sce, sep=""), vector(length=T))
  assign(paste("I_K_ratio_mean",sce, sep=""), vector(length=T))
  assign(paste("ur_per_mean",sce, sep=""), vector(length=T))
  assign(paste("POP_mean",sce, sep=""), vector(length=T))
  assign(paste("beta_S1_mean",sce, sep=""), vector(length=T))
  assign(paste("greenI_mean",sce, sep=""), vector(length=T))
  assign(paste("kappa_mean",sce, sep=""), vector(length=T))
  assign(paste("theta_per_mean",sce, sep=""), vector(length=T))
  assign(paste("r_total_mean",sce, sep=""), vector(length=T))
  assign(paste("dsr_mean",sce, sep=""), vector(length=T))
  assign(paste("L_K_mean",sce, sep=""), vector(length=T))
  assign(paste("L_mean",sce, sep=""), vector(length=T))
  assign(paste("E_mean",sce, sep=""), vector(length=T))
  assign(paste("SEC_mean",sce, sep=""), vector(length=T))
  assign(paste("CO_GOV_mean",sce, sep=""), vector(length=T))
  assign(paste("I_mean",sce, sep=""), vector(length=T))
  assign(paste("def_per_mean",sce, sep=""), vector(length=T))
  assign(paste("lev_B_mean",sce, sep=""), vector(length=T))
  assign(paste("PORT_BCH_mean",sce, sep=""), vector(length=T))
  assign(paste("PORT_BGH_mean",sce, sep=""), vector(length=T))
  assign(paste("p_C_mean",sce, sep=""), vector(length=T))
  assign(paste("p_G_mean",sce, sep=""), vector(length=T))
  assign(paste("h_mean",sce, sep=""), vector(length=T))
  assign(paste("EMIS_mean",sce, sep=""), vector(length=T))
  assign(paste("T_AT_mean",sce, sep=""), vector(length=T))
  assign(paste("dep_E_mean",sce, sep=""), vector(length=T))
  assign(paste("dep_M_mean",sce, sep=""), vector(length=T))
  assign(paste("CO2_CUM_mean",sce, sep=""), vector(length=T))
  assign(paste("rho_mean",sce, sep=""), vector(length=T))
  assign(paste("mu_mean",sce, sep=""), vector(length=T))
  assign(paste("epsilon_mean",sce, sep=""), vector(length=T))
  assign(paste("W_mean",sce, sep=""), vector(length=T))
  assign(paste("W_POP_ratio_mean",sce, sep=""), vector(length=T))
  assign(paste("hazratio_mean",sce, sep=""), vector(length=T))
  assign(paste("CAR_per_mean",sce, sep=""), vector(length=T))
  assign(paste("CR_mean",sce, sep=""), vector(length=T))
  assign(paste("CR_C_S1_mean",sce, sep=""), vector(length=T))
  assign(paste("CR_G_mean",sce, sep=""), vector(length=T))
  assign(paste("spr_mean",sce, sep=""), vector(length=T))
  assign(paste("spr_C_S1_mean",sce, sep=""), vector(length=T))
  assign(paste("spr_G_mean",sce, sep=""), vector(length=T))
  assign(paste("SEC_Y_per_mean",sce, sep=""), vector(length=T))
  assign(paste("fiscal_balance_per_mean",sce, sep=""), vector(length=T))
  assign(paste("fiscal_balance_help_per_mean",sce, sep=""), vector(length=T))
  assign(paste("firms_balance_per_mean",sce, sep=""), vector(length=T))
  assign(paste("households_balance_per_mean",sce, sep=""), vector(length=T))
  assign(paste("banks_balance_per_mean",sce, sep=""), vector(length=T))
  assign(paste("illiq_mean",sce, sep=""), vector(length=T))
  assign(paste("LF_mean",sce, sep=""), vector(length=T))
  assign(paste("lf_mean",sce, sep=""), vector(length=T))
  assign(paste("yield_C_mean",sce, sep=""), vector(length=T))
  assign(paste("yield_G_mean",sce, sep=""), vector(length=T))
  assign(paste("yield_C_per_mean",sce, sep=""), vector(length=T))
  assign(paste("yield_G_per_mean",sce, sep=""), vector(length=T))
  
  ################################################
  #2. VALUES FOR PARAMETERS & VARIABLES 
  ################################################
  
  ###################################
  #Parameters necessary for the scenarios
  ###################################
  
  for (k in 1:monte) {
    GVA_S1<-4 #trillion US$ 
    GVA_S2<-17.37 #trillion US$ 
    GVA_S3<-6.8 #trillion US$ 
    GVA_S4<-47.7 #trillion US$ 
    GVA<- GVA_S1+GVA_S2+GVA_S3+GVA_S4# trillion US$
    carbon_S1<-15.30 #GtCO2 
    carbon_S2<-6.11#GtCO2
    carbon_S3<-8.03 #GtCO2 
    carbon_S4<-3.38 #GtCO2
    carbon<- carbon_S1+carbon_S2+carbon_S3+carbon_S4# GtCO2
    
    sh_GVA_S1<-GVA_S1/(GVA_S1+GVA_S2+GVA_S3+GVA_S4)
    sh_GVA_S2<-GVA_S2/(GVA_S1+GVA_S2+GVA_S3+GVA_S4)
    sh_GVA_S3<-GVA_S3/(GVA_S1+GVA_S2+GVA_S3+GVA_S4)
    sh_GVA_S4<-1-sh_GVA_S1-sh_GVA_S2-sh_GVA_S3
    sh_EMIS_IN_S1<-carbon_S1/(carbon_S1+carbon_S2+carbon_S3+carbon_S4)
    sh_EMIS_IN_S2<- carbon_S2/(carbon_S1+carbon_S2+carbon_S3+carbon_S4)
    sh_EMIS_IN_S3<- carbon_S3/(carbon_S1+carbon_S2+carbon_S3+carbon_S4)
    sh_EMIS_IN_S4<- 1-sh_EMIS_IN_S1-sh_EMIS_IN_S2-sh_EMIS_IN_S3
    dd_S1<- (carbon_S1/GVA_S1)/(carbon/GVA)
    dd_S2<-(carbon_S2/GVA_S2)/(carbon/GVA)
    dd_S3<- (carbon_S3/GVA_S3)/(carbon/GVA)
    dd_S4<- (carbon_S4/GVA_S4)/(carbon/GVA)
    g_baseline<-0.029 # auxiliary parameter 
    
    B_G_initial<-0.4#(Bi) trillion US$
    B_initial<-12.6 #(Bi) trillion US$
    B_C_initial<-B_initial-B_G_initial #auxiliary exogenous variable
    B_GCB_initial<- 0.1*0.25 #trillion US$; 0.1 is the proportion of green bonds in the total corporate banks hedl by central banks, $0.25tr is the estimated amount of corporate sector holdings
    s_G<-B_GCB_initial/(B_G_initial/(1+g_baseline)) #(Ci)
    B_CCB_initial<-0.25-B_GCB_initial # $0.25tr is the estimated amount of corporate sector holdings 
    s_C<-B_CCB_initial/(B_C_initial/(1+g_baseline)) #(Ci)
    prop<-0.72 
    IG_initial<- 0.7#energy-related green investment was US$0.5 trillion in 2018 
    IG_PRI_initial<- prop*IG_initial #trillion US$
    IG_GOV_initial<-(1-prop)*IG_initial #trillion US$
    Y_initial<-85.93# trillion US 
    w_G_initial<-1  #(Biii)
    w_C_S1_initial<-1 #(Biii)
    w_C_S2_initial<-1 #(Biii)
    w_C_S3_initial<-1 #(Biii)
    w_C_S4_initial<-1 #(Biii)
    w_G_2022<-1  #(Biii)
    w_C_S1_2022<-1 #(Biii)
    w_C_S2_2022<-1 #(Biii)
    w_C_S3_2022<-1 #(Biii)
    w_C_S4_2022<-1 #(Biii)
    r_3_dCR_dCAR_max<- -1.40683353# econometrically estimated coefficient
    l_1<-1 
    spr_1<- 0.03122215 #econometrically estimated coefficient
    spr_3<-1
    beta_2<-1 #(Biii)
    gov_IG<-IG_GOV_initial/(Y_initial/(1+g_baseline))
    con_change<-0 #auxiliary parameter
    tau_C_dummy_H<-0
    tau_C_dummy_M<-0
    CR_dummy<-0
    damage_dummy<-1 #auxiliary parameter
    random_dummy<-0 #auxiliary parameter
    stop_recycling_dummy<-0
    subsidy_increase_dummy<-0
    subsidy_ratio<-1
    
    for (i in 1:T) {
      ##########################
      #Parameters per scenario
      ##########################
      
      #########
      #Baseline
      #########
      if (j==1){} 
      ###########################################
      #Baseline scenario with randomness (Appendix B)
      ############################################
      if (j==2){
        random_dummy<-1 
      } 
      ####################################################################
      #Baseline without damage feedback effects and with randomness (Appendix B)
      ####################################################################
      if (j==3) {
        damage_dummy<-0 
        random_dummy<-1 
      }
      #####
      #GSF
      #####
      if (j==4 & i<5) {}
      if (j==4 & i>=5) {
        w_G_2022<-1-0.25
      }
      #####
      #DPF
      #####
      if (j==5 & i<5) {}
      if (j==5 & i>=5) {
        w_C_S1_2022<-1+min(dd_S1,1)*0.25  
        w_C_S2_2022<-1+min(dd_S2,1)*0.25  
        w_C_S3_2022<-1+ min(dd_S3,1)*0.25  
        w_C_S4_2022<-1+ min(dd_S4,1)*0.25  
      }
      #########
      #GSF+DPF
      #########
      if (j==6 & i<5) {}
      if (j==6 & i>=5) {
        w_G_2022<-1-0.25
        w_C_S1_2022<-1+min(dd_S1,1)*0.25
        w_C_S2_2022<-1+min(dd_S2,1)*0.25  
        w_C_S3_2022<-1+ min(dd_S3,1)*0.25  
        w_C_S4_2022<-1+ min(dd_S4,1)*0.25  
      }
      ########################################
      #High carbon tax + green subsidy 
      ########################################
      if (j==7 & i<5) {}
      if (j==7 & i>=5) {
        tau_C_dummy_H<-1
      }
      ###########################################
      #High carbon tax+green subsidy+GSF 
      ###########################################
      if (j==8 & i<5) {}
      if (j==8 & i>=5) {
        w_G_2022<-1-0.25
        tau_C_dummy_H<-1
      }
      ################################
      #High Carbon tax+green subsidy+DPF 
      ################################
      if (j==9 & i<5) {}
      if (j==9 & i>=5) {
        w_C_S1_2022<-1+min(dd_S1,1)*0.25  
        w_C_S2_2022<-1+min(dd_S2,1)*0.25  
        w_C_S3_2022<-1+ min(dd_S3,1)*0.25  
        w_C_S4_2022<-1+ min(dd_S4,1)*0.25  
        tau_C_dummy_H<-1
      }
      ################################
      #Medium carbon tax + green subsidy
      ################################
      if (j==10 & i<5) {}
      if (j==10 & i>=5) {
        tau_C_dummy_M<-1
      }
      ###################################
      #Medium carbon tax+green subsidy+GSF 
      ###################################
      if (j==11 & i<5) {}
      if (j==11 & i>=5) {
        w_G_2022<-1-0.25
        tau_C_dummy_M<-1
      }
      ###################################
      #Medium carbon tax+green subsidy+DPF 
      ###################################
      if (j==12 & i<5) {}
      if (j==12 & i>=5) {
        w_C_S1_2022<-1+min(dd_S1,1)*0.25  
        w_C_S2_2022<-1+min(dd_S2,1)*0.25  
        w_C_S3_2022<-1+ min(dd_S3,1)*0.25  
        w_C_S4_2022<-1+ min(dd_S4,1)*0.25  
        tau_C_dummy_M<-1
      }
      #############################
      #Sensitivity [Baseline: minimum]
      #############################
      if (j==13 & i<5){
        r_3_dCR_dCAR_max<- -0.88360892 
        l_1<-1*(1-0.5)
        spr_1<-  0.0166877
        spr_3<-1*(1-0.5)
        beta_2<-1*(1-0.5)
      }
      if (j==13 & i>=5){
        r_3_dCR_dCAR_max<- -0.88360892 
        l_1<-1*(1-0.5)
        spr_1<-  0.0166877
        spr_3<-1*(1-0.5)
        beta_2<-1*(1-0.5)
      }
      #############################
      #Sensitivity [Baseline: maximum]
      ##############################
      if (j==14 & i<5){
        r_3_dCR_dCAR_max<- -2.13561067 
        l_1<-1*(1+0.5) 
        spr_1<- 0.04306504
        spr_3<-1*(1+0.5)
        beta_2<-1*(1+0.5) 
      }
      if (j==14 & i>=5){
        r_3_dCR_dCAR_max<- -2.13561067 
        l_1<-1*(1+0.5) 
        spr_1<- 0.04306504
        spr_3<-1*(1+0.5)
        beta_2<-1*(1+0.5)
      }
      #########################
      #Sensitivity [GSF: minimum]
      #########################
      if (j==15 & i<5){
        r_3_dCR_dCAR_max<- -0.88360892 
        l_1<-1*(1-0.5)
        spr_1<-  0.0166877
        spr_3<-1*(1-0.5)
        beta_2<-1*(1-0.5)
      }
      if (j==15 & i>=5){
        w_G_2022<-1-0.25  
        r_3_dCR_dCAR_max<- -0.88360892 
        l_1<-1*(1-0.5)
        spr_1<-  0.0166877
        spr_3<-1*(1-0.5)
        beta_2<-1*(1-0.5)
      }
      #########################
      #Sensitivity [GSF: maximum]
      #########################
      if (j==16 & i<5){
        r_3_dCR_dCAR_max<- -2.13561067 
        l_1<-1*(1+0.5) 
        spr_1<- 0.04306504
        spr_3<-1*(1+0.5)
        beta_2<-1*(1+0.5) 
      }
      if (j==16 & i>=5){
        w_G_2022<-1-0.25  
        r_3_dCR_dCAR_max<- -2.13561067 
        l_1<-1*(1+0.5) 
        spr_1<- 0.04306504
        spr_3<-1*(1+0.5)
        beta_2<-1*(1+0.5)
      }
      #########################
      #Sensitivity [DPF: minimum]
      #########################
      if (j==17 & i<5){
        r_3_dCR_dCAR_max<- -0.88360892 
        l_1<-1*(1-0.5)
        spr_1<-  0.0166877
        spr_3<-1*(1-0.5)
        beta_2<-1*(1-0.5)
      }
      if (j==17 & i>=5){
        w_C_S1_2022<-1+min(dd_S1,1)*0.25  
        w_C_S2_2022<-1+min(dd_S2,1)*0.25  
        w_C_S3_2022<-1+ min(dd_S3,1)*0.25  
        w_C_S4_2022<-1+ min(dd_S4,1)*0.25  
        r_3_dCR_dCAR_max<- -0.88360892 
        l_1<-1*(1-0.5)
        spr_1<-0.0166877
        spr_3<-1*(1-0.5)
        beta_2<-1*(1-0.5)
      }
      #########################
      #Sensitivity [DPF: maximum]
      #########################
      if (j==18 & i<5){
        r_3_dCR_dCAR_max<- -2.13561067  
        l_1<-1*(1+0.5) 
        spr_1<- 0.04306504
        spr_3<-1*(1+0.5)
        beta_2<-1*(1+0.5)
      }
      if (j==18 & i>=5){
        w_C_S1_2022<-1+min(dd_S1,1)*0.25  
        w_C_S2_2022<-1+min(dd_S2,1)*0.25  
        w_C_S3_2022<-1+ min(dd_S3,1)*0.25  
        w_C_S4_2022<-1+ min(dd_S4,1)*0.25  
        r_3_dCR_dCAR_max<- -2.13561067 
        l_1<-1*(1+0.5) 
        spr_1<- 0.04306504
        spr_3<-1*(1+0.5)
        beta_2<-1*(1+0.5)
      }
      #############################
      #Sensitivity [GSF+DPF: minimum]
      ##############################
      if (j==19 & i<5){
        r_3_dCR_dCAR_max<- -0.88360892 
        l_1<-1*(1-0.5)
        spr_1<-  0.0166877
        spr_3<-1*(1-0.5)
        beta_2<-1*(1-0.5)
      }
      if (j==19 & i>=5){
        w_G_2022<-1-0.25  
        w_C_S1_2022<-1+min(dd_S1,1)*0.25  
        w_C_S2_2022<-1+min(dd_S2,1)*0.25  
        w_C_S3_2022<-1+ min(dd_S3,1)*0.25  
        w_C_S4_2022<-1+ min(dd_S4,1)*0.25  
        r_3_dCR_dCAR_max<- -0.88360892 
        l_1<-1*(1-0.5)
        spr_1<-  0.0166877
        spr_3<-1*(1-0.5)
        beta_2<-1*(1-0.5) 
      }
      ##############################
      #Sensitivity [GSF+DPF: maximum]
      ##############################
      if (j==20 & i<5){
        r_3_dCR_dCAR_max<- -2.13561067 
        l_1<-1*(1+0.5) 
        spr_1<- 0.04306504
        spr_3<-1*(1+0.5)
        beta_2<-1*(1+0.5)
      }
      if (j==20 & i>=5){
        w_G_2022<-1-0.25  
        w_C_S1_2022<-1+min(dd_S1,1)*0.25  
        w_C_S2_2022<-1+min(dd_S2,1)*0.25  
        w_C_S3_2022<-1+ min(dd_S3,1)*0.25  
        w_C_S4_2022<-1+ min(dd_S4,1)*0.25  
        r_3_dCR_dCAR_max<- -2.13561067 
        l_1<-1*(1+0.5) 
        spr_1<- 0.04306504
        spr_3<-1*(1+0.5)
        beta_2<-1*(1+0.5)
      }
      ################################################
      #Sensitivity [High carbon tax+Green subsidy: minimum]
      ################################################
      if (j==21 & i<5){
        r_3_dCR_dCAR_max<- -0.88360892 
        l_1<-1*(1-0.5)
        spr_1<-  0.0166877
        spr_3<-1*(1-0.5)
        beta_2<-1*(1-0.5)
      }
      if (j==21 & i>=5){
        r_3_dCR_dCAR_max<- -0.88360892 
        l_1<-1*(1-0.5)
        spr_1<-  0.0166877
        spr_3<-1*(1-0.5)
        beta_2<-1*(1-0.5) 
        tau_C_dummy_H<-1
      }
      ################################################
      #Sensitivity [High carbon tax+Green subsidy: maximum]
      ################################################
      if (j==22 & i<5){
        r_3_dCR_dCAR_max<- -2.13561067 
        l_1<-1*(1+0.5) 
        spr_1<- 0.04306504
        spr_3<-1*(1+0.5)
        beta_2<-1*(1+0.5)
      }
      if (j==22 & i>=5){
        r_3_dCR_dCAR_max<- -2.13561067 
        l_1<-1*(1+0.5) 
        spr_1<- 0.04306504
        spr_3<-1*(1+0.5)
        beta_2<-1*(1+0.5)
        tau_C_dummy_H<-1
      }
      ####################################################
      #Sensitivity [High carbon tax+Green subsidy+GFS: minimum]
      ####################################################
      if (j==23 & i<5){
        r_3_dCR_dCAR_max<- -0.88360892 
        l_1<-1*(1-0.5)
        spr_1<-  0.0166877
        spr_3<-1*(1-0.5)
        beta_2<-1*(1-0.5)
      }
      if (j==23 & i>=5){
        w_G_2022<-1-0.25  
        r_3_dCR_dCAR_max<- -0.88360892 
        l_1<-1*(1-0.5)
        spr_1<-  0.0166877
        spr_3<-1*(1-0.5)
        beta_2<-1*(1-0.5)
        tau_C_dummy_H<-1
      }
      ####################################################
      #Sensitivity [High carbon tax+Green subsidy+GFS: maximum]
      ####################################################
      if (j==24 & i<5){
        r_3_dCR_dCAR_max<- -2.13561067 
        l_1<-1*(1+0.5) 
        spr_1<- 0.04306504
        spr_3<-1*(1+0.5)
        beta_2<-1*(1+0.5)
      }
      if (j==24 & i>=5){
        w_G_2022<-1-0.25  
        r_3_dCR_dCAR_max<- -2.13561067 
        l_1<-1*(1+0.5) 
        spr_1<- 0.04306504
        spr_3<-1*(1+0.5)
        beta_2<-1*(1+0.5)
        tau_C_dummy_H<-1
      }
      ####################################################
      #Sensitivity [High carbon tax+DPF+Green subsidy: minimum]
      ####################################################
      if (j==25 & i<5){
        r_3_dCR_dCAR_max<- -0.88360892 
        l_1<-1*(1-0.5)
        spr_1<-  0.0166877
        spr_3<-1*(1-0.5)
        beta_2<-1*(1-0.5)
      }
      if (j==25 & i>=5){
        w_C_S1_2022<-1+min(dd_S1,1)*0.25  
        w_C_S2_2022<-1+min(dd_S2,1)*0.25  
        w_C_S3_2022<-1+ min(dd_S3,1)*0.25  
        w_C_S4_2022<-1+ min(dd_S4,1)*0.25  
        r_3_dCR_dCAR_max<- -0.88360892 
        l_1<-1*(1-0.5)
        spr_1<-  0.0166877
        spr_3<-1*(1-0.5)
        beta_2<-1*(1-0.5)
        tau_C_dummy_H<-1
      }
      ################################################
      #Sensitivity [High carbon tax+Green subsidy: maximum]
      ################################################
      if (j==26 & i<5){
        r_3_dCR_dCAR_max<- -2.13561067 
        l_1<-1*(1+0.5) 
        spr_1<- 0.04306504
        spr_3<-1*(1+0.5)
        beta_2<-1*(1+0.5)
      }
      if (j==26 & i>=5){
        w_C_S1_2022<-1+min(dd_S1,1)*0.25  
        w_C_S2_2022<-1+min(dd_S2,1)*0.25  
        w_C_S3_2022<-1+ min(dd_S3,1)*0.25  
        w_C_S4_2022<-1+ min(dd_S4,1)*0.25  
        r_3_dCR_dCAR_max<- -2.13561067 
        l_1<-1*(1+0.5) 
        spr_1<- 0.04306504
        spr_3<-1*(1+0.5)
        beta_2<-1*(1+0.5)
        tau_C_dummy_H<-1
      }
      #####################################################
      #Medium carbon tax without recycling 0f carbon tax revenues
      #####################################################
      if (j==27 & i<5){}
      if (j==27 & i>=5){
        stop_recycling_dummy<-1
        tau_C_dummy_M<-1
      }
      ######################################
      #Green subsidy increase (without recycling)
      ######################################
      if (j==28 & i<5){}
      if (j==28 & i>=5){
        stop_recycling_dummy<-1
        subsidy_increase_dummy<-1
        subsidy_ratio<-20
      }
      
      
      if (i == 1) {
        for (iterations in 1:15){
          
          #######################
          #Non-scenario parameters
          #######################
          #(A) Econometrically estimated parameters
          #(B) Directly calibrated parameters using (i) data, (ii) previous studies or (iii) a reasonable range of values
          #(C) Parameters that have been indirectly calibrated such that the model (i) matches the initial values obtained from the data or (ii) generates the baseline scenario 
          
          #Auxiliary parameters for the simulation of the model
          cr_rationing_dummy<-1 #auxiliary parameter
          stagnation_dummy<-1 #auxiliary parameter
          lambda_is_optimal<-0 #auxiliary parameter
          beta_is_endogenous_different<-1 #auxiliary parameter
          ucr_dummy<-1
          cr_adjust<-0.0138*cr_rationing_dummy #auxiliary parameter
          dam_adjust<-0.0015*damage_dummy #auxiliary parameter
          stag_adjust<-0.5*stagnation_dummy #auxiliary parameter
          ############
          beta_dummy_0<-0#auxiliary parameter
          h_is_optimal<-0 #auxiliary parameter
          
          #2.1.1 Matter, recycling and waste
          haz<-0.04 #(Bi) 
          con_M<-0.0015#(Biii)
          xi<-(g_baseline*(DEM[i]-mu[i]*delta[i]*(K[i]/(1+g_baseline))))/(CO_PRI[i]*mu[i]-DEM[i]+mu[i]*delta[i]*(K[i]/(1+g_baseline))) #ensures that the growth rate of DC is equal to g_baseline; (Cii)
          
          #2.1.2 Energy
          con_E<-0.003 #(Biii)
          
          #2.1.3 Emissions and climate change
          t_1<-0.5 #(Bii)
          t_2<-1.1 #(Bii)
          car<-3.667 #auxiliary parameter; (Bii)
          phi<-1.72/(car*1000) # 0C/GtCO2; 1.72 is TCRE measured in 0C/TtC; (Bii)
          
          #2.1.4 Ecological efficiency and technology
          zeta_1<-0.0001 #(Cii)
          g_omega_initial<--0.004 #auxiliary parameter; (Cii)
          mu_max<-1.5 #(Biii)
          mu_min<-0.3 #(Biii)
          mu_rat<-0.9 #auxiliary parameter
          mu_2049<-mu_rat*mu[i] #auxiliary parameter; (Cii)
          rho_max<-0.8 #(Biii)
          rho_rat<-1.4 #auxiliary parameter
          rho_2049<-rho_rat*rho[i] #auxiliary parameter; (Cii)
          epsilon_max<-12 #(Biii)
          epsilon_min<-2 #(Biii)
          epsilon_rat<-0.71 #auxiliary parameter; (Cii)
          epsilon_2049<-epsilon_rat*epsilon[i] #auxiliary parameter
          theta_2049<-0.23 #auxiliary parameter
          seq_2049<-0.009 #auxiliary parameter 
          KGE_KCE_2049<-0.0929#auxiliary parameter
          KGNE_KCNE_2049<-0.0453#auxiliary parameter
          KSEQ_KCE_2049<-0.0017 #auxiliary parameter
          pi_2<-(log((mu[i]-mu_min)/(mu_max-mu[i]))-log((mu_2049-mu_min)/(mu_max-mu_2049)))/(KGNE_KCNE_2049-1*(KGNE[i])/(KCNE[i])) #(Cii)
          pi_1<-((mu[i]-mu_min)/(mu_max-mu[i]))*exp(pi_2*1*(KGNE[i])/(KCNE[i])) #(Cii)
          pi_4<-(log(rho_max-rho[i])-log(rho_max-rho_2049)-log(rho[i])+log(rho_2049))/(KGNE_KCNE_2049-1*(KGNE[i])/(KCNE[i])) #(Cii)
          pi_3<-((rho_max-rho[i])*exp(pi_4*1*(KGNE[i])/(KCNE[i])))/rho[i] #(Cii)
          pi_6<-(log((epsilon[i]-epsilon_min)/(epsilon_max-epsilon[i]))-log((epsilon_2049-epsilon_min)/(epsilon_max-epsilon_2049)))/(KGE_KCE_2049-1*(KGE[i])/(KCE[i])) #(Cii)
          pi_5<-((epsilon[i]-epsilon_min)/(epsilon_max-epsilon[i]))*exp(pi_6*1*(KGE[i])/(KCE[i])) #(Cii)
          pi_8<-(log(1-theta[i])-log(1-theta_2049)-log(theta[i])+log(theta_2049))/(KGE_KCE_2049-1*(KGE[i])/(KCE[i]))
          pi_7<-((1-theta[i])*exp(pi_8*1*(KGE[i])/(KCE[i])))/theta[i] #(Cii)
          pi_10<-(log(1-seq[i])-log(1-seq_2049)-log(seq[i])+log(seq_2049))/(KSEQ_KCE_2049-1*(KSEQ[i])/(KCE_PRI_S1[i]+KCE_PRI_S2[i]))
          pi_9<-((1-seq[i])*exp(pi_10*1*(KSEQ[i])/(KCE_PRI_S1[i]+KCE_PRI_S2[i])))/seq[i] #(Cii)
          
          #2.2.1 Output determination and climate damages
          eta_1<-0 #(Bii)
          eta_2<-0.00284 #(Bii)
          eta_3<-0.0000819#(Bii)
          ad_K<-0.8 #(Biii)
          ad_LF<-0.95 #(Biii)
          ad_P<-0.95 #(Biii)
          p<-0.1 #(Biii)
          
          #2.2.2 Firms
          alpha_00<-(g_baseline+cr_adjust+dam_adjust+delta[1])*2 #stabilises I/K; (Cii)
          alpha_1_dI_du_max<-0.073#econometrically estimated coefficient
          alpha_1<-(alpha_1_dI_du_max*4)/(alpha_00+2*delta[1]) #(A)
          alpha_2_dI_dr_max<-0.084#econometrically estimated coefficient
          alpha_2<-(alpha_2_dI_dr_max*4)/(alpha_00+2*delta[1]) #(A)
          alpha_31_dI_dur_max<-0.006#econometrically estimated coefficient
          alpha_31<-(alpha_31_dI_dur_max*4)/(alpha_00+2*delta[1]) #(A)
          alpha_32<-0.5 #(Biii)
          alpha_41<-0.1 #(Biii)
          alpha_42<-0.99 #(Biii)
          alpha_51<-0.1 #(Biii)
          alpha_52<-0.99 #(Biii)
          delta_0<-0.048 #(Bi) 
          K_Y_ratio<-(0.24*(1+g_baseline))/(g_baseline+delta_0) #auxiliary parameter; (Bi) 
          s_W<-0.54 #(Bi) 
          hours<-1900 #auxiliary parameter; (Bi) 
          beta_1<-20 #(Cii)
          
          tau_F<-0.15 #(Biii)
          x_11<-0.25 #(Biii)
          x_10<-x_1[i]+x_11*yield_C[i]  #(Ci)
          x_21<-0.25 #(Biii)
          x_20[i]<-x_2[i]+x_21*yield_G[i] #(Ci)
          coupon_C[i]<-yield_C[i]*p_C[i] #(Ci)
          coupon_G[i]<-yield_G[i]*p_G[i] #(Ci)
          p_C_bar<-100 #(Biii)
          p_G_bar<-100 #(Biii)
          g_beta_initial<-0.001 #auxiliary parameter; (Cii)
          zeta_2<-0.1 #(Cii)
          beta_initial<-IG_PRI_D[i]/I_PRI_D[i] #auxiliary parameter; (Ci); stabilizes KG(PRI)/K(PRI)
          beta_initial_S1<-IG_PRI_D_S1[i]/I_PRI_D_S1[i] #auxiliary parameter; (Ci)
          beta_initial_S2<-IG_PRI_D_S2[i]/I_PRI_D_S2[i] #auxiliary parameter; (Ci)
          beta_initial_S3<-IG_PRI_D_S3[i]/I_PRI_D_S3[i] #auxiliary parameter; (Ci)
          beta_initial_S4<-IG_PRI_D_S4[i]/I_PRI_D_S4[i] #auxiliary parameter; (Ci)
          zeta_3<-0.006 #(Cii)
          sigma_0_initial<- -0.01393*(1/(1-zeta_3)) #auxiliary parameter; (Cii) 
          sigma_1<-0.0095 #(Cii)
          sigma_2<-0.825  #econometrically estimated; (A)  
          g_x_20[i]<- 0.02
          zeta_4<-0.4 #(Cii)
          endog[i]<-1
          zeta_7<- 0.1
          zeta_8<- 0.05
          
          #$0.328tr was the renewable energy, transmission and waste investment in 2018
          #$0.122tr was the low-carbon transport investment in 2018
          #$0.009tr was the agriculture energy-related in 2018
          #$0.5tr was the green energy-related investment in 2018
          #$0.5-0.328-0.122-0.009=0.041tr was the residual green energy-related investment in 2018 which is allocated among the sectors based on their GVA
          #$0.2tr is the assumed non-energy related investment in 2018 which is allocated among the sectors based on their GVA
          
          sh_GREEN_S1<-(0.328+sh_GVA_S1*(0.041+0.2))/IG_initial 
          sh_GREEN_S2<- sh_GVA_S2*(0.041+0.2)/IG_initial 
          sh_GREEN_S3<-(0.122+sh_GVA_S3*(0.041+0.2))/IG_initial 
          sh_GREEN_S4<-(0.009+sh_GVA_S4*(0.041+0.2))/IG_initial 
          gamma_E1<-(0.328+sh_GVA_S1*0.041)*prop/IG_PRI_S1[i]
          gamma_E2<- (sh_GVA_S2*0.041)*prop/IG_PRI_S2[i] 
          gamma_E3<-(0.122+sh_GVA_S3*0.041)*prop/IG_PRI_S3[i] 
          gamma_E4<-(0.009+sh_GVA_S4*0.041)*prop/IG_PRI_S4[i] 
          
          gamma_E<-0.5/IG_initial
          gamma_SEQ1<-0.000833/IG_PRI_S1[i] # 2.5/(3*1000) 2.5 billion is the cumulative sequestration investment in the power sector; we divide by 3 to approximate 2018 investment
          gamma_SEQ2<-0.0007/IG_PRI_S2[i] # 2.1/(3*1000)  2.1 billion is the cumulative sequestration investment in the power sector; we divide by 3 to approximate 2018 investment
          
          #2.2.3 Households 
          c_2<-0.0498 #econometrically estimated; (A)
          pr<-0.99 #(Biii)
          lambda_12<- -0.01 #(Biii)
          lambda_13<- -0.01 #(Biii)
          lambda_14<- -0.01 #(Biii)
          lambda_15<- -0.01 #(Biii)
          lambda_23<- -0.01 #(Biii)
          lambda_24<- -0.01 #(Biii)
          lambda_25<- -0.01 #(Biii)
          lambda_34<- -0.01 #(Biii)
          lambda_35<- -0.01 #(Biii)
          portfolio_dam<-0.2 #auxiliary parameter; (Biii)
          lambda_10prime<- 0.5*portfolio_dam #(Biii)
          lambda_20prime<- -1*portfolio_dam #(Biii)
          lambda_30prime<- 0*portfolio_dam #(Biii)
          lambda_40prime<- -lambda_10prime-lambda_20prime-lambda_30prime #(Biii)
          lambda_21<- lambda_12 #(Biii)
          lambda_41<- lambda_14 #(Biii)
          lambda_42<- lambda_24 #(Biii)
          lambda_43<- lambda_34 #(Biii)
          lambda_31<- lambda_13 #(Biii)
          lambda_32<- lambda_23 #(Biii)
          lambda_11<- -lambda_21- lambda_31-lambda_41 #(Biii)
          lambda_22<- -lambda_12-lambda_32-lambda_42 #(Biii)
          lambda_33<- -lambda_13-lambda_23-lambda_43 #(Biii)
          lambda_44<- -lambda_14-lambda_24-lambda_34 #(Biii)
          lambda_45<- -lambda_15-lambda_25-lambda_35 #(Biii)
          int_A<-0.03 #(Bi)  
          int_D<-0.025 #(Bi) 
          int_S<-0.025 #(Bi) 
          zeta_5<-0.016 #(Cii)
          zeta_6<-0.0003 #(Cii)
          lf_2<-0.001 #(Biii)
          lf_1[i]<-LF[i]/(POP[i]*(1-(1-ad_LF)*D_TF[i]))+lf_2*hazratio[i]#(Ci)
          PORT_BCH[i]<-B_CH[i]/ (V_HF[i]/(1+g_baseline)) #auxiliary parameter; (Ci)
          PORT_BGH[i]<-B_GH[i]/(V_HF[i]/(1+g_baseline)) #auxiliary parameter; (Ci)
          PORT_SECH[i]<-SEC_H[i]/(V_HF[i]/(1+g_baseline)) #auxiliary parameter; (Ci)
          PORT_D[i]<-1- PORT_BCH[i]-PORT_BGH[i]- PORT_SECH[i]  #auxiliary parameter; (Ci)
          lambda_10<-(SEC_H[i]/(V_HF[i]/(1+g_baseline)))-lambda_10prime*D_T[i]-lambda_11*int_S-lambda_12*yield_C[i]-lambda_13*yield_G[i]-lambda_14*int_D-lambda_15*(Y_H[i]/V_HF[i]) #(Ci)
          lambda_20<-(B_CH[i]/(V_HF[i]/(1+g_baseline)))-lambda_20prime*D_T[i]-lambda_21*int_S-lambda_22*yield_C[i]-lambda_23*yield_G[i]-lambda_24*int_D-lambda_25*(Y_H[i]/V_HF[i]) #(Ci)
          lambda_30[i]<-(B_GH[i]/(V_HF[i]/(1+g_baseline)))-lambda_30prime*D_T[i]-lambda_31*int_S-lambda_32*yield_C[i]-lambda_33*yield_G[i]-lambda_34*int_D-lambda_35*(Y_H[i]/V_HF[i]) #(Ci)
          lambda_40[i]<-1- lambda_10- lambda_20- lambda_30[i]#(Ci) 
          zeta_10<-0.18 #(Cii)
          g_lambda_30[i]<-zeta_10*g_bG[i]
          
          #2.2.4 Commercial banks
          rep<-0.1 #(Biii)
          h_1<-0.18 #(Bi) 
          def_max<-0.2 #(Biii)
          def_0<-(def_max/def[i])-1 #(Ci)
          def_2_ddef_dilliq_max<-0.25 #auxiliary parameter 
          def_2<-(def_2_ddef_dilliq_max)*((1+def_0)^2)/(def_max*def_0) #(Biii)
          def_1<-def_2*illiq[i] #(Ci)
          CR_max<-0.5 #(Biii)
          CAR_min<-0.08#(Biii)
          lev_B_max<-1/0.03#(Biii)
          CR_initial<-0.2 #auxiliary parameter; (Biii)
          h_2<-SEC_B[i]/D[i] #(Ci)
          w_H<-0 #(Biii)
          w_S<-0 #(Biii)
          
          w_G[i]<-w_G_initial 
          w_C_S1[i]<- w_C_S1_initial 
          w_C_S2[i]<- w_C_S2_initial
          w_C_S3[i]<- w_C_S3_initial
          w_C_S4[i]<- w_C_S4_initial
          sh_LG[i]<-LG[i]/L[i]
          sh_LC_S1[i]<-LC_S1[i]/L[i]
          sh_LC_S2[i]<-LC_S2[i]/L[i]
          sh_LC_S3[i]<-LC_S3[i]/L[i]
          sh_LC_S4[i]<-1-sh_LG[i]-sh_LC_S1[i]-sh_LC_S2[i]-sh_LC_S3[i]
          spr_initial<- 0.05
          spr_2<- 0.0243 #econometrically estimated coefficient
          spr_0<-spr_initial+spr_1*(CAR[i]-CAR_min)-spr_2*dsr[i]
          int_G[i]<-spr_initial+int_A
          int_C_S1[i]<-spr_initial+int_A
          int_C_S2[i]<-spr_initial+int_A
          int_C_S3[i]<-spr_initial+int_A
          int_C_S4[i]<-spr_initial+int_A
          spr[i]<-spr_initial
          spr_G[i]<-(1+spr_3*(w_G[i]-w_LT[i]))*spr[i]
          spr_C_S1[i]<-(1+spr_3*(w_C_S1[i]-w_LT[i]))*spr[i]
          spr_C_S2[i]<- (1+spr_3*(w_C_S2[i]-w_LT[i]))*spr[i]
          spr_C_S3[i]<- (1+spr_3*(w_C_S3[i]-w_LT[i]))*spr[i]
          spr_C_S4[i]<-(spr[i]-sh_LG[i]*spr_G[i]-sh_LC_S1[i]*spr_C_S1[i]-sh_LC_S2[i]*spr_C_S2[i]-sh_LC_S3[i]*spr_C_S3[i])/sh_LC_S4[i]
          
          #2.2.5 Government sector
          gov_C<-0.17 #(Bi)  
          gov_IC<- IC_GOV[i]/(Y[i]/(1+g_baseline))
          
          ###############################################
          #Parameters/endogenous variables for the steady state
          ###############################################
          
          c_1<-((Y[i]*g_baseline)/K_PRI[i]-(c_2*V_HF[i])/K_PRI[i]-g_baseline-delta[i]-(Y[i]*gov_C)/K_PRI[i] -(Y[i]*gov_IC)/K_PRI[i] -(Y[i]*gov_IG)/K_PRI[i]+Y[i]/K_PRI[i])/(Y_H[i]/K_PRI[i]) #stabilises Y/K; (Cii)
          
          s_F<-(I_PRI_D[i]-CR_C[i]*IC_PRI_D[i]-CR_G[i]*IG_PRI_D[i]-delta[i]*(K_PRI[i]-CR_C[i]*KC_PRI[i]-CR_G[i]*KG_PRI[i])/(1+g_baseline)-CR_C[i]*rep*LC[i]/(1+g_baseline)-CR_G[i]*rep*LG[i]/(1+g_baseline)-def[i]*L[i]/(1+g_baseline)-(1-CR_C[i])*x_1[i]*IC_PRI_D[i]-(1-CR_G[i])*x_2[i]*IG_PRI_D[i]-L[i]*g_baseline/(1+g_baseline))/(TP[i]*(1-CR_C[i]*(1-beta[i])-CR_G[i]*beta[i])/(1+g_baseline)) #stabilises L/K; (Cii)
          
          tau_H<-(gov_C+gov_IC+gov_IG-(tau_F*TP_G[i])/Y[i] -(tau_C[i]*EMIS_IN[i])/Y[i]+(int_S*SEC[i])/Y[i]-(coupon_C[i]*b_CCB[i]+coupon_G[i]*b_GCB[i]+int_A*A[i]+int_S*SEC_CB[i])/Y[i]+SUB[i]*(1+g_baseline)/Y[i]-(SEC[i]*g_baseline)/Y[i])/(Y_HG[i]/Y[i]) #stabilises SEC/Y; (Cii)
          
          s_B<-(def[i]*L[i]+CAP[i]*g_baseline)/BP[i] #stabilises CAP/K; (Cii)
          
          V_CB[i]<-0 #stabilises V_CB/K; (Cii)
          
          g_bC[i]<-(g_baseline-g_pC[i])/(1+ g_pC[i]) #stabilises the growth rate of B_C; (Cii)
          g_bG[i]<-(g_baseline-g_pG[i])/(1+ g_pG[i]) #stabilises the growth rate of B_G; (Cii)
          g_pC[i]<- 0.0 #(Cii)
          g_pG[i]<-0.0 #(Cii)
          
          x_1[i]<-(B_C[i]*g_bC[i])/(IC_PRI_D[i]*(1+g_bC[i])) #stabilises the growth rate of b_C; (Cii)
          x_2[i]<-(B_G[i]*g_bG[i])/(IG_PRI_D[i]*(1+g_bG[i])) #stabilises the growth rate of b_G; (Cii) 
          
          I_PRI[i]<-prop*0.24*Y[i]#trillion US$ 
          
          KG_GOV[i]<- gov_IG*Y[i]/(delta[i]+g_baseline) #stabilises KG_GOV/Y; trillion US$
          
          KC_GOV[i]<-K_GOV[i]-KG_GOV[i] #stabilises KC_GOV/Y; trillion US$ 
          
          alpha_01<-alpha_1*u[i]+alpha_2*r[i]-alpha_31*ur[i]^(-alpha_32)-alpha_41*(1-ue[i])^(-alpha_42)-alpha_51*(1-um[i])^(-alpha_52) #stabilises I/K; (Cii)
          
          ##################################
          #Initial values for endogenous variables 
          ##################################
          
          #2.1.1 Matter, waste and recycling
          MY[i]<- 52.22 #Gt 
          REC[i]<- 4.80 #Gt 
          M[i]<- MY[i]-REC[i]#Gt 
          DEM[i]<- 17.67 #Gt 
          SES[i]<-mu[i]*(K[i]+DC[i]) #Gt
          W[i]<-DEM[i]-REC[i]#Gt
          CEN[i]<-(1/car)*EMIS_IN[i] #Gt
          O2[i]<-EMIS_IN[i]-CEN[i] #Gt
          HW_CUM[i]<- 14.64 #Gt 
          hazratio[i]<-HW_CUM[i]/POP[i]
          REV_M[i]<-M[i]/dep_M[i] #Gt
          CON_M[i]<-con_M*RES_M[i] #Gt
          RES_M[i]<-(35/0.54-1)*REV_M[i] #Gt
          dep_M[i]<-0.02
          
          #2.1.2 Energy
          E[i]<-590 # EJ 
          E_NF[i]<-theta[i]*E[i] #EJ
          E_F[i]<-E[i]-E_NF[i] #EJ
          ED[i]<- E_F[i]+E_NF[i] #EJ
          REV_E[i]<- 40237 #EJ 
          CON_E[i]<-con_E*RES_E[i] #Gt
          RES_E[i]<- 550183 # EJ
          dep_E[i]<-E_F[i]/REV_E[i]
          
          #2.1.3 Emissions and climate change 
          EMIS_IN[i]<-36.60# GtCO2
          zeta_9<-0.014 #(Cii)
          g_EMIS_L[i]<-0.016
          EMIS_L[i]<- 5.53 # GtCO2
          EMIS[i]<- EMIS_IN[i]+EMIS_L[i] # GtCO2
          T_AT[i]<-1.14 #0C 
          CO2_CUM[i]<-T_AT[i]/(t_2*phi) # GtCO2
          
          #2.1.4 Ecological efficiency and technology
          omega[i]<-EMIS_IN[i]/(E_F[i]*(1-seq[i])) #GtCO2/EJ
          g_omega[i]<-g_omega_initial
          mu[i]<- MY[i]/(Y[i]-CO_GOV[i]) #Gt/trillion $ or kg/$
          rho[i]<- REC[i]/DEM[i]  
          epsilon[i]<-E[i]/Y[i] #EJ/trillion $
          theta[i]<-0.15
          seq[i]<- 80/(1000*EMIS_IN[i])# 80Mt is the value of sequestrated emissions in 2018 
          
          #2.2.1 Output determination
          Y_star_M[i]<-(REV_M[i]+REC[i])/mu[i] #trillion US$
          Y_star_E[i]<-REV_E[i]/((1-theta[i])*epsilon[i]) #trillion US$
          Y_star_K[i]<-v[i]*K_PRI[i] #trillion US$
          Y_star_N[i]<-lambda[i]*LF[i]*h[i] #trillion US$
          Y_star[i]<-min(min(Y_star_M[i], Y_star_E[i]),min(Y_star_N[i], Y_star_K[i])) #trillion US$
          Y[i]<-Y_initial#trillion US$
          um[i]<- (Y[i]-CO_GOV[i])/Y_star_M[i]
          ue[i]<-Y[i]/Y_star_E[i] 
          u[i]<-0.72
          re[i]<-1- ur[i] 
          D_T[i]<-(1-(1/(1+eta_1*T_AT[i]+eta_2*T_AT[i]^2+eta_3*T_AT[i]^6.754)))*damage_dummy
          D_TP[i]<-p*D_T[i]
          D_TF[i]<-1-(1-D_T[i])/(1-D_TP[i])
          
          #2.2.2 Firms
          TP_G[i]<-Y[i]-wage[i]*N[i]-int_C_S1[i]*LC_S1[i]/(1+g_baseline)-int_C_S2[i]*LC_S2[i]/(1+g_baseline)-int_C_S3[i]*LC_S3[i]/(1+g_baseline)-int_C_S4[i]*LC_S4[i]/(1+g_baseline)-int_G[i]*(LG_S1[i]/(1+g_baseline)+LG_S2[i]/(1+g_baseline)+LG_S3[i]/(1+g_baseline)+LG_S4[i]/(1+g_baseline))-delta[i]*(K_PRI[i]/(1+g_baseline))-coupon_C[i]*(b_C[i]/(1+g_bC[i]))-coupon_G[i]*(b_G[i]/(1+g_bG[i])) #trillion US$
          TP[i]<-TP_G[i]-TAX_F[i]-TAX_C[i]+SUB[i] #trillion US$
          RP[i]<-s_F*(TP[i]/(1+g_baseline)) #trillion US$
          DP[i]<-TP[i]-RP[i] #trillion US$
          r[i]<-TP[i]/K_PRI[i]
          random1[i]<-0
          I_PRI_D[i]<-(g_baseline+delta[i]+cr_adjust*cr_rationing_dummy)*(K_PRI[i]/(1+g_baseline)) #stabilises I_PRI_D/K; trillion US$
          I_PRI_D_S1[i]<-sh_GVA_S1*I_PRI_D[i] #trillion US$
          I_PRI_D_S2[i]<-sh_GVA_S2*I_PRI_D[i] #trillion US$
          I_PRI_D_S3[i]<-sh_GVA_S3*I_PRI_D[i] #trillion US$
          I_PRI_D_S4[i]<-sh_GVA_S4*I_PRI_D[i] #trillion US$
          
          IG_PRI_D_S1[i]<-(1+CR_G[i]*cr_rationing_dummy)*IG_PRI_S1[i] #trillion US$
          IG_PRI_D_S2[i]<- (1+CR_G[i]*cr_rationing_dummy)*IG_PRI_S2[i]#trillion US$
          IG_PRI_D_S3[i]<- (1+CR_G[i]*cr_rationing_dummy)*IG_PRI_S3[i]#trillion US$
          IG_PRI_D_S4[i]<- (1+CR_G[i]*cr_rationing_dummy)*IG_PRI_S4[i] #trillion US$
          IC_PRI_D_S1[i]<-I_PRI_D_S1[i]-IG_PRI_D_S1[i] #trillion US$
          IC_PRI_D_S2[i]<-I_PRI_D_S2[i]-IG_PRI_D_S2[i] #trillion US$
          IC_PRI_D_S3[i]<-I_PRI_D_S3[i]-IG_PRI_D_S3[i] #trillion US$
          IC_PRI_D_S4[i]<-I_PRI_D_S4[i]-IG_PRI_D_S4[i] #trillion US$
          
          tucr[i]<- ucr[i]*(1-gov_SUB[i])
          tucn[i]<- ucn[i]+ tau_C[i]*omega[i]*(1-seq[i])
          ucr[i]<- 1.2*ucn[i]
          g_ucr[i]<- 0.01
          ucn[i]<- 0.028
          g_ucn[i]<- 0.005
          
          beta_0_S1[i]<-beta_S1[i]+beta_1*sh_EMIS_IN_S1*(tucr[i]/(1-g_ucr[i])-tucn[i]/(1+g_ucn[i]))+beta_2*(sh_L[i]*(int_G[i]-int_C_S1[i])+(1-sh_L[i])*(yield_G[i]-yield_C[i]))
          beta_0_S2[i]<- beta_S2[i]+beta_1*sh_EMIS_IN_S2*( tucr[i]/(1-g_ucr[i])-tucn[i]/(1+g_ucn[i]))+beta_2*(sh_L[i]*(int_G[i]-int_C_S2[i])+(1-sh_L[i])*(yield_G[i]-yield_C[i]))
          beta_0_S3[i]<-beta_S3[i]+beta_1*sh_EMIS_IN_S3*( tucr[i]/(1-g_ucr[i])-tucn[i]/(1+g_ucn[i]))+beta_2*(sh_L[i]*(int_G[i]-int_C_S3[i])+(1-sh_L[i])*(yield_G[i]-yield_C[i])) 
          beta_0_S4[i]<-beta_S4[i]+beta_1*sh_EMIS_IN_S4*( tucr[i]/(1-g_ucr[i])-tucn[i]/(1+g_ucn[i]))+beta_2*(sh_L[i]*(int_G[i]-int_C_S4[i])+(1-sh_L[i])*(yield_G[i]-yield_C[i]))
          
          beta_S1[i]<-beta_initial*(1-beta_is_endogenous_different)+beta_initial_S1* beta_is_endogenous_different
          beta_S2[i]<-beta_initial*(1-beta_is_endogenous_different)+beta_initial_S2* beta_is_endogenous_different
          beta_S3[i]<-beta_initial*(1-beta_is_endogenous_different)+beta_initial_S3* beta_is_endogenous_different
          beta_S4[i]<-beta_initial*(1-beta_is_endogenous_different)+beta_initial_S4* beta_is_endogenous_different
          g_beta_0[i]<-g_beta_initial 
          NLG_D_S1[i]<-IG_PRI_D_S1[i]-sh_GVA_S1*beta_S1[i]*RP[i]+rep*LG_S1[i]/(1+g_baseline)-delta[i]*KG_PRI_S1[i]/(1+g_baseline)-sh_GVA_S1*p_G_bar*(b_G[i]-b_G[i]/(1+g_bG[i])) #trillion US$
          NLG_D_S2[i]<- IG_PRI_D_S2[i]-sh_GVA_S2*beta_S2[i]*RP[i]+rep*LG_S2[i]/(1+g_baseline)-delta[i]*KG_PRI_S2[i]/(1+g_baseline)-sh_GVA_S2*p_G_bar*(b_G[i]-b_G[i]/(1+g_bG[i])) #trillion US$
          NLG_D_S3[i]<- IG_PRI_D_S3[i]-sh_GVA_S3*beta_S3[i]*RP[i]+rep*LG_S3[i]/(1+g_baseline)-delta[i]*KG_PRI_S3[i]/(1+g_baseline)-sh_GVA_S3*p_G_bar*(b_G[i]-b_G[i]/(1+g_bG[i])) #trillion US$
          NLG_D_S4[i]<- IG_PRI_D_S4[i]-sh_GVA_S4*beta_S4[i]*RP[i]+rep*LG_S4[i]/(1+g_baseline)-delta[i]*KG_PRI_S4[i]/(1+g_baseline)-sh_GVA_S4*p_G_bar*(b_G[i]-b_G[i]/(1+g_bG[i])) #trillion US$
          NLC_D_S1[i]<-IC_PRI_D_S1[i]-sh_GVA_S1*(1-beta_S1[i])*RP[i]+rep*LC_S1[i]/(1+g_baseline)-delta[i]*KC_PRI_S1[i]/(1+g_baseline)-sh_GVA_S1*p_C_bar*(b_C[i]-b_C[i]/(1+g_bC[i]))   #trillion US$
          NLC_D_S2[i]<-IC_PRI_D_S2[i]-sh_GVA_S2*(1-beta_S2[i])*RP[i]+rep*LC_S2[i]/(1+g_baseline)-delta[i]*KC_PRI_S2[i]/(1+g_baseline)-sh_GVA_S2*p_C_bar*(b_C[i]-b_C[i]/(1+g_bC[i]))   #trillion US$
          NLC_D_S3[i]<-IC_PRI_D_S3[i]-sh_GVA_S3*(1-beta_S3[i])*RP[i]+rep*LC_S3[i]/(1+g_baseline)-delta[i]*KC_PRI_S3[i]/(1+g_baseline)-sh_GVA_S3*p_C_bar*(b_C[i]-b_C[i]/(1+g_bC[i]))   #trillion US$
          NLC_D_S4[i]<-IC_PRI_D_S4[i]-sh_GVA_S4*(1-beta_S4[i])*RP[i]+rep*LC_S4[i]/(1+g_baseline)-delta[i]*KC_PRI_S4[i]/(1+g_baseline)-sh_GVA_S4*p_C_bar*(b_C[i]-b_C[i]/(1+g_bC[i]))   #trillion US$
          IG_PRI_S1[i]<-sh_GREEN_S1*IG_PRI[i] #trillion US$
          IG_PRI_S2[i]<-sh_GREEN_S2*IG_PRI[i] #trillion US$
          IG_PRI_S3[i]<-sh_GREEN_S3*IG_PRI[i] #trillion US$
          IG_PRI_S4[i]<-sh_GREEN_S4*IG_PRI[i] #trillion US$
          IC_PRI_S1[i]<-I_PRI_S1[i]-IG_PRI_S1[i] #trillion US$
          IC_PRI_S2[i]<-I_PRI_S2[i]-IG_PRI_S2[i] #trillion US$
          IC_PRI_S3[i]<-I_PRI_S3[i]-IG_PRI_S3[i] #trillion US$
          IC_PRI_S4[i]<-I_PRI_S4[i]-IG_PRI_S4[i] #trillion US$
          IG_PRI[i]<- IG_PRI_initial #trillion US$
          IC_PRI[i]<-I_PRI[i]-IG_PRI[i] #trillion US$
          L[i]<-(0.914-B[i]/Y[i])*Y[i] #0.914 is the credit to non-financial corporations/GDP ratio; trillion US$ 
          kappa[i]<- IG_PRI[i]/I_PRI[i]# stabilises KG(PRI)/K(PRI)
          KG_PRI[i]<-kappa[i]*K_PRI[i] #trillion US$; stabilizes KG(PRI)/K(PRI)
          KG_PRI_S1[i]<- sh_GVA_S1*KG_PRI[i]*(1-beta_is_endogenous_different)+sh_GREEN_S1*KG_PRI[i]* beta_is_endogenous_different #trillion US$
          KG_PRI_S2[i]<- sh_GVA_S2*KG_PRI[i]*(1-beta_is_endogenous_different)+sh_GREEN_S2*KG_PRI[i]* beta_is_endogenous_different #trillion US$
          KG_PRI_S3[i]<- sh_GVA_S3*KG_PRI[i]*(1-beta_is_endogenous_different)+sh_GREEN_S3*KG_PRI[i]* beta_is_endogenous_different #trillion US$
          KG_PRI_S4[i]<- sh_GVA_S4*KG_PRI[i]*(1-beta_is_endogenous_different)+sh_GREEN_S4*KG_PRI[i]* beta_is_endogenous_different #trillion US$
          KC_PRI_S1[i]<-sh_GVA_S1*KC_PRI[i] #trillion US$
          KC_PRI_S2[i]<-sh_GVA_S2*KC_PRI[i] #trillion US$
          KC_PRI_S3[i]<-sh_GVA_S3*KC_PRI[i] #trillion US$
          KC_PRI_S4[i]<-sh_GVA_S4*KC_PRI[i] #trillion US$
          KC_PRI[i]<-K_PRI[i]-KG_PRI[i] #trillion US$
          
          KGE_PRI_S1[i]<-gamma_E1*KG_PRI_S1[i] #trillion US$
          KGE_PRI_S2[i]<-gamma_E2*KG_PRI_S2[i] #trillion US$
          KGE_PRI_S3[i]<-gamma_E3*KG_PRI_S3[i] #trillion US$
          KGE_PRI_S4[i]<-gamma_E4*KG_PRI_S4[i] #trillion US$
          KGNE_PRI_S1[i]<-(1-gamma_E1)*KG_PRI_S1[i]#trillion US$
          KGNE_PRI_S2[i]<-(1-gamma_E2)*KG_PRI_S2[i]#trillion US$
          KGNE_PRI_S3[i]<-(1-gamma_E3)*KG_PRI_S3[i]#trillion US$
          KGNE_PRI_S4[i]<-(1-gamma_E4)*KG_PRI_S4[i]#trillion US$
          KSEQ_PRI_S1[i]<-gamma_SEQ1*KGE_PRI_S1[i] #trillion US$
          KSEQ_PRI_S2[i]<-gamma_SEQ2*KGE_PRI_S2[i] #trillion US$
          KSEQ[i]<- KSEQ_PRI_S1[i]+KSEQ_PRI_S2[i] #trillion US$
          KCE_PRI_S1[i]<-gamma_E1*KC_PRI_S1[i] #trillion US$
          KCE_PRI_S2[i]<-gamma_E2*KC_PRI_S2[i] #trillion US$
          KCE_PRI_S3[i]<-gamma_E3*KC_PRI_S3[i] #trillion US$
          KCE_PRI_S4[i]<-gamma_E4*KC_PRI_S4[i] #trillion US$
          KCNE_PRI_S1[i]<-(1-gamma_E1)*KC_PRI_S1[i] #trillion US$
          KCNE_PRI_S2[i]<-(1-gamma_E2)*KC_PRI_S2[i] #trillion US$
          KCNE_PRI_S3[i]<-(1-gamma_E3)*KC_PRI_S3[i] #trillion US$
          KCNE_PRI_S4[i]<-(1-gamma_E4)*KC_PRI_S4[i] #trillion US$
          KGE[i]<- KGE_PRI_S1[i]+KGE_PRI_S2[i]+KGE_PRI_S3[i]+KGE_PRI_S4[i]+gamma_E*KG_GOV[i] #trillion US$
          KGNE[i]<- KGNE_PRI_S1[i]+KGNE_PRI_S2[i]+KGNE_PRI_S3[i]+KGNE_PRI_S4[i]+ (1-gamma_E)*KG_GOV[i] #trillion US$
          KCE[i]<- KCE_PRI_S1[i]+KCE_PRI_S2[i]+KCE_PRI_S3[i]+KCE_PRI_S4[i]+gamma_E*KC_GOV[i] #trillion US$
          KCNE[i]<- KCNE_PRI_S1[i]+KCNE_PRI_S2[i]+KCNE_PRI_S3[i]+KCNE_PRI_S4[i]+ (1-gamma_E)*KC_GOV[i] #trillion US$
          
          K_PRI[i]<- prop*K[i] #trillion US$ 
          delta[i]<-delta_0+(1-delta_0)*(1-ad_K)*D_TF[i]
          v[i]<-Y[i]/(K_PRI[i]*u[i])
          g_lambda[i]<-sigma_0[i]+sigma_1+sigma_2*g_Y[i]
          sigma_0[i]<-sigma_0_initial
          lambda[i]<-(Y[i]/N[i])/h[i] #trillion US$/(millions of empoyees*annual hours worked per employee)
          wage[i]<-s_W*lambda[i]*h[i] #trillion US$/millions of employees
          N[i]<-re[i]*LF[i] #billion people
          ur[i]<- 0.054 
          b_C[i]<-B_C[i]/p_C[i]#trillion US$
          b_G[i]<-B_G[i]/p_G[i]#trillion US$
          B[i]<-B_initial #trillion US$
          B_G[i]<-B_G_initial #trillion US$
          B_C[i]<-B_C_initial  #trillion US$
          p_C[i]<- p_C_bar #index
          p_G[i]<- p_G_bar #index
          
          yield_C[i]<-0.05
          yield_G[i]<-0.05
          DL[i]<-def[i]*(L[i]/(1+g_baseline)) #trillion US$
          def[i]<- 0.037 
          
          illiq[i]<-((int_C_S1[i]+rep)*LC_S1[i]/(1+g_baseline)+(int_C_S2[i]+rep)*LC_S2[i]/(1+g_baseline)+(int_C_S3[i]+rep)*LC_S3[i]/(1+g_baseline)+(int_C_S4[i]+rep)*LC_S4[i]/(1+g_baseline)+(int_G[i]+rep)*LG[i]/(1+g_baseline)+coupon_C[i]*b_C[i]/(1+g_bC[i])+coupon_G[i]*b_G[i]/(1+g_bG[i])+wage[i]*N[i]+TAX_F[i]+TAX_C[i]-SUB[i]+delta[i]*K_PRI[i]/(1+g_baseline))/(Y[i]+(1-CR_C_S1[i])*NLC_D_S1[i]+(1-CR_C_S2[i])*NLC_D_S2[i]+(1-CR_C_S3[i])*NLC_D_S3[i]+(1-CR_C_S4[i])*NLC_D_S4[i]+ (1-CR_G[i])*(NLG_D_S1[i]+ NLG_D_S2[i]+ NLG_D_S3[i]+NLG_D_S4[i])+(b_C[i]-b_C[i]/(1+g_bC[i]))*p_C_bar+ (b_G[i]-b_G[i]/(1+g_bG[i]))*p_G_bar) 
          
          dsr[i]<-((int_C_S1[i]+rep)*LC_S1[i]/(1+g_baseline)+(int_C_S2[i]+rep)*LC_S2[i]/(1+g_baseline)+(int_C_S3[i]+rep)*LC_S3[i]/(1+g_baseline)+(int_C_S4[i]+rep)*LC_S4[i]/(1+g_baseline)+(int_G[i]+rep)*(LG[i]/(1+g_baseline))+coupon_C[i]*(b_C[i]/(1+g_bC[i]))+coupon_G[i]*(b_G[i]/(1+g_bG[i])))/(TP[i]+int_C_S1[i]*LC_S1[i]/(1+g_baseline)+int_C_S2[i]*LC_S2[i]/(1+g_baseline)+int_C_S3[i]*LC_S3[i]/(1+g_baseline)+int_C_S4[i]*LC_S4[i]/(1+g_baseline)+(int_G[i])*(LG[i]/(1+g_baseline))+coupon_C[i]*(b_C[i]/(1+g_bC[i]))+coupon_G[i]*(b_G[i]/(1+g_bG[i])))
          
          #2.2.3 Households 
          Y_HG[i]<-wage[i]*N[i]+DP[i]+BP_D[i]+int_D*(D[i]/(1+g_baseline))+int_S*(SEC_H[i]/(1+g_baseline))+coupon_C[i]*(b_CH[i]/(1+g_bC[i]))+coupon_G[i]*(b_GH[i]/(1+g_bG[i]))  #trillion US$
          Y_H[i]<-Y_HG[i]-TAX_H[i] #trillion US$
          CO_PRI_N[i]<-Y[i]-I_PRI[i]-CO_GOV[i]-I_GOV[i] #trillion US$
          CO_PRI[i]<-CO_PRI_N[i] #trillion US$
          
          V_HF[i]<-D[i]+p_C[i]*b_CH[i]+p_G[i]*b_GH[i]+SEC_H[i]   #trillion US$
          SEC_H[i]<-SEC[i]-SEC_CB[i]-SEC_B[i] #trillion US$
          B_CH[i]<-B_C[i]-B_CCB[i] #trillion US$
          B_GH[i]<-B_G[i]-B_GCB[i] #trillion US$
          D[i]<-70 #trillion US$ 
          b_CH[i]<-B_CH[i]/p_C[i] #trillion US$
          b_GH[i]<-B_GH[i]/p_G[i] #trillion US$
          DC[i]<-((DEM[i]-mu[i]*delta[i]*(K_PRI[i]+K_GOV[i])/(1+g_baseline))/(mu[i]*xi))*(1+g_baseline)  #trillion US$
          LF[i]<- 3.43 #billion people  
          POP[i]<- 7.63 #billion people 
          g_POP[i]<- 0.011 
          g_POP_stag[i]<-0 
          
          #2.2.4 Banks
          BP[i]<-int_C_S1[i]*LC_S1[i]/(1+g_baseline)+int_C_S2[i]*LC_S2[i]/(1+g_baseline) +int_C_S3[i]*LC_S3[i]/(1+g_baseline)+int_C_S4[i]*LC_S4[i]/(1+g_baseline)+int_G[i]*(LG_S1[i]/(1+g_baseline)+LG_S2[i]/(1+g_baseline)+LG_S3[i]/(1+g_baseline)+LG_S4[i]/(1+g_baseline))+int_S*(SEC_B[i]/(1+g_baseline))-int_D*(D[i]/(1+g_baseline))-int_A*(A[i]/(1+g_baseline)) #trillion US$
          HPM[i]<-h_1*D[i] #trillion US$
          CAP[i]<-(HPM[i]+LG[i]+LC[i]+SEC_B[i])/lev_B[i]
          CAP_before[i]<-CAP[i]
          BP_U[i]<-s_B*(BP[i]/(1+g_baseline)) #trillion US$
          BP_D[i]<-BP[i]-BP_U[i] #trillion US$
          SEC_B[i]<-0.15*SEC[i] #trillion US$
          A[i]<-HPM[i]+LG[i]+LC[i]+ SEC_B[i]-D[i]-CAP[i]#trillion US$
          A_N[i]<-A[i] #trillion US$
          random2[i]<-0
          LG[i]<- kappa[i]*L[i] #trillion US$ 
          LG_S1[i]<-sh_GVA_S1*LG[i] #trillion US$
          LG_S2[i]<-sh_GVA_S2*LG[i] #trillion US$
          LG_S3[i]<-sh_GVA_S3*LG[i] #trillion US$
          LG_S4[i]<-sh_GVA_S4*LG[i] #trillion US$
          LC[i]<-L[i]-LG[i] #trillion US$ 
          LC_S1[i]<-sh_GVA_S1*LC[i] #trillion US$
          LC_S2[i]<-sh_GVA_S2*LC[i] #trillion US$
          LC_S3[i]<-sh_GVA_S3*LC[i] #trillion US$
          LC_S4[i]<-sh_GVA_S4*LC[i] #trillion US$
          lev_B[i]<-1/0.1067 
          BAILOUT[i]<-0 #trillion US$
          BAILOUT_levB[i]<-0 #trillion US$
          BAILOUT_CAR[i]<-0 #trillion US$
          
          r_0<-(CR_max/CR_initial)-1 #(Ci)
          r_2_dCR_ddsr_max<- 1.18362275#econometrically estimated coefficient
          r_2<-(r_2_dCR_ddsr_max *((1+r_0)^2))/(CR_max*r_0) #(Biii)
          r_3<- -(r_3_dCR_dCAR_max *((1+r_0)^2))/(CR_max*r_0) #(Biii)
          r_1<-r_2*dsr[i]-r_3*(CAR[i]-CAR_min)#(Ci)
          sh_NLG[i]<- 0.03401567#NLG_D[i]/(NLG_D[i]+NLC_D_S1[i]+NLC_D_S2[i] +NLC_D_S3[i]+NLC_D_S4[i])
          sh_NLC_S1[i]<- 0.05092839#NLC_D_S1[i]/(NLG_D_S1[i]+NLG_D_S2[i]+NLG_D_S3[i]+NLG_D_S4[i]+NLC_D_S1[i]+NLC_D_S2[i] +NLC_D_S3[i]+NLC_D_S4[i])
          sh_NLC_S2[i]<- 0.22115655#NLC_D_S2[i]/(NLG_D_S1[i]+NLG_D_S2[i]+NLG_D_S3[i]+NLG_D_S4[i]+NLC_D_S1[i]+NLC_D_S2[i] +NLC_D_S3[i]+NLC_D_S4[i])
          sh_NLC_S3[i]<- 0.08657827#NLC_D_S3[i]/(NLG_D_S1[i]+NLG_D_S2[i]+NLG_D_S3[i]+NLG_D_S4[i]+NLC_D_S1[i]+NLC_D_S2[i] +NLC_D_S3[i]+NLC_D_S4[i])
          sh_NLC_S4[i]<-1-sh_NLG[i]-sh_NLC_S1[i]-sh_NLC_S2[i]-sh_NLC_S3[i]
          w_LT[i]<-sh_LG[i]*w_G[i]+sh_LC_S1[i]*w_C_S1[i]+sh_LC_S2[i]*w_C_S2[i]+sh_LC_S3[i]*w_C_S3[i]+sh_LC_S4[i]*w_C_S4[i]
          CR_G[i]<-(1+l_1*(w_G[i]-w_LT[i]))*CR[i]
          CR_C_S1[i]<-(1+l_1*(w_C_S1[i]-w_LT[i]))*CR[i]
          CR_C_S2[i]<-(1+l_1*(w_C_S2[i]-w_LT[i]))*CR[i]
          CR_C_S3[i]<-(1+l_1*(w_C_S3[i]-w_LT[i]))*CR[i]
          CR_C_S4[i]<-(CR[i]-sh_NLG[i]*CR_G[i]-sh_NLC_S1[i]*CR_C_S1[i]-sh_NLC_S2[i]*CR_C_S2[i]-sh_NLC_S3[i]*CR_C_S3[i])/sh_NLC_S4[i]
          
          CR_C[i]<-0.2*cr_rationing_dummy
          CR[i]<-CR_max/(1+r_0*exp(random2[i]*random_dummy+r_1-r_2*dsr[i]+r_3*(CAR[i]-CAR_min)))*cr_rationing_dummy  
          CAR[i]<-CAP[i]/(w_G[i]*LG[i]+w_C_S1[i]*LC_S1[i]+w_C_S2[i]*LC_S2[i]+w_C_S3[i]*LC_S3[i] +w_C_S4[i]*LC_S4[i]+w_S*SEC_B[i]+w_H*HPM[i]) 
          lev_B_before[i]<-lev_B[i]
          CAR_before[i]<-CAR[i]
          
          #2.2.5 Government
          GNS[i]<- -CO_GOV[i]+TAX[i]-int_S*SEC[i]/(1+g_baseline)+CBP[i]-SUB[i]-delta[i]*K_GOV[i]/(1+g_baseline) #trillion US$
          SEC[i]<-0.815*Y[i] #trillion US$, 81.5 general government debt to GDP 
          IG_GOV[i]<-IG_GOV_initial#trillion US$
          IC_GOV[i]<-I_GOV[i]-IG_GOV[i] #trillion US$ 
          I_GOV[i]<- (1-prop)*0.24*Y[i] #trillion US$ 
          K_GOV[i]<-(1-prop)*K[i] #trillion US$
          K[i]<-K_Y_ratio*Y[i]#trillion US$ 
          KG[i]<-KG_PRI[i]+KG_GOV[i] #trillion US$
          KC[i]<-KC_PRI[i]+KC_GOV[i] #trillion US$
          CO_GOV[i]<-gov_C*(Y[i]/(1+g_baseline)) #trillion US$
          SUB[i]<- TAX_C[i]#trillion US$ 
          gov_SUB[i]<- min(SUB[i]/(E_NF[i]*ucr[i]),0.999) #trillion US$
          TAX_H[i]<-tau_H*(Y_HG[i]/(1+g_baseline))   #trillion US$
          TAX_F[i]<-tau_F*(TP_G[i]/(1+g_baseline)) #trillion US$
          TAX_C[i]<- 0.044# trillion US$ 
          tau_C[i]<- TAX_C[i]/EMIS_IN[i]
          tau_C_baseline[i]<-Data[1,c(1)]/1000
          tau_C_policy_Medium[i]<- Data[1,c(2)]/1000
          tau_C_policy_High[i]<- Data[1,c(3)]/1000
          TAX[i]<-TAX_H[i]+TAX_F[i]+TAX_C[i]  #trillion US$
          
          #2.2.6 Central bank
          CBP[i]<-coupon_C[i]*(b_CCB[i]/(1+g_bC[i])) +coupon_G[i]*(b_GCB[i]/(1+g_bG[i]))+int_A*(A[i]/(1+g_baseline))+ int_S*(SEC_CB[i]/(1+g_baseline))    #trillion US$
          B_GCB[i]<-B_GCB_initial #trillion US$
          B_CCB[i]<-B_CCB_initial #trillion US$
          b_CCB[i]<-B_CCB[i]/p_C[i] #trillion US$
          b_GCB[i]<-B_GCB[i]/p_G[i] #trillion US$
          SEC_CB[i]<-HPM[i]+V_CB[i]-p_C[i]*b_CCB[i]-p_G[i]*b_GCB[i]-A[i] #trillion US$
          
          #Auxiliary variables
          g_Y[i]<- 0.0304#g_baseline
          IG_E[i]<-0.5 #trillion US$
          IG_E_Y[i]<- IG_E[i]/Y[i]*100
          IG_E_Y_cum[i]<- IG_E[i]/Y[i]*100
          g_K[i]<-g_baseline
          ur_per[i]<-ur[i]*100
          g_Y_per[i]<-g_Y[i]*100
          theta_per[i]<-theta[i]*100
          greenI[i]<-IG_PRI[i]/I_PRI[i]
          r_total[i]<-(TP[i]/K_PRI[i]) *100
          I_Y_ratio[i]<-I_PRI[i]/Y[i]
          C_Y_ratio[i]<-CO_PRI[i]/Y[i] 
          C_K_ratio[i]<-CO_PRI[i]/K_PRI[i]
          Y_K_ratio[i]<-Y[i]/K_PRI[i]
          Y_H_Y_ratio[i]<-Y_H[i]/Y[i]
          Y_H_K_ratio[i]<-Y_H[i]/K_PRI[i]
          A_K_ratio[i]<-A[i]/K_PRI[i]
          h[i]<-hours
          h_optimal[i]<-h[i]
          Wbill[i]<-wage[i]*N[i]
          Interest[i]<- int_C_S1[i]*LC_S1[i]/(1+g_baseline)+int_C_S2[i]*LC_S2[i]/(1+g_baseline)+int_C_S3[i]*LC_S3[i]/(1+g_baseline)+int_C_S4[i]*LC_S4[i]/(1+g_baseline)+int_G[i]*LG[i]/(1+g_baseline)
          Depreciation[i]<-delta[i]*K_PRI[i]
          SEC_Y[i]<-SEC[i]/Y[i]
          SEC_Y_per[i]<-SEC[i]/Y[i]*100
          fiscal_balance_per[i]<-fiscal_balance[i]*100
          L_K[i]<-L[i]/K_PRI[i]
          I[i]<-I_PRI[i]+I_GOV[i]
          sh_L[i]<-L[i]/(L[i]+B[i])
          B_C_issue[i]<-x_1[i]*IC_PRI_D[i] #2 is the actual bond issuance
          B_G_issue[i]<-x_2[i]*IG_PRI_D[i] #0.2 is the actual bond issuance
          IG_cum[i]<-IG_PRI[i]+IG_GOV[i]
          IC_cum[i]<-IC_PRI[i]+IC_GOV[i]
          D_K_ratio[i]<-D[i]/K_PRI[i]
          HPM_K_ratio[i]<-HPM[i]/K_PRI[i]
          fiscal_balance[i]<-(TAX[i]-CO_GOV[i]-IC_GOV[i]-IG_GOV[i]-SUB[i]-BAILOUT[i]-int_S*SEC[i]/(1+g_baseline))/Y[i]
          I_K_ratio[i]<-I_PRI[i]/(K_PRI[i]/(1+g_baseline))
          ID_K_ratio[i]<-I_PRI_D[i]/(K_PRI[i]/(1+g_baseline))
          Y_POP_ratio[i]<-(Y[i]/POP[i])
          W_POP_ratio[i]<-W[i]/POP[i]
          E_N_ratio[i]<-E[i]/N[i]
          Y_E_ratio[i]<-Y[i]/E[i]
          lambda_perworker[i]<-lambda[i]*h[i]
          lf[i]<-LF[i]/POP[i]
          haz_flow[i]<-haz*W[i]
          Y_HD[i]<-wage[i]*N[i]-TAX_H[i]+DP[i]+BP_D[i]+int_D*(D[i]/(1+g_baseline))+int_S*(SEC_H[i]/(1+g_baseline))+coupon_C[i]*(b_CH[i]/(1+g_bC[i]))+coupon_G[i]*(b_GH[i]/(1+g_bG[i]))- xi*DC[i]/(1+g_baseline)   #trillion US$
          V_H[i]<-DC[i]+D[i]+p_C_bar*b_CH[i]+p_G_bar*b_GH[i]+SEC_H[i]   #trillion US$
          KG_GOV_Y_ratio[i]<- KG_GOV[i]/Y[i]
          KC_GOV_Y_ratio[i]<-KC_GOV[i]/Y[i] 
          IG_GOV_Y_ratio[i]<- IG_GOV[i]/Y[i]
          IC_GOV_Y_ratio[i]<-IC_GOV[i]/Y[i] 
          KG_Y_ratio[i]<- KG_PRI[i]/Y[i]
          KC_Y_ratio[i]<-KC_PRI[i]/Y[i] 
          KGE_KCE[i]<-KGE[i]/KCE[i]
          KGNE_KCNE[i]<-KGNE[i]/KCNE[i]
          KSEQ_KCE[i]<-KSEQ[i]/(KCE_PRI_S1[i]+KCE_PRI_S2[i])
          KG_K_ratio[i]<- KG_PRI[i]/K_PRI[i]
          LG_L_ratio_pseudo[i]<- (beta[i]*(g_baseline*K_PRI[i]/(1+g_baseline)-RP[i]-x_2[i]*I_PRI_D[i])/L[i])/(def[i]/(1+g_baseline)+g_baseline/(1+g_baseline))
          LG_L_ratio[i]<- LG[i]/L[i]
          g_Y_cum[i]<-g_Y[i]
          def_per[i]<-def[i]*100
          CAR_per[i]<-CAR[i]*100
          IG_PRI_D[i]<-IG_PRI_D_S1[i]+ IG_PRI_D_S2[i]+ IG_PRI_D_S3[i]+ IG_PRI_D_S4[i] #trillion US$
          IC_PRI_D[i]<-I_PRI_D[i]-IG_PRI_D[i]  #trillion US$
          beta[i]<-beta_initial
          I_PRI_S1[i]<-sh_GVA_S1*I_PRI[i] 
          I_PRI_S2[i]<-sh_GVA_S2*I_PRI[i] 
          I_PRI_S3[i]<-sh_GVA_S3*I_PRI[i] 
          I_PRI_S4[i]<-sh_GVA_S4*I_PRI[i] 
          TAX_C_change[i]<- TAX_C[i]
          SUB_change[i]<- SUB[i]
          IG_GOV_change[i]<- IG_GOV[i]
          omega_ratio[i]<-omega[i]/omega[1]
          mu_ratio[i]<-mu[i]/mu[1]
          rho_ratio[i]<-rho[i]/rho[1]
          epsilon_ratio[i]<-epsilon[i]/epsilon[1]
          yield_C_per[i]<- yield_C[i]*100
          yield_G_per[i]<- yield_G[i]*100
          private_balance[i]<-Y[i]-CO_PRI[i]-I_PRI[i]-TAX[i]+int_S*SEC_H[i]/(1+g_baseline)+int_S*SEC_B[i]/(1+g_baseline)+SUB[i]+BAILOUT[i]-coupon_C[i]*b_CCB[i]/(1+g_bC[i])-coupon_G[i]*b_GCB[i]/(1+g_bG[i])-int_A*A[i]/(1+g_baseline)
          firms_balance[i]<-RP[i]+s_F*delta[i]*K_PRI[i]/(1+g_baseline)-I_PRI[i]
          households_balance[i]<- Y_H[i]+(1-s_F)*delta[i]*K_PRI[i]/(1+g_baseline)-CO_PRI[i]
          banks_balance[i]<- BP_U[i]+BAILOUT[i]
          fiscal_balance_help[i]<-TAX[i]-CO_GOV[i]-IC_GOV[i]-IG_GOV[i]-SUB[i]-BAILOUT[i]-int_S*SEC[i]/(1+g_baseline)+CBP[i]
          private_balance_per[i]<- private_balance[i]/Y[i]*100
          firms_balance_per[i]<- firms_balance[i]/Y[i]*100
          households_balance_per[i]<- households_balance[i]/Y[i]*100
          banks_balance_per[i]<- banks_balance[i]/Y[i]*100
          fiscal_balance_help_per[i]<-fiscal_balance_help[i]/Y[i]*100
          tau_C_baseline_per[i]<-Data[1,c(1)]
          tau_C_policy_Medium_per[i]<- Data[1,c(2)]
          tau_C_policy_High_per[i]<- Data[1,c(3)]
          
          r_2_dCR_ddsr_max_pseudo<- -(-1.326*0.18)*(L[i]/(NLC_D_S1[i]+ NLC_D_S2[i]+ NLC_D_S3[i]+ NLC_D_S4[i]+ NLG_D_S1[i]+ NLG_D_S2[i]+ NLG_D_S3[i]+ NLG_D_S4[i]))# -1.752 is taken from the econometric estimates; 0.18 is the sum of the interest rate and repayment ratio
          
          r_3_dCR_dCAR_max_pseudo<- -(0.527*(0.0836/0.1553))*(L[i]/(NLC_D_S1[i]+ NLC_D_S2[i]+ NLC_D_S3[i]+ NLC_D_S4[i]+ NLG_D_S1[i]+ NLG_D_S2[i]+ NLG_D_S3[i]+ NLG_D_S4[i]))# 0.630 is taken from the econometric estimates; 0.1553 is taken from IMF and is the average value of CAR for the period 2005-2018; 0.0836 is taken from OECD and is the average value of Capital/assets for the same period
          
          r_3_dCR_dCAR_max_pseudo_case_I<- -((0.331)*(0.0836/0.1553))*(L[i]/(NLC_D_S1[i]+ NLC_D_S2[i]+ NLC_D_S3[i]+ NLC_D_S4[i]+ NLG_D_S1[i]+ NLG_D_S2[i]+ NLG_D_S3[i]+ NLG_D_S4[i]))# 0.421 is taken from the econometric estimates; 0.1553 is taken from IMF and is the average value of CAR for the period 2005-2018; 0.0836 is taken from OECD and is the average value of Capital/assets for the same period
          
          r_3_dCR_dCAR_max_pseudo_case_II<- -(0.8*(0.0836/0.1553))*(L[i]/(NLC_D_S1[i]+ NLC_D_S2[i]+ NLC_D_S3[i]+ NLC_D_S4[i]+ NLG_D_S1[i]+ NLG_D_S2[i]+ NLG_D_S3[i]+ NLG_D_S4[i]))# 1.000 is taken from the econometric estimates; 0.1553 is taken from IMF and is the average value of CAR for the period 2005-2018; 0.0836 is taken from OECD and is the average value of Capital/assets for the same period
          
          spr_1_pseudo<- (-0.058*(0.0836/0.1553))# -0.058 is taken from the econometric estimates; 0.1553 is taken from IMF and is the average value of CAR for the period 2005-2018; 0.0836 is taken from OECD and is the average value of Capital/assets for the same period
          
          spr_1_pseudo_case_I<- (-0.031*(0.0836/0.1553))# -0.031 is taken from the econometric estimates; 0.1553 is taken from IMF and is the average value of CAR for the period 2005-2018; 0.0836 is taken from OECD and is the average value of Capital/assets for the same period
          
          spr_1_pseudo_case_II<- (-0.080*(0.0836/0.1553))# -0.080 is taken from the econometric estimates; 0.1553 is taken from IMF and is the average value of CAR for the period 2005-2018; 0.0836 is taken from OECD and is the average value of Capital/assets for the same period
          
          spr_2_pseudo <- (0.135*0.18)# 0.135 is taken from the econometric estimates; 0.18 is the sum of the interest rate and repayment ratio
          
          EMIS_IN_Y[i]<- EMIS_IN[i]/Y[i]
          EMIS_IN_Y_SSP360[i]<- Data[1,c(4)]
        }
      }
      
      ####################
      #3. MODEL EQUATIONS
      ####################
      
      else {
        
        if (i==3) {covid<-1} else {covid<-0} 
        random1[i]<-0.9*random1[i-1]+rnorm(1,0, 0.1)
        random2[i]<-0.9*random2[i-1]+rnorm(1,0, 0.003)
        
        for (iterations in 1:15){
          
          #2.1.1 Matter, recycling and waste 
          MY[i]<- mu[i]*(Y[i]-CO_GOV[i]) #Eq. (1)
          M[i]<-MY[i]-REC[i] #Eq. (2)
          REC[i]<-rho[i]*DEM[i] #Eq. (3)
          DEM[i]<-mu[i]*(delta[i]*(K_PRI[i-1]+ K_GOV[i-1])+xi*DC[i-1]) #Eq. (4)
          SES[i]<-SES[i-1]+MY[i]-DEM[i] #Eq. (5)
          W[i]<-M[i]+CEN[i]+O2[i]-EMIS_IN[i]-(SES[i]- SES[i-1])  #Eq. (6)
          CEN[i]<-(1/car)*EMIS_IN[i] #Eq. (7)
          O2[i]<-EMIS_IN[i]-CEN[i] #Eq. (8)
          HW_CUM[i]<-HW_CUM[i-1]+haz*W[i] #Eq. (9)
          hazratio[i]<-HW_CUM[i]/POP[i] #Eq. (10)
          REV_M[i]<-REV_M[i-1]+CON_M[i]-M[i] #Eq. (11)
          CON_M[i]<-con_M*RES_M[i-1] #Eq. (12)
          RES_M[i]<-RES_M[i-1]-CON_M[i] #Eq. (13)
          dep_M[i]<-M[i]/REV_M[i-1] #Eq. (14)
          
          #2.1.2 Energy
          E[i]<-epsilon[i]*Y[i] #Eq. (15)
          E_NF[i]<-theta[i]*E[i] #Eq. (16)
          E_F[i]<-E[i]-E_NF[i] #Eq. (17)
          ED[i]<-E_F[i]+E_NF[i] #Eq. (18)
          REV_E[i]<-REV_E[i-1]+CON_E[i]-E_F[i] #Eq. (19)
          CON_E[i]<-con_E*RES_E[i-1] #Eq. (20)
          RES_E[i]<-RES_E[i-1]-CON_E[i] #Eq. (21)
          dep_E[i]<-E_F[i]/REV_E[i-1] #Eq. (22)
          
          #2.1.3 Emissions and climate change
          EMIS_IN[i]<-omega[i]*(1-seq[i])*E_F[i] #Eq. (23)
          g_EMIS_L[i]<-g_EMIS_L[i-1]*(1-zeta_9) #Eq. (24)
          EMIS_L[i]<- EMIS_L[i-1]*(1-g_EMIS_L[i]) #Eq. (25)
          EMIS[i]<-EMIS_IN[i]+EMIS_L[i] #Eq. (26)
          CO2_CUM[i]<-CO2_CUM[i-1]+EMIS[i] #Eq. (27)
          T_AT[i]<-T_AT[i-1]+t_1*(t_2*phi*CO2_CUM[i-1]-T_AT[i-1]) #Eq. (28)
          
          #2.1.4 Ecological efficiency and technology
          omega[i]<-omega[i-1]*(1+g_omega[i]) #Eq. (29)
          g_omega[i]<-g_omega[i-1]*(1-zeta_1) #Eq. (30)
          mu[i]<-mu_max-((mu_max-mu_min)/(1+pi_1*exp(-pi_2*((KGNE[i-1])/(KCNE[i-1]))))) #Eq. (31)
          rho[i]<-rho_max/(1+pi_3*exp(-pi_4*((KGNE[i-1])/(KCNE[i-1])))) #Eq. (32)
          epsilon[i]<-epsilon_max-((epsilon_max-epsilon_min)/(1+pi_5*exp(-pi_6*((KGE[i-1])/(KCE[i-1]))))) #Eq. (33)
          theta[i]<-1/(1+pi_7*exp(-pi_8*((KGE[i-1])/(KCE[i-1])))) #Eq. (34)
          seq[i]<-1/(1+pi_9*exp(-pi_10*((KSEQ[i-1])/(KCE_PRI_S1[i-1]+KCE_PRI_S2[i-1])))) #Eq. (35) 
          
          #2.2.1 Output determination
          Y_star_M[i]<-(REV_M[i-1]+REC[i])/mu[i] #Eq. (36)
          Y_star_E[i]<-REV_E[i-1]/((1-theta[i])*epsilon[i]) #Eq. (37)
          Y_star_K[i]<-v[i]*K_PRI[i] #Eq. (38)
          Y_star_N[i]<-lambda[i]*LF[i]*h[i] #Eq. (39)
          Y_star[i]<-min(min(Y_star_M[i], Y_star_E[i]),min(Y_star_N[i], Y_star_K[i])) #Eq. (40)
          Y[i]<-CO_PRI[i]+I_PRI[i]+CO_GOV[i]+I_GOV[i] #Eq. (41)
          Y_N[i]<-CO_PRI_N[i]+I_PRI[i]+CO_GOV[i]+I_GOV[i]
          um[i]<- (Y[i]-CO_GOV[i])/Y_star_M[i] #Eq. (42)
          ue[i]<-Y[i]/Y_star_E[i] #Eq. (43)
          u[i]<-Y[i]/Y_star_K[i] #Eq. (44)
          re[i]<-Y[i]/Y_star_N[i] #Eq. (45)
          D_T[i]<-(1-(1/(1+eta_1*T_AT[i]+eta_2*T_AT[i]^2+eta_3*T_AT[i]^6.754)))*damage_dummy #Eq. (46)
          D_TP[i]<-p*D_T[i] #Eq. (47)
          D_TF[i]<-1-(1-D_T[i])/(1-D_TP[i]) #Eq. (48)
          
          #2.2.2 Firms
          TP_G[i]<-Y[i]-wage[i]*N[i]-int_C_S1[i-1]*LC_S1[i-1]-int_C_S2[i-1]*LC_S2[i-1]-int_C_S3[i-1]*LC_S3[i-1]-int_C_S4[i-1]*LC_S4[i-1]-int_G[i-1]*(LG_S1[i-1]+LG_S2[i-1]+LG_S3[i-1]+LG_S4[i-1])-delta[i]*K_PRI[i-1]-coupon_C[i-1]*b_C[i-1]-coupon_G[i-1]*b_G[i-1] #Eq. (49)
          TP[i]<-TP_G[i]-TAX_F[i] -TAX_C[i]+SUB[i]    #Eq. (50)
          RP[i]<-s_F*TP[i-1]#Eq. (51)
          DP[i]<-TP[i]-RP[i] #Eq. (52)
          r[i]<-TP[i]/K_PRI[i] #Eq. (53)
          I_PRI_D[i]<-(((alpha_00)/(1+exp(random1[i]*random_dummy+alpha_01-alpha_1*u[i-1]-alpha_2*r[i-1]+(alpha_31*ur[i-1]^(-alpha_32))+(alpha_41*(1-ue[i-1])^(-alpha_42))+(alpha_51*(1-um[i-1])^(-alpha_52))))-delta[1])*K_PRI[i-1])*(1-D_T[i-1])*(1-stag_adjust*g_POP_stag[i]-covid*0.5) +delta[i]*K_PRI[i-1] #Eq. (54)
          I_PRI_D_S1[i]<-sh_GVA_S1*I_PRI_D[i] #Eq. (55)
          I_PRI_D_S2[i]<-sh_GVA_S2*I_PRI_D[i] #Eq. (55)
          I_PRI_D_S3[i]<-sh_GVA_S3*I_PRI_D[i] #Eq. (55)
          I_PRI_D_S4[i]<-sh_GVA_S4*I_PRI_D[i] #Eq. (55)
          IG_PRI_D_S1[i]<-beta_S1[i]*I_PRI_D_S1[i] #Eq. (56)
          IG_PRI_D_S2[i]<-beta_S2[i]*I_PRI_D_S2[i] #Eq. (56)
          IG_PRI_D_S3[i]<-beta_S3[i]*I_PRI_D_S3[i] #Eq. (56)
          IG_PRI_D_S4[i]<-beta_S4[i]*I_PRI_D_S4[i] #Eq. (56)
          
          beta_S1[i]<-min(beta[i-1]*(1-beta_is_endogenous_different)+(beta_0_S1[i]-beta_1*sh_EMIS_IN_S1* (tucr[i-1]-tucn[i-1])-beta_2*(sh_L[i-1]*(int_G[i-1]-int_C_S1[i-1])+(1-sh_L[i-1])*(yield_G[i-1]-yield_C[i-1])))* beta_is_endogenous_different,1) #Eq. (57)
          beta_S2[i]<- min(beta[i-1]* (1-beta_is_endogenous_different)+(beta_0_S2[i]-beta_1*sh_EMIS_IN_S2*( tucr[i-1]-tucn[i-1])-beta_2*(sh_L[i-1]*(int_G[i-1]-int_C_S2[i-1])+(1-sh_L[i-1])*(yield_G[i-1]-yield_C[i-1])))*beta_is_endogenous_different ,1)#Eq. (57)
          beta_S3[i]<-min(beta[i-1]*(1-beta_is_endogenous_different)+(beta_0_S3[i]-beta_1*sh_EMIS_IN_S3*( tucr[i-1]-tucn[i-1])-beta_2*(sh_L[i-1]*(int_G[i-1]-int_C_S3[i-1])+(1-sh_L[i-1])*(yield_G[i-1]-yield_C[i-1])))*beta_is_endogenous_different,1)#Eq. (57)
          beta_S4[i]<- min(beta[i-1]*(1-beta_is_endogenous_different)+(beta_0_S4[i]-beta_1*sh_EMIS_IN_S4*( tucr[i-1]-tucn[i-1])-beta_2*(sh_L[i-1]*(int_G[i-1]-int_C_S4[i-1])+(1-sh_L[i-1])*(yield_G[i-1]-yield_C[i-1])))*beta_is_endogenous_different ,1)#Eq. (57)
          
          if (i==6) {beta_0_S1[i]<-beta_0_S1[i-1]*(1+g_beta_0[i])*(1+beta_dummy_0*1) } else {beta_0_S1[i]<-beta_0_S1[i-1]*(1+g_beta_0[i])} #Eq. (58)
          if (i==6) {beta_0_S2[i]<-beta_0_S2[i-1]*(1+g_beta_0[i])*(1+beta_dummy_0*4) } else {beta_0_S2[i]<-beta_0_S2[i-1]*(1+g_beta_0[i])} #Eq. (58)
          if (i==6) {beta_0_S3[i]<-beta_0_S3[i-1]*(1+g_beta_0[i])*(1+beta_dummy_0*6) } else {beta_0_S3[i]<-beta_0_S3[i-1]*(1+g_beta_0[i])} #Eq. (58)
          if (i==6) {beta_0_S4[i]<-beta_0_S4[i-1]*(1+g_beta_0[i])*(1+beta_dummy_0*7) } else {beta_0_S4[i]<-beta_0_S4[i-1]*(1+g_beta_0[i])} #Eq. (58)
          g_beta_0[i]<-g_beta_0[i-1]*(1-zeta_2) #Eq. (59)
          
          tucr[i]<- ucr[i]*(1-gov_SUB[i]) #Eq. (60)
          tucn[i]<- ucn[i]+ tau_C[i]*omega[i]*(1-seq[i]) #Eq. (61)
          ucn[i]<- ucn[i-1]*(1+g_ucn[i]) #Eq. (62)
          g_ucn[i]<- g_ucn[i-1]*(1-zeta_8) #Eq. (63)
          endog[i]<-(1-ucr_dummy)+ucr_dummy*(1-theta[i])/(1-theta[i-1]) 
          ucr[i]<- max(ucr[i-1]*(1-g_ucr[i])*endog[i],0.4*ucr[1]) #Eq. (64)
          g_ucr[i]<- g_ucr[i-1]*(1-zeta_7) #Eq. (65)
          IC_PRI_D_S1[i]<-I_PRI_D_S1[i]-IG_PRI_D_S1[i]  #Eq. (66)
          IC_PRI_D_S2[i]<-I_PRI_D_S2[i]-IG_PRI_D_S2[i]#Eq. (66)
          IC_PRI_D_S3[i]<-I_PRI_D_S3[i]-IG_PRI_D_S3[i]  #Eq. (66)
          IC_PRI_D_S4[i]<-I_PRI_D_S4[i]-IG_PRI_D_S4[i]#Eq. (66)
          
          NLG_DN_S1[i]<-IG_PRI_D_S1[i]-sh_GVA_S1*beta_S1[i]*RP[i]+rep*LG_S1[i-1]-delta[i]*KG_PRI_S1[i-1]-sh_GVA_S1*p_G_bar*(b_G[i]-b_G[i-1]) #Eq. (67)
          NLG_DN_S2[i]<- IG_PRI_D_S2[i]-sh_GVA_S2*beta_S2[i]*RP[i]+rep*LG_S2[i-1]-delta[i]*KG_PRI_S2[i-1]-sh_GVA_S2*p_G_bar*(b_G[i]-b_G[i-1]) #Eq. (67)
          NLG_DN_S3[i]<- IG_PRI_D_S3[i]-sh_GVA_S3*beta_S3[i]*RP[i]+rep*LG_S3[i-1]-delta[i]*KG_PRI_S3[i-1]-sh_GVA_S3*p_G_bar*(b_G[i]-b_G[i-1]) #Eq. (67)
          NLG_DN_S4[i]<- IG_PRI_D_S4[i]-sh_GVA_S4*beta_S4[i]*RP[i]+rep*LG_S4[i-1]-delta[i]*KG_PRI_S4[i-1]-sh_GVA_S4*p_G_bar*(b_G[i]-b_G[i-1]) #Eq. (67)
          if (NLC_DN_S1[i]>=0) {NLG_D_S1[i]<-NLG_DN_S1[i]}  else {NLG_D_S1[i]<-NLG_DN_S1[i]+NLC_DN_S1[i]}  
          if (NLC_DN_S2[i]>=0) {NLG_D_S2[i]<-NLG_DN_S2[i]}  else {NLG_D_S2[i]<-NLG_DN_S2[i]+NLC_DN_S2[i]}  
          if (NLC_DN_S3[i]>=0) {NLG_D_S3[i]<-NLG_DN_S3[i]}  else {NLG_D_S3[i]<-NLG_DN_S3[i]+NLC_DN_S3[i]}  
          if (NLC_DN_S4[i]>=0) {NLG_D_S4[i]<-NLG_DN_S4[i]}  else {NLG_D_S4[i]<-NLG_DN_S4[i]+NLC_DN_S4[i]}  
          
          NLC_DN_S1[i]<-IC_PRI_D_S1[i]-sh_GVA_S1*(1-beta_S1[i])*RP[i]+rep*LC_S1[i-1]-delta[i]*KC_PRI_S1[i-1]-sh_GVA_S1*p_C_bar*(b_C[i]-b_C[i-1])   #Eq. (68)
          NLC_DN_S2[i]<-IC_PRI_D_S2[i]-sh_GVA_S2*(1-beta_S2[i])*RP[i]+rep*LC_S2[i-1]-delta[i]*KC_PRI_S2[i-1]-sh_GVA_S2*p_C_bar*(b_C[i]-b_C[i-1])   #Eq. (68)
          NLC_DN_S3[i]<-IC_PRI_D_S3[i]-sh_GVA_S3*(1-beta_S3[i])*RP[i]+rep*LC_S3[i-1]-delta[i]*KC_PRI_S3[i-1]-sh_GVA_S3*p_C_bar*(b_C[i]-b_C[i-1])   #Eq. (68)
          NLC_DN_S4[i]<-IC_PRI_D_S4[i]-sh_GVA_S4*(1-beta_S4[i])*RP[i]+rep*LC_S4[i-1]-delta[i]*KC_PRI_S4[i-1]-sh_GVA_S4*p_C_bar*(b_C[i]-b_C[i-1])   #Eq. (68)
          if (NLC_DN_S1[i]>=0) {NLC_D_S1[i]<-NLC_DN_S1[i]}  else { NLC_D_S1[i]<-0} 
          if (NLC_DN_S2[i]>=0) {NLC_D_S2[i]<-NLC_DN_S2[i]}  else { NLC_D_S2[i]<-0} 
          if (NLC_DN_S3[i]>=0) {NLC_D_S3[i]<-NLC_DN_S3[i]}  else { NLC_D_S3[i]<-0} 
          if (NLC_DN_S4[i]>=0) {NLC_D_S4[i]<-NLC_DN_S4[i]}  else { NLC_D_S4[i]<-0} 
          IG_PRI_S1[i]<-sh_GVA_S1*beta_S1[i]*RP[i]+(LG_S1[i]-LG_S1[i-1])+delta[i]*KG_PRI_S1[i-1]+sh_GVA_S1*p_G_bar*(b_G[i]-b_G[i-1])+def[i]*LG_S1[i-1]#Eq. (69)
          IG_PRI_S2[i]<-sh_GVA_S2*beta_S2[i]*RP[i]+(LG_S2[i]-LG_S2[i-1])+delta[i]*KG_PRI_S2[i-1]+sh_GVA_S2*p_G_bar*(b_G[i]-b_G[i-1])+def[i]*LG_S2[i-1]#Eq. (69)
          IG_PRI_S3[i]<-sh_GVA_S3*beta_S3[i]*RP[i]+(LG_S3[i]-LG_S3[i-1])+delta[i]*KG_PRI_S3[i-1]+sh_GVA_S3*p_G_bar*(b_G[i]-b_G[i-1])+def[i]*LG_S3[i-1]#Eq. (69)
          IG_PRI_S4[i]<-sh_GVA_S4*beta_S4[i]*RP[i]+(LG_S4[i]-LG_S4[i-1])+delta[i]*KG_PRI_S4[i-1]+sh_GVA_S4*p_G_bar*(b_G[i]-b_G[i-1])+def[i]*LG_S4[i-1]#Eq. (69)
          IC_PRI_S1[i]<-sh_GVA_S1*(1-beta_S1[i])*RP[i]+(LC_S1[i]-LC_S1[i-1])+delta[i]*KC_PRI_S1[i-1]+sh_GVA_S1*p_C_bar*(b_C[i]-b_C[i-1])+def[i]*LC_S1[i-1]#Eq. (70)
          IC_PRI_S2[i]<-sh_GVA_S2*(1-beta_S2[i])*RP[i]+(LC_S2[i]-LC_S2[i-1])+delta[i]*KC_PRI_S2[i-1]+sh_GVA_S2*p_C_bar*(b_C[i]-b_C[i-1])+def[i]*LC_S2[i-1]#Eq. (70)
          IC_PRI_S3[i]<-sh_GVA_S3*(1-beta_S3[i])*RP[i]+(LC_S3[i]-LC_S3[i-1])+delta[i]*KC_PRI_S3[i-1]+sh_GVA_S3*p_C_bar*(b_C[i]-b_C[i-1])+def[i]*LC_S3[i-1]#Eq. (70)
          IC_PRI_S4[i]<- RP[i]+(LC[i]-LC[i-1])+(LG[i]-LG[i-1])+delta[i]*K_PRI[i-1]-IG_PRI_S1[i]-IG_PRI_S2[i]-IG_PRI_S3[i] -IG_PRI_S4[i]-IC_PRI_S1[i]-IC_PRI_S2[i]-IC_PRI_S3[i]+p_G_bar*(b_G[i]-b_G[i-1])+p_C_bar*(b_C[i]-b_C[i-1])+DL[i]    #Eq. (71)
          IG_PRI[i]<-IG_PRI_S1[i]+IG_PRI_S2[i]+IG_PRI_S3[i]+IG_PRI_S4[i]# Eq. (72)
          IC_PRI[i]<- IC_PRI_S1[i]+IC_PRI_S2[i]+IC_PRI_S3[i]+IC_PRI_S4[i]#Eq. (73) 
          I_PRI[i]<-IC_PRI[i]+IG_PRI[i] #Eq. (74)
          kappa[i]<- IG_PRI[i]/I_PRI[i] # Eq. (75)
          L[i]<-LC[i]+LG[i] #Eq. (76)
          KG_PRI_S1[i]<-KG_PRI_S1[i-1]+IG_PRI_S1[i]-delta[i]*KG_PRI_S1[i-1] #Eq. (77)
          KG_PRI_S2[i]<-KG_PRI_S2[i-1]+IG_PRI_S2[i]-delta[i]*KG_PRI_S2[i-1] #Eq. (77)
          KG_PRI_S3[i]<-KG_PRI_S3[i-1]+IG_PRI_S3[i]-delta[i]*KG_PRI_S3[i-1] #Eq. (77)
          KG_PRI_S4[i]<-KG_PRI_S4[i-1]+IG_PRI_S4[i]-delta[i]*KG_PRI_S4[i-1] #Eq. (77)
          KC_PRI_S1[i]<-KC_PRI_S1[i-1]+IC_PRI_S1[i]-delta[i]*KC_PRI_S1[i-1] #Eq. (78)
          KC_PRI_S2[i]<-KC_PRI_S2[i-1]+IC_PRI_S2[i]-delta[i]*KC_PRI_S2[i-1] #Eq. (78)
          KC_PRI_S3[i]<-KC_PRI_S3[i-1]+IC_PRI_S3[i]-delta[i]*KC_PRI_S3[i-1] #Eq. (78)
          KC_PRI_S4[i]<-KC_PRI_S4[i-1]+IC_PRI_S4[i]-delta[i]*KC_PRI_S4[i-1] #Eq. (78)
          KG_PRI[i]<- KG_PRI_S1[i]+KG_PRI_S2[i]+KG_PRI_S3[i]+KG_PRI_S4[i] #Eq. (79)
          KC_PRI[i]<-KC_PRI_S1[i]+KC_PRI_S2[i]+KC_PRI_S3[i]+KC_PRI_S4[i] #Eq. (80)
          K_PRI[i]<-KC_PRI[i]+KG_PRI[i] #Eq. (81)
          KGE_PRI_S1[i]<-gamma_E1*KG_PRI_S1[i] #Eq. (82)
          KGE_PRI_S2[i]<-gamma_E2*KG_PRI_S2[i] #Eq. (82)
          KGE_PRI_S3[i]<-gamma_E3*KG_PRI_S3[i] #Eq. (82)
          KGE_PRI_S4[i]<-gamma_E4*KG_PRI_S4[i] #Eq. (82)
          KGNE_PRI_S1[i]<-(1-gamma_E1)*KG_PRI_S1[i] #Eq. (83)
          KGNE_PRI_S2[i]<-(1-gamma_E2)*KG_PRI_S2[i] #Eq. (83)
          KGNE_PRI_S3[i]<-(1-gamma_E3)*KG_PRI_S3[i] #Eq. (83)
          KGNE_PRI_S4[i]<-(1-gamma_E4)*KG_PRI_S4[i] #Eq. (83)
          KCE_PRI_S1[i]<-gamma_E1*KC_PRI_S1[i] #Eq. (84)
          KCE_PRI_S2[i]<-gamma_E2*KC_PRI_S2[i] #Eq. (84)
          KCE_PRI_S3[i]<-gamma_E3*KC_PRI_S3[i] #Eq. (84)
          KCE_PRI_S4[i]<-gamma_E4*KC_PRI_S4[i] #Eq. (84)
          KCNE_PRI_S1[i]<-(1-gamma_E1)*KC_PRI_S1[i] #Eq. (85)
          KCNE_PRI_S2[i]<-(1-gamma_E2)*KC_PRI_S2[i] #Eq. (85)
          KCNE_PRI_S3[i]<-(1-gamma_E3)*KC_PRI_S3[i] #Eq. (85)
          KCNE_PRI_S4[i]<-(1-gamma_E4)*KC_PRI_S4[i] #Eq. (85)
          KSEQ_PRI_S1[i]<-gamma_SEQ1*KGE_PRI_S1[i] #Eq. (86)
          KSEQ_PRI_S2[i]<-gamma_SEQ2*KGE_PRI_S2[i] #Eq. (86)
          KGE[i]<- KGE_PRI_S1[i]+KGE_PRI_S2[i]+KGE_PRI_S3[i]+KGE_PRI_S4[i]+gamma_E*KG_GOV[i]#Eq. (87)
          KGNE[i]<- KGNE_PRI_S1[i]+KGNE_PRI_S2[i]+KGNE_PRI_S3[i]+KGNE_PRI_S4[i]+ (1-gamma_E)*KG_GOV[i] #Eq. (88)
          KCE[i]<- KCE_PRI_S1[i]+KCE_PRI_S2[i]+KCE_PRI_S3[i]+KCE_PRI_S4[i]+gamma_E*KC_GOV[i] #Eq. (89)
          KCNE[i]<- KCNE_PRI_S1[i]+KCNE_PRI_S2[i]+KCNE_PRI_S3[i]+KCNE_PRI_S4[i]+ (1-gamma_E)*KC_GOV[i] #Eq. (90)
          KSEQ[i]<- KSEQ_PRI_S1[i]+KSEQ_PRI_S2[i] #Eq. (91)
          
          delta[i]<-delta_0+(1-delta_0)*(1-ad_K)*D_TF[i-1]  #Eq. (92)
          v[i]<-v[i-1]*(1-(1-ad_P)*D_TP[i-1]) #Eq. (93)
          g_lambda[i]<-sigma_0[i]+sigma_1+sigma_2*g_Y[i-1]#Eq. (94)
          if (lambda_is_optimal==0) { sigma_0[i]<-sigma_0[i-1]*(1-zeta_3)} else { sigma_0[i]<-sigma_0_optimal[i]} #Eq. (95)
          lambda[i]<-lambda[i-1]*(1+g_lambda[i])*(1-(1-ad_P)*D_TP[i-1]) #Eq. (96)
          wage[i]<-s_W*lambda[i]*h[i] #Eq. (97)
          N[i]<-Y[i]/(lambda[i]*h[i])  #Eq. (98)	
          ur[i]<-1-re[i] #Eq. (99)
          b_C[i]<-b_C[i-1]+(x_1[i]*IC_PRI_D[i])/p_C_bar #Eq. (100)
          b_G[i]<-b_G[i-1]+(x_2[i]*IG_PRI_D[i])/p_G_bar #Eq. (101) 
          x_1[i]<-max(x_10-x_11*yield_C[i-1],0) #Eq. (102)  
          x_2[i]<-x_20[i]-x_21*yield_G[i-1]    #Eq. (103)
          x_20[i]<-x_20[i-1]*(1+g_x_20[i]) #Eq. (104)
          g_x_20[i]<-g_x_20[i-1]*(1-zeta_4) #Eq. (105)
          yield_C[i]<-coupon_C[i]/p_C[i] #Eq. (106)
          yield_G[i]<-coupon_G[i]/p_G[i] #Eq. (107)
          coupon_C[i]<-yield_C[i-1]*p_C_bar #Eq. (108)
          coupon_G[i]<-yield_G[i-1]*p_G_bar #Eq. (109)
          B_C[i]<-B_CH[i]+B_CCB[i]  #Eq. (110)
          B_G[i]<-B_GH[i]+B_GCB[i]  #Eq. (111)
          p_C[i]<-B_C[i]/b_C[i] #Eq. (112)
          p_G[i]<-B_G[i]/b_G[i]  #Eq. (113)
          B[i]<-B_C[i]+B_G[i] #Eq. (114)
          DL[i]<-def[i]*L[i-1] #Eq. (115)
          def[i]<-def_max/(1+def_0*exp(def_1-def_2*illiq[i-1])) #Eq. (116)
          
          illiq[i]<-((int_C_S1[i-1]+rep)*LC_S1[i-1]+(int_C_S2[i-1]+rep)*LC_S2[i-1]+(int_C_S3[i-1]+rep)*LC_S3[i-1]+(int_C_S4[i-1]+rep)*LC_S4[i-1]+(int_G[i-1]+rep)*LG[i-1]+coupon_C[i-1]*b_C[i-1]+coupon_G[i-1]*b_G[i-1]+wage[i]*N[i]+TAX_F[i]+TAX_C[i]-SUB[i]+delta[i]*K_PRI[i-1])/(Y[i]+(1-CR_C_S1[i])*NLC_D_S1[i]+(1-CR_C_S2[i])*NLC_D_S2[i]+(1-CR_C_S3[i])*NLC_D_S3[i]+(1-CR_C_S4[i])*NLC_D_S4[i]+ (1-CR_G[i])*(NLG_D_S1[i]+ NLG_D_S2[i]+ NLG_D_S3[i]+NLG_D_S4[i])+(b_C[i]-b_C[i-1])*p_C_bar+ (b_G[i]-b_G[i-1])*p_G_bar)  #Eq. (117)
          
          dsr[i]<-((int_C_S1[i-1]+rep)*LC_S1[i-1]+(int_C_S2[i-1]+rep)*LC_S2[i-1]+(int_C_S3[i-1]+rep)*LC_S3[i-1]+(int_C_S4[i-1]+rep)*LC_S4[i-1]+(int_G[i-1]+rep)*LG[i-1]+coupon_C[i-1]*b_C[i-1]+coupon_G[i-1]*b_G[i-1])/(TP[i]+int_C_S1[i-1]*LC_S1[i-1]+int_C_S2[i-1]*LC_S2[i-1]+int_C_S3[i-1]*LC_S3[i-1]+int_C_S4[i-1]*LC_S4[i-1]+(int_G[i-1])*LG[i-1]+coupon_C[i-1]*b_C[i-1]+coupon_G[i-1]*b_G[i-1]) #Eq. (118)
          
          #2.2.3 Households
          Y_HG[i]<-wage[i]*N[i]+DP[i]+BP_D[i]+int_D*D[i-1]+int_S*SEC_H[i-1]+coupon_C[i-1]*b_CH[i-1]+coupon_G[i-1]*b_GH[i-1]  #Eq. (119)
          Y_H[i]<-Y_HG[i]-TAX_H[i] #Eq. (120)
          CO_PRI_N[i]<-((c_1*(1-covid*0.2)-con_change*c_1)*Y_H[i-1]+(c_2-con_change*c_2)*V_HF[i-1])*(1-D_T[i-1]) #Eq. (121)
          if (Y_N[i]<Y_star[i]) {CO_PRI[i]<-CO_PRI_N[i]} else {CO_PRI[i]<- pr*(Y_star[i]-CO_GOV[i]-I_GOV[i]-I_PRI[i])} #Eq. (122)
          V_HF[i]<-V_HF[i-1]+Y_H[i]-CO_PRI[i]+b_CH[i-1]*(p_C[i]-p_C[i-1])+b_GH[i-1]*(p_G[i]-p_G[i-1])#Eq. (123)
          SEC_H[i]<-(lambda_10+lambda_10prime*D_T[i-1]+lambda_11*int_S+lambda_12*yield_C[i-1]+lambda_13*yield_G[i-1]+lambda_14*int_D+lambda_15*(Y_H[i-1]/V_HF[i-1]))*V_HF[i-1] #Eq. (124)
          B_CH[i]<-(lambda_20+lambda_20prime*D_T[i-1]+lambda_21*int_S+lambda_22*yield_C[i-1]+lambda_23*yield_G[i-1]+lambda_24*int_D+lambda_25*(Y_H[i-1]/V_HF[i-1]))*V_HF[i-1] #Eq. (125)
          B_GH[i]<-(lambda_30[i]+lambda_30prime*D_T[i-1]+lambda_31*int_S+lambda_32*yield_C[i-1]+lambda_33*yield_G[i-1]+lambda_34*int_D+lambda_35*(Y_H[i-1]/V_HF[i-1]))*V_HF[i-1] #Eq. (126)
          D_N[i]<-(lambda_40[i]+lambda_40prime*D_T[i-1]+lambda_41*int_S+lambda_42*yield_C[i-1]+lambda_43*yield_G[i-1]+lambda_44*int_D+lambda_45*(Y_H[i-1]/V_HF[i-1]))*V_HF[i-1] #Eq. (127n)
          D[i]<-D[i-1]+Y_H[i]-CO_PRI[i]-(SEC_H[i]-SEC_H[i-1])-(b_CH[i]-b_CH[i-1])*p_C_bar-(b_GH[i]-b_GH[i-1])*p_G_bar  #Eq. (127)
          lambda_30[i]<-lambda_30[i-1]*(1+g_lambda_30[i]) #Eq. (128)
          g_lambda_30[i]<-zeta_10*g_bG[i-1] #Eq. (129)
          lambda_40[i]<-1- lambda_10- lambda_20- lambda_30[i] #(Ci) 
          b_CH[i]<-B_CH[i]/p_C[i] #Eq. (130)
          b_GH[i]<-B_GH[i]/p_G[i] #Eq. (131)
          DC[i]<-DC[i-1]+CO_PRI[i]-xi*DC[i-1] #Eq. (132)
          g_POP[i]<-g_POP[i-1]*(1-zeta_5) #Eq. (133)
          POP[i]<-POP[i-1]*(1+g_POP[i]) #Eq. (134)
          LF[i]<-(lf_1[i]-lf_2*hazratio[i-1])*(1-(1-ad_LF)*D_TF[i-1])*POP[i] #Eq. (135)
          lf_1[i]<-lf_1[i-1]*(1-zeta_6) #Eq. (136)
          
          #2.2.4 Banks
          BP[i]<-int_C_S1[i-1]*LC_S1[i-1]+int_C_S2[i-1]*LC_S2[i-1] +int_C_S3[i-1]*LC_S3[i-1] +int_C_S4[i-1]*LC_S4[i-1]+int_G[i-1]*(LG_S1[i-1]+LG_S2[i-1]+LG_S3[i-1]+LG_S4[i-1])+int_S*SEC_B[i-1]-int_D*D[i-1]-int_A*A[i-1] #Eq. (137)
          CAP[i]<-CAP[i-1]+BP_U[i]-DL[i]+BAILOUT[i] #Eq. (138)
          CAP_before[i]<-CAP[i-1]+BP_U[i]-DL[i]
          BP_U[i]<-s_B*BP[i-1] #Eq. (139)
          BP_D[i]<-BP[i]-BP_U[i] #Eq. (140)
          HPM[i]<-h_1*D[i] #Eq. (141)
          SEC_BN[i]<-SEC_B[i-1]+(A[i]-A[i-1])+(D[i]-D[i-1])+BP_U[i]+BAILOUT[i] -(LG[i]-LG[i-1])-(LC[i]-LC[i-1])-(HPM[i]-HPM[i-1])-DL[i] 
          if (iterations>10 & A_N[i]>0) {SEC_B[i]<-h_2*D[i]} else {SEC_B[i]<-SEC_BN[i]} #Eq. (142)
          A_N[i]<-A[i-1]+(h_2*D[i]-SEC_B[i-1])+(LG[i]-LG[i-1])+(LC[i]-LC[i-1])+(HPM[i]-HPM[i-1])+DL[i]-(D[i]-D[i-1])-BP_U[i]-BAILOUT[i]   
          if (iterations>10 & A_N[i]>0) {A[i]<-A_N[i]} else {A[i]<-0} #Eq. (143)
          
          if (i<5) {w_G[i]<-w_G_initial} else {w_G[i]<-w_G_2022} 
          if (i<5) {w_C_S1[i]<-w_C_S1_initial} else {w_C_S1[i]<-w_C_S1_2022} 
          if (i<5) {w_C_S2[i]<-w_C_S2_initial} else {w_C_S2[i]<-w_C_S2_2022} 
          if (i<5) {w_C_S3[i]<-w_C_S3_initial} else {w_C_S3[i]<-w_C_S3_2022} 
          if (i<5) {w_C_S4[i]<-w_C_S4_initial} else {w_C_S4[i]<-w_C_S4_2022} 
          
          if (i<=6) { CR[i]<-CR_max/(1+r_0*exp(r_1-r_2*dsr[i-1]+r_3*(CAR[i-1]-CAR_min)))*cr_rationing_dummy } else { CR[i]<- (1-CR_dummy)*CR_max/(1+r_0*exp(r_1-r_2*dsr[i-1]+r_3*(CAR[i-1]-CAR_min)))*cr_rationing_dummy +CR_dummy*CR[i-1] } #Eq. (144)
          if (i<=6) { CR_G[i]<- (1+l_1*(w_G[i-1]-w_LT[i-1]))*CR[i] } else { CR_G[i]<- (1-CR_dummy)* (1+l_1*(w_G[i-1]-w_LT[i-1]))*CR[i] +CR_dummy*CR_G[i-1] } #Eq. (145)
          if (i<=6) {CR_C_S1[i]<-(1+l_1*(w_C_S1[i-1]-w_LT[i-1]))*CR[i]  } else {CR_C_S1[i]<- (1-CR_dummy)* (1+l_1*(w_C_S1[i-1]-w_LT[i-1]))*CR[i] +CR_dummy*CR_C_S1[i-1] } #Eq. (146)
          if (i<=6) {CR_C_S2[i]<-(1+l_1*(w_C_S2[i-1]-w_LT[i-1]))*CR[i]  } else {CR_C_S2[i]<- (1-CR_dummy)* (1+l_1*(w_C_S2[i-1]-w_LT[i-1]))*CR[i]+CR_dummy*CR_C_S2[i-1] } #Eq. (146)
          if (i<=6) {CR_C_S3[i]<-(1+l_1*(w_C_S3[i-1]-w_LT[i-1]))*CR[i]  } else {CR_C_S3[i]<- (1-CR_dummy)* (1+l_1*(w_C_S3[i-1]-w_LT[i-1]))*CR[i] +CR_dummy*CR_C_S3[i-1] } #Eq. (146)
          CR_C_S4[i]<-(CR[i]-sh_NLG[i-1]*CR_G[i]-sh_NLC_S1[i-1]*CR_C_S1[i]-sh_NLC_S2[i-1]*CR_C_S2[i]-sh_NLC_S3[i-1]*CR_C_S3[i])/(sh_NLC_S4[i-1]+0.001)  #0.001 has been used in order to allow for CR_C_S4[i[ to be defined even if sh_NLC_4 is zero; #Eq. (147) 
          LC_S1[i]<-LC_S1[i-1]+(1-CR_C_S1[i])*NLC_D_S1[i]-rep*LC_S1[i-1]-def[i]*LC_S1[i-1] #Eq. (148)
          LC_S2[i]<-LC_S2[i-1]+(1-CR_C_S2[i])*NLC_D_S2[i]-rep*LC_S2[i-1]-def[i]*LC_S2[i-1] #Eq. (148)
          LC_S3[i]<-LC_S3[i-1]+(1-CR_C_S3[i])*NLC_D_S3[i]-rep*LC_S3[i-1]-def[i]*LC_S3[i-1] #Eq. (148)
          LC_S4[i]<-LC_S4[i-1]+(1-CR_C_S4[i])*NLC_D_S4[i]-rep*LC_S4[i-1]-def[i]*LC_S4[i-1] #Eq. (148)
          LG_S1[i]<-LG_S1[i-1]+(1-CR_G[i])*NLG_D_S1[i]-rep*LG_S1[i-1]-def[i]*LG_S1[i-1]  #Eq. (149)
          LG_S2[i]<-LG_S2[i-1]+(1-CR_G[i])*NLG_D_S2[i]-rep*LG_S2[i-1]-def[i]*LG_S2[i-1]  #Eq. (149)
          LG_S3[i]<-LG_S3[i-1]+(1-CR_G[i])*NLG_D_S3[i]-rep*LG_S3[i-1]-def[i]*LG_S3[i-1]  #Eq. (149)
          LG_S4[i]<-LG_S4[i-1]+(1-CR_G[i])*NLG_D_S4[i]-rep*LG_S4[i-1]-def[i]*LG_S4[i-1]  #Eq. (149)
          LC[i]<-LC_S1[i]+LC_S2[i]+LC_S3[i]+LC_S4[i] #Eq. (150)
          LG[i]<-LG_S1[i]+LG_S2[i]+LG_S3[i]+LG_S4[i] #Eq. (151)
          lev_B[i]<-(HPM[i]+LG[i]+LC[i]+SEC_B[i])/CAP[i] #Eq. (152)
          CAR[i]<-CAP[i]/(w_G[i]*LG[i]+w_C_S1[i]*LC_S1[i]+w_C_S2[i]*LC_S2[i]+w_C_S3[i]*LC_S3[i] +w_C_S4[i]*LC_S4[i]+w_S*SEC_B[i]+w_H*HPM[i]) #Eq. (153)
          lev_B_before[i]<- (HPM[i-1]+LG[i-1]+LC[i-1]+SEC_B[i-1])*(1+g_Y[i-1])/CAP_before[i]
          CAR_before[i]<-CAP_before[i]/((w_G[i]*LG[i-1]+w_C_S1[i]*LC_S1[i-1]+w_C_S2[i]*LC_S2[i-1]+w_C_S3[i]*LC_S3[i-1]+w_C_S4[i]*LC_S4[i-1]+w_S*SEC_B[i-1]+w_H*HPM[i-1])*(1+g_Y[i-1]))
          w_LT[i]<-sh_LG[i-1]*w_G[i]+sh_LC_S1[i-1]*w_C_S1[i]+sh_LC_S2[i-1]*w_C_S2[i]+sh_LC_S3[i-1]*w_C_S3[i]+sh_LC_S4[i-1]*w_C_S4[i]#Eq. (154)
          int_G[i]<-spr_G[i]+int_A #Eq. (155)
          int_C_S1[i]<-spr_C_S1[i]+int_A #Eq. (156)
          int_C_S2[i]<-spr_C_S2[i]+int_A #Eq. (156)
          int_C_S3[i]<-spr_C_S3[i]+int_A #Eq. (156)
          int_C_S4[i]<-spr_C_S4[i]+int_A #Eq. (156)
          
          if (i<=7) { spr[i]<-spr_0-spr_1*(CAR[i-1]-CAR_min)+spr_2*dsr[i-1] } else { spr[i]<- (1-CR_dummy)*(spr_0-spr_1*(CAR[i-1]-CAR_min)+spr_2*dsr[i-1])+CR_dummy*spr[i-1] } #Eq. (157)
          if (i<=7) { spr_G[i]<-(1+spr_3*(w_G[i-1]-w_LT[i-1]))*spr[i]  } else { spr_G[i]<- (1-CR_dummy)*(1+spr_3*(w_G[i-1]-w_LT[i-1]))*spr[i] +CR_dummy*spr_G[i-1] } #Eq. (158)
          if (i<=7) { spr_C_S1[i]<-(1+spr_3*(w_C_S1[i-1]-w_LT[i-1]))*spr[i]  } else { spr_C_S1[i]<- (1-CR_dummy)*(1+spr_3*(w_C_S1[i-1]-w_LT[i-1]))*spr[i]  +CR_dummy*spr_C_S1[i-1] } #Eq. (159)
          if (i<=7) { spr_C_S2[i]<- (1+spr_3*(w_C_S2[i-1]-w_LT[i-1]))*spr[i]  } else { spr_C_S2[i]<- (1-CR_dummy)*(1+spr_3*(w_C_S2[i-1]-w_LT[i-1]))*spr[i] +CR_dummy*spr_C_S2[i-1] } #Eq. (159)
          if (i<=7) { spr_C_S3[i]<- (1+spr_3*(w_C_S3[i-1]-w_LT[i-1]))*spr[i]  } else { spr_C_S3[i]<- (1-CR_dummy)*(1+spr_3*(w_C_S3[i-1]-w_LT[i-1]))*spr[i]  +CR_dummy*spr_C_S3[i-1] } #Eq. (159)
          spr_C_S4[i]<-(spr[i]-sh_LG[i-1]*spr_G[i]-sh_LC_S1[i-1]*spr_C_S1[i]-sh_LC_S2[i-1]*spr_C_S2[i]-sh_LC_S3[i-1]*spr_C_S3[i])/sh_LC_S4[i-1] #Eq. (160)
          
          sh_NLG[i]<-(NLG_D_S1[i]+NLG_D_S2[i]+NLG_D_S3[i]+NLG_D_S4[i])/(NLG_D_S1[i]+NLG_D_S2[i]+NLG_D_S3[i]+NLG_D_S4[i]+NLC_D_S1[i]+NLC_D_S2[i]+NLC_D_S3[i]+NLC_D_S4[i])
          sh_NLC_S1[i]<-NLC_D_S1[i]/(NLG_D_S1[i]+NLG_D_S2[i]+NLG_D_S3[i]+NLG_D_S4[i]+NLC_D_S1[i]+NLC_D_S2[i]+NLC_D_S3[i]+NLC_D_S4[i])
          sh_NLC_S2[i]<-NLC_D_S2[i]/(NLG_D_S1[i]+NLG_D_S2[i]+NLG_D_S3[i]+NLG_D_S4[i]+NLC_D_S1[i]+NLC_D_S2[i]+NLC_D_S3[i]+NLC_D_S4[i])
          sh_NLC_S3[i]<-NLC_D_S3[i]/(NLG_D_S1[i]+NLG_D_S2[i]+NLG_D_S3[i]+NLG_D_S4[i]+NLC_D_S1[i]+NLC_D_S2[i]+NLC_D_S3[i]+NLC_D_S4[i])
          sh_NLC_S4[i]<-1-sh_NLG[i]-sh_NLC_S1[i]-sh_NLC_S2[i]-sh_NLC_S3[i]
          
          sh_LG[i]<-LG[i]/L[i]
          sh_LC_S1[i]<-LC_S1[i]/L[i]
          sh_LC_S2[i]<- LC_S2[i]/L[i]
          sh_LC_S3[i]<- LC_S3[i]/L[i]
          sh_LC_S4[i]<- 1-sh_LG[i]-sh_LC_S1[i] -sh_LC_S2[i]-sh_LC_S3[i]
          
          #2.2.5 Government
          GNS[i]<- -CO_GOV[i]+TAX[i]-int_S*SEC[i-1]+CBP[i]-SUB[i]-delta[i]*K_GOV[i-1] #Eq. (161)
          SEC[i]<-SEC[i-1]+I_GOV[i]-GNS[i]-delta[i]*K_GOV[i-1]+BAILOUT[i]#Eq. (162)
          IG_GOV[i]<-gov_IG*Y[i-1] #Eq. (163)
          IC_GOV[i]<-gov_IC*Y[i-1] #Eq. (164)
          I_GOV[i]<- IC_GOV[i]+IG_GOV[i] #Eq. (165)
          KG_GOV[i]<-KG_GOV[i-1]+IG_GOV[i]-delta[i]*KG_GOV[i-1] #Eq. (166)
          KC_GOV[i]<-KC_GOV[i-1]+IC_GOV[i]-delta[i]*KC_GOV[i-1] #Eq. (167)
          K_GOV[i]<-KC_GOV[i]+KG_GOV[i] #Eq. (168)
          K[i]<-K_PRI[i]+K_GOV[i] #Eq. (169)
          KG[i]<-KG_PRI[i]+KG_GOV[i] #Eq. (170)
          KC[i]<-KC_PRI[i]+KC_GOV[i] #Eq. (171)
          BAILOUT_CAR[i]<-if (CAR_before[i]>CAR_min) {BAILOUT_CAR[i]<-0}  else {BAILOUT_CAR[i]<-CAR_min*(w_G[i]*LG[i-1]+w_C_S1[i]*LC_S1[i-1]+w_C_S2[i]*LC_S2[i-1]+w_C_S3[i]*LC_S3[i-1]+w_C_S4[i]*LC_S4[i-1]+w_S*SEC_B[i-1])*(1+g_Y[i-1])-CAP_before[i]}
          BAILOUT_levB[i]<-if (lev_B_before[i]<lev_B_max) {BAILOUT_levB[i]<-0}  else {BAILOUT_levB[i]<-(HPM[i-1]+LG[i-1]+LC[i-1]+SEC_B[i-1])*(1+g_Y[i-1])/lev_B_max-CAP_before[i]}
          BAILOUT[i]<-max(BAILOUT_CAR[i],BAILOUT_levB[i])
          CO_GOV[i]<-gov_C*(1+covid*0.4)*Y[i-1] #Eq. (172)
          
          if (i<5) {SUB[i]<- tau_C[i]*EMIS_IN[i-1]} else {SUB[i]<- (1-stop_recycling_dummy)* tau_C[i]*EMIS_IN[i-1] + stop_recycling_dummy*((1-subsidy_increase_dummy)*tau_C[5]*EMIS_IN[4]+subsidy_increase_dummy*gov_SUB[i]*E_NF[i-1]*ucr[i-1])} # Eq. (173)
          
          if (i<5) {gov_SUB[i]<- min(SUB[i]/(E_NF[i-1]*ucr[i-1]),0.999)} else {gov_SUB[i]<-(1-subsidy_increase_dummy)*min(SUB[i]/(E_NF[i-1]*ucr[i-1]),0.999)+subsidy_increase_dummy*subsidy_ratio*gov_SUB[4]} # Eq. (174)
          
          TAX_H[i]<-tau_H*Y_HG[i-1] #Eq. (175)
          TAX_F[i]<-tau_F*TP_G[i-1] #Eq. (176)
          TAX_C[i]<-tau_C[i]*EMIS_IN[i-1] #Eq. (177)
          
          if (i<5) {tau_C[i]<- tau_C_baseline[i]} else {tau_C[i]<- (1- tau_C_dummy_H - tau_C_dummy_M)*tau_C_baseline[i]+tau_C_dummy_H*tau_C_policy_High[i] +tau_C_dummy_M*tau_C_policy_Medium[i]} 
          tau_C_baseline[i]<-Data[i,c(1)]/1000
          tau_C_policy_Medium[i]<-Data[i,c(2)]/1000
          tau_C_policy_High[i]<-Data[i,c(3)]/1000
          TAX[i]<-TAX_H[i]+TAX_F[i]+TAX_C[i]    #Eq. (178)
          
          #2.2.6 Central bank
          CBP[i]<-coupon_C[i-1]*b_CCB[i-1]+coupon_G[i-1]*b_GCB[i-1]+int_A*A[i-1]+int_S*SEC_CB[i-1]  #Eq. (179)
          B_GCB[i]<-s_G*B_G[i-1] #Eq. (180)
          B_CCB[i]<-s_C*B_C[i-1] #Eq. (181)
          b_CCB[i]<-B_CCB[i]/p_C[i] #Eq. (182)
          b_GCB[i]<-B_GCB[i]/p_G[i] #Eq. (183)
          SEC_CB[i]<-SEC[i]-SEC_H[i]-SEC_B[i]  #Eq. (184)
          SEC_CBred[i]<-SEC_CB[i-1]+(HPM[i]-HPM[i-1])-(A[i]-A[i-1])-p_C_bar*(b_CCB[i]-b_CCB[i-1])-p_G_bar*(b_GCB[i]-b_GCB[i-1]) #Eq. (185-red)
          
          #Auxiliary equations
          g_Y[i]<-(Y[i]-Y[i-1])/Y[i-1] 
          gI[i]<-(I_PRI[i]-I_PRI[i-1])/I_PRI[i-1]  
          gC[i]<-(CO_PRI[i]-CO_PRI[i-1])/CO_PRI[i-1]  
          g_K[i]<-(K_PRI[i]-K_PRI[i-1])/K_PRI[i-1]
          ur_per[i]<-ur[i]*100
          g_Y_per[i]<-g_Y[i]*100
          theta_per[i]<-theta[i]*100
          greenI[i]<-IG_PRI[i]/I_PRI[i]
          r_total[i]<-(TP[i]/K_PRI[i]) *100
          Y_POP_ratio[i]<-(Y[i]/POP[i])
          I_Y_ratio[i]<-I_PRI[i]/Y[i]
          C_Y_ratio[i]<-CO_PRI[i]/Y[i] 
          C_K_ratio[i]<-CO_PRI[i]/K_PRI[i]
          Y_K_ratio[i]<-Y[i]/K_PRI[i]
          Y_H_Y_ratio[i]<-Y_H[i]/Y[i]
          Y_H_K_ratio[i]<-Y_H[i]/K_PRI[i]
          A_K_ratio[i]<-A[i]/K_PRI[i]
          sigma_0_optimal[i]<-(((g_Y[i-1]-g_lf[i])/(1+ g_lf[i])*(1-(1-ad_P)*D_TP[i-1]))+(1-ad_P)*D_TP[i-1])/(1-(1-ad_P)*D_TP[i-1]) -sigma_1-sigma_2*g_Y[i-1]
          h_optimal[i]<-h_optimal[i-1]*(1+g_h[i])
          g_h[i]<-((((g_Y[i]-g_lf[i])/(1+ g_lf[i]))+(1-ad_P)*D_TP[i-1])/(1-(1-ad_P)*D_TP[i-1]) -g_lambda[i])/(1+g_lambda[i])
          if (h_is_optimal==0) {h[i]<-hours} else {h[i]<-h_optimal[i]}
          omega_ratio[i]<-omega[i]/omega[1]
          mu_ratio[i]<-mu[i]/mu[1]
          rho_ratio[i]<-rho[i]/rho[1]
          epsilon_ratio[i]<-epsilon[i]/epsilon[1]
          Wbill[i]<-wage[i]*N[i]
          Interest[i]<- int_C_S1[i]*LC_S1[i-1]+int_C_S2[i]*LC_S2[i-1]+int_C_S3[i]*LC_S3[i-1]+int_C_S4[i]*LC_S4[i-1]+int_G[i]*LG[i-1]
          Depreciation[i]<-delta[i]*K_PRI[i-1]
          SEC_Y[i]<-SEC[i]/Y[i]
          SEC_Y_per[i]<-SEC[i]/Y[i]*100
          fiscal_balance_per[i]<-fiscal_balance[i]*100
          PORT_BCH[i]<-B_CH[i]/V_HF[i-1]
          PORT_BGH[i]<-B_GH[i]/V_HF[i-1]
          PORT_SECH[i]<-SEC_H[i]/V_HF[i-1]
          PORT_D[i]<-D_N[i]/V_HF[i-1]
          L_K[i]<-L[i]/K_PRI[i]
          I[i]<-I_PRI[i]+I_GOV[i]
          sh_L[i]<-L[i]/(L[i]+B[i])
          g_bC[i]<-(b_C[i]-b_C[i-1])/b_C[i-1]
          g_bG[i]<-(b_G[i]-b_G[i-1])/b_G[i-1]
          g_pG[i]<-(p_G[i]-p_G[i-1])/p_G[i-1]
          g_pC[i]<-(p_C[i]-p_C[i-1])/p_C[i-1]
          B_C_issue[i]<-x_1[i]*IC_PRI_D[i]  
          B_G_issue[i]<-x_2[i]*IG_PRI_D[i] 
          IG_cum[i]<-IG_cum[i-1]+IG_PRI[i]+IG_GOV[i]
          IC_cum[i]<-IC_cum[i-1]+IC_PRI[i]+IC_GOV[i]
          IG_E_Y_cum[i]<-IG_E_Y_cum[i-1]+IG_E_Y[i]
          IG_E[i]<- gamma_E1*IG_PRI_S1[i]+ gamma_E2*IG_PRI_S2[i]+ gamma_E3*IG_PRI_S3[i]+ gamma_E4*IG_PRI_S4[i]+ gamma_E*IG_GOV[i]
          IG_E_Y[i]<- IG_E[i]/Y[i]*100
          yield_C_per[i]<- yield_C[i]*100
          yield_G_per[i]<- yield_G[i]*100
          
          V_CB[i]<- p_C_bar*b_CCB[i]+ p_G_bar*b_GCB[i]+A[i]+SEC_CB[i]-HPM[i] 
          gGOVCO[i]<-(CO_GOV[i]-CO_GOV[i-1])/CO_GOV[i-1]
          D_K_ratio[i]<-D[i]/K_PRI[i]
          HPM_K_ratio[i]<-HPM[i]/K_PRI[i]
          I_K_ratio[i]<-I_PRI[i]/K_PRI[i-1]
          ID_K_ratio[i]<-I_PRI_D[i]/K_PRI[i-1]
          fiscal_balance[i]<-(TAX[i]-CO_GOV[i]-IC_GOV[i]-IG_GOV[i]-SUB[i]-BAILOUT[i]-int_S*SEC[i-1])/Y[i]
          g_lf[i]<-(LF[i]-LF[i-1])/LF[i-1]
          W_POP_ratio[i]<-W[i]/POP[i]
          E_N_ratio[i]<-E[i]/N[i]
          g_EN_ratio[i]<-(E_N_ratio[i]-E_N_ratio[i-1])/E_N_ratio[i-1]
          Y_E_ratio[i]<-Y[i]/E[i]
          g_YE_ratio[i]<-(Y_E_ratio[i]-Y_E_ratio[i-1])/Y_E_ratio[i-1]
          lambda_perworker[i]<-lambda[i]*h[i]
          E_ratio[i]<-E[i]/E[1]
          CO2_ratio[i]<-EMIS_IN[i]/EMIS_IN[1]
          lf[i]<-LF[i]/POP[i]
          haz_flow[i]<-haz*W[i]
          Y_HD[i]<-wage[i]*N[i]-TAX_H[i]+DP[i]+BP_D[i]+int_D*D[i-1]+int_S*SEC_H[i-1]+coupon_C[i]*b_CH[i-1]+coupon_G[i]*b_GH[i-1]- xi*DC[i-1]
          V_H[i]<- DC[i]+D[i]+p_C_bar*b_CH[i]+p_G_bar*b_GH[i]+SEC_H[i]   
          KG_GOV_Y_ratio[i]<- KG_GOV[i]/Y[i]
          KC_GOV_Y_ratio[i]<-KC_GOV[i]/Y[i] 
          IG_GOV_Y_ratio[i]<- IG_GOV[i]/Y[i]
          IC_GOV_Y_ratio[i]<-IC_GOV[i]/Y[i] 
          KG_Y_ratio[i]<- KG_PRI[i]/Y[i]
          KC_Y_ratio[i]<-KC_PRI[i]/Y[i] 
          KG_K_ratio[i]<- KG_PRI[i]/K_PRI[i]
          LG_L_ratio_pseudo[i]<- (beta[i]*(g_baseline*K_PRI[i]/(1+g_baseline)-RP[i]-x_2[i]*I_PRI_D[i])/L[i])/(def[i]/(1+g_baseline)+g_baseline/(1+g_baseline))
          LG_L_ratio[i]<- LG[i]/L[i]
          KGE_KCE[i]<-KGE[i]/KCE[i]
          KGNE_KCNE[i]<-KGNE[i]/KCNE[i]
          KSEQ_KCE[i]<-KSEQ[i]/(KCE_PRI_S1[i]+KCE_PRI_S2[i])
          g_Y_cum[i]<-g_Y_cum[i-1]+g_Y[i]
          def_per[i]<-def[i]*100
          CAR_per[i]<-CAR[i]*100
          beta[i]<-beta[i-1] 
          I_PRI_S1[i]<-IC_PRI_S1[i]+IG_PRI_S1[i] #Eq. (72)
          I_PRI_S2[i]<- IC_PRI_S2[i]+IG_PRI_S2[i] #Eq. (72)
          I_PRI_S3[i]<- IC_PRI_S3[i]+IG_PRI_S3[i] #Eq. (72)
          I_PRI_S4[i]<- IC_PRI_S4[i]+IG_PRI_S4[i] #Eq. (72)
          IG_PRI_D[i]<- IG_PRI_D_S1[i]+ IG_PRI_D_S2[i]+ IG_PRI_D_S3[i]+ IG_PRI_D_S4[i]  
          IC_PRI_D[i]<-I_PRI_D[i]-IG_PRI_D[i] 
          TAX_C_change[i]<- TAX_C[i]-TAX_C[i-1]
          SUB_change[i]<- SUB[i]-SUB[i-1]
          IG_GOV_change[i]<- IG_GOV[i]-IG_GOV[i-1]
          EMIS_IN_Y[i]<- EMIS_IN[i]/Y[i]
          EMIS_IN_Y_SSP360[i]<- Data[i,c(4)]
          
          private_balance[i]<-Y[i]-CO_PRI[i]-I_PRI[i]-TAX[i]+int_S*SEC_H[i-1]+int_S*SEC_B[i-1]+SUB[i]+BAILOUT[i]-coupon_C[i]*b_CCB[i-1]-coupon_G[i]*b_GCB[i-1]-int_A*A[i-1]
          firms_balance[i]<-RP[i]+s_F*delta[i]*K_PRI[i-1]-I_PRI[i]
          households_balance[i]<- Y_H[i]+(1-s_F)*delta[i]*K_PRI[i-1]-CO_PRI[i]
          banks_balance[i]<- BP_U[i]+BAILOUT[i]
          fiscal_balance_help[i]<- TAX[i]-CO_GOV[i]-IC_GOV[i]-IG_GOV[i]-SUB[i]-BAILOUT[i]-int_S*SEC[i-1]+CBP[i]
          private_balance_per[i]<- private_balance[i]/Y[i]*100
          firms_balance_per[i]<- firms_balance[i]/Y[i]*100
          households_balance_per[i]<- households_balance[i]/Y[i]*100
          banks_balance_per[i]<- banks_balance[i]/Y[i]*100
          fiscal_balance_help_per[i]<-fiscal_balance_help[i]/Y[i]*100
          g_POP_stag[i]<-(POP[i]-POP[1])/POP[1]
          tau_C_baseline_per[i]<-Data[i,c(1)]
          tau_C_policy_Medium_per[i]<- Data[i,c(2)]
          tau_C_policy_High_per[i]<- Data[i,c(3)]
          
          
        }
      }
    }
    
    #####################################
    #4. FILLING IN THE MONTE-CARLO VARIABLES
    #####################################
    
    eval(parse(text=(paste("Monte_g_Y_per",sce,"[,",k,"]","<-", "g_Y_per", sep=""))))
    eval(parse(text=(paste("Monte_epsilon_ratio",sce,"[,",k,"]","<-", "epsilon_ratio", sep=""))))
    eval(parse(text=(paste("Monte_mu_ratio",sce,"[,",k,"]","<-", "mu_ratio", sep=""))))
    eval(parse(text=(paste("Monte_IG_PRI",sce,"[,",k,"]","<-", "IG_PRI", sep=""))))
    eval(parse(text=(paste("Monte_IG_GOV",sce,"[,",k,"]","<-", "IG_GOV", sep=""))))
    eval(parse(text=(paste("Monte_IG_E_Y",sce,"[,",k,"]","<-", "IG_E_Y", sep=""))))
    eval(parse(text=(paste("Monte_SUB",sce,"[,",k,"]","<-", "SUB", sep=""))))
    eval(parse(text=(paste("Monte_ur_per",sce,"[,",k,"]","<-", "ur_per", sep=""))))
    eval(parse(text=(paste("Monte_Y",sce,"[,",k,"]","<-", "Y", sep=""))))
    eval(parse(text=(paste("Monte_CO_PRI",sce,"[,",k,"]","<-", "CO_PRI", sep=""))))
    eval(parse(text=(paste("Monte_I_PRI",sce,"[,",k,"]","<-", "I_PRI", sep=""))))
    eval(parse(text=(paste("Monte_lambda",sce,"[,",k,"]","<-", "lambda", sep=""))))
    eval(parse(text=(paste("Monte_N",sce,"[,",k,"]","<-", "N", sep=""))))
    eval(parse(text=(paste("Monte_Y_POP_ratio",sce,"[,",k,"]","<-", " Y_POP_ratio", sep=""))))
    eval(parse(text=(paste("Monte_lambda_perworker",sce,"[,",k,"]","<-","lambda_perworker", sep=""))))
    eval(parse(text=(paste("Monte_EMIS_IN_Y",sce,"[,",k,"]","<-","EMIS_IN_Y", sep=""))))
    eval(parse(text=(paste("Monte_EMIS_IN",sce,"[,",k,"]","<-","EMIS_IN", sep=""))))
    eval(parse(text=(paste("Monte_E_N_ratio",sce,"[,",k,"]","<-","E_N_ratio", sep=""))))
    eval(parse(text=(paste("Monte_I_K_ratio",sce,"[,",k,"]","<-","I_K_ratio", sep=""))))
    eval(parse(text=(paste("Monte_ur_per",sce,"[,",k,"]","<-","ur_per", sep=""))))
    eval(parse(text=(paste("Monte_POP",sce,"[,",k,"]","<-","POP", sep=""))))
    eval(parse(text=(paste("Monte_beta_S1",sce,"[,",k,"]","<-","beta_S1", sep=""))))
    eval(parse(text=(paste("Monte_greenI",sce,"[,",k,"]","<-","greenI", sep=""))))
    eval(parse(text=(paste("Monte_kappa",sce,"[,",k,"]","<-","kappa", sep=""))))
    eval(parse(text=(paste("Monte_theta_per",sce,"[,",k,"]","<-","theta_per", sep=""))))
    eval(parse(text=(paste("Monte_r_total",sce,"[,",k,"]","<-","r_total", sep=""))))
    eval(parse(text=(paste("Monte_dsr",sce,"[,",k,"]","<-","dsr", sep=""))))
    eval(parse(text=(paste("Monte_L_K",sce,"[,",k,"]","<-","L_K", sep=""))))
    eval(parse(text=(paste("Monte_L",sce,"[,",k,"]","<-","L", sep=""))))
    eval(parse(text=(paste("Monte_E",sce,"[,",k,"]","<-","E", sep=""))))
    eval(parse(text=(paste("Monte_SEC",sce,"[,",k,"]","<-","SEC", sep=""))))
    eval(parse(text=(paste("Monte_CO_GOV",sce,"[,",k,"]","<-","CO_GOV", sep=""))))
    eval(parse(text=(paste("Monte_I",sce,"[,",k,"]","<-","I", sep=""))))
    eval(parse(text=(paste("Monte_def_per",sce,"[,",k,"]","<-","def_per", sep=""))))
    eval(parse(text=(paste("Monte_lev_B",sce,"[,",k,"]","<-","lev_B", sep=""))))
    eval(parse(text=(paste("Monte_PORT_BCH",sce,"[,",k,"]","<-","PORT_BCH", sep=""))))
    eval(parse(text=(paste("Monte_PORT_BGH",sce,"[,",k,"]","<-","PORT_BGH", sep=""))))
    eval(parse(text=(paste("Monte_p_C",sce,"[,",k,"]","<-","p_C", sep=""))))
    eval(parse(text=(paste("Monte_p_G",sce,"[,",k,"]","<-","p_G", sep=""))))
    eval(parse(text=(paste("Monte_h",sce,"[,",k,"]","<-","h", sep=""))))
    eval(parse(text=(paste("Monte_EMIS",sce,"[,",k,"]","<-","EMIS", sep=""))))
    eval(parse(text=(paste("Monte_T_AT",sce,"[,",k,"]","<-","T_AT", sep=""))))
    eval(parse(text=(paste("Monte_dep_E",sce,"[,",k,"]","<-","dep_E", sep=""))))
    eval(parse(text=(paste("Monte_dep_M",sce,"[,",k,"]","<-","dep_M", sep=""))))
    eval(parse(text=(paste("Monte_CO2_CUM",sce,"[,",k,"]","<-","CO2_CUM", sep=""))))
    eval(parse(text=(paste("Monte_rho",sce,"[,",k,"]","<-","rho", sep=""))))
    eval(parse(text=(paste("Monte_mu",sce,"[,",k,"]","<-","mu", sep=""))))
    eval(parse(text=(paste("Monte_epsilon",sce,"[,",k,"]","<-","epsilon", sep=""))))
    eval(parse(text=(paste("Monte_W",sce,"[,",k,"]","<-","W", sep=""))))
    eval(parse(text=(paste("Monte_W_POP_ratio",sce,"[,",k,"]","<-","W_POP_ratio", sep=""))))
    eval(parse(text=(paste("Monte_hazratio",sce,"[,",k,"]","<-","hazratio", sep=""))))
    eval(parse(text=(paste("Monte_CAR_per",sce,"[,",k,"]","<-","CAR_per", sep=""))))
    eval(parse(text=(paste("Monte_CR",sce,"[,",k,"]","<-","CR", sep=""))))
    eval(parse(text=(paste("Monte_CR_C_S1",sce,"[,",k,"]","<-","CR_C_S1", sep=""))))
    eval(parse(text=(paste("Monte_CR_G",sce,"[,",k,"]","<-","CR_G", sep=""))))
    eval(parse(text=(paste("Monte_spr",sce,"[,",k,"]","<-","spr", sep=""))))
    eval(parse(text=(paste("Monte_spr_C_S1",sce,"[,",k,"]","<-","spr_C_S1", sep=""))))
    eval(parse(text=(paste("Monte_spr_G",sce,"[,",k,"]","<-","spr_G", sep=""))))
    eval(parse(text=(paste("Monte_SEC_Y_per",sce,"[,",k,"]","<-","SEC_Y_per", sep=""))))
    eval(parse(text=(paste("Monte_fiscal_balance_per",sce,"[,",k,"]","<-","fiscal_balance_per", sep=""))))
    eval(parse(text=(paste("Monte_fiscal_balance_help_per",sce,"[,",k,"]","<-","fiscal_balance_help_per", sep=""))))
    eval(parse(text=(paste("Monte_firms_balance_per",sce,"[,",k,"]","<-","firms_balance_per", sep=""))))
    eval(parse(text=(paste("Monte_households_balance_per",sce,"[,",k,"]","<-","households_balance_per", sep=""))))
    eval(parse(text=(paste("Monte_banks_balance_per",sce,"[,",k,"]","<-","banks_balance_per", sep=""))))
    eval(parse(text=(paste("Monte_illiq",sce,"[,",k,"]","<-","illiq", sep=""))))
    eval(parse(text=(paste("Monte_LF",sce,"[,",k,"]","<-","LF", sep=""))))
    eval(parse(text=(paste("Monte_lf",sce,"[,",k,"]","<-","lf", sep=""))))
    eval(parse(text=(paste("Monte_yield_C",sce,"[,",k,"]","<-", "yield_C", sep=""))))
    eval(parse(text=(paste("Monte_yield_G",sce,"[,",k,"]","<-", "yield_G", sep=""))))
    eval(parse(text=(paste("Monte_yield_C_per",sce,"[,",k,"]","<-", "yield_C_per", sep=""))))
    eval(parse(text=(paste("Monte_yield_G_per",sce,"[,",k,"]","<-", "yield_G_per", sep=""))))
    
    for (i in 1:T){
      
      eval(parse(text=(paste("g_Y_per_mean", sce, "[",i,"]", "<-", "mean(Monte_g_Y_per", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("epsilon_ratio_mean", sce, "[",i,"]", "<-", "mean(Monte_epsilon_ratio", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("mu_ratio_mean", sce, "[",i,"]", "<-", "mean(Monte_mu_ratio", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("IG_PRI_mean", sce, "[",i,"]", "<-", "mean(Monte_IG_PRI", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("IG_GOV_mean", sce, "[",i,"]", "<-", "mean(Monte_IG_GOV", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("IG_E_Y_mean", sce, "[",i,"]", "<-", "mean(Monte_IG_E_Y", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("SUB_mean", sce, "[",i,"]", "<-", "mean(Monte_SUB", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("Y_sdplus", sce, "[",i,"]", "<-", "mean(Monte_Y", sce, "[",i,",])", "+","sd(Monte_Y", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("Y_sdminus", sce, "[",i,"]", "<-", "mean(Monte_Y", sce, "[",i,",])", "-","sd(Monte_Y", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("CR_C_S1_sdplus", sce, "[",i,"]", "<-", "mean(Monte_CR_C_S1", sce, "[",i,",])", "+","sd(Monte_CR_C_S1", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("CR_C_S1_sdminus", sce, "[",i,"]", "<-", "mean(Monte_CR_C_S1", sce, "[",i,",])", "-","sd(Monte_CR_C_S1", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("EMIS_sdplus", sce, "[",i,"]", "<-", "mean(Monte_EMIS", sce, "[",i,",])", "+","sd(Monte_EMIS", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("EMIS_sdminus", sce, "[",i,"]", "<-", "mean(Monte_EMIS", sce, "[",i,",])", "-","sd(Monte_EMIS", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("r_total_sdplus", sce, "[",i,"]", "<-", "mean(Monte_r_total", sce, "[",i,",])", "+","sd(Monte_r_total", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("r_total_sdminus", sce, "[",i,"]", "<-", "mean(Monte_r_total", sce, "[",i,",])", "-","sd(Monte_r_total", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("def_per_sdplus", sce, "[",i,"]", "<-", "mean(Monte_def_per", sce, "[",i,",])", "+","sd(Monte_def_per", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("def_per_sdminus", sce, "[",i,"]", "<-", "mean(Monte_def_per", sce, "[",i,",])", "-","sd(Monte_def_per", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("lev_B_sdplus", sce, "[",i,"]", "<-", "mean(Monte_lev_B", sce, "[",i,",])", "+","sd(Monte_lev_B", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("lev_B_sdminus", sce, "[",i,"]", "<-", "mean(Monte_lev_B", sce, "[",i,",])", "-","sd(Monte_lev_B", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("ur_per_mean", sce, "[",i,"]", "<-", "mean(Monte_ur_per", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("Y_mean", sce, "[",i,"]", "<-", "mean(Monte_Y", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("CO_PRI_mean", sce, "[",i,"]", "<-", "mean(Monte_CO_PRI", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("I_PRI_mean", sce, "[",i,"]", "<-", "mean(Monte_I_PRI", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("lambda_mean", sce, "[",i,"]", "<-", "mean(Monte_lambda", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("N_mean", sce, "[",i,"]", "<-", "mean(Monte_N", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("Y_POP_ratio_mean", sce, "[",i,"]", "<-", "mean(Monte_Y_POP_ratio", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("lambda_perworker_mean", sce, "[",i,"]", "<-", "mean(Monte_lambda_perworker", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("EMIS_IN_Y_mean", sce, "[",i,"]", "<-", "mean(Monte_EMIS_IN_Y", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("EMIS_IN_mean", sce, "[",i,"]", "<-", "mean(Monte_EMIS_IN", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("E_N_ratio_mean", sce, "[",i,"]", "<-", "mean(Monte_E_N_ratio", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("I_K_ratio_mean", sce, "[",i,"]", "<-", "mean(Monte_I_K_ratio", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("ur_per_mean", sce, "[",i,"]", "<-", "mean(Monte_ur_per", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("POP_mean", sce, "[",i,"]", "<-", "mean(Monte_POP", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("beta_S1_mean", sce, "[",i,"]", "<-", "mean(Monte_beta_S1", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("greenI_mean", sce, "[",i,"]", "<-", "mean(Monte_greenI", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("kappa_mean", sce, "[",i,"]", "<-", "mean(Monte_kappa", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("theta_per_mean", sce, "[",i,"]", "<-", "mean(Monte_theta_per", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("r_total_mean", sce, "[",i,"]", "<-", "mean(Monte_r_total", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("dsr_mean", sce, "[",i,"]", "<-", "mean(Monte_dsr", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("L_K_mean", sce, "[",i,"]", "<-", "mean(Monte_L_K", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("L_mean", sce, "[",i,"]", "<-", "mean(Monte_L", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("E_mean", sce, "[",i,"]", "<-", "mean(Monte_E", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("SEC_mean", sce, "[",i,"]", "<-", "mean(Monte_SEC", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("CO_GOV_mean", sce, "[",i,"]", "<-", "mean(Monte_CO_GOV", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("I_mean", sce, "[",i,"]", "<-", "mean(Monte_I", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("def_per_mean", sce, "[",i,"]", "<-", "mean(Monte_def_per", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("lev_B_mean", sce, "[",i,"]", "<-", "mean(Monte_lev_B", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("PORT_BCH_mean", sce, "[",i,"]", "<-", "mean(Monte_PORT_BCH", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("PORT_BGH_mean", sce, "[",i,"]", "<-", "mean(Monte_PORT_BGH", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("p_C_mean", sce, "[",i,"]", "<-", "mean(Monte_p_C", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("p_G_mean", sce, "[",i,"]", "<-", "mean(Monte_p_G", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("h_mean", sce, "[",i,"]", "<-", "mean(Monte_h", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("EMIS_mean", sce, "[",i,"]", "<-", "mean(Monte_EMIS", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("T_AT_mean", sce, "[",i,"]", "<-", "mean(Monte_T_AT", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("dep_E_mean", sce, "[",i,"]", "<-", "mean(Monte_dep_E", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("dep_M_mean", sce, "[",i,"]", "<-", "mean(Monte_dep_M", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("CO2_CUM_mean", sce, "[",i,"]", "<-", "mean(Monte_CO2_CUM", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("rho_mean", sce, "[",i,"]", "<-", "mean(Monte_rho", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("mu_mean", sce, "[",i,"]", "<-", "mean(Monte_mu", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("epsilon_mean", sce, "[",i,"]", "<-", "mean(Monte_epsilon", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("W_mean", sce, "[",i,"]", "<-", "mean(Monte_W", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("W_POP_ratio_mean", sce, "[",i,"]", "<-", "mean(Monte_W_POP_ratio", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("hazratio_mean", sce, "[",i,"]", "<-", "mean(Monte_hazratio", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("CAR_per_mean", sce, "[",i,"]", "<-", "mean(Monte_CAR_per", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("CR_mean", sce, "[",i,"]", "<-", "mean(Monte_CR", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("CR_C_S1_mean", sce, "[",i,"]", "<-", "mean(Monte_CR_C_S1", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("CR_G_mean", sce, "[",i,"]", "<-", "mean(Monte_CR_G", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("spr_mean", sce, "[",i,"]", "<-", "mean(Monte_spr", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("spr_C_S1_mean", sce, "[",i,"]", "<-", "mean(Monte_spr_C_S1", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("spr_G_mean", sce, "[",i,"]", "<-", "mean(Monte_spr_G", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("SEC_Y_per_mean", sce, "[",i,"]", "<-", "mean(Monte_SEC_Y_per", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("fiscal_balance_per_mean", sce, "[",i,"]", "<-", "mean(Monte_fiscal_balance_per", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("fiscal_balance_help_per_mean", sce, "[",i,"]", "<-", "mean(Monte_fiscal_balance_help_per", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("firms_balance_per_mean", sce, "[",i,"]", "<-", "mean(Monte_firms_balance_per", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("households_balance_per_mean", sce, "[",i,"]", "<-", "mean(Monte_households_balance_per", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("banks_balance_per_mean", sce, "[",i,"]", "<-", "mean(Monte_banks_balance_per", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("illiq_mean", sce, "[",i,"]", "<-", "mean(Monte_illiq", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("LF_mean", sce, "[",i,"]", "<-", "mean(Monte_LF", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("lf_mean", sce, "[",i,"]", "<-", "mean(Monte_lf", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("yield_C_mean", sce, "[",i,"]", "<-", "mean(Monte_yield_C", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("yield_G_mean", sce, "[",i,"]", "<-", "mean(Monte_yield_G", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("yield_C_per_mean", sce, "[",i,"]", "<-", "mean(Monte_yield_C_per", sce, "[",i,",])", sep=""))))
      eval(parse(text=(paste("yield_G_per_mean", sce, "[",i,"]", "<-", "mean(Monte_yield_G_per", sce, "[",i,",])", sep=""))))
      
      
    }
    
  }
  
  #-##########################
  #5. PRINTING THE RESULTS IN R
  ############################
  
  if (sce<(scenarios_print+1)) {
    
    var1<-sprintf ("%.3f", g_Y[2])
    var2<-sprintf ("%.3f", g_Y[3])
    var3<-sprintf ("%.3f", g_Y[33])
    var4<-sprintf ("%.3f", ur[33])
    var5<-sprintf ("%.2f", POP[33])
    var6<-sprintf ("%.2f", lf[33])
    var7<-sprintf ("%.2f", omega_ratio[33])
    var8<-sprintf ("%.2f", mu_ratio[33])
    var9<-sprintf ("%.2f", epsilon_ratio[33])
    var10<-sprintf ("%.2f", rho_ratio[33])
    var11<-sprintf ("%.2f", theta[33])
    var12<-sprintf ("%.2f", theta[83])
    var13<-sprintf ("%.4f", KGE_KCE[32])
    var14<-sprintf ("%.4f", KGNE_KCNE[32])
    var15<-sprintf ("%.4f", KSEQ_KCE[32])
    var16<-sprintf ("%.5f", sigma_0_optimal[2])
    var17<-sprintf ("%.2f", IG_cum[33]/33)
    var18<-sprintf ("%.2f", IC_cum[33]/33) 
    var19<-sprintf ("%.2f", B_GCB[5])
    var20<-sprintf ("%.2f", E_ratio[33])
    var21<-sprintf ("%.3f", CO2_ratio[33])
    var22<-sprintf ("%.2f", T_AT[33])
    var23<-sprintf ("%.3f", (g_Y_cum[33])/33)
    var24<-sprintf ("%.3f", (TAX_C_change[6]))
    var25<-sprintf ("%.3f", (SUB_change[6]))
    var26<-sprintf ("%.3f", (IG_GOV_change[5]))
    var27<-sprintf ("%.3f", (seq[33]))
    var28<-sprintf ("%.3f", ((g_Y_cum[83])-g_Y_cum[33])/50)
    var29<-sprintf ("%.2f", epsilon_ratio[83])
    var30<-sprintf ("%.3f", (seq[83]))
    var31<-sprintf ("%.2f", E_ratio[83])
    var32<-sprintf ("%.3f", CO2_ratio[83])
    var33<-sprintf ("%.2f", T_AT[83])
    var34<-sprintf ("%.1f", POP[83])
    var35<-sprintf ("%.2f", omega_ratio[83])
    var36<-sprintf ("%.1f", EMIS_L[33])
    var37<-sprintf ("%.1f", EMIS_L[83])
    var38<-sprintf ("%.2f", IG_E_Y_cum[33]/33) 
    
    cat ("SCENARIO", sce, " ")
    cat ("\n")
    cat ("g_Y(2019):", var1," ")
    cat ("g_Y(2020):", var2," ")
    cat ("g_Y(2050):", var3," ")
    cat ("g_Y_mean(2018-2050):", var23," ")
    cat ("g_Y_mean(2051-2100):", var28," ")
    cat ("\n")
    cat ("ur(2050):", var4," ")
    cat ("POP(2050):", var5," ")
    cat ("POP(2100):", var34," ")
    cat ("lf(2050):", var6," ")
    cat ("sigma_0_optimal:", var16," ")
    cat ("\n")
    cat ("omega_rat(2050):", var7," ")
    cat ("omega_rat(2100):", var35," ")
    cat ("mu_rat(2050):", var8," ")
    cat ("rho_rat(2050):", var10," ")
    cat ("epsilon_rat(2050):", var9," ")
    cat ("epsilon_rat(2100):", var29," ")
    cat ("theta(2050):", var11," ")
    cat ("theta(2100):", var12," ")
    cat ("seq(2050):", var27," ")
    cat ("seq(2100):", var30," ")
    cat ("\n")
    cat ("KGE_KCE(2049):", var13," ")
    cat ("KGNE_KCNE(2049):", var14," ")
    cat ("KSEQ_KCE(2049):", var15," ")
    cat ("IG_annual(2018-2050):", var17," ")
    cat ("IC_annual(2018-2050):", var18," ")
    cat ("IG_E_Y_annual(2018-2050):", var38," ")
    cat ("\n")
    cat ("E_rat(2050):", var20," ")
    cat ("E_rat(2100):", var31," ")
    cat ("CO2_rat(2050):", var21," ")
    cat ("CO2_rat(2100):", var32," ")
    cat ("EMIS_L(2050):", var36," ")
    cat ("EMIS_L(2100):", var37," ")
    cat ("T_AT(2050):", var22," ")
    cat ("T_AT(2100):", var33," ")
    cat ("\n")
    cat ("B_GCB(2050):", var19," ")
    cat ("TAX_C_change(2022):", var24," ")
    cat ("SUB_change(2022):", var25," ")
    cat ("IG_GOV_change(2022):", var26," ")
    
    cat ("\n")
    cat ("\n")
    
  }
  
  ####################################################
  #6. EXCEL FILES WITH VARIABLES AND PARAMETER VALUES
  ####################################################       
  
  matrixname<-paste("Macroeconomy",sce, sep="")
  assign (matrixname, (round(cbind(SEC_CB, SEC_CBred, A_K_ratio, HPM_K_ratio, D_K_ratio, A, A_N, HPM, D, V_HF, B_GCB, B_CCB, g_Y, gI, gC, g_K, gGOVCO, ID_K_ratio, I_K_ratio, I_PRI_D, I_PRI_D_S1, I_PRI_D_S2, I_PRI_D_S3, I_PRI_D_S4, I_PRI, I_PRI_S1, I_PRI_S2, I_PRI_S3, I_PRI_S4, SEC, SEC_H, SEC_B, SEC_Y, x_1, x_2, PORT_BCH, PORT_BGH, PORT_SECH, PORT_D, B, B_G, B_GH, B_C, B_CH, b_C, b_G, p_C, p_G, HPM, D, A, LC, LC_S1, LC_S2, LC_S3, LC_S4, LG, LG_S1, LG_S2, LG_S3, LG_S4, CAP, sh_NLG, sh_NLC_S1, sh_NLC_S2, sh_NLC_S3, sh_NLC_S4, CR, CR_C_S1, CR_C_S2, CR_C_S3, CR_C_S4, CR_G, dsr, def, illiq, L_K, I, lev_B, CAR, fiscal_balance, TAX, TAX_H, TAX_F, TAX_C, CO_GOV, D_T, D_TP, D_TF,  delta, LF, lf, POP, g_POP, um, ue, u, re, v, Y_POP_ratio, g_lambda, g_EN_ratio, g_YE_ratio, lambda, lambda_perworker, E_N_ratio, Y_E_ratio, CO_PRI, I_PRI, gI, theta, kappa, beta, beta_S1, beta_S2, beta_S3, beta_S4, Y, Wbill, Interest, Depreciation,  Y_star, Y_star_M, Y_star_E, Y_star_N, Y_star_K, Y_H, TP, RP, DP, L, K_PRI, KC_PRI, KC_PRI_S1, KC_PRI_S2, KC_PRI_S3, KC_PRI_S4, KG_PRI, KG_PRI_S1, KG_PRI_S2, KG_PRI_S3, KG_PRI_S4, L_K, r, r_total, IC_PRI_D, IC_PRI_D_S1, IC_PRI_D_S2, IC_PRI_D_S3, IC_PRI_D_S4, IC_PRI, IC_PRI_S1, IC_PRI_S2, IC_PRI_S3, IC_PRI_S4, IG_PRI_D, IG_PRI_D_S1, IG_PRI_D_S2, IG_PRI_D_S3, IG_PRI_D_S4,IG_PRI, IG_PRI_S1, IG_PRI_S2, IG_PRI_S3, IG_PRI_S4, NLC_D_S1, NLC_D_S2, NLC_D_S3, NLC_D_S4, NLG_D_S1, NLG_D_S2, NLG_D_S3, NLG_D_S4,  BP, N, Y_K_ratio, C_K_ratio, I_Y_ratio, Y_H_Y_ratio, Y_H_K_ratio, C_Y_ratio, sigma_0, sigma_0_optimal, h_optimal, h, g_h, g_Y_per, ur_per, greenI, p_C, yield_C,V_CB,  g_pC, g_bC, B_C_issue, p_G, yield_G, g_pG, g_bG, B_G_issue, IG_cum, IC_cum, w_G, w_C_S1, w_C_S2, w_C_S3, w_C_S4,BAILOUT, BAILOUT_CAR, BAILOUT_levB, int_G, int_C_S1, int_C_S2, int_C_S3, int_C_S4, spr_G, spr_C_S1, spr_C_S2, spr_C_S3, spr_C_S4, w_LT, sh_LG, sh_LC_S1, sh_LC_S2, sh_LC_S3, sh_LC_S4, KG_GOV_Y_ratio, KC_GOV_Y_ratio, KG_GOV, KC_GOV , IG_GOV, IC_GOV, IG_GOV_Y_ratio, IC_GOV_Y_ratio, KC_Y_ratio, KG_Y_ratio, KG_K_ratio, GVA_S1, GVA_S2, GVA_S3, GVA_S4, sh_EMIS_IN_S1, sh_EMIS_IN_S2, sh_EMIS_IN_S3, sh_EMIS_IN_S4, dd_S1, dd_S2, dd_S3, dd_S4, LG_L_ratio, LG_L_ratio_pseudo, tau_C, spr_1_pseudo, spr_1_pseudo_case_I, spr_1_pseudo_case_II, spr_2_pseudo, r_2_dCR_ddsr_max_pseudo, r_3_dCR_dCAR_max_pseudo, r_3_dCR_dCAR_max_pseudo_case_I, r_3_dCR_dCAR_max_pseudo_case_II, sh_GREEN_S1, sh_GREEN_S2, sh_GREEN_S3, sh_GREEN_S4, seq, EMIS_IN_Y, SUB, gov_SUB, tau_C,tau_C_baseline, tau_C_policy_High, tau_C_policy_Medium, omega, EMIS_IN_Y_SSP360, endog, KGE_PRI_S1, KGE_PRI_S2, KGE_PRI_S3, KGE_PRI_S4, KGE, KGNE_PRI_S1, KGNE_PRI_S2, KGNE_PRI_S3, KGNE_PRI_S4, KGNE, KCE_PRI_S1, KCE_PRI_S2, KCE_PRI_S3, KCE_PRI_S4, KCE, KCNE_PRI_S1, KCNE_PRI_S2, KCNE_PRI_S3, KCNE_PRI_S4, KCNE, KSEQ, KSEQ_PRI_S1, KSEQ_PRI_S2, seq, fiscal_balance_help, private_balance, firms_balance, households_balance, banks_balance, fiscal_balance_help_per, private_balance_per, firms_balance_per, households_balance_per, banks_balance_per, g_POP_stag, g_ucn, g_ucr, tucn, tucr, ucn, ucr, IG_E, IG_E_Y, IG_E_Y_cum), digits=8)))
  
  matrixname<-paste("Ecosystem",sce, sep="")
  assign (matrixname, (round(cbind(T_AT, CO2_CUM, dep_E, dep_M, REV_E, REV_M, D_T, hazratio, HW_CUM, EMIS, theta, M, REC, CEN, O2, MY, W, W_POP_ratio, DEM, EMIS_IN, EMIS_L, g_EMIS_L, E, E_NF, E_F, ED, mu, rho, omega, epsilon, omega_ratio, mu_ratio, rho_ratio, epsilon_ratio, theta_per, haz_flow), digits=4)))
  
  matrixname<-paste("Variables",sce, sep="")
  assign (matrixname, (round(cbind(A, B, BAILOUT, B_C, b_C, B_CCB, b_CCB, B_CH, b_CH, B_G, b_G, B_GCB, b_GCB, B_GH, b_GH, BP, BP_D, BP_U, CO_GOV, CO_PRI, CO_PRI_N, CAP, CAR, CBP, CEN, CO2_CUM, CON_E, CON_M, coupon_C, coupon_G, CR, CR_C_S1, CR_C_S2, CR_C_S3, CR_C_S4, CR_G, D, DC, def, DEM, dep_E, dep_M, DL, DP, dsr, D_T, D_TF, D_TP, E, ED, EMIS, EMIS_IN, EMIS_L, zero, E_F, E_NF, GNS, g_EMIS_L, g_POP, gov_SUB, g_ucn, g_ucr, g_x_20, g_Y, g_beta_0, g_lambda, g_lambda_30, g_omega, hazratio, HPM, HW_CUM, I_GOV, I_PRI, IC_GOV,IC_PRI, IC_PRI_S1, IC_PRI_S2, IC_PRI_S3, IC_PRI_S4, IG_GOV, IG_PRI, IG_PRI_S1, IG_PRI_S2, IG_PRI_S3, IG_PRI_S4, I_PRI_D, I_PRI_D_S1, I_PRI_D_S2, I_PRI_D_S3, I_PRI_D_S4, IC_PRI_D_S1, IC_PRI_D_S2, IC_PRI_D_S3, IC_PRI_D_S4, IG_PRI_D_S1, IG_PRI_D_S2, IG_PRI_D_S3, IG_PRI_D_S4,  illiq, int_C_S1, zero, int_C_S2, int_C_S3, int_C_S4, int_G, K, K_GOV, K_PRI, KC, KC_GOV, KC_PRI, KC_PRI_S1, KC_PRI_S2, KC_PRI_S3, KC_PRI_S4, KCE, KCE_PRI_S1, KCE_PRI_S2, KCE_PRI_S3, KCE_PRI_S4, KCNE, KCNE_PRI_S1, KCNE_PRI_S2, KCNE_PRI_S3, KCNE_PRI_S4, KG, KG_GOV, KG_PRI, KG_PRI_S1, KG_PRI_S2, KG_PRI_S3, KG_PRI_S4, KGE, KGE_PRI_S1, KGE_PRI_S2, KGE_PRI_S3, KGE_PRI_S4, KGNE, KGNE_PRI_S1, KGNE_PRI_S2, KGNE_PRI_S3, KGNE_PRI_S4, KSEQ, KSEQ_PRI_S1, KSEQ_PRI_S2, L, LC, LC_S1, LC_S2, LC_S3, LC_S4, LG, LG_S1, LG_S2, LG_S3, LG_S4, lev_B, LF, lf_1, zero,  M, MY, N, NLC_D_S1, NLC_D_S2, NLC_D_S3, NLC_D_S4, NLG_D_S1, NLG_D_S2, NLG_D_S3, NLG_D_S4, O2, p_C, p_G, POP, r, re, REC, RES_E, RES_M, REV_E, REV_M, RP, SEC, SEC_B, SEC_CB, SEC_H, seq, SES, sh_NLG, sh_NLC_S1, sh_NLC_S2, sh_NLC_S3, sh_NLC_S4, sh_L, sh_LC_S1, sh_LC_S2, sh_LC_S3, sh_LC_S4, sh_LG, spr, spr_G, spr_C_S1, spr_C_S2, spr_C_S3, spr_C_S4, SUB, TAX, T_AT, TAX_C, TAX_F, zero,TAX_H, TP, TP_G, tucn, tucr, u, ucn, ucr, ue, um, ur, v, V_CB, V_H, V_HF, wage, W, w_C_S1, w_C_S2, w_C_S3, w_C_S4, w_G, w_LT, x_1, x_2, x_20, Y, Y_star, Y_star_E, Y_H, Y_HD, Y_HG,  yield_C, yield_G, Y_star_K, Y_star_M, Y_star_N, beta_S1, beta_S2, beta_S3, beta_S4, beta_0_S1, beta_0_S2, beta_0_S3, beta_0_S4, delta, epsilon, theta, kappa, lambda, lambda_30, mu, rho, sigma_0, tau_C, omega), digits=6)))
  
  matrixname<-paste("Parameters",sce, sep="")
  assign (matrixname, (round(cbind(ad_K, ad_LF, ad_P, c_1, c_2, CAR_min, con_E, con_M, CR_max, dd_S1, dd_S2, dd_S3, dd_S4, def_max, def_0, def_1, def_2, gov_C, gov_IC, gov_IG, h, h_1, h_2, haz, int_A, int_D, int_S, l_1, lev_B_max, lf_2, p, p_C_bar, p_G_bar, pr, r_0, r_1, r_2, r_3, zero, rep, s_B, s_C, s_F, s_G, s_W, sh_EMIS_IN_S1, sh_EMIS_IN_S2, sh_EMIS_IN_S3, sh_EMIS_IN_S4, sh_GREEN_S1, sh_GREEN_S2, sh_GREEN_S3, sh_GREEN_S4, sh_GVA_S1, sh_GVA_S2, sh_GVA_S3, sh_GVA_S4, spr_0, spr_1, spr_2, spr_3, t_1, t_2, w_H, w_S, x_10, x_11, x_21, alpha_00, alpha_01, alpha_1, alpha_2, alpha_31, alpha_32, alpha_41, alpha_42, alpha_51, alpha_52, beta_1, beta_2, zero, gamma_E, gamma_E1, gamma_E2, gamma_E3, gamma_E4, gamma_SEQ1, gamma_SEQ2, delta_0, epsilon_max, epsilon_min, zeta_1, zeta_2, zeta_3, zeta_4, zeta_5, zeta_6, zeta_7, zeta_8, zeta_9, zeta_10, eta_1, eta_2, eta_3, lambda_10, lambda_10prime, lambda_11, lambda_12, lambda_13, lambda_14, lambda_15, lambda_20, lambda_20prime, lambda_21, lambda_22, lambda_23, lambda_24, lambda_25, lambda_30prime, lambda_31, lambda_32, lambda_33, lambda_34, lambda_35, lambda_40, lambda_40prime, lambda_41, lambda_42, lambda_43, lambda_44, lambda_45, mu_max, mu_min, xi, pi_1, pi_2, pi_3, pi_4, zero, pi_5, pi_6, pi_7, pi_8, pi_9, pi_10, rho_max, sigma_1, sigma_2, tau_F, tau_H, phi), digits=8)))
  
  write.csv(get(paste0("Macroeconomy",sce)), paste0("DEFINE_Macroeconomy",sce,".csv"))
  write.csv(get(paste0("Ecosystem",sce)), paste0("DEFINE_Ecosystem",sce,".csv"))
  write.csv(get(paste0("Variables",sce)), paste0("DEFINE_Variables",sce,".csv"))
  write.csv(get(paste0("Parameters",sce)), paste0("DEFINE_Parameters",sce,".csv"))
  
  sce=sce+1
  
}


#######################
#7. FIGURES
#######################
#Font
font="Helvetica"
windowsFonts(Selectedfont=windowsFont(font))
savefont<-par(family="Selectedfont")
#Positions and labels of tick marks
xlabel1<-c(1, 13, 23, 33, 43, 53, 63, 73, 83)  
xlabel2<-c("2018", "2030","2040","2050","2060", "2070", "2080", "2090", "2100") 
#Names/colours/lines for scenarios
height_graph=800
width_graph=1.625*height_graph
res_graph=0.3* height_graph
Colour1="black"    
Colour2="darkgreen"     
Colour3= "darkred"   
Colour4a="darkblue"
Colour4b="transparent"
Colour_range1="#C0C0C0" 
Colour_range2="#04b976"
Line1=1
Line2=1
Line3=1
Line4=1
Line_range1=1
Line_range2=1
Width1=1
Width2=0.5
Width3=0.5
Width4=0.5
pch1=NA
pch2=16
pch3=4
pch4=6
Width_range1=3.5
Width_range2=3.5
Scenario1="Baseline"
Scenario2="Baseline without damages (mean)" 
Range1= "Baseline (+/- 1 st. dev.)"
Range2= "Baseline without damages (+/- 1 st. dev.)"
Mean="Baseline (mean)"
Balance1="Government balance"
Balance2="Firm balance"
Balance3="Household balance"
Balance4="Bank balance"
Carbon_Baseline="Baseline (SSP3 6.0 W/m^2)"
Carbon_M="SSP3 4.5 W/m^2"


library(Cairo)
dev.off()
Figures<-function(Figure,x2,x3,x4, Sce2, Sce3, Sce4, Colour4){
  
  #Figure
  Cairo(file=paste("DEFINE-Fig",Figure,"-1-Growth_rate_of_output.jpeg", sep=""), width=width_graph, height=height_graph, res=res_graph)
  par(mar = c(3,3.3,1,1)+0.1,mgp=c(1.8,0.5,0), family=font)
  plot(g_Y_per_mean1, type="l", col=Colour1, lwd=Width1, ylab="Growth rate of output (%)",  ylim=c(0,5), yaxp=c(0,5,5), xlab=" ", xaxt="n", cex.lab=0.8, cex.axis=0.7)
  axis(side=1, at=xlabel1, labels=xlabel2, cex.lab=0.8, cex.axis=0.7)
  lines(eval(parse(text=(paste("g_Y_per_mean",x2, sep="")))), col=Colour2, pch=pch2, type="o", lwd=Width2, cex=0.4) 
  lines(eval(parse(text=(paste("g_Y_per_mean",x3, sep="")))), col=Colour3, pch=pch3, type="o", lwd=Width3, cex=0.4) 
  lines(eval(parse(text=(paste("g_Y_per_mean",x4, sep="")))), col=Colour4, pch=pch4, type="o", lwd=Width4, cex=0.4) 
  lines(g_Y_per_mean1, lty=Line1, col=Colour1, lwd=Width1)
  legend("topright", legend=c(Scenario1, Sce2, Sce3, Sce4), lty=c(Line1, Line2, Line3, Line4), col=c(Colour1, Colour2, Colour3, Colour4), pch=c(pch1, pch2, pch3, pch4), bty="n", cex=0.7) 
  dev.off() 
  
  #Figure
  Cairo(file=paste("DEFINE-Fig",Figure,"-2-Share_of_non-fossil.jpeg", sep=""), width=width_graph, height=height_graph, res=res_graph)
  par(mar = c(3,3.3,1,1)+0.1,mgp=c(1.8,0.5,0), family=font)
  plot(theta_per_mean1, type="l", col=Colour1, lwd=Width1, ylab=" Share of non-fossil energy (%)", ylim=c(10, 50), yaxp=c(10,50,4), xlab=" ", xaxt="n", cex.lab=0.8, cex.axis=0.7)
  axis(side=1,at=xlabel1, labels=xlabel2, cex.lab=0.8, cex.axis=0.7)
  lines(eval(parse(text=(paste("theta_per_mean",x2, sep="")))), col=Colour2, pch=pch2, type="o", lwd=Width2, cex=0.4) 
  lines(eval(parse(text=(paste("theta_per_mean",x3, sep="")))), col=Colour3, pch=pch3, type="o", lwd=Width3,cex=0.4) 
  lines(eval(parse(text=(paste("theta_per_mean",x4, sep="")))), col=Colour4, pch=pch4, type="o", lwd=Width4, cex=0.4) 
  lines(theta_per_mean1, lty=Line1, col=Colour1, lwd=Width1)
  legend("topleft", legend=c(Scenario1, Sce2, Sce3, Sce4), lty=c(Line1, Line2, Line3, Line4), col=c(Colour1, Colour2, Colour3, Colour4), pch=c(pch1, pch2, pch3, pch4), bty="n", cex=0.7) 
  dev.off() 
  
  
  #Figure
  Cairo(file=paste("DEFINE-Fig",Figure,"-3-Emissions.jpeg", sep=""), width=width_graph, height=height_graph, res=res_graph)
  par(mar = c(3,3.3,1,1)+0.1,mgp=c(1.8,0.5,0), family=font)
  plot(EMIS_mean1, type="l", col=Colour1, lwd=Width1, ylab=expression("CO" [2]*" emissions (GtCO" [2]*"/year)"), ylim=c(0,60), yaxp=c(0,60,6), xlab=" ", xaxt="n", cex.lab=0.8, cex.axis=0.7)
  axis(side=1,at=xlabel1, labels=xlabel2, cex.lab=0.8, cex.axis=0.7)
  lines(eval(parse(text=(paste("EMIS_mean",x2, sep="")))), col=Colour2, pch=pch2, type="o", lwd=Width2, cex=0.4) 
  lines(eval(parse(text=(paste("EMIS_mean",x3, sep="")))), col=Colour3, pch=pch3, type="o", lwd=Width3,cex=0.4) 
  lines(eval(parse(text=(paste("EMIS_mean",x4, sep="")))), col=Colour4, pch=pch4, type="o", lwd=Width4, cex=0.4) 
  lines(EMIS_mean1, lty=Line1, col=Colour1, lwd=Width1)
  legend("bottomleft", legend=c(Scenario1, Sce2, Sce3, Sce4), lty=c(Line1, Line2, Line3, Line4), col=c(Colour1, Colour2, Colour3, Colour4), pch=c(pch1, pch2, pch3, pch4), bty="n", cex=0.7) 
  dev.off() 
  
  #Figure 
  Cairo(file=paste("DEFINE-Fig",Figure,"-4-Temperature.jpeg", sep=""), width=width_graph, height=height_graph, res=res_graph)
  par(mar = c(3,3.3,1,1)+0.1,mgp=c(1.8,0.5,0), family=font)
  plot(T_AT_mean1, type="l", col=Colour1, lwd=Width1, ylab=expression("Temperature change ("^{o}*"C)"), ylim=c(1,4), yaxp=c(1,4,6),  cex.lab=1.2, xlab=" ", xaxt="n", cex.lab=0.8, cex.axis=0.7)
  axis(side=1, at=xlabel1, labels=xlabel2, cex.lab=0.8, cex.axis=0.7)
  lines(eval(parse(text=(paste("T_AT_mean",x2, sep="")))), col=Colour2, pch=pch2, type="o", lwd=Width2, cex=0.4) 
  lines(eval(parse(text=(paste("T_AT_mean",x3, sep="")))), col=Colour3,  pch=pch3, type="o", lwd=Width3, cex=0.4) 
  lines(eval(parse(text=(paste("T_AT_mean",x4, sep="")))), col=Colour4,  pch=pch4, type="o", lwd=Width4, cex=0.4) 
  lines(T_AT_mean1, lty=Line1, col=Colour1, lwd=Width1)
  legend("topleft", legend=c(Scenario1, Sce2, Sce3, Sce4), lty=c(Line1, Line2, Line3, Line4), col=c(Colour1, Colour2, Colour3, Colour4), pch=c(pch1, pch2, pch3, pch4), bty="n", cex=0.7) 
  dev.off() 
  
  #Figure
  Cairo(file=paste("DEFINE-Fig",Figure,"-5-Rate_of_profit.jpeg", sep=""), width=width_graph, height=height_graph, res=res_graph)
  par(mar = c(3,3.3,1,1)+0.1,mgp=c(1.8,0.5,0), family=font)
  plot(r_total_mean1, type="l", col=Colour1, lwd=Width1, ylab="Firms' profit rate (%)",  ylim=c(4, 12), yaxp=c(4,12,4),  xlab=" ", xaxt="n", cex.lab=0.8, cex.axis=0.7)
  axis(side=1,at=xlabel1, labels=xlabel2, cex.lab=0.8, cex.axis=0.7)
  lines(eval(parse(text=(paste("r_total_mean",x2, sep="")))), col=Colour2, pch=pch2, type="o", lwd=Width2, cex=0.4) 
  lines(eval(parse(text=(paste("r_total_mean",x3, sep="")))), col=Colour3, pch=pch3, type="o", lwd=Width3,cex=0.4) 
  lines(eval(parse(text=(paste("r_total_mean",x4, sep="")))), col=Colour4, pch=pch4, type="o", lwd=Width4,  cex=0.4) 
  lines(r_total_mean1, lty=Line1, col=Colour1, lwd=Width1)
  legend("bottomleft", legend=c(Scenario1, Sce2, Sce3, Sce4), lty=c(Line1, Line2, Line3, Line4), col=c(Colour1, Colour2, Colour3, Colour4), pch=c(pch1, pch2, pch3, pch4), bty="n", cex=0.7) 
  dev.off() 
  
  #Figure
  Cairo(file=paste("DEFINE-Fig",Figure,"-6-Default_rate.jpeg", sep=""), width=width_graph, height=height_graph, res=res_graph)
  par(mar = c(3,3.3,1,1)+0.1,mgp=c(1.8,0.5,0), family=font)
  plot(def_per_mean1, type="l", col=Colour1, lwd=Width1, ylab= "Firms' default rate (%)", ylim=c(2,8), xlab=" ", xaxt="n", cex.lab=0.8, cex.axis=0.7)
  axis(side=1,at=xlabel1, labels=xlabel2, cex.lab=0.8, cex.axis=0.7)
  lines(eval(parse(text=(paste("def_per_mean",x2, sep="")))), col=Colour2, pch=pch2, type="o", lwd=Width2,cex=0.4) 
  lines(eval(parse(text=(paste("def_per_mean",x3, sep="")))), col=Colour3, pch=pch3, type="o", lwd=Width3,cex=0.4) 
  lines(eval(parse(text=(paste("def_per_mean",x4, sep="")))), col=Colour4, pch=pch4, type="o", lwd=Width4, cex=0.4) 
  lines(def_per_mean1, lty=Line1, col=Colour1, lwd=Width1)
  legend("topleft", legend=c(Scenario1, Sce2, Sce3, Sce4), lty=c(Line1, Line2, Line3, Line4), col=c(Colour1, Colour2, Colour3, Colour4), pch=c(pch1, pch2, pch3, pch4), bty="n", cex=0.7) 
  dev.off() 
  
  #Figure
  Cairo(file=paste("DEFINE-Fig",Figure,"-7-Bank_leverage_ratio.jpeg", sep=""), width=width_graph, height=height_graph, res=res_graph)
  par(mar = c(3,3.3,1,1)+0.1,mgp=c(1.8,0.5,0), family=font)
  plot(lev_B_mean1, type="l", col=Colour1, lwd=Width1, ylab= "Banks' leverage ratio", ylim=c(0,40), xlab=" ", xaxt="n", cex.lab=0.8, cex.axis=0.7)
  axis(side=1,at=xlabel1, labels=xlabel2, cex.lab=0.8, cex.axis=0.7)
  lines(eval(parse(text=(paste("lev_B_mean",x2, sep="")))), col=Colour2, pch=pch2, type="o", lwd=Width2, cex=0.4) 
  lines(eval(parse(text=(paste("lev_B_mean",x3, sep="")))), col=Colour3, pch=pch3, type="o", lwd=Width3,cex=0.4) 
  lines(eval(parse(text=(paste("lev_B_mean",x4, sep="")))), col=Colour4, pch=pch4, type="o", lwd=Width4, cex=0.4) 
  lines(lev_B_mean1, lty=Line1, col=Colour1, lwd=Width1)
  legend("topleft", legend=c(Scenario1, Sce2, Sce3, Sce4), lty=c(Line1, Line2, Line3, Line4), col=c(Colour1, Colour2, Colour3, Colour4), pch=c(pch1, pch2, pch3, pch4), bty="n", cex=0.7) 
  dev.off() 
  
  
  #Figure
  Cairo(file=paste("DEFINE-Fig",Figure,"-8-Capital_adequacy_ratio.jpeg", sep=""), width=width_graph, height=height_graph, res=res_graph)
  par(mar = c(3,3.3,1,1)+0.1,mgp=c(1.8,0.5,0), family=font)
  plot(CAR_per_mean1, type="l", col=Colour1, lwd=Width1, ylab="Capital adequacy ratio (%)",  ylim=c(5,20), yaxp=c(5,20,6), xlab=" ", xaxt="n", cex.lab=0.8, cex.axis=0.7)
  axis(side=1,at=xlabel1, labels=xlabel2, cex.lab=0.8, cex.axis=0.7)
  lines(eval(parse(text=(paste("CAR_per_mean",x2, sep="")))), col=Colour2, pch=pch2, type="o", lwd=Width2, cex=0.4) 
  lines(eval(parse(text=(paste("CAR_per_mean",x3, sep="")))), col=Colour3, pch=pch3, type="o", lwd=Width3, cex=0.4) 
  lines(eval(parse(text=(paste("CAR_per_mean",x4, sep="")))), col=Colour4, pch=pch4, type="o", lwd=Width4, cex=0.4) 
  lines(CAR_per_mean1, lty=Line1, col=Colour1, lwd=Width1)
  legend("bottomleft", legend=c(Scenario1, Sce2, Sce3, Sce4), lty=c(Line1, Line2, Line3, Line4), col=c(Colour1, Colour2, Colour3, Colour4), pch=c(pch1, pch2, pch3, pch4), bty="n", cex=0.7) 
  dev.off() 
  
  #Figure
  Cairo(file=paste("DEFINE-Fig",Figure,"-9-Green_credit_rationing.jpeg", sep=""), width=width_graph, height=height_graph, res=res_graph)
  par(mar = c(3,3.3,1,1)+0.1,mgp=c(1.8,0.5,0), family=font)
  plot(CR_G_mean1, type="l", col=Colour1, lwd=Width1, ylab= "Credit rat. on green loans", ylim=c(0,0.4), yaxp=c(0,0.4,4), xlab=" ", xaxt="n", cex.lab=0.8, cex.axis=0.7)
  axis(side=1,at=xlabel1, labels=xlabel2, cex.lab=0.8, cex.axis=0.7)
  lines(eval(parse(text=(paste("CR_G_mean",x2, sep="")))), col=Colour2, pch=pch2, type="o", lwd=Width2, cex=0.4) 
  lines(eval(parse(text=(paste("CR_G_mean",x3, sep="")))), col=Colour3, pch=pch3, type="o", lwd=Width3,cex=0.4) 
  lines(eval(parse(text=(paste("CR_G_mean",x4, sep="")))), col=Colour4, pch=pch4, type="o", lwd=Width4, cex=0.4) 
  lines(CR_G_mean1, lty=Line1, col=Colour1, lwd=Width1)
  legend("topleft", legend=c(Scenario1, Sce2, Sce3, Sce4), lty=c(Line1, Line2, Line3, Line4), col=c(Colour1, Colour2, Colour3, Colour4), pch=c(pch1, pch2, pch3, pch4), bty="n", cex=0.7) 
  dev.off() 
  
  #Figure
  Cairo(file=paste("DEFINE-Fig",Figure,"-10-Conventional_credit_rationing.jpeg", sep=""), width=width_graph, height=height_graph, res=res_graph)
  par(mar = c(3,3.3,1,1)+0.1,mgp=c(1.8,0.5,0), family=font)
  plot(CR_C_S1_mean1, type="l", col=Colour1, lwd=Width1, ylab= "Credit rat. on conv. loans, mining & utilities", ylim=c(0,0.4), yaxp=c(0,0.4,4), xlab=" ", xaxt="n", cex.lab=0.8, cex.axis=0.7)
  axis(side=1,at=xlabel1, labels=xlabel2, cex.lab=0.8, cex.axis=0.7)
  lines(eval(parse(text=(paste("CR_C_S1_mean",x2, sep="")))), col=Colour2, pch=pch2, type="o", lwd=Width2, cex=0.4) 
  lines(eval(parse(text=(paste("CR_C_S1_mean",x3, sep="")))), col=Colour3, pch=pch3, type="o", lwd=Width3,cex=0.4) 
  lines(eval(parse(text=(paste("CR_C_S1_mean",x4, sep="")))), col=Colour4, pch=pch4, type="o", lwd=Width4, cex=0.4) 
  lines(CR_C_S1_mean1, lty=Line1, col=Colour1, lwd=Width1)
  legend("topleft", legend=c(Scenario1, Sce2, Sce3, Sce4), lty=c(Line1, Line2, Line3, Line4), col=c(Colour1, Colour2, Colour3, Colour4), pch=c(pch1, pch2, pch3, pch4), bty="n", cex=0.7) 
  dev.off() 
  
  #Figure
  Cairo(file=paste("DEFINE-Fig",Figure,"-11-Green_spread.jpeg", sep=""), width=width_graph, height=height_graph, res=res_graph)
  par(mar = c(3,3.3,1,1)+0.1,mgp=c(1.8,0.5,0), family=font)
  plot(spr_G_mean1, type="l", col=Colour1, lwd=Width1, ylab= "Spread on green loans", ylim=c(0,0.08), yaxp=c(0,0.08,4), xlab=" ", xaxt="n", cex.lab=0.8, cex.axis=0.7)
  axis(side=1,at=xlabel1, labels=xlabel2, cex.lab=0.8, cex.axis=0.7)
  lines(eval(parse(text=(paste("spr_G_mean",x2, sep="")))), col=Colour2, pch=pch2, type="o", lwd=Width2, cex=0.4) 
  lines(eval(parse(text=(paste("spr_G_mean",x3, sep="")))), col=Colour3, pch=pch3, type="o", lwd=Width3, cex=0.4) 
  lines(eval(parse(text=(paste("spr_G_mean",x4, sep="")))), col=Colour4, pch=pch4, type="o", lwd=Width4, cex=0.4) 
  lines(spr_G_mean1, lty=Line1, col=Colour1, lwd=Width1)
  legend("bottomleft", legend=c(Scenario1, Sce2, Sce3, Sce4), lty=c(Line1, Line2, Line3, Line4), col=c(Colour1, Colour2, Colour3, Colour4), pch=c(pch1, pch2, pch3, pch4), bty="n", cex=0.7) 
  dev.off() 
  
  #Figure
  Cairo(file=paste("DEFINE-Fig",Figure,"-12-Conventional_spread.jpeg", sep=""), width=width_graph, height=height_graph, res=res_graph)
  par(mar = c(3,3.3,1,1)+0.1,mgp=c(1.8,0.5,0), family=font)
  plot(spr_C_S1_mean1, type="l", col=Colour1, lwd=Width1, ylab= "Spread on conv. loans, mining & utilities", ylim=c(0,0.08), yaxp=c(0,0.1,5), xlab=" ", xaxt="n", cex.lab=0.8, cex.axis=0.7)
  axis(side=1,at=xlabel1, labels=xlabel2, cex.lab=0.8, cex.axis=0.7)
  lines(eval(parse(text=(paste("spr_C_S1_mean",x2, sep="")))), col=Colour2, pch=pch2, type="o", lwd=Width2, cex=0.4) 
  lines(eval(parse(text=(paste("spr_C_S1_mean",x3, sep="")))), col=Colour3, pch=pch3, type="o", lwd=Width3, cex=0.4) 
  lines(eval(parse(text=(paste("spr_C_S1_mean",x4, sep="")))), col=Colour4, pch=pch4, type="o", lwd=Width4, cex=0.4) 
  lines(spr_C_S1_mean1, lty=Line1, col=Colour1, lwd=Width1)
  legend("bottomleft", legend=c(Scenario1, Sce2, Sce3, Sce4), lty=c(Line1, Line2, Line3, Line4), col=c(Colour1, Colour2, Colour3, Colour4), pch=c(pch1, pch2, pch3, pch4), bty="n", cex=0.7) 
  dev.off() 
  
  
  #Figure
  Cairo(file="DEFINE-FigB-1-Output.jpeg", width=width_graph, height=height_graph, res=res_graph)
  par(mar = c(3,3.3,1,1)+0.1,mgp=c(1.8,0.5,0), family=font)
  plot(Y_mean2, type="l", col=Colour1, lwd=Width1, ylab= "Output (US$ trillion)", ylim=c(0, 800), yaxp=c(0,800,4),  xlab=" ", xaxt="n", cex.lab=0.8, cex.axis=0.7)
  axis(side=1,at=xlabel1, labels=xlabel2, cex.lab=0.8, cex.axis=0.7)
  polygon(c(1:T, T:1), c(Y_sdplus2, rev(Y_sdminus2)), col= c(Colour_range1), border=NA)
  polygon(c(1:T, T:1), c(Y_sdplus3, rev(Y_sdminus3)), col= c(Colour_range2), border=NA)
  lines(Y_mean3, col=Colour2, pch=pch2, type="o", lwd=Width2, cex=0.4) 
  lines(Y_mean2, lty=Line1, col=Colour1, lwd=Width1) 
  legend("topleft", legend=c(Mean, Range1, Scenario2, Range2), lty=c(Line1, Line_range1, Line2, Line_range2), col=c(Colour1, Colour_range1, Colour2, Colour_range2), pch=c(pch1, pch1, pch2, pch1),  lwd=c(Width1, Width_range1, Width2, Width_range2), bty="n", cex=0.7) 
  dev.off() 
  
  #Figure
  Cairo(file="DEFINE-FigB-2-Leverage_of_banks.jpeg", width=width_graph, height=height_graph, res=res_graph)
  par(mar = c(3,3.3,1,1)+0.1,mgp=c(1.8,0.5,0), family=font)
  plot(lev_B_mean2, type="l", col=Colour1, lwd=Width1, ylab= "Banks' leverage ratio", ylim=c(0,40), xlab=" ", xaxt="n", cex.lab=0.8, cex.axis=0.7)
  axis(side=1,at=xlabel1, labels=xlabel2, cex.lab=0.8, cex.axis=0.7)
  polygon(c(1:T, T:1), c(lev_B_sdplus2, rev(lev_B_sdminus2)), col= c(Colour_range1), border=NA)
  polygon(c(1:T, T:1), c(lev_B_sdplus3, rev(lev_B_sdminus3)), col= c(Colour_range2), border=NA)
  lines(lev_B_mean3, col=Colour2, pch=pch2, type="o", lwd=Width2,cex=0.4) 
  lines(lev_B_mean2, lty=Line1, col=Colour1, lwd=Width1) 
  legend("topleft", legend=c(Mean, Range1, Scenario2, Range2), lty=c(Line1, Line_range1, Line2, Line_range2), col=c(Colour1, Colour_range1, Colour2, Colour_range2), pch=c(pch1, pch1, pch2, pch1),  lwd=c(Width1, Width_range1, Width2, Width_range2), bty="n", cex=0.7) 
  dev.off() 
  
  
  #Figure
  Cairo(file="DEFINE-FigB-3-Emissions.jpeg", width=width_graph, height=height_graph, res=res_graph)
  par(mar = c(3,3.3,1,1)+0.1,mgp=c(1.8,0.5,0), family=font)
  plot(EMIS_mean2, type="l", col=Colour1, lwd=Width1, ylab=expression("CO" [2]*" emissions (GtCO" [2]*"/year)"), ylim=c(0,120), yaxp=c(0,120,4), xlab=" ", xaxt="n", cex.lab=0.8, cex.axis=0.7)
  axis(side=1,at=xlabel1, labels=xlabel2, cex.lab=0.8, cex.axis=0.7)
  polygon(c(1:T, T:1), c(EMIS_sdplus2, rev(EMIS_sdminus2)), col= c(Colour_range1), border=NA)
  polygon(c(1:T, T:1), c(EMIS_sdplus3, rev(EMIS_sdminus3)), col= c(Colour_range2), border=NA)
  lines(EMIS_mean3, col=Colour2, pch=pch2, type="o", lwd=Width2,cex=0.4) 
  lines(EMIS_mean2, lty=Line1, col=Colour1, lwd=Width1) 
  legend("topleft", legend=c(Mean, Range1, Scenario2, Range2), lty=c(Line1, Line_range1, Line2, Line_range2), col=c(Colour1, Colour_range1, Colour2, Colour_range2), pch=c(pch1, pch1, pch2, pch1),  lwd=c(Width1, Width_range1, Width2, Width_range2), bty="n", cex=0.7) 
  dev.off() 
  
  #Figure
  Cairo(file="DEFINE-FigB-4-Conventional_credit_rationing.jpeg", width=width_graph, height=height_graph, res=res_graph)
  par(mar = c(3,3.3,1,1)+0.1,mgp=c(1.8,0.5,0), family=font)
  plot(CR_C_S1_mean2, type="l", col=Colour1, lwd=Width1, ylab= "Credit rat. on conv. loans, mining & utilities", ylim=c(0,0.4), yaxp=c(0,0.4,8), xlab=" ", xaxt="n", cex.lab=0.8, cex.axis=0.7)
  axis(side=1,at=xlabel1, labels=xlabel2, cex.lab=0.8, cex.axis=0.7)
  polygon(c(1:T, T:1), c(CR_C_S1_sdplus2, rev(CR_C_S1_sdminus2)), col= c(Colour_range1), border=NA)
  polygon(c(1:T, T:1), c(CR_C_S1_sdplus3, rev(CR_C_S1_sdminus3)), col= c(Colour_range2), border=NA)
  lines(CR_C_S1_mean3, col=Colour2, pch=pch2, type="o", lwd=Width2, cex=0.4)
  lines(CR_C_S1_mean2, lty=Line1, col=Colour1, lwd=Width1)
  legend("bottomleft", legend=c(Mean, Range1, Scenario2, Range2), lty=c(Line1, Line_range1, Line2, Line_range2), col=c(Colour1, Colour_range1, Colour2, Colour_range2), pch=c(pch1, pch1, pch2, pch1),  lwd=c(Width1, Width_range1, Width2, Width_range2), bty="n", cex=0.7) 
  dev.off() 
  
  #Figure
  Cairo(file="DEFINE-FigB-5-Rate_of_profit.jpeg", width=width_graph, height=height_graph, res=res_graph)
  par(mar = c(3,3.3,1,1)+0.1,mgp=c(1.8,0.5,0), family=font)
  plot(r_total_mean2, type="l", col=Colour1, lwd=Width1, ylab="Firms' profit rate (%)",  ylim=c(4, 12), yaxp=c(4,12,4),  xlab=" ", xaxt="n", cex.lab=0.8, cex.axis=0.7)
  axis(side=1,at=xlabel1, labels=xlabel2, cex.lab=0.8, cex.axis=0.7)
  polygon(c(1:T, T:1), c(r_total_sdplus2, rev(r_total_sdminus2)), col= c(Colour_range1), border=NA)
  polygon(c(1:T, T:1), c(r_total_sdplus3, rev(r_total_sdminus3)), col= c(Colour_range2), border=NA)
  lines(r_total_mean3, col=Colour2, pch=pch2, type="o", lwd=Width2, cex=0.4) 
  lines(r_total_mean2, lty=Line1, col=Colour1, lwd=Width1) 
  legend("bottomleft", legend=c(Mean, Range1, Scenario2, Range2), lty=c(Line1, Line_range1, Line2, Line_range2), col=c(Colour1, Colour_range1, Colour2, Colour_range2), pch=c(pch1, pch1, pch2, pch1),  lwd=c(Width1, Width_range1, Width2, Width_range2), bty="n", cex=0.7) 
  dev.off() 
  
  #Figure
  Cairo(file="DEFINE-FigB-6-Default_rate.jpeg", width=width_graph, height=height_graph, res=res_graph)
  par(mar = c(3,3.3,1,1)+0.1,mgp=c(1.8,0.5,0), family=font)
  plot(def_per_mean2, type="l", col=Colour1, lwd=Width1, ylab= "Firms' default rate (%)", ylim=c(2,8), xlab=" ", xaxt="n", cex.lab=0.8, cex.axis=0.7)
  axis(side=1,at=xlabel1, labels=xlabel2, cex.lab=0.8, cex.axis=0.7)
  polygon(c(1:T, T:1), c(def_per_sdplus2, rev(def_per_sdminus2)), col= c(Colour_range1), border=NA)
  polygon(c(1:T, T:1), c(def_per_sdplus3, rev(def_per_sdminus3)), col= c(Colour_range2), border=NA)
  lines(def_per_mean3, col=Colour2, pch=pch2, type="o", lwd=Width2, cex=0.4) 
  lines(def_per_mean2, lty=Line1, col=Colour1, lwd=Width1) 
  legend("topleft", legend=c(Mean, Range1, Scenario2, Range2), lty=c(Line1, Line_range1, Line2, Line_range2), col=c(Colour1, Colour_range1, Colour2, Colour_range2), pch=c(pch1, pch1, pch2, pch1),  lwd=c(Width1, Width_range1, Width2, Width_range2), bty="n", cex=0.7) 
  dev.off() 
  
}

Figures(2,4,5,6, "GSF","DPF", "GSF+DPF", Colour4a)

if (figures_dummy==1) {
  Figures(3,7,4,8, "CT+GS","GSF", "GSF+CT+GS", Colour4a)
  Figures(4,7,5,9, "CT+GS","DPF", "DPF+CT+GS", Colour4a)
  Figures("D1",10,4,11, "CT+GS","GSF", "GSF+CT+GS", Colour4a)
  Figures("D2",10,5,12, "CT+GS","DPF", "DPF+CT+GS", Colour4a)
  Figures("D3",10, 27, 27, "CT+GS","Only CT", "", Colour4b)
}


#######################
#8. TABLES
#######################


matrixname<-paste("Table_1")
assign (matrixname, (round(rbind(c(g_Y_per_mean1[1], g_Y_per_mean1[33], mean(g_Y_per_mean1[1:33]), sd(g_Y_per_mean1[1:33])), c(ur_per_mean1[1], ur_per_mean1[33], mean(ur_per_mean1[1:33]), sd(ur_per_mean1[1:33])), c(POP_mean1[1], POP_mean1[33], mean(POP_mean1[1:33]), sd(POP_mean1[1:33])), c(theta_per_mean1[1], theta_per_mean1[33], mean(theta_per_mean1[1:33]), sd(theta_per_mean1[1:33])) , c(epsilon_ratio_mean1[1], epsilon_ratio_mean1[33], mean(epsilon_ratio_mean1[1:33]), sd(epsilon_ratio_mean1[1:33])) ,c(mu_ratio_mean1[1], mu_ratio_mean1[33], mean(mu_ratio_mean1[1:33]), sd(mu_ratio_mean1[1:33])) , c(EMIS_mean1[1], EMIS_mean1[33], mean(EMIS_mean1[1:33]), sd(EMIS_mean1[1:33])) ,c(tau_C_baseline[1]*1000, tau_C_baseline[33]*1000, mean(tau_C_baseline[1:33])*1000, sd(tau_C_baseline[1:33])*1000) , c(IG_E_Y_mean1[1], IG_E_Y_mean1[33], mean(IG_E_Y_mean1[1:33]), sd(IG_E_Y_mean1[1:33])) , c(def_per_mean1[1], def_per_mean1[33], mean(def_per_mean1[1:33]), sd(def_per_mean1[1:33])) , c(yield_C_per_mean1[1], yield_C_per_mean1[33], mean(yield_C_per_mean1[1:33]), sd(yield_C_per_mean1[1:33])) ,c(yield_G_per_mean1[1], yield_G_per_mean1[33], mean(yield_G_per_mean1[1:33]), sd(yield_G_per_mean1[1:33]))), digits=8)))

write.csv(Table_1,"DEFINE_Table_1.csv")

if (tables_dummy==1) {
  
  Tables<-function(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, y15, y16, y17, y18, y19, y20, y21){
    
    matrixname<-paste("Table_3")
    assign (matrixname, (round(rbind(c(eval(parse(text=(paste("g_Y_per_mean",y8,"[8]",sep="")))),
                                       eval(parse(text=(paste("g_Y_per_mean",y1,"[8]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y9,"[8]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("g_Y_per_mean",y8,"[13]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y1,"[13]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y9,"[13]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("g_Y_per_mean",y8,"[33]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y1,"[33]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y9,"[33]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("g_Y_per_mean",y8,"[58]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y1,"[58]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y9,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("g_Y_per_mean",y8,"[83]",sep="")))),
                                       eval(parse(text=(paste("g_Y_per_mean",y1,"[83]",sep="")))),
                                       eval(parse(text=(paste("g_Y_per_mean",y9,"[83]",sep=""))))),
                                     c(eval(parse(text=(paste("g_Y_per_mean",y10,"[8]",sep="")))),
                                       eval(parse(text=(paste("g_Y_per_mean",y2,"[8]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y11,"[8]",sep="")))), 
                                       0,  
                                       eval(parse(text=(paste("g_Y_per_mean",y10,"[13]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y2,"[13]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y11,"[13]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("g_Y_per_mean",y10,"[33]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y2,"[33]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y11,"[33]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("g_Y_per_mean",y10,"[58]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y2,"[58]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y11,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("g_Y_per_mean",y10,"[83]",sep="")))),
                                       eval(parse(text=(paste("g_Y_per_mean",y2,"[83]",sep="")))),
                                       eval(parse(text=(paste("g_Y_per_mean",y11,"[83]",sep=""))))),
                                     c(eval(parse(text=(paste("g_Y_per_mean",y12,"[8]",sep="")))),
                                       eval(parse(text=(paste("g_Y_per_mean",y3,"[8]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y13,"[8]",sep="")))), 
                                       0,  
                                       eval(parse(text=(paste("g_Y_per_mean",y12,"[13]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y3,"[13]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y13,"[13]",sep="")))), 
                                       0,  
                                       eval(parse(text=(paste("g_Y_per_mean",y12,"[33]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y3,"[33]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y13,"[33]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("g_Y_per_mean",y12,"[58]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y3,"[58]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y13,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("g_Y_per_mean",y12,"[83]",sep="")))),
                                       eval(parse(text=(paste("g_Y_per_mean",y3,"[83]",sep="")))),
                                       eval(parse(text=(paste("g_Y_per_mean",y13,"[83]",sep=""))))),
                                     c(eval(parse(text=(paste("g_Y_per_mean",y14,"[8]",sep="")))),
                                       eval(parse(text=(paste("g_Y_per_mean",y4,"[8]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y15,"[8]",sep="")))), 
                                       0,  
                                       eval(parse(text=(paste("g_Y_per_mean",y14,"[13]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y4,"[13]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y15,"[13]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("g_Y_per_mean",y14,"[33]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y4,"[33]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y15,"[33]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("g_Y_per_mean",y14,"[58]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y4,"[58]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y15,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("g_Y_per_mean",y14,"[83]",sep="")))),
                                       eval(parse(text=(paste("g_Y_per_mean",y4,"[83]",sep="")))),
                                       eval(parse(text=(paste("g_Y_per_mean",y15,"[83]",sep=""))))),
                                     c(eval(parse(text=(paste("g_Y_per_mean",y16,"[8]",sep="")))),
                                       eval(parse(text=(paste("g_Y_per_mean",y5,"[8]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y17,"[8]",sep="")))), 
                                       0,  
                                       eval(parse(text=(paste("g_Y_per_mean",y16,"[13]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y5,"[13]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y17,"[13]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("g_Y_per_mean",y16,"[33]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y5,"[33]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y17,"[33]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("g_Y_per_mean",y16,"[58]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y5,"[58]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y17,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("g_Y_per_mean",y16,"[83]",sep="")))),
                                       eval(parse(text=(paste("g_Y_per_mean",y5,"[83]",sep="")))),
                                       eval(parse(text=(paste("g_Y_per_mean",y17,"[83]",sep=""))))),
                                     c(eval(parse(text=(paste("g_Y_per_mean",y18,"[8]",sep="")))),
                                       eval(parse(text=(paste("g_Y_per_mean",y6,"[8]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y19,"[8]",sep="")))), 
                                       0,  
                                       eval(parse(text=(paste("g_Y_per_mean",y18,"[13]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y6,"[13]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y19,"[13]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("g_Y_per_mean",y18,"[33]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y6,"[33]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y19,"[33]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("g_Y_per_mean",y18,"[58]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y6,"[58]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y19,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("g_Y_per_mean",y18,"[83]",sep="")))),
                                       eval(parse(text=(paste("g_Y_per_mean",y6,"[83]",sep="")))),
                                       eval(parse(text=(paste("g_Y_per_mean",y19,"[83]",sep=""))))),
                                     c(eval(parse(text=(paste("g_Y_per_mean",y20,"[8]",sep="")))),
                                       eval(parse(text=(paste("g_Y_per_mean",y7,"[8]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y21,"[8]",sep="")))), 
                                       0,  
                                       eval(parse(text=(paste("g_Y_per_mean",y20,"[13]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y7,"[13]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y21,"[13]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("g_Y_per_mean",y20,"[33]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y7,"[33]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y21,"[33]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("g_Y_per_mean",y20,"[58]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y7,"[58]",sep="")))), 
                                       eval(parse(text=(paste("g_Y_per_mean",y21,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("g_Y_per_mean",y20,"[83]",sep="")))),
                                       eval(parse(text=(paste("g_Y_per_mean",y7,"[83]",sep="")))),
                                       eval(parse(text=(paste("g_Y_per_mean",y21,"[83]",sep=""))))),
                                     c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                     c(eval(parse(text=(paste("lev_B_mean",y8,"[8]",sep="")))),
                                       eval(parse(text=(paste("lev_B_mean",y1,"[8]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y9,"[8]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("lev_B_mean",y8,"[13]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y1,"[13]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y9,"[13]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("lev_B_mean",y8,"[33]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y1,"[33]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y9,"[33]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("lev_B_mean",y8,"[58]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y1,"[58]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y9,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("lev_B_mean",y8,"[83]",sep="")))),
                                       eval(parse(text=(paste("lev_B_mean",y1,"[83]",sep="")))),
                                       eval(parse(text=(paste("lev_B_mean",y9,"[83]",sep=""))))),
                                     c(eval(parse(text=(paste("lev_B_mean",y10,"[8]",sep="")))),
                                       eval(parse(text=(paste("lev_B_mean",y2,"[8]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y11,"[8]",sep="")))), 
                                       0,  
                                       eval(parse(text=(paste("lev_B_mean",y10,"[13]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y2,"[13]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y11,"[13]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("lev_B_mean",y10,"[33]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y2,"[33]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y11,"[33]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("lev_B_mean",y10,"[58]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y2,"[58]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y11,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("lev_B_mean",y10,"[83]",sep="")))),
                                       eval(parse(text=(paste("lev_B_mean",y2,"[83]",sep="")))),
                                       eval(parse(text=(paste("lev_B_mean",y11,"[83]",sep=""))))),
                                     c(eval(parse(text=(paste("lev_B_mean",y12,"[8]",sep="")))),
                                       eval(parse(text=(paste("lev_B_mean",y3,"[8]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y13,"[8]",sep="")))), 
                                       0,  
                                       eval(parse(text=(paste("lev_B_mean",y12,"[13]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y3,"[13]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y13,"[13]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("lev_B_mean",y12,"[33]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y3,"[33]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y13,"[33]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("lev_B_mean",y12,"[58]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y3,"[58]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y13,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("lev_B_mean",y12,"[83]",sep="")))),
                                       eval(parse(text=(paste("lev_B_mean",y3,"[83]",sep="")))),
                                       eval(parse(text=(paste("lev_B_mean",y13,"[83]",sep=""))))),
                                     c(eval(parse(text=(paste("lev_B_mean",y14,"[8]",sep="")))),
                                       eval(parse(text=(paste("lev_B_mean",y4,"[8]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y15,"[8]",sep="")))), 
                                       0,  
                                       eval(parse(text=(paste("lev_B_mean",y14,"[13]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y4,"[13]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y15,"[13]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("lev_B_mean",y14,"[33]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y4,"[33]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y15,"[33]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("lev_B_mean",y14,"[58]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y4,"[58]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y15,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("lev_B_mean",y14,"[83]",sep="")))),
                                       eval(parse(text=(paste("lev_B_mean",y4,"[83]",sep="")))),
                                       eval(parse(text=(paste("lev_B_mean",y15,"[83]",sep=""))))),
                                     c(eval(parse(text=(paste("lev_B_mean",y16,"[8]",sep="")))),
                                       eval(parse(text=(paste("lev_B_mean",y5,"[8]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y17,"[8]",sep="")))), 
                                       0,  
                                       eval(parse(text=(paste("lev_B_mean",y16,"[13]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y5,"[13]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y17,"[13]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("lev_B_mean",y16,"[33]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y5,"[33]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y17,"[33]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("lev_B_mean",y16,"[58]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y5,"[58]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y17,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("lev_B_mean",y16,"[83]",sep="")))),
                                       eval(parse(text=(paste("lev_B_mean",y5,"[83]",sep="")))),
                                       eval(parse(text=(paste("lev_B_mean",y17,"[83]",sep=""))))),
                                     c(eval(parse(text=(paste("lev_B_mean",y18,"[8]",sep="")))),
                                       eval(parse(text=(paste("lev_B_mean",y6,"[8]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y19,"[8]",sep="")))), 
                                       0,  
                                       eval(parse(text=(paste("lev_B_mean",y18,"[13]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y6,"[13]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y19,"[13]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("lev_B_mean",y18,"[33]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y6,"[33]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y19,"[33]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("lev_B_mean",y18,"[58]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y6,"[58]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y19,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("lev_B_mean",y18,"[83]",sep="")))),
                                       eval(parse(text=(paste("lev_B_mean",y6,"[83]",sep="")))),
                                       eval(parse(text=(paste("lev_B_mean",y19,"[83]",sep=""))))),
                                     c(eval(parse(text=(paste("lev_B_mean",y20,"[8]",sep="")))),
                                       eval(parse(text=(paste("lev_B_mean",y7,"[8]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y21,"[8]",sep="")))), 
                                       0,  
                                       eval(parse(text=(paste("lev_B_mean",y20,"[13]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y7,"[13]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y21,"[13]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("lev_B_mean",y20,"[33]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y7,"[33]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y21,"[33]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("lev_B_mean",y20,"[58]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y7,"[58]",sep="")))), 
                                       eval(parse(text=(paste("lev_B_mean",y21,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("lev_B_mean",y20,"[83]",sep="")))),
                                       eval(parse(text=(paste("lev_B_mean",y7,"[83]",sep="")))),
                                       eval(parse(text=(paste("lev_B_mean",y21,"[83]",sep=""))))),
                                     c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                     c(eval(parse(text=(paste("def_per_mean",y8,"[8]",sep="")))),
                                       eval(parse(text=(paste("def_per_mean",y1,"[8]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y9,"[8]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("def_per_mean",y8,"[13]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y1,"[13]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y9,"[13]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("def_per_mean",y8,"[33]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y1,"[33]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y9,"[33]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("def_per_mean",y8,"[58]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y1,"[58]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y9,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("def_per_mean",y8,"[83]",sep="")))),
                                       eval(parse(text=(paste("def_per_mean",y1,"[83]",sep="")))),
                                       eval(parse(text=(paste("def_per_mean",y9,"[83]",sep=""))))),
                                     c(eval(parse(text=(paste("def_per_mean",y10,"[8]",sep="")))),
                                       eval(parse(text=(paste("def_per_mean",y2,"[8]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y11,"[8]",sep="")))), 
                                       0,  
                                       eval(parse(text=(paste("def_per_mean",y10,"[13]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y2,"[13]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y11,"[13]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("def_per_mean",y10,"[33]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y2,"[33]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y11,"[33]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("def_per_mean",y10,"[58]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y2,"[58]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y11,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("def_per_mean",y10,"[83]",sep="")))),
                                       eval(parse(text=(paste("def_per_mean",y2,"[83]",sep="")))),
                                       eval(parse(text=(paste("def_per_mean",y11,"[83]",sep=""))))),
                                     c(eval(parse(text=(paste("def_per_mean",y12,"[8]",sep="")))),
                                       eval(parse(text=(paste("def_per_mean",y3,"[8]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y13,"[8]",sep="")))), 
                                       0,  
                                       eval(parse(text=(paste("def_per_mean",y12,"[13]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y3,"[13]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y13,"[13]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("def_per_mean",y12,"[33]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y3,"[33]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y13,"[33]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("def_per_mean",y12,"[58]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y3,"[58]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y13,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("def_per_mean",y12,"[83]",sep="")))),
                                       eval(parse(text=(paste("def_per_mean",y3,"[83]",sep="")))),
                                       eval(parse(text=(paste("def_per_mean",y13,"[83]",sep=""))))),
                                     c(eval(parse(text=(paste("def_per_mean",y14,"[8]",sep="")))),
                                       eval(parse(text=(paste("def_per_mean",y4,"[8]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y15,"[8]",sep="")))), 
                                       0,  
                                       eval(parse(text=(paste("def_per_mean",y14,"[13]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y4,"[13]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y15,"[13]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("def_per_mean",y14,"[33]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y4,"[33]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y15,"[33]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("def_per_mean",y14,"[58]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y4,"[58]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y15,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("def_per_mean",y14,"[83]",sep="")))),
                                       eval(parse(text=(paste("def_per_mean",y4,"[83]",sep="")))),
                                       eval(parse(text=(paste("def_per_mean",y15,"[83]",sep=""))))),
                                     c(eval(parse(text=(paste("def_per_mean",y16,"[8]",sep="")))),
                                       eval(parse(text=(paste("def_per_mean",y5,"[8]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y17,"[8]",sep="")))), 
                                       0,  
                                       eval(parse(text=(paste("def_per_mean",y16,"[13]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y5,"[13]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y17,"[13]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("def_per_mean",y16,"[33]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y5,"[33]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y17,"[33]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("def_per_mean",y16,"[58]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y5,"[58]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y17,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("def_per_mean",y16,"[83]",sep="")))),
                                       eval(parse(text=(paste("def_per_mean",y5,"[83]",sep="")))),
                                       eval(parse(text=(paste("def_per_mean",y17,"[83]",sep=""))))),
                                     c(eval(parse(text=(paste("def_per_mean",y18,"[8]",sep="")))),
                                       eval(parse(text=(paste("def_per_mean",y6,"[8]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y19,"[8]",sep="")))), 
                                       0,  
                                       eval(parse(text=(paste("def_per_mean",y18,"[13]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y6,"[13]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y19,"[13]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("def_per_mean",y18,"[33]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y6,"[33]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y19,"[33]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("def_per_mean",y18,"[58]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y6,"[58]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y19,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("def_per_mean",y18,"[83]",sep="")))),
                                       eval(parse(text=(paste("def_per_mean",y6,"[83]",sep="")))),
                                       eval(parse(text=(paste("def_per_mean",y19,"[83]",sep=""))))),
                                     c(eval(parse(text=(paste("def_per_mean",y20,"[8]",sep="")))),
                                       eval(parse(text=(paste("def_per_mean",y7,"[8]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y21,"[8]",sep="")))), 
                                       0,  
                                       eval(parse(text=(paste("def_per_mean",y20,"[13]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y7,"[13]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y21,"[13]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("def_per_mean",y20,"[33]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y7,"[33]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y21,"[33]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("def_per_mean",y20,"[58]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y7,"[58]",sep="")))), 
                                       eval(parse(text=(paste("def_per_mean",y21,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("def_per_mean",y20,"[83]",sep="")))),
                                       eval(parse(text=(paste("def_per_mean",y7,"[83]",sep="")))),
                                       eval(parse(text=(paste("def_per_mean",y21,"[83]",sep=""))))),
                                     c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                     c(eval(parse(text=(paste("T_AT_mean",y8,"[8]",sep="")))),
                                       eval(parse(text=(paste("T_AT_mean",y1,"[8]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y9,"[8]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("T_AT_mean",y8,"[13]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y1,"[13]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y9,"[13]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("T_AT_mean",y8,"[33]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y1,"[33]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y9,"[33]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("T_AT_mean",y8,"[58]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y1,"[58]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y9,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("T_AT_mean",y8,"[83]",sep="")))),
                                       eval(parse(text=(paste("T_AT_mean",y1,"[83]",sep="")))),
                                       eval(parse(text=(paste("T_AT_mean",y9,"[83]",sep=""))))),
                                     c(eval(parse(text=(paste("T_AT_mean",y10,"[8]",sep="")))),
                                       eval(parse(text=(paste("T_AT_mean",y2,"[8]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y11,"[8]",sep="")))), 
                                       0,  
                                       eval(parse(text=(paste("T_AT_mean",y10,"[13]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y2,"[13]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y11,"[13]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("T_AT_mean",y10,"[33]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y2,"[33]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y11,"[33]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("T_AT_mean",y10,"[58]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y2,"[58]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y11,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("T_AT_mean",y10,"[83]",sep="")))),
                                       eval(parse(text=(paste("T_AT_mean",y2,"[83]",sep="")))),
                                       eval(parse(text=(paste("T_AT_mean",y11,"[83]",sep=""))))),
                                     c(eval(parse(text=(paste("T_AT_mean",y12,"[8]",sep="")))),
                                       eval(parse(text=(paste("T_AT_mean",y3,"[8]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y13,"[8]",sep="")))), 
                                       0,  
                                       eval(parse(text=(paste("T_AT_mean",y12,"[13]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y3,"[13]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y13,"[13]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("T_AT_mean",y12,"[33]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y3,"[33]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y13,"[33]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("T_AT_mean",y12,"[58]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y3,"[58]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y13,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("T_AT_mean",y12,"[83]",sep="")))),
                                       eval(parse(text=(paste("T_AT_mean",y3,"[83]",sep="")))),
                                       eval(parse(text=(paste("T_AT_mean",y13,"[83]",sep=""))))),
                                     c(eval(parse(text=(paste("T_AT_mean",y14,"[8]",sep="")))),
                                       eval(parse(text=(paste("T_AT_mean",y4,"[8]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y15,"[8]",sep="")))), 
                                       0,  
                                       eval(parse(text=(paste("T_AT_mean",y14,"[13]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y4,"[13]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y15,"[13]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("T_AT_mean",y14,"[33]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y4,"[33]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y15,"[33]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("T_AT_mean",y14,"[58]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y4,"[58]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y15,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("T_AT_mean",y14,"[83]",sep="")))),
                                       eval(parse(text=(paste("T_AT_mean",y4,"[83]",sep="")))),
                                       eval(parse(text=(paste("T_AT_mean",y15,"[83]",sep=""))))),
                                     c(eval(parse(text=(paste("T_AT_mean",y16,"[8]",sep="")))),
                                       eval(parse(text=(paste("T_AT_mean",y5,"[8]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y17,"[8]",sep="")))), 
                                       0,  
                                       eval(parse(text=(paste("T_AT_mean",y16,"[13]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y5,"[13]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y17,"[13]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("T_AT_mean",y16,"[33]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y5,"[33]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y17,"[33]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("T_AT_mean",y16,"[58]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y5,"[58]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y17,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("T_AT_mean",y16,"[83]",sep="")))),
                                       eval(parse(text=(paste("T_AT_mean",y5,"[83]",sep="")))),
                                       eval(parse(text=(paste("T_AT_mean",y17,"[83]",sep=""))))),
                                     c(eval(parse(text=(paste("T_AT_mean",y18,"[8]",sep="")))),
                                       eval(parse(text=(paste("T_AT_mean",y6,"[8]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y19,"[8]",sep="")))), 
                                       0,  
                                       eval(parse(text=(paste("T_AT_mean",y18,"[13]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y6,"[13]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y19,"[13]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("T_AT_mean",y18,"[33]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y6,"[33]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y19,"[33]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("T_AT_mean",y18,"[58]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y6,"[58]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y19,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("T_AT_mean",y18,"[83]",sep="")))),
                                       eval(parse(text=(paste("T_AT_mean",y6,"[83]",sep="")))),
                                       eval(parse(text=(paste("T_AT_mean",y19,"[83]",sep=""))))),
                                     c(eval(parse(text=(paste("T_AT_mean",y20,"[8]",sep="")))),
                                       eval(parse(text=(paste("T_AT_mean",y7,"[8]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y21,"[8]",sep="")))), 
                                       0,  
                                       eval(parse(text=(paste("T_AT_mean",y20,"[13]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y7,"[13]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y21,"[13]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("T_AT_mean",y20,"[33]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y7,"[33]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y21,"[33]",sep="")))), 
                                       0, 
                                       eval(parse(text=(paste("T_AT_mean",y20,"[58]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y7,"[58]",sep="")))), 
                                       eval(parse(text=(paste("T_AT_mean",y21,"[58]",sep="")))), 
                                       0,
                                       eval(parse(text=(paste("T_AT_mean",y20,"[83]",sep="")))),
                                       eval(parse(text=(paste("T_AT_mean",y7,"[83]",sep="")))),
                                       eval(parse(text=(paste("T_AT_mean",y21,"[83]",sep="")))))
    ), digits=8)))
    
    write.csv(Table_3,"DEFINE_Table_3.csv") 
  }  
  Tables(1, 4, 5, 6, 7, 8, 9, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26)
  
}


