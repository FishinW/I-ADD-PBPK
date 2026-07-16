## Revised version
## Delete fu in other tissues (except brain)
## Change the parameters values to Mouse
## I-ADD device is set in the artery of brain
## Three layer brain model


DexPBPK<-
  "
$PROB
## Minimal DEX PBPK (rat): Brain membrane-limited; Liver/Yes Kidney/Rest flow-limited
## QRC and VRC are derived:
##   QRC = 1 - (QLC + QKC + QBRC)
##   VRC = 1 - (VLC + VKC + VBRC + VBC)

$PARAM @annotated
// --- Body weight & cardiac output ---
BW     : 0.02     : kg,        body weight (rat)
QCC    : 16.5     : L/h/kg^0.75, cardiac output scaling

// --- Flow fractions (of QC) ---
QLC    : 0.161    : -, fraction of QC to liver
QKC    : 0.091    : -, fraction of QC to kidney
QBRC   : 0.033    : -, fraction of QC to brain
// QRC is derived: QRC = 1 - (QLC+QKC+QBRC)

// --- Tissue volume fractions (of BW) ---
VLC    : 0.065    : -, fraction of BW for liver volume
VKC    : 0.017    : -, fraction of BW for kidney volume
VBRC   : 0.0048    : -, fraction of BW for brain volume
VBC    : 0.085    : -, fraction of BW for total blood volume
// VRC is derived: VRC = 1 - (VLC+VKC+VBRC+VBC)


// --- Brain capillary fraction & arterial split ---
BVBR_VAS   : 0.03     : -, fraction of VBR that is vascular      (brain)
BVBR_INT   : 0.074    : -, fraction of VBR that is interstitial  (brain)
BVBR_CEL   : 0.896    : -, fraction of VBR that is intracellular  (brain) BVBR_VAS + BCBR_INT + BVBR_CEL = 1.0
BVB_A      : 0.25     : -, fraction of VB that is arterial blood (rest venous)


// --- Partition coefficients (placeholders) ---
PL     : 6.76      : -, liver tissue:plasma partition
PK     : 1.51      : -, kidney tissue:plasma partition
PR     : 1.0      : -, rest tissue:plasma partition
PBR    : 1.2      : -, brain tissue:plasma partition

// --- Brain permeability scaling (only brain is membrane-limited) ---
PABRC_VAS   : 0.6   : -, Brain vascular to interstitial permeability coefficient. (PABR_VAS = PABRC_VAS*QBR)
PABRC_INT   : 0.6   : -, Brain interstitial to vascular permeability coefficient. (PABR_INT = PABRC_INT*QBR)
PABRC_OUT   : 1.0   : -, Brain intracellular to interstitial permeability coefficient. (PABR_OUT = PABRC_OUT*QBR)
PABRC_IN   : 0.24   : -, Brain interstitial to intracellular permeability coefficient. (PABR_IN = PABRC_IN*QBR)

// --- Binding / blood-plasma & clearance ---
fucel_brain  : 0.24       : -, unbound fraction in bain cellular
fuint_brain  : 0.29       : -, unbound fraction in brain intestitial
fu     : 0.175      : -, unbound fraction in plasma
BP     : 0.725      : -, blood:plasma ratio
CLint  : 0.40     : L/h, intrinsic hepatic clearance (applied to liver tissue)


$INIT @annotated
AA     : 0 : mg, arterial blood amount
AV     : 0 : mg, venous blood amount
AL     : 0 : mg, liver tissue (flow-limited)
AK     : 0 : mg, kidney
AR     : 0 : mg, rest tissue (flow-limited)
ABR_VAS: 0 : mg, brain vascular blood
ABR_INT: 0 : mg, brain interstital
ABR_CEL: 0 : mg, brain intracellular
ADose  : 0 : mg, cumulative administered dose (duplicate every input here)
AMet   : 0 : mg, cumulative hepatic elimination


$ODE
// --- Cardiac output & primary organ flows (L/h) ---
double QC   = QCC*pow(BW, 0.75);  
double QL   = QC*QLC;
double QK   = QC*QKC;
double QBR  = QC*QBRC;

// --- Derived rest flow fraction and flow ---
double QRC_frac = 1.0 - (QLC + QKC + QBRC);   // should be >= 0
double QR   = QC * QRC_frac;

// --- Volumes (L) ---
double VL   = BW*VLC;
double VK   = BW*VKC;
double VBR  = BW*VBRC;
double VB   = BW*VBC;

// --- Derived rest volume fraction and volume ---
double VRC_frac = 1.0 - (VLC + VKC + VBRC + VBC);  // should be >= 0
double VR   = BW*VRC_frac;

// --- Sub-volumes (L) ---
double VBR_VAS = VBR*BVBR_VAS;
double VBR_INT = VBR*BVBR_INT;
double VBR_CEL = VBR - VBR_VAS - VBR_INT ;

double VA   = VB*BVB_A;
double VV   = VB - VA;

// --- Concentrations (mg/L) ---
double CA    = AA/VA;
double CV    = AV/VV;
double CL    = (VL>0.0)  ? AL/VL    : 0.0;
double CK    = (VK>0.0) ?  AK/VK    : 0.0;


double CR    = (VR>0.0)  ? AR/VR    : 0.0;
double CBR_VAS  = (VBR_VAS>0.) ? ABR_VAS/VBR_VAS: 0.0;
double CBR_INT  = (VBR_INT>0.) ? ABR_INT/VBR_INT: 0.0;
double CBR_CEL  = (VBR_CEL>0.) ? ABR_CEL/VBR_CEL: 0.0;


// --- Brain permeability (L/h) ---
double PABR_VAS  = PABRC_VAS*QBR;                  // Brain vascular to interstitial permeability 
double PABR_INT  = PABRC_INT*QBR;                  // Brain interstitial to vascular permeability
double PABR_OUT  = PABRC_OUT*QBR;                            // Brain intracellular to interstitial permeability (L/h) starts as “how much blood can leave/enter per time” for each organ
double PABR_IN   = PABRC_IN*QBR;                      // Brain interstitial to intracellular  permeability (L/h) 


// --- Blood recirculation (lungs implicit) ---
double RAA = QC*CV - QC*CA;
double RAV = QL*(CL*BP/PL) + QK*(CK*BP/PK) + QR*(CR*BP/PR) + QBR*CBR_VAS - QC*CV;

// --- Flow-limited tissues ---
double RAL = QL*(CA - (CL*BP/PL)) - CLint*CL*fu/PL;      // liver
double RAK = QK*(CA - (CK*BP/PK));                      // kidney
double RAR = QR*(CA - (CR*BP/PR));                       // rest


// --- Brain (membrane-limited: capillary tissue) ---
// --- Three layer model ---
double CBRvas_u = (CBR_VAS/BP)*fu;
double CBRint_u = (CBR_INT/PBR)*fuint_brain;
double CBRcel_u = (CBR_CEL/PBR)*fucel_brain;                                // if you use same PBR; ideally separate P for cell

double JBR_vas_int = PABR_VAS * CBRvas_u - PABR_INT*CBRint_u ;          // mg/L * L/h = mg/h
double JBR_int_cel = PABR_IN * (CBRint_u) - PABR_OUT * (CBRcel_u);      // 

double RABR_VAS = QBR*(CA - CBR_VAS) - JBR_vas_int;
double RABR_INT = +JBR_vas_int - JBR_int_cel;
double RABR_CEL = +JBR_int_cel;



// --- ODEs ---
dxdt_AA      = RAA;
dxdt_AV      = RAV;
dxdt_AL      = RAL;
dxdt_AK      = RAK;
dxdt_AR      = RAR;
dxdt_ABR_VAS = RABR_VAS;
dxdt_ABR_INT = RABR_INT;
dxdt_ABR_CEL = RABR_CEL;


// --- Accumulators ---
dxdt_ADose = 0.0;               // track cumulative released
dxdt_AMet  = CLint*fu*CL/PL;       // hepatic elimination (mg/h)


$TABLE
// Totals (mg)
double MASS_TOTAL = AA+AV+AL+AK+AR+ABR_VAS+ABR_INT+ABR_CEL;
double MASS_IN    = ADose;
double MASS_ELIM  = AMet;
double MASS_RESID = MASS_IN - (MASS_TOTAL + MASS_ELIM);   // ~0 with tight tolerances

// Selected concentrations (mg/L)
double CART   = AA/(VB*BVB_A);
double CVEN   = AV/(VB*(1.0-BVB_A));
double CLIVER = CL;
double CKID   = CK;
double CREST  = CR;
double CBRvas = (VBR_VAS>0.0)? ABR_VAS/VBR_VAS : 0.0;
double CBRint = (VBR_INT>0.0)? ABR_INT/VBR_INT : 0.0;
double CBRcel = (VBR_CEL>0.0)? ABR_CEL/VBR_CEL : 0.0;
double CBR    = (ABR_CEL + ABR_INT + ABR_VAS)/(VBR_CEL + VBR_INT + VBR_VAS);


capture CART_out        = CART;
capture CVEN_out        = CVEN;
capture CPLASMA_out	    = CVEN/BP;
capture CLIVER_out      = CLIVER;
capture CKID_out        = CKID;
capture CREST_out       = CREST;
capture MASS_TOTAL_out  = MASS_TOTAL;
capture MASS_IN_out     = MASS_IN;
capture MASS_ELIM_out   = MASS_ELIM;
capture MASS_RESID_out  = MASS_RESID;
// capture CBR_VAS_out     = CBRvas;
// capture CBR_INT_out     = CBRint;
capture CBR_out         = CBR;

"

# 
# library(mrgsolve)
# mod <- mcode("DexPBPK",DexPBPK) %>%
#   param(fu=0.1, BP=1)
# 
# param(mod)
# 
# # User-defined schedule
# DOSE_IADD   <- 1        # mg per dose
# dur_h       <- 0.05     # h infusion duration per dose
# tinterval   <- 0.5      # h between starts of each dose
# NDOSE       <- 6        # total number of doses
# 
# rate_iadd <- DOSE_IADD / dur_h  # mg/h
# 
# # Dose into brain capillary
# e_iadd <- ev(
#   ID = 1,
#   cmt = "ABR_VAS",
#   amt = DOSE_IADD,
#   rate = rate_iadd,       # infusion (omit rate for bolus)
#   ii = tinterval,
#   addl = NDOSE - 1,
#   replicate = FALSE
# )
# 
# # Duplicate into ADose for mass-balance purpose
# e_in <- ev(
#   ID = 1,
#   cmt = "ADose",
#   amt = DOSE_IADD,
#   rate = rate_iadd,
#   ii = tinterval,
#   addl = NDOSE - 1,
#   replicate = FALSE
# )
# 
# e <- e_iadd + e_in
# out <- mrgsim(mod, ev = e, end = NDOSE*tinterval + 4, delta = 0.02)
# head(out)
# 
# # # The model do not need any dose event as it already has dose event
# # # Simulate for 10000 h with 0.5 h step
# # out <- mod %>%
# #   mrgsim(end = 1000, delta = 0.5) %>%
# #   as.data.frame()
# #
# # head(out)
# 
# 
# library(ggplot2)
# 
# ggplot(as.data.frame(out), aes(x = time, y = MASS_RESID_out)) +
#   geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) +
#   geom_line(linewidth = 0.7) +
#   labs(
#     x = "Time (h)",
#     y = "Mass residual (mg)",
#     title = "Mass balance residual: MASS_RESID_out = MASS_IN_out − (MASS_TOTAL_out + MASS_ELIM_out)"
#   ) +
#   theme_bw() +
#   theme(plot.title = element_text(size = 11))
# 
