# =============================================================================
# File: 03_model_rat_brain_iadd_30day.R
# Purpose: mrgsolve model definition for 03 model rat brain iadd 30day.
# Note: Model equations and parameter values are preserved from the validated script.
# The object created by this file is `DexPBPK`.
# =============================================================================

DexPBPK<-
  "
$PROB
## Minimal DEX PBPK (rat): Brain membrane-limited; Liver/Yes Kidney/Rest flow-limited
## QRC and VRC are derived:
##   QRC = 1 - (QLC + QKC + QBRC)
##   VRC = 1 - (VLC + VKC + VBRC + VBC)

$PARAM @annotated
// --- Body weight & cardiac output ---
BW     : 0.25     : kg,        body weight (rat)
QCC    : 14     : L/h/kg^0.75, cardiac output scaling

// --- Flow fractions (of QC) ---
QLC    : 0.183    : -, fraction of QC to liver
QKC    : 0.141    : -, fraction of QC to kidney
QBRC   : 0.02    : -, fraction of QC to brain
// QRC is derived: QRC = 1 - (QLC+QKC+QBRC)

// --- Tissue volume fractions (of BW) ---
VLC    : 0.0784    : -, fraction of BW for liver volume
VKC    : 0.0148    : -, fraction of BW for kidney volume
VBRC   : 0.0048    : -, fraction of BW for brain volume
VBC    : 0.054    : -, fraction of BW for total blood volume
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

// --- Dose release ---
DOSE_IADD   : 0.16      : mg, total amount loaded in implant
FMAX        : 1.0       : Max release propotion of 60-days release profile
k           : 0.272     : -, weibull shape
lambda      : 9.736     : h, weibull scale


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
IADD   : 0.16 : mg, I-ADD compartment (releasable)   // average based on mass of device



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



// Weibull total releasable dose
double t  = SOLVERTIME;
double tt = (t > 0.01) ? t : 0.01;           // Use a practical lower bound to avoid singularity when k<1
double z = tt / lambda;
double RATE_IADD = DOSE_IADD * FMAX * (k / lambda) *
                   pow(z, k - 1.0) * exp(-pow(z, k)); 

// stop release if device is depleted (prevents negative IADD)
if(IADD <= 0.0) RATE_IADD = 0.0;


// --- Brain (membrane-limited: capillary tissue) ---
// --- Three layer model ---
double CBRvas_u = (CBR_VAS/BP)*fu;
double CBRint_u = (CBR_INT/PBR)*fuint_brain;
double CBRcel_u = (CBR_CEL/PBR)*fucel_brain;                                // if you use same PBR; ideally separate P for cell

double JBR_vas_int = PABR_VAS * CBRvas_u - PABR_INT*CBRint_u ;          // mg/L * L/h = mg/h
double JBR_int_cel = PABR_IN * (CBRint_u) - PABR_OUT * (CBRcel_u);      // 

double RABR_VAS = QBR*(CA - CBR_VAS) - JBR_vas_int + RATE_IADD;
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
dxdt_IADD    = -RATE_IADD;



// --- Accumulators ---
dxdt_ADose = RATE_IADD;               // track cumulative released
dxdt_AMet  = CLint*fu*CL/PL;       // hepatic elimination (mg/h)



$TABLE
// Totals (mg)
double MASS_TOTAL = AA+AV+AL+AK+AR+ABR_VAS+ABR_INT+ABR_CEL+IADD;
double MASS_IN    = DOSE_IADD;   // ADose
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
// capture CKID_VAS_out     = CKIDvas;
// capture CKID_INT_out     = CKIDint;
// capture CKID_CEL_out     = CKIDint + CKIDcel;
capture CKID_out         = CKID;
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
# 
# library(mrgsolve)
# mod <- mcode("DexPBPK",DexPBPK) %>%
#   param(fu=0.1, BP=1, CLint=0.30)
# 
# param(mod)
# 
# 
# # The model do not need any dose event as it already has dose event
# # Simulate for 10000 h with 0.5 h step
# out <- mod %>%
#   mrgsim(end = 1000, delta = 0.5) %>%
#   as.data.frame()
# 
# head(out)
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
