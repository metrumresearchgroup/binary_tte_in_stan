$PARAM  
VC = 2.44, 
KA = 0.92, 
CL=1, 
ECOG1=1, 
ECOG2=0, 
cAGE=0, 
LAMBDA0 = 0.0217, 
EMAX = 0.692, 
EC50=4.956,
beta1 = .095, // 10% increase in hazard ECOG=1 vs ECOG=0
beta2 = .223, // 25% increase in hazard ECOG>1 vs ECOG=0
beta3 = .095, // 10% increase in hazard 



$SET delta=0.1, end=180

$CMT GUT CENT CHAZARD

$ODE 
double CP = CENT/VC;

dxdt_GUT = -KA*GUT;
dxdt_CENT = KA*GUT - (CL/VC)*CENT;
double NU = exp(beta1*ECOG1+beta2*ECOG2+beta3*cAGE);
dxdt_CHAZARD = LAMBDA0 * (1 + EMAX * CP / (EC50 + CP))*NU;

$TABLE
double EFFECT = LAMBDA0 * (1 + EMAX * CP / (EC50 + CP))*NU;
CP = CENT/VC;

$CAPTURE CP EFFECT
