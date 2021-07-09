$PROB
# Model: `pk2cmt`
  - Two-compartment PK model
      - Dual first-order absorption
      - Optional nonlinear clearance from `CENT`
  - Source: `mrgsolve` internal library
  - Date: `r Sys.Date()`
  - Version: `r packageVersion("mrgsolve")`
  
$PARAM @annotated
CL   :  1  : Clearance (dL/hr)
VC   : 20  : Central volume (dL)
Q    :  2  : Inter-compartmental clearance (dL/hr)
VP   : 10  : Peripheral volume of distribution (dL)
KA1  :  1  : Absorption rate constant 1 (1/time)        
KA2  :  1  : Absorption rate constant 2 (1/time)       
VMAX :  0  : Maximum velocity (mass/time)         
KM   :  2  : Michaelis Constant (mass/volume)
lambda : 2.96 : baseline hazard constant (1/year)
IF50   : 10.2 : IC50 (IU/dL)
gamma  : -0.566 : rate (1/year)
  
$CMT  @annotated
EV1    : First extravascular compartment (mass)
CENT   : Central compartment (mass)
PERIPH : Peripheral compartment (mass) 
EV2    : Second extravascular compartment (mass)
hazard : hazard
    
$GLOBAL 
#define CP (CENT/VC)
#define CT (PERIPH/VP)
#define CLNL (VMAX/(KM+CP))

$ODE
dxdt_EV1 = -KA1*EV1;
dxdt_EV2 = -KA2*EV2;
dxdt_CENT = KA1*EV1 + KA2*EV2 - (CL+CLNL+Q)*CP  + Q*CT;
dxdt_PERIPH = Q*CP - Q*CT;

double haz = (lambda/24/365) * exp(gamma * (SOLVERTIME/24/365 - 1)) * (1 - CP / (IF50 + CP));
dxdt_hazard = haz;

$CAPTURE  @annotated
CP : Plasma concentration (mass/time)
  
