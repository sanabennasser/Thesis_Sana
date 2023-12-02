//*********************************************************************
//Geometrical data
//*********************************************************************
mm=1e-3;
NbWires =1;
dc = 14.3*mm; //Diameter of conductor
txlpe = 8*mm; //Thickness of XLPE insulation
tcu=0.1*mm; //thickness of copper tape
topvc= 2.2*mm; //thickness of outer sheath (pvc)
dtot = 38*mm; //Outer Diameter
dinf=5*dtot; //Electromagnetc domain

//*********************************************************************
//material properties
//*********************************************************************
// relative permittivity
espr_xlpe = 2.5;
espr_pvc= 4;
// 1 for Cu, Al, soil(dry)

//relative permeability
mur_pvc=1;
mur_cu = 1;
// 1 for Cu, Al, polyethylene, semiconductor, XLPE, soil(dry)

//electrical conductivity [S/m]
sigma_cu = 5.99e7;
sigma_pvc = 10e-6;
sigma_xlpe = 1.0e-18;
sigma_soil = 28;
V0 = 36000;

// //*********************************************************************
// //mesh properties
// //*********************************************************************
// DefineConstant[ s= {1., Name "{00Parameters/Global mesh size factor"}];
// //DefineConstant[ Freq= {1., Name "Frequency"}];
//*********************************************************************
//Physical numbers (none of these numbers can onverlap)
//*********************************************************************
WIRE=100;
XLPE=200;
CU=300;
OPVC=400;
SOIL_EM=500;
OUTBND_EM=511;
