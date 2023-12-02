// Include "cable_data.geo"; // Optional file but recommended as data are often common to geometry and finite-element description
Include "cable_data.geo";
DefineConstant[
  Flag_AnalysisType = {0,
    Choices{
      0="Electric",
      1="Magnetic"
    },
    Name "{00Parameters/00Type of analysis", Highlight "Blue",
    ServerAction Str["Reset","GetDP/1ResolutionChoices"]
  }
  nb_iter = 20, // Maximum number of nonlinear iterations (You may adapt)
  relaxation_factor = 1, // value in [0,1]; if 1, there is no relaxation; if <1, you used the solution of previous iteration for helping convergence
  stop_criterion = 1e-6, // prescribed tolerance, iterative process stops when the difference between two consecutive iterations is smaller than this value

  // You can predefine the default Resolution, ComputeCommand and Operation (-solve, -pos)
  // In these files:
  // Resolution depends on Flag_AnalysisType => the Default depends on the formulation file that is included
  // PostOperations are called in the Resolution => No need to indicate anything else
  // r_ = {"", Name "GetDP/1ResolutionChoices", Visible 1}
  c_ = {"-solve -v2", Name "GetDP/9ComputeCommand", Visible 1},
  p_ = {"", Name "GetDP/2PostOperationChoices", Visible 1}
];

Group {
  XLPE = Region[{XLPE}];
  CU = Region[{CU}];
  OuterPVC = Region[{OPVC}];
  SoilEM = Region[{SOIL_EM}];
  Soil = Region[{SoilEM}];
//  Wire=Region[{WIRE}];
  For k In {1:NbWires}
    Ind~{k} = Region[{(WIRE+k-1)}];
    Inds += Region[{(WIRE+k-1)}];
  EndFor
  // electrodynamics
  Cable = Region[{Inds, XLPE,CU,OuterPVC}]; // All the regions in the cable, for convenience and postprocessing
  SurfaceGe0 = Region[{OUTBND_EM}]; //n.b=0 on this boundary
  Ele=Region[{Inds,XLPE,CU,OuterPVC,SoilEM}];

//  Ele += Region[{}];
  Domain_Ele= Region[{Ele}];
  //magnetodynamics
  DomainS_Mag= Region[{Inds}]; // If using Current_2D, it allows accounting for the dependance of sigma with T
  DomainCC_Mag = Region[{SoilEM,XLPE,OuterPVC,Inds}];// Non-conduction domain
  DomainC_Mag = Region [{CU}];// Conducting domain //the aluminum shealth, the steel armour and the steel pipe.
  DomainS0_Mag= Region[{}]; //if imposing source with js0[]

  DomainCWtihI_Mag = Region[{}];
  Domain_Mag= Region[{DomainCC_Mag,DomainC_Mag}];
  DomainDummy = Region[{12345}];
}

Function {
  mu0 = 4.e-7 * Pi;
  eps0 = 8.854187818e-12;

  // TO DEFINE FOR ALL MATERIALS
  // nu[], sigma[], epsilon[]... are piecewise defined functions
  // - if no Region is indicated, the program assumes the same value all over: e.g. nu[] = 1/mu0;
  // - if we have only one region for a particular value, you can put write: nu[MyRegion]=  1/(mu0*mur_steel);
  // - if you have more than one region with the same characteristic: nu[Region[{MyRegion1, MyRegion2}]] = 1/(mu0*mur_steel);
  // ATTENTION: You can't define the function twice for a Physical Region, you will get a runtime error
  nu[Region[{Inds}]]=1./mu0;
  nu[Region[{XLPE,Soil,CU}]]=1./mu0;
  nu[Region[{OuterPVC}]]=1./(mu0*mur_pvc);

  sigma[OuterPVC]=sigma_pvc;
  sigma[CU]=sigma_cu;
  sigma[XLPE]= sigma_xlpe;
  sigma[Soil]=sigma_soil;
  sigma[Inds]= sigma_cu;


  epsilon[Region[{Inds,Soil,CU}]]=eps0;
  epsilon[Region[{XLPE}]]=eps0*espr_xlpe;
  epsilon[Region[{OuterPVC}]]=eps0*espr_pvc;
 
  DefineConstant[
    Freq = {50, Min 0, Max 1000000, Step 50000, Name "{00Parameters/Frequency", Highlight "Red" }
  ];
  Omega = 2*Pi*Freq;


  Pa = 0.;
  I = 410; // maximum value current in data sheet
  DefineFunction [js0];

  Ns[]= 1;
  Sc[]= SurfaceArea[];
  // Using second order hierarchical basis functions if set to 2
 Flag_Degree_a = 1;
  Flag_Degree_v = 1;

}

Constraint {
  // All the constraint hereafter must be adapted to your problem. Commented definitions are kept as example.

  // Electrical constraints
  { Name ElectricScalarPotential;
    Case {
      {Region Ind_1; Value V0; TimeFunction F_Cos_wt_p[]{2*Pi*Freq,Pa};}
      {Region SurfaceGe0; Value 0;}
    }
  }

  { Name ZeroElectricScalarPotential; // Only if second order basis functions
    Case {
      For k In {1:NbWires}
      {Region Ind~{k}; Value 0;}
      EndFor
      {Region SurfaceGe0; Value 0;}

      }
  }

  //Magnetic constraints
  { Name MagneticVectorPotential_2D;
    Case {
      {Region SurfaceGe0; Value 0.;}
    }
  }
  { Name Voltage_2D;
    Case {}}
  { Name Current_2D;
  Case{
//contraint used if Inds in DomainS_Mag
    {Region Ind_1; Value I; TimeFunction F_Cos_wt_p[]{2*Pi*Freq,Pa};}

    //Usefeul in DomainCWtihI_Mag
    // {Region Steel; Value 0.;}
    // {Region APLSheath; Value 0.;}
    }}
  }

Include "Jacobian_Integration.pro"; // Normally no modification is needed

// The following files contain: basis functions, formulations, resolution, post-processing, post-operation
// Some adaptations may be needed
If (Flag_AnalysisType ==0)
  Include "electrostatic_formulation.pro";
EndIf
If (Flag_AnalysisType ==1)
  Include "darwin_formulation.pro";
EndIf
