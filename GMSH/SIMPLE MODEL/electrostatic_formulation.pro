FunctionSpace {

  { Name Hgrad_v_Ele; Type Form0;
    BasisFunction {
      // v = \sum_n v_n  s_n,  for all nodes
      { Name sn; NameOfCoef vn; Function BF_Node;
        Support Domain_Ele; Entity NodesOf[ All ]; }

    }

    Constraint {
      { NameOfCoef vn; EntityType NodesOf;
        NameOfConstraint ElectricScalarPotential; }

    }

  }

}

Formulation {
//   { Name Electrostatic_v; Type FemEquation;
//     Quantity {
//       { Name v; Type Local; NameOfSpace Hgrad_v_Ele; }
//     }
//     Equation {
//      Galerkin { [ sigma[] * Dof{d v} , {d v} ] ; //commenting out this line gives a wrong capacitance (0.03 instead of 0.18uF/km)
//         In Domain_Ele; Jacobian Vol ; Integration I1 ; } 
//       Galerkin { [ Dt[epsilon[] * Dof{d v} , {d v}] ] ;
//         In Domain_Ele; Jacobian Vol ; Integration I1 ; }

//     }
//   }

// }
{ Name Electrodynamics_v; Type FemEquation;
    Quantity {
      { Name v; Type Local; NameOfSpace Hgrad_v_Ele; }
    }
    Equation {
      Galerkin { [ sigma[] * Dof{d v} , {d v} ] ;
        In Domain_Ele; Jacobian Vol ; Integration I1 ; }
      Galerkin { DtDof[ epsilon[] * Dof{d v} , {d v} ];
        In Domain_Ele; Jacobian Vol; Integration I1; }
    }
  }
}

Resolution {
  // { Name Electrostatics;
  //   System {
  //     { Name Sys_Ele; NameOfFormulation Electrostatic_v;
  //     }
  //   }
  //   Operation {
  //     CreateDir["res"];

  //     Generate[Sys_Ele]; Solve[Sys_Ele]; SaveSolution[Sys_Ele];
  //     PostOperation[Ele_Maps];
  //   }
  // // }
  { Name Electrodynamics;
    System {
      { Name Sys_Ele; NameOfFormulation Electrodynamics_v;
        Type Complex; Frequency Freq; }
    }
    Operation {
      CreateDir["res"];

      Generate[Sys_Ele]; Solve[Sys_Ele]; SaveSolution[Sys_Ele];
      PostOperation[Ele_Maps];
     // PostOperation[Ele_Cuts]; // To adapt for your cable
    }
  }

}


PostProcessing {

            { Name EleSta_v; NameOfFormulation Electrodynamics_v;
              Quantity {

               { Name v; Value { Term { [ {v} ]; In Domain_Ele; Jacobian Vol; } } }
                { Name e; Value { Term { [ -{d v} ]; In Domain_Ele; Jacobian Vol; } } }
                { Name em; Value { Term { [ Norm[-{d v}] ]; In Domain_Ele; Jacobian Vol; } } }

                 { Name d; Value { Term { [ -epsilon[] * {d v} ]; In Domain_Ele; Jacobian Vol; } } }
                 { Name dm; Value { Term { [ Norm[-epsilon[] * {d v}] ]; In Domain_Ele; Jacobian Vol; } } }

                { Name ElectricEnergy; Value {
                    Integral {
                      [ 0.5 * epsilon[] * SquNorm[{d v}] ];
                      In Domain_Ele; Jacobian Vol; Integration I1;
          }
	}
      }

      { Name V0 ; Value {
          // For recovering the imposed voltage in post-pro
          // Most likely you will need to adapt for your cable
          // The default hereafter is for a three-phase cable
          Term { Type Global ; [ V0 * F_Cos_wt_p[]{2*Pi*Freq, Pa}] ; In Ind_1 ; }
        } }

      { Name C_from_Energy ; Value { Term { Type Global; [1e9* 2*$We/SquNorm[$voltage] ] ; In DomainDummy ; } } }
    }
  }}



PostOperation{

  // Electric
  //-------------------------------

  po0 = "{01Capacitance/"; // Useful only for the GUI

  { Name Ele_Maps; NameOfPostProcessing EleSta_v;
    Operation {
      Print[ v,  OnElementsOf Domain_Ele, File "res/v.pos" ];
      Print[ em, OnElementsOf Domain_Ele, Name "|E| [V/m]",  File "res/em.pos" ]; // Name is not compulsory, it may be adapted
      Print[ dm, OnElementsOf Domain_Ele, Name "|D| [A/m²]", File "res/dm.pos" ];
      Print[ ElectricEnergy[Domain_Ele], OnGlobal, Format Table, StoreInVariable $We,
        SendToServer StrCat[po0,"0Electric energy"], File "res/energy.dat" ];

      Print[ V0, OnRegion Ind_1, Format Table, StoreInVariable $voltage,
        SendToServer StrCat[po0,"0U1"], Units "V", File "res/U.dat" ];
      Print[ C_from_Energy, OnRegion DomainDummy, Format Table, StoreInVariable $C1,
        SendToServer StrCat[po0,"1Cpha"], Units "uF/km", File "res/C.dat" ];
    }
  }

}
