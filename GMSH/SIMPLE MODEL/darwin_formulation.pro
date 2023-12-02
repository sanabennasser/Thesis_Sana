FunctionSpace {

  { Name Hcurl_a_Mag_2D; Type Form1P;
    BasisFunction {
      { Name se; NameOfCoef ae; Function BF_PerpendicularEdge;
        Support Domain_Mag; Entity NodesOf[ All ]; }
    }
    Constraint {
      { NameOfCoef ae;
        EntityType NodesOf; NameOfConstraint MagneticVectorPotential_2D; }
    }
  }

  { Name Hregion_i_2D ; Type Vector ;
    BasisFunction {
      { Name sr ; NameOfCoef ir ; Function BF_RegionZ ;
        Support DomainS_Mag ; Entity DomainS_Mag ; }
    }
    GlobalQuantity {
      { Name Is ; Type AliasOf        ; NameOfCoef ir ; }
      { Name Us ; Type AssociatedWith ; NameOfCoef ir ; }
    }
    Constraint {
      { NameOfCoef Us ; EntityType Region ; NameOfConstraint Voltage_2D ; }
      { NameOfCoef Is ; EntityType Region ; NameOfConstraint Current_2D ; }
    }
  }

}

Formulation {

  { Name Darwin_a_2D; Type FemEquation; // Magnetodynamics + displacement current, no coupling
    Quantity {
      { Name a;  Type Local; NameOfSpace Hcurl_a_Mag_2D; }

      { Name ir ; Type Local  ; NameOfSpace Hregion_i_2D ; }
      { Name Us ; Type Global ; NameOfSpace Hregion_i_2D[Us] ; }
      { Name Is ; Type Global ; NameOfSpace Hregion_i_2D[Is] ; }
    }

    Equation {
      Galerkin { [ nu[] * Dof{d a} , {d a} ];
        In Domain_Mag; Jacobian Vol; Integration I1; }
      Galerkin { DtDof [ sigma[] * Dof{a} , {a} ];
        In DomainC_Mag; Jacobian Vol; Integration I1; }
      Galerkin { DtDtDof [ epsilon[] * Dof{a} , {a} ]; // Added term => Darwin approximation
        In DomainC_Mag; Jacobian Vol; Integration I1; }

      // Either you impose directly the function js0[]
      Galerkin { [ -js0[] , {a} ];
        In DomainS0_Mag; Jacobian Vol; Integration I1; }

      // or you use the constraints => allows accounting for sigma[{T}]
      Galerkin { [ -Ns[]/Sc[] * Dof{ir}, {a} ] ;
        In DomainS_Mag ; Jacobian Vol ; Integration I1 ; }
      Galerkin { DtDof [ Ns[]/Sc[] * Dof{a}, {ir} ] ;
        In DomainS_Mag ; Jacobian Vol ; Integration I1 ; }

      Galerkin { [ Ns[]/Sc[] / sigma[] * Ns[]/Sc[]* Dof{ir} , {ir} ] ; // resistance term
      In DomainS_Mag ; Jacobian Vol ; Integration I1 ; }
    //  GlobalTerm { [ 0.124 * Dof{Is} , {Is} ] ; In DomainS_Mag ; } // OR this resitance term
      GlobalTerm { [ Dof{Us}, {Is} ] ; In DomainS_Mag ; }

    }
  }

}


Resolution {

  { Name Darwin;
    System {
      { Name Sys_Mag; NameOfFormulation Darwin_a_2D;
        Type Complex; Frequency Freq; 
      }
    }
    Operation {
      CreateDir["res"];

      InitSolution[Sys_Mag];
      Generate[Sys_Mag]; Solve[Sys_Mag]; SaveSolution[Sys_Mag];
      PostOperation[Mag_Maps];
      PostOperation[Mag_Global];
    }
  }

}


PostProcessing {

  { Name Darwin_a_2D; NameOfFormulation Darwin_a_2D;
    PostQuantity {
      { Name a; Value { Term { [ {a} ]; In Domain_Mag; Jacobian Vol; } } }
      { Name az; Value { Term { [ CompZ[{a}] ]; In Domain_Mag; Jacobian Vol; } } }
      { Name b; Value { Term { [ {d a} ]; In Domain_Mag; Jacobian Vol; } } }
      
      { Name bm; Value { Term { [ Norm[{d a}] ]; In Domain_Mag; Jacobian Vol; } } }
      { Name j; Value { Term { [ -sigma[]*Dt[{a}] ]; In DomainC_Mag; Jacobian Vol; } } }
      { Name jz; Value { Term { [ CompZ[-sigma[]*Dt[{a}]] ]; In DomainC_Mag; Jacobian Vol; } } }
      { Name jm; Value { Term { [ Norm[-sigma[]*Dt[{a}]] ]; In DomainC_Mag; Jacobian Vol; } } }

      { Name d; Value { Term { [ epsilon[] * Dt[Dt[{a}]] ]; In DomainC_Mag; Jacobian Vol; } } }
      { Name dz; Value { Term { [ CompZ[epsilon[] *  Dt[Dt[{a}]] ] ]; In DomainC_Mag; Jacobian Vol; } } }
      { Name dm; Value { Term { [ Norm[epsilon[]  *  Dt[Dt[{a}]] ] ]; In DomainC_Mag; Jacobian Vol; } } }

      { Name JouleLosses; Value {
          Integral { [ 0.5*sigma[]*SquNorm[Dt[{a}]] ]   ; In DomainC_Mag  ; Jacobian Vol ; Integration I1 ; }
          Integral { [ 0.5/sigma[]*SquNorm[js0[]] ]          ; In DomainS0_Mag ; Jacobian Vol ; Integration I1 ; }
          Integral { [ 0.5/sigma[]*SquNorm[Ns[]/Sc[]*{ir}] ] ; In DomainS_Mag  ; Jacobian Vol ; Integration I1 ; }
        }
      }

      { Name U ; Value {
          Term { [ {Us} ] ; In DomainS_Mag ; }
        }
      }

      { Name I ; Value {
          Term { [ {Is} ] ; In DomainS_Mag ; }
        }
      }

      { Name S ; Value {
          Term { [ {Us}*Conj[{Is}] ] ; In DomainS_Mag ; }
        }
      }

      { Name R ; Value {
          Term { [ -Re[{Us}/{Is}*1e3] ] ; In DomainS_Mag ; }
        }
      }

      { Name L ; Value {
          Term { [ -Im[{Us}/{Is}]/(2*Pi*Freq) *1e6] ; In DomainS_Mag ; }
        }
      }
        //     { Name R_from_losses ; 
        // Value { Term { Type Global; [ $losses/SquNorm[$current]/0.5 ] ; In DomainDummy ; } } }

      // { Name MagneticEnergy; Value { 
      //     Integral { [ 0.5* nu[] * SquNorm[{d a}] ]; 
      //     In DomainS_Mag; Jacobian Vol; Integration I1;}
      //   }
      // }

      // { Name Inductance_from_MagEnergy ; // three phase system => 1/3
      //   Value { Term { Type Global; [  1e6*2*$magenergy/SquNorm[$current] ] ; In DomainDummy ; } } }

    }
  }

}


PostOperation{
  // Magnetic
  //-------------------------------

  { Name Mag_Maps; NameOfPostProcessing Darwin_a_2D;
    Operation {
      // You may want to see the result only in part of the domain => adapt domain that follows OnElementsOf
      // Name is not compulsory, it can be modified
      Print[ bm , OnElementsOf Domain_Mag,
        Name "|B| [T]", File "res/bm.pos" ];
      Print[ b , OnElementsOf Domain_Mag,
        Name "B [T]", File "res/b.pos" ];  

         Print[ az , OnElementsOf Domain_Mag,
        Name " Az [T]", File "res/az.pos" ];     

      Print[ jm , OnElementsOf DomainC_Mag,
        Name "|j| [A/m^2] Al sheath", File "res/jm.pos" ];
      Print[ dm , OnElementsOf DomainC_Mag,
        Name "|D| [A/m²]", File "res/dm.pos" ];
    }
  }

  // This two definitions refer to the GUI => organisation of the results
  po = "{01Losses/";
  po2 = "{02PU-parameters/";
  { Name Mag_Global; NameOfPostProcessing Darwin_a_2D;
    Operation {
      // You may restrict DomainC_Mag to part of it
    Print[ JouleLosses[DomainC_Mag], OnGlobal, Format Table,
        SendToServer StrCat[po,"All conducting domain"],
        Units "W/m", File "res/losses_total.dat" ];
        Print[ JouleLosses[Inds], OnGlobal, Format Table,
          SendToServer StrCat[po,"Losses in conductor"],
          Units "W/m", File "res/losses_inds.dat" ];
          	Print[ JouleLosses[Ind_1], OnGlobal, Format Table, StoreInVariable $losses,
        	SendToServer "Output/10Joule_loss_ind1", File "res/losses_ind_1.dat" ];
      //        	Print[ R_from_losses, OnRegion DomainDummy, Format Table,
     	//File "res/RRinds.dat",
    //	SendToServer StrCat["po/","20R1 from losses"], Color "LightYellow" ];
        Print[ R, OnRegion Ind_1, Format Table,
          SendToServer StrCat[po2,"0R"],
          Units "Ω/km", File "res/Rinds.dat" ];
        Print[ L, OnRegion Ind_1, Format Table,
          SendToServer StrCat[po2,"1L"],
          Units "mH/km", File "res/Linds.dat" ];
        Print[ U, OnRegion Ind_1, Format Table,
          SendToServer StrCat[po2,"2U_1"],
          Units "V", File "res/Uinds.dat" ];
        Print[ I, OnRegion Ind_1, Format Table,
          SendToServer StrCat[po2,"4I_2"],
          Units "A", File "res/Iinds.dat" ];
          	Print[ I, OnRegion Ind_1, Format Table, StoreInVariable $current, Color "Pink", 
        	SendToServer "Output/001Re(I_ind1)"{0}, File"res/Iinds.dat"]; // real part in GUI and complex value stored in $current
        Print[ S, OnRegion Ind_1, Format Table,
          SendToServer StrCat[po2,"6S"],
          Units "VA", File "res/Sinds.dat" ];
	// Print[ MagneticEnergy[Region[{Domain_Mag}]], OnGlobal, Format Table, StoreInVariable $magenergy,
  //       	SendToServer "Output/60MagEnergy", File "me.pos", Color "Orange" ];
  //     	Print[ Inductance_from_MagEnergy, OnRegion DomainDummy, Format Table,
  //       	File "L_from_me.dat",
  //       	SendToServer StrCat["Output/","61L from ME"], Color "Yellow" ];
    }
  }
}
