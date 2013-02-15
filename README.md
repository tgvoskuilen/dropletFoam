dropletFoam
===========

Compressible, multiphase OpenFOAM combustion solver


Solution Flow
-----------------------
## Problem Setup

  *  Create mixtureThermo
     *  create vapor phase
        *   create fields
        *   create subspecies
            *   link with global mass fraction Y
            *   create local mass fraction Yp
     *  create liquid phase
        *   create fields
        *   create subspecies
            *   link with global mass fraction Y
            *   create local mass fraction Yp
            *   create viscosity model
     *  create evaporation models in liquid phase subspecies
        
     *  calculate thermo density, set to old time slot
        
     *  set hs_ based on T_
        
     *  calculate()
        *   Updates T_ based on hs_
        *   updates alpha_ (not used?)
        *   liquid.correct(p,T)
            *   updates rhoAlpha_
        *   vapor.correct(p,T)
            *   updates rhoAlpha_
        *   updates mu_
        *   updates rho_
        *   updates psi_
        
     *  set vapor specie relative Ys
        *   Yp = Y / sum(Y in phase)
            
     *  set liquid specie relative Ys
        *   Yp = Y / sum(Y in phase)
            
## Solution Routine

    mixture.solve
        update alphaVaporSharp_
        combustion->correct
        calc evaporation
        
        solve alphas
            updates rhoPhi_ and rhoPhiAlpha_ in phases
        update alphaVaporSharp_
        calculate surface tension
        update rhoE to satisfy continuity with rhoPhi_
        
        solve subspecies
---phase rho-rhoPhi pair only needed valid in here-----
            update rho_ to satisfy rhoPhiAlpha_ ?
                really should update alpha*rho_ to satisfy...
            solve Yp flux using rho_ and rhoPhiAlpha_ pair?


---end validity zone-----------------------------------
    
    UEqn
        needs rhoE-rhoPhi pair
    TEqn
        needs rhoE-rhoPhi pair
    mixture.setHs(Te)
    thermo.correct
        updates rho_ in thermo and phases
        updates mu_
        updates psi_
        
        
    
    
        
