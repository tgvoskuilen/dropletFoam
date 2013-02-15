dropletFoam
===========

Compressible, multiphase OpenFOAM combustion solver


Solution Flow
-----------------------
## Problem Setup

 * Create mixtureThermo
     * Create vapor phase
        * Create fields
        * Create subspecies
            * Link with global mass fraction Y
            * Create local mass fraction Yp
     * Create liquid phase
        * Create fields
        * Create subspecies
            * Link with global mass fraction Y
            * Create local mass fraction Yp
            * Create viscosity model
            
     * Create evaporation models in liquid phase subspecies
     * Calculate thermo density, set to old time slot
     * Set hs_ based on T_
     
     * Call calculate()
        * Updates T_ based on hs_
        * Updates alpha_ (not used?)
        * Liquid.correct(p,T)
            * Updates rhoAlpha_
        * Vapor.correct(p,T)
            * Updates rhoAlpha_
        * Updates mu_
        * Updates rho_
        * Updates psi_
        
     * Set vapor specie relative Ys
        * Yp = Y / sum(Y in phase)
            
     * Set liquid specie relative Ys
        * Yp = Y / sum(Y in phase)
            
## Solution Routine
 * Update rho to satisfy continuity with rhoPhi (not needed?)
 
 * Call `mixture.solve( rho )`
   * Call `combustion->correct()`
     * Uses current `thermo.rho` and Ys
   * Calculate evaporation rates
   * Solve VOF field advection (alpha)
     * Updates `rhoPhi_` and `rhoPhiAlpha_` in phases
     * Uses volume source due to evaporation
   * Update `alphaVaporSharp_`
   * Calculate surface tension fields (`kappa_` and `sigma_`)
   * Update `rho` to satisfy continuity with `rhoPhi_`
     * Solves `fvm::ddt(rho) + fvc::div(rhoPhi_) == 0`
     
   * Set phase masks (`faceMask` and `cellMask`)
   * Calculate diffusion coefficients using Schmidt number
     * Set `Df = Sc * muf * faceMask`
   
   * Solve subspecies
     * Apply `faceMask` to `rhoPhiAlpha_` to prevent advection into other phase
     * Update `rhoAlpha_` to satisfy continuity with masked `rhoPhiAlpha_`
       * `fvm::ddt(rhoAlpha_) + fvc::div(rhoPhiAlpha_*faceMask) == S_evap`
     * Solve Yp flux using `rhoAlpha_` and masked `rhoPhiAlpha_` pair
   * Set global mass fractions
     * Set `Y = Yp * rhoAlpha_ / (rhoAlpha_ + otherRhoAlpha_)`

 * UEqn
   * Needs rho-rhoPhi pair that satisfy continuity
 * TEqn
   * Needs rho-rhoPhi pair that satisfy continuity
   * Call `mixture.setHs(Te)`
   * Call `thermo.correct()`
     * Updates `rho_ = alphaL*rhoL(p,T) + alphaV*rhoV(p,T)`
     * Updates `mu_ = alphaV*muV + alphaL*muL(p,T)`
       * muV is calculated from Sutherland mixture of Ys
     * Updates `psi_ = alphaV*alphaV.psi(T)`
     * Updates `T_` using `hs_`
 * pEqn (inside pressure corrector loop)
   * `rho = thermo.rho`
   * `thermo.rho -= p*psi`
   * `solve for updated p`
   * `thermo.rho += p*psi`

 * Update external density field from thermo density field
   * rho = thermo.rho
   
    
    
        
