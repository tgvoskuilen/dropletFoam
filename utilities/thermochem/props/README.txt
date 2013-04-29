Raw Data Sources
    transport.dat    <-  viscosity and conductivity polynomials
    lennardjones.csv <-  Lennard jones coefficients, can be used for Dij and viscosity (not used)
    thermo.dat       <-  JANAF polynomials for H, Cp, and S, as well as atomic composition and molar mass
    
    
"processed" data sources
    thermo.pydb  <-  Python-readable database for creating "thermo". Contains molar mass, 
                     JANAF polynomials, and Sutherland coefficients. Create by calling
                     updateThermoDb.py
