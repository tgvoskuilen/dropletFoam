Properties Database
=====================================
This folder contains the database of thermodynamic and transport properties
used in the simulations.

# Raw Data Sources
The raw data is contained in the three files described below.

### `transport.dat`

Viscosity and conductivity NASA polynomials


### `lennardjones.csv`

Lennard jones coefficients, can be used for Dij and viscosity (not used)


### `thermo.dat`

JANAF polynomials for H, Cp, and S, as well as atomic composition and molar mass


# Processed Data

To combine the raw data into a single readable format, use the `updateThermoDb`
function to generate the `thermo.pydb` file. If any of the inputs change (any of
the raw data listed above) re-run the script with `./updateThermoDb`. This
processed data file contains a Python-readable database with molar mass, JANAF
polynomials, and Sutherland coefficients (fit to the NASA polynomials).

To run the processing script and see plots of all the transport property 
fitting, use `./updateThermoDb -v`.
