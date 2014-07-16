Cases Folder README
=====================

## Using the `makeRunCase` Script

This folder contains templates for various useful cases to run with dropletFoam.
The cases do not have generated files in them. To keep the case folder clean,
when you want to run one of the cases, use the `makeRunCase` Python script. You
may have to make the script executable with `chmod +x makeRunCase` first.
For example, to make a run folder for `vialTest`, use:

    ./makeRunCase vialTest
  
which will create a `../run/vialTest` folder (relative to the case
directory). If there is already a `vialTest` case in the `run` folder, it will
give an error. To force it to overwrite the existing case, use

    ./makeRunCase -f vialTest
  
To view the help menu for the makeRunCase script, use

    ./makeRunCase -h
  
If you want to rename the run case to something different from the template, you
can do so with

    ./makeRunCase.py -o vialTestB vialTest

or, combining options,

    ./makeRunCase.py -f -o vialTestB vialTest
    
which can be nice for setting up some parametric runs. After creating all your
run cases, you can go edit the relevant parameter in each case in the `run`
folder without altering the base template.
  

## List of Template Cases

There are several case templates currently in the case folder. As the solver
has evolved, some interfaces have changes so they may not all work. Cases which
have been tested to work with the current interface are listed below with a
short description.

1. `burningDrop`

This is a 2D axisymmetric case with a liquid MMH drop hanging by surface
tension force on a small cylinder. There is a hot mixture of N2/N2O flowing
over it and it burns. The case does not use dynamic mesh adaptation, and
does not require liquid reactions since the oxidizer is in the gas phase.

2. `dropCaseGels`

This is a 3D case with two droplets of 2 mm diameter, one of MMH and one of
RFNA, with liquid and gas reactions. The droplets are shear-thinning gels 
of Aerosil (intert species SI added).

3. `dropCaseLiquids`

This is a 3D case with the same conditions as `dropCaseGels` except that the
droplets are Newtonian liquids rather than gels.

4. `impingingJets`

This is a 3D case focused on the impingement region of the impinging jets
setup, with one MMH jet and one RFNA jet.

5. `vialTest`

This is a 2D planar case to mimic the drop test experiments, with a 2mm diameter
drop of RFNA falling into a 1 mm deep pool of MMH. The domain is 1 cm wide
and 2 cm tall, with the top open to atmosphere and the sides as adiabatic
walls.


