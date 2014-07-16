Cases Folder README
=====================

## Using the `makeRunCase` Script

This folder contains templates for various useful cases to run with dropletFoam.
The cases do not have generated files in them. To keep the case folder clean,
when you want to run one of the cases, use the `makeRunCase.py` Python script.
For example, to make a run folder for `vialTest`, use:

    python makeRunCase.py vialTest
  
which will create a `../run/vialTest` folder (relative to the case
directory). If there is already a `vialTest` case in the `run` folder, it will
give an error. To force it to overwrite the existing case, use

    python makeRunCase.py vialTest -f
  
To view the help menu for the makeRunCase script, use

    python makeRunCase.py -h
  
If you want to rename the run case to something different from the template, you
can do so with

    python makeRunCase.py vialTest vialTestB

which can be nice for setting up some parametric runs. After creating all your
run cases, you can go edit the relevant parameter in each case in the `run`
folder without altering the base template.
  

## List of Template Cases

There are several case templates currently in the case folder. As the solver
has evolved, some interfaces have changes so they may not all work. Cases which
have been tested to work with the current interface are listed below with a
short description.

1. `burningDrop`

2. `dropCaseGels`

3. `dropCaseLiquids`

4. `impingingJets`

5. `vialTest`

6. `dropPoolRFNA`



