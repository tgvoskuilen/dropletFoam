
Chemistry Database
===========================
 
This module builds liquid/gas chemistry input files from a selection of 
reaction sets. Details about the individual sets are included below. To build a
multiphase chemistry set in a python script, use:

    import chemBuilder
    reaction_set = chemBuilder.build_set('GasChem2','LiquidChem2')
    reaction_set.write_file('chem.inp')

Gas Reaction Sets
=========================
The following reaction sets are for gaseous reaction with MMH and various 
oxidizers

## GasChem0

| Species | Reactions | References   |
| ------- | --------- | ------------ |
| 2       | 1         | N/A          |

    
Empty gas mechanism with just N and N2 species and HCON elements as basis 
for non-reacting cases. The N2<=>N decomposition reaction is included since
the chemkinToFoam converter crashes with 0 reactions.
    

## GasChem1

| Species | Reactions | References   |
| ------- | --------- | ------------ |
| 80      | 512       | 1            |

This is the original ARL mechanism, but without the 'NAMMH' reaction. Note, this
requires you to have implemented the T&H reaction type and requires CH2(S), which
can sometimes cause problems because of the parentheses in its name.


## GasChem2

| Species | Reactions | References   |
| ------- | --------- | ------------ |
| 25      | 98        | 2, 3         |

This is Labbe's first reduced version of the ARL mechanism.


## GasChem3

| Species | Reactions | References   |
| ------- | --------- | ------------ |
| 29      | 120       | 4            |

This is Labbe's second reduced version of the ARL mechanism. This requires
CH2(S), which can sometimes cause problems because of the parentheses in 
its name.


## GasChem4

| Species | Reactions | References   |
| ------- | --------- | ------------ |
| 41      | 198       | 3, 5         |

This is the reduced version of Labbe's new mechanism.


## GasChem5

| Species | Reactions | References   |
| ------- | --------- | ------------ |
| 41      | 194       | 3            |

This is the reduced and stabilized version of Labbe's new mechanism, which is
the same as GasChem4 but with 4 stiff reactions removed


## GasChem6

| Species | Reactions | References   |
| ------- | --------- | ------------ |
| 35      | 49        | 6            |


This is a reduced version of the original ARL mechanism provided by Mike
McQuaid from ARL. The only modification made internally is the removal of 
the NAMMH reaction. Note, this requires you to have implemented the T&H 
reaction type.

   
    
Liquid reaction sets
===================================
The following liquid reaction sets are predominantly to enable reaction pathways
between MMH and nitric acid in liquid phase.

    
## LiquidChem1

| Species | Reactions | References   |
| ------- | --------- | ------------ |
| 9       | 1         | 3            |

This is a global one-step reaction between MMH and nitric acid producing
heat and intermediates based on some observations by Thynell (ref). It also
includes liquid NO2 so it can evaporate and work with RFNA, but there is no
reaction between liquid MMH and liquid NO2
    
## LiquidChem2

| Species | Reactions | References   |
| ------- | --------- | ------------ |
| 11      | 2         | 3            |

This includes the global 1-step reaction from LiquidChem1 and the forward
aerosol reaction.


References
=====================================
1. W. Anderson, M. McQuaid, M. Nusca, A. Koltar, "A detailed,
   finite-rate, chemical kinetics mechanism for monomethylhydrazine
   red fuming nitric acid systems," Army Research Laboratory,
   ARL-TR-5088, February 2010
   
2. N. Labbe, Y. Kim, P. Westmoreland, "Computational mechanism
   development for hypergolic propellant systems: MMH and DMAZ,"
   Proceedings of the 2010 AIChE Annual Meeting, Salt Lake City, UT
               
3. N. Labbe, C. Needham, P. Westmoreland, T. Voskuilen, T. Pourpoint,
   "Predictive Modeling of Hypergolic-Propellant Performance. Part I:
   Reaction-Kinetics Model for Monomethylhydrazine and Nitric Acid Liquids
   and Vapors," 2014, In Preparation

4. N. Labbe, P. Westmoreland, Internal communications, MURI Task 3.2

5. N. Labbe, "Determining detailed reaction kinetics for nitrogen
   and oxygen containing fuels," Ph.D. Thesis, University of
   Massachusetts Amherst, February 2013

6. M. McQuaid, Internal communication, Spring 2014

