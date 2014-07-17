
Chemistry Database
===========================
 
This module builds liquid/gas chemistry input files from a selection of 
reaction sets. Details about the individual sets are included below. To build a
multiphase chemistry set in a python script, use:

    import chemBuilder
    reaction_set = chemBuilder.build_set('GasChem2','LiquidChem2')
    reaction_set.write_file('chem.inp')

The inclusion of liquid chemistry is optional, and can be omitted. Additional
species can also be added. A few example uses are given below:

    reaction_set = chemBuilder.build_set('GasChem1')
    reaction_set = chemBuilder.build_set(gas='GasChem2',extra_species=['CH3NHNH2L'])
    reaction_set = chemBuilder.build_set(gas='GasChem3',liquid='LiquidChem1',extra_species=['SI'])
    

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

Adding Tsang & Herron Reactions
======================================
Several of the mechanisms above require the use of Tsang & Herron reactions,
which are not a part of OpenFOAM by default. The following section provides instructions
for how to add them to OpenFOAM so you can use the mechanisms which require them.

The Tsang & Herron reaction form is a simplified version of Troe, which is already
a part of OpenFOAM. The difference between these lies in the calculation of 
F for the pressure-dependence. The line numbers and file locations are from
OpenFOAM-2.1.x, so they may differ slightly from other versions. For most
modifications you will just be searching for `Troe` and replacing it with
`TsangHerron`. 

## Part 1: Adding a TsangHerron fallOffFunction to OpenFOAM

1. Locate the `fallOffFunctions` folder. In my version, it is at
   `thermophysicalModels/specie/reaction/reactionRate/fallOffFunctions`
    
2. Copy the `TroeFallOffFunction` folder and rename the new folder 
   `TsangHerronFallOffFunction`
    
3. Rename `TroeFallOffFunction.H` and `TroeFallOffFunctionI.H` in the new 
   folder to `TsangHerronFallOffFunction.H` and `TsangHerronFallOffFunctionI.H`
  
4. Open `TsangHerronFallOffFunction.H` and make the following edits:

  1. Replace all `Troe` with `TsangHerron`

  2. In the `// Private data` section, remove the 
     `scalar alpha_; scalar Tsss_, Tss_, Ts_;` lines and add the line 
     `scalar a0_, a1_`;

  3. In the `inline TsangHerronFallOffFunction` input list, change the inputs to
     ```
     const scalar a0,
     const scalar a1
     ```


5. Open `TsangHerronFallOffFunctionI.H` and make the following edits:

  1. Replace all `Troe` with `TsangHerron`
   
  2. Change the three constructor functions to be:
     ```
        inline Foam::TsangHerronFallOffFunction::TsangHerronFallOffFunction
        (
            const scalar a0,
            const scalar a1
        )
        :
        a0_(a0),
        a1_(a1)
        {}
     ```
     ```
     inline Foam::TsangHerronFallOffFunction::TsangHerronFallOffFunction(Istream& is)
     :
         a0_(readScalar(is.readBegin("TsangHerronFallOffFunction(Istream&)"))),
         a1_(readScalar(is))
     {
         is.readEnd("TsangHerronFallOffFunction(Istream&)");
     }
     ```
     ```
     inline Foam::TsangHerronFallOffFunction::TsangHerronFallOffFunction(const dictionary& dict)
     :
     a0_(readScalar(dict.lookup("a0"))),
     a1_(readScalar(dict.lookup("a1")))
     {}
     ```
  3. Change the `operator()` function to use the Tsang and Herron approach
```   
    inline Foam::scalar Foam::TsangHerronFallOffFunction::operator()
    (
        const scalar T,
        const scalar Pr
    ) const
    {
        scalar logFcent = log10(max(a0_ + a1_*T, SMALL));

        scalar logPr = log10(max(Pr, SMALL));
        return pow(10.0, logFcent/(1.0 + sqr(logPr)));
    }
```
    4. Change the next two writer functions to be
``` 
    inline void Foam::TsangHerronFallOffFunction::write(Ostream& os) const
    {
        os.writeKeyword("a0") << a0_ << token::END_STATEMENT << nl;
        os.writeKeyword("a1") << a1_ << token::END_STATEMENT << nl;
    }

    inline Foam::Ostream& Foam::operator<<
    (
        Foam::Ostream& os,
        const Foam::TsangHerronFallOffFunction& tfof
    )
    {
        os  << token::BEGIN_LIST
            << tfof.a0_
            << token::SPACE << tfof.a1_
            << token::END_LIST;

        return os;
    }
```

5. Edit the reaction macros
    1. Go to `thermophysicalModels/specie/reaction/reactions`
    2. Open `makeChemkinReactions.C` and make the following edits:
        1. Add `#include "TsangHerronFallOffFunction.H"` after the line with the 
           Troe function
        2. Add the following around line 80 (copy the Troe version and replace 
           Troe with TsangHerron)
```
    makePressureDependentReactions
    (
        gasThermoPhysics,
        ArrheniusReactionRate,
        TsangHerronFallOffFunction
    )
```
  * Open `makeReactionThermoReactions.C` and do the following edits
    * Add `#include "TsangHerronFallOffFunction.H"` after the line with the 
      Troe function
    * Add the following lines to the macro after the Troe part, remembering 
      the `\` at the end of each line
```
    makePressureDependentReactions                                             \
    (                                                                          \
       Thermo,                                                                 \
       ArrheniusReactionRate,                                                  \
       TsangHerronFallOffFunction                                              \
    )                                                                          \
```
* Recompile the specie library by going to `thermophysicalModels/specie` and
  running `wclean`, `rmdepall`, and `wmake libso`
    
Congratulations, you have added a new reaction type to openfoam. Next up is to 
modify the Chemkin reader so it can read the T&H inputs from a chemkin 
formatted input and create the newly defined reaction.


## Part 2: Adding the T&H Form to the Chemkin Reader

* Go to `thermophysicalModels/reactionThermo/chemistryReaders/chemkinReader`
* Open `chemkinReader.H` and make the following edits:
  * Add `TsangHerronReactionType` after the `TroeReactionType` in the `enum reactionKeyword`
  * Add `TsangHerron` after `Troe` in the `enum fallOffFunctionType`
  * Change the size of `fallOffFunctionNames` to 5

* Open `chemkinReader.C` and make the following edits:
  * Add `#include "TsangHerronFallOffFunction.H"` after the Troe include around line 37
  * Add `TsangHerron` after `Troe` in the `fallOffFunctionNames` around 
    line 78 and change the `4` in the size to a `5`
  * In the `initReactionKeywordTable` after the line with `TROE`, add the line 
    to look for the `TH` keyword in the chemkin file (note that using `T&H` 
    here results in errors, since the parser hangs on the `&` character for some reason)
```
reactionKeywordTable_.insert("TH", TsangHerronReactionType);
```
  * Locate the case structure with a `Troe` entry around line 287 (search 
    for `Troe`) and copy the entire `Troe` case. Pay particular attention to 
    the numbers in the section below. The T&H form needs 1 or 2 
    coefficients, while the Troe needs more. Edit it to be:
```
        case TsangHerron:
        {
            scalarList TsangHerronCoeffs
            (
                reactionCoeffsTable[fallOffFunctionNames[fofType]]
            );

            if (TsangHerronCoeffs.size() != 2 && TsangHerronCoeffs.size() != 1)
            {
                FatalErrorIn("chemkinReader::addPressureDependentReaction")
                    << "Wrong number of coefficients for T&H rate expression"
                       " on line " << lineNo_-1 << ", should be 1 or 2 but "
                    << TsangHerronCoeffs.size() << " supplied." << nl
                    << "Coefficients are "
                    << TsangHerronCoeffs << nl
                    << exit(FatalError);
            }

            if (TsangHerronCoeffs.size() == 1)
            {
                TsangHerronCoeffs.setSize(2);
                TsangHerronCoeffs[1] = 0.0;
            }

            addReactionType
            (
                rType,
                lhs, rhs,
                PressureDependencyType
                    <ArrheniusReactionRate, TsangHerronFallOffFunction>
                (
                    ArrheniusReactionRate
                    (
                        Afactor0*k0Coeffs[0],
                        k0Coeffs[1],
                        k0Coeffs[2]/RR
                    ),
                    ArrheniusReactionRate
                    (
                        AfactorInf*kInfCoeffs[0],
                        kInfCoeffs[1],
                        kInfCoeffs[2]/RR
                    ),
                    TsangHerronFallOffFunction
                    (
                        TsangHerronCoeffs[0],
                        TsangHerronCoeffs[1]
                    ),
                    thirdBodyEfficiencies(speciesTable_, efficiencies)
                )
            );
        }
        break;
```

* Open `chemkinLexer.L` and make the following edits:
  * Locate the `TroeReactionType` case entry around line 926 (search for Troe)
    and copy the entire Troe case entry, modifying it to be:
```
        case TsangHerronReactionType:
        {
            if (!pDependentSpecieName.size())
            {
                FatalErrorIn("chemkinReader::lex()")
                    << "T&H keyword given for a"
                       " reaction which does not contain a pressure"
                       " dependent specie" << " on line " << lineNo_
                    << exit(FatalError);
            }

            if
            (
                fofType == unknownFallOffFunctionType
             || fofType == Lindemann
            )
            {
                fofType = TsangHerron;
            }
            else
            {
                FatalErrorIn("chemkinReader::lex()")
                    << "Attempt to set fall-off function type to T&H"
                       " when it is already set to "
                    << fallOffFunctionNames[fofType]
                    << " on line " << lineNo_
                    << exit(FatalError);
            }

            reactionCoeffsName = fallOffFunctionNames[fofType];
            BEGIN(readReactionCoeffs);
            break;
        }
```

* Recompile using the following steps
  * Go up to the `reactionThermo` folder and run `wclean`, `rmdepall`, 
    then `wmake libso` to recompile the library. (There are always a load of 
    warnings from `chemkinLexer.L` about old-style casts. You can ignore those.)
    
  * Go to the `applications/utilities/thermophysical/chemkinToFoam` folder 
    and recompile with �wmake�
       




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

