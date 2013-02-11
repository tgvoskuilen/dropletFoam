/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::PLICInterface

Description


SourceFiles
    PLICInterface.C

\*---------------------------------------------------------------------------*/

#ifndef PLICInterface_H
#define PLICInterface_H

#include "fvCFD.H"
#include "MULES.H"
#include "volFields.H"
#include "dictionaryEntry.H"
#include "dimensionedScalar.H"
#include "dimensionedScalarFwd.H"
#include "CPCCellToCellStencil.H"
#include "cuttableCell.H"
#include "cuttableFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                           Class phase Declaration
\*---------------------------------------------------------------------------*/

class PLICInterface
{
    // Private data

        //- Host mesh
        const fvMesh& mesh_;
        
        //- Parameter dictionary
        IOdictionary plicDict_;

        //- Volume fraction of Vapor
        volScalarField& alphaVapor_;
        
        //- Volume fraction of liquid
        volScalarField& alphaLiquid_;
        
        const volVectorField& U_;
        
        //- Face area fraction of vapor
        surfaceScalarField alphaf_;
        
        //- Volume flux of gas phase
        surfaceScalarField phiAlpha_;
        
        volScalarField sumalphaf_;
        
        //- Interface area
        volScalarField iArea_;
        
        //- Interface normal vector
        volVectorField iNormal_;
        
        //- Interface plane point
        volVectorField iPoint_;
        
        //- Centroid of gas portion of a cell (zero in liquid)
        volVectorField gasC_;
        
        //- Centroid of liquid portion of a cell (zero in gas)
        volVectorField liquidC_;
        
        volScalarField deltaV_;
        volScalarField deltaL_;
        surfaceScalarField wL_;
        surfaceScalarField wV_;
        
        volScalarField intermeds_;
        volScalarField alphaLsmooth_;
        volVectorField smoothN_;
        
        //TEMPORARY
        volScalarField liquidCells_;
        volScalarField smallLiquidCells_;
        volScalarField noLiquidCells_;
        volScalarField gasCells_;
        volScalarField smallGasCells_;
        volScalarField noGasCells_;
        
        //- Interface curvature
        volScalarField kappaI_;
        volScalarField kappaTest_;
        //- Alpha cutoff for small cell treatment
        scalar alphaMin_;
        
        //- Alpha interface reconstruction cutoff
        scalar reconstructTol_;

    // Private member functions        
        //- Calculate the interface normals
        void calculateInterfaceNormal();
        
        //- Get the outward facing normal vector from faceI relative to cellI
        vector outwardNormal(label faceI, label cellI) const;
        
public:

    // Constructors

        //- Construct on host mesh
        PLICInterface
        (
            volScalarField& alphaLiquid,
            volScalarField& alphaVapor,
            const volVectorField& U
        );

        //- Return clone
        autoPtr<PLICInterface> clone() const;

    // Public Member Functions

        scalar alphaMin() const
        {
            return alphaMin_;
        }
        
        const surfaceScalarField& alphaLiquidf() const
        {
            return alphaf_;
        }
        
        tmp<surfaceScalarField> alphaVaporf() const
        {
            return 1.0 - alphaf_;
        }

        const volScalarField& alphaVapor() const
        {
            return alphaVapor_;
        }
        
        const volScalarField& alphaLiquid() const
        {
            return alphaLiquid_;
        }
        
        const surfaceScalarField& wL() const
        {
            return wL_;
        }
        
        const surfaceScalarField& wV() const
        {
            return wV_;
        }

        const volScalarField& area() const
        {
            return iArea_;
        }
        
        const volVectorField& normal() const
        {
            return iNormal_;
        }
        
        const surfaceScalarField& phiAlpha() const
        {
            return phiAlpha_;
        }
        
        const volScalarField& kappa() const
        {
            return kappaI_;
        }
        
        tmp<volScalarField> AbydeltaLV() const;
      
       // tmp<volScalarField> AbydeltaV() const;
        
        tmp<volVectorField> shearVec
        (
            word region, 
            const volVectorField& uL, 
            const volVectorField& uV,
            const volScalarField& muL,
            const volScalarField& muV
        ) const;
        
        scalar calcCurvature
        (
            const vector& n1,
            const point& p1,
            const vector& n2,
            const point& p2
        ) const;
        
        //- Get the "used" alpha (zero in small cells)
        tmp<volScalarField> alphaVaporCorr() const;

        //- Evolve the interface given a solid consumption rate
        void advect
        (
            const volScalarField& liquidVolSource
        );
        
        //- Update the stored curvature field
        void updateKappa(bool changeNormals);
        void updateKappaTest(bool changeNormals);
        tmp<surfaceVectorField> liquidUf() const;
        
        //- Calculate the interface parameters given alpha
        void correct();
        
        //- Update the interface after a mesh update
        void update();
        
        void calcTransferWeights();
        
        //- Get the mesh refinement criteria field
        tmp<volScalarField> getRefinementField() const;
        
        //- Calculate the small cell transfer weights
        tmp<surfaceScalarField> scTransferWeights();
        
        
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh> >
        scAverage
        (
            word region,
            const GeometricField<Type, fvPatchField, volMesh>& src
        ) const;

        //- Transfer a field from small cells to neighbour cells
        template<class Type>
        void transfer
        (
            const surfaceScalarField& w,
            const GeometricField<Type, fvPatchField, volMesh>& src,
            GeometricField<Type, fvPatchField, volMesh>& dest,
            Type zero
        );
        
        void solveCholesky
        (
            scalarList& A,
            scalarList& b,
            scalarList& x
        ) const;
        
        void addPoint
        (
            const vector& r,
            const Tensor<scalar>& T,
            scalarList& A,
            scalarList& b,
            scalar w,
            label dims
        ) const;
        
        void calcInterfacePlanes();
        
        void calcAlphaf();
        //- Set source terms in small cells to set values in governing eqns
        /*template<class Type>
        void setScValue
        (
            const surfaceScalarField& w,
            GeometricField<Type, fvPatchField, volMesh>& Su,
            GeometricField<scalar, fvPatchField, volMesh>& Sp,
            const GeometricField<Type, fvPatchField, volMesh>& src,
            const dimensioned<Type>& solidValue,
            const dimensionedScalar& diag,
            const word& mode
        );*/
        
        //-
        tmp<surfaceScalarField> phiAlphaFrac() const;
        
        
        //- Return a mask for cells with no gas
        tmp<volScalarField> noGasCells() const;
        
        //- Return a mask for cells with no liquid
        tmp<volScalarField> noLiquidCells() const;
        
        //- Return a mask for small cells
        tmp<volScalarField> smallGasCells() const;
        
        //- Return a mask for small and liquid cells
        tmp<volScalarField> smallLiquidCells() const;
        
        //- Return a mask for liquid cells
        tmp<volScalarField> liquidCells() const;
        
        //- Return a mask for gas cells
        tmp<volScalarField> gasCells() const;
        
        //- Return a mask for no cells
        tmp<volScalarField> noCells() const;
        
        //- IB-aware interpolation
        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > ibInterpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;
        
       
        template<class Type>
        tmp
        <
            GeometricField
            <
                typename outerProduct<vector,Type>::type, fvPatchField, volMesh
            >
        >
        reconstruct
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
        ) const;

        template<class Type>
        tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
        > reconstruct
        (
            const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >&
        ) const;
        
        // Overloading

        inline tmp<surfaceScalarField> interpolate
        (
            const volScalarField& vf
        ) const
        {
            return ibInterpolate<scalar>(vf);
        }

        inline tmp<surfaceScalarField> interpolate
        (
            tmp<volScalarField> tvf
        ) const
        {
            return ibInterpolate<scalar>(tvf());
        }

        inline tmp<surfaceVectorField> interpolate
        (
            const volVectorField& vf
        ) const
        {
            return ibInterpolate<vector>(vf);
        }

        inline tmp<surfaceVectorField> interpolate
        (
            tmp<volVectorField> tvf
        ) const
        {
            return ibInterpolate<vector>(tvf());
        }
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PLICInterfaceTemplates.H"

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
