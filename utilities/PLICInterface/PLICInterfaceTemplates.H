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

#ifndef PLICInterfaceTemplates_H
#define PLICInterfaceTemplates_H

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > PLICInterface::scAverage
(
    word region,
    const GeometricField<Type, fvPatchField, volMesh>& src
) const
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tavg = src*0.0;
    GeometricField<Type, fvPatchField, volMesh>& avg = tavg();
    
    const surfaceScalarField& w = (region=="Vapor") ? wV_ : wL_;
    const volScalarField& alpha = (region=="Vapor") ? alphaVapor_ : alphaLiquid_;
    
    volScalarField wtot = fvc::surfaceSum(w);
    
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbor = mesh_.neighbour();
    
    forAll(w, faceI)
    {
        if (w[faceI] > 0.0)
        {
            label own = owner[faceI];
            label nei = neighbor[faceI];

            label sc = (alpha[own] < alphaMin_) ? own : nei; 
            label nbc = (sc == own) ? nei : own;

            avg[sc] += w[faceI]/wtot[sc] * src[nbc];
        }
    }
    
    
    // Now do averaging across parallel patches
    const typename GeometricField<Type, fvPatchField, volMesh>::GeometricBoundaryField&
        srcBf = src.boundaryField();
                         
    const surfaceScalarField::GeometricBoundaryField& 
        wBf = w.boundaryField();
    
    forAll(srcBf, patchI)
    {
        const fvPatchField<Type>& srcPf = srcBf[patchI];
        const scalarField& wPf = wBf[patchI];
        
        const labelList& pFaceCells = mesh_.boundary()[patchI].faceCells();
    
        if (srcPf.coupled())
        {
            // Get values across parallel patch
            const Field<Type> srcPNf(srcPf.patchNeighbourField());

            forAll(srcPf, pFaceI)
            {
                label pfCellI = pFaceCells[pFaceI];
                
                if ( wPf[pFaceI] > 0.0)
                {
                    if (alpha[pfCellI] < alphaMin_)
                    {//I own the small cell (pfCellI is the small cell)
                        avg[pfCellI] += wPf[pFaceI]/wtot[pfCellI]
                                          * srcPNf[pFaceI];
                    }
                }
            }
        }
    }


    return tavg;
}



template<class Type>
void PLICInterface::transfer
(
    const surfaceScalarField& w,
    const GeometricField<Type, fvPatchField, volMesh>& src,
    GeometricField<Type, fvPatchField, volMesh>& dest,
    Type zero
)
{
    volScalarField wtot = fvc::surfaceSum(w);
    
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbor = mesh_.neighbour();
    
    forAll(w, faceI)
    {
        if (w[faceI] > 0.0)
        {
            label own = owner[faceI];
            label nei = neighbor[faceI];

            label sc = (alphaVapor_[own] < alphaMin_) ? own : nei; 
            label rc = (sc == own) ? nei : own;

            dest[rc] += w[faceI]/wtot[sc] * src[sc];
            dest[sc] = zero;
        }
    }
    
    // Now do transfer at parallel patches
    const typename GeometricField<Type, fvPatchField, volMesh>::GeometricBoundaryField&
        srcBf = src.boundaryField();
        
    const volScalarField::GeometricBoundaryField& 
        wtotBf = wtot.boundaryField();
                  
    const surfaceScalarField::GeometricBoundaryField& 
        wBf = w.boundaryField();
    
    forAll(srcBf, patchI)
    {
        const fvPatchField<Type>& srcPf = srcBf[patchI];
        const fvPatchScalarField& wtotPf = wtotBf[patchI];
        const scalarField& wPf = wBf[patchI];
        
        const labelList& pFaceCells = mesh_.boundary()[patchI].faceCells();
    
        if (srcPf.coupled())
        {
            // Get values across parallel patch
            const Field<Type> srcPNf(srcPf.patchNeighbourField());
            const scalarField wtotPNf(wtotPf.patchNeighbourField());
            
            forAll(srcPf, pFaceI)
            {
                label pfCellI = pFaceCells[pFaceI];
                
                if ( wPf[pFaceI] > 0.0)
                {
                    if (alphaVapor_[pfCellI] < alphaMin_)
                    {//I own the small cell
                        dest[pfCellI] = zero;
                    }
                    else
                    {//I own the receiving cell
                        dest[pfCellI] += wPf[pFaceI]/wtotPNf[pFaceI]
                                          * srcPNf[pFaceI];
                    }
                }
            }
        }
    }
}

/*
template<class Type>
void PLICInterface::setScValue
(
    const surfaceScalarField& w,
    GeometricField<Type, fvPatchField, volMesh>& Su,
    GeometricField<scalar, fvPatchField, volMesh>& Sp,
    const GeometricField<Type, fvPatchField, volMesh>& src,
    const dimensioned<Type>& solidValue,
    const dimensionedScalar& diag,
    const word& mode
)
{
    Su = liquidCells()*solidValue*diag;
    Sp = smallAndLiquidCells()*diag;
    
    if (mode == "fix")
    {
        Su += smallCells()*src*diag;
    }
    else if (mode == "avg")
    {
        volScalarField wtot = fvc::surfaceSum(w);
        
        const labelUList& owner = mesh_.owner();
        const labelUList& neighbor = mesh_.neighbour();
        
        forAll(w, faceI)
        {
            if (w[faceI] > 0.0)
            {
                label own = owner[faceI];
                label nei = neighbor[faceI];

                label sc = (alphaVapor_[own] < alphaMin_) ? own : nei; 
                label rc = (sc == own) ? nei : own;

                Su[sc] += w[faceI]/wtot[sc] * src[rc] * diag.value();
            }
        }
        
        // Now do transfer at parallel patches
        const typename GeometricField
            <Type, fvPatchField, volMesh>::GeometricBoundaryField&
            srcBf = src.boundaryField();
                      
        const surfaceScalarField::GeometricBoundaryField& 
            wBf = w.boundaryField();
        
        forAll(srcBf, patchI)
        {
            const fvPatchField<Type>& srcPf = srcBf[patchI];
            const scalarField& wPf = wBf[patchI];
            
            const labelList& pFaceCells = mesh_.boundary()[patchI].faceCells();
        
            if (srcPf.coupled())
            {
                // Get values across parallel patch
                const Field<Type> srcPNf(srcPf.patchNeighbourField());
                
                forAll(srcPf, pFaceI)
                {
                    label pfCellI = pFaceCells[pFaceI];
                    
                    if ( wPf[pFaceI] > 0.0)
                    {                        
                        if (alphaVapor_[pfCellI] < alphaMin_)
                        {//I own the small cell
                            Su[pfCellI] += wPf[pFaceI]/wtot[pfCellI] 
                                          * srcPNf[pFaceI] * diag.value();
                        }
                    }
                }
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "PLICInterface::setScValue()"
        )   << "Setting mode  '" << mode
            << "' not recognized"
            << exit(FatalError);
    }
}
*/

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
PLICInterface::ibInterpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return fvc::interpolate(vf);
}


template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector,Type>::type, fvPatchField, volMesh
    >
>
PLICInterface::reconstruct
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    tmp<GeometricField<GradType, fvPatchField, volMesh> > treconField
    (
        new GeometricField<GradType, fvPatchField, volMesh>
        (
            IOobject
            (
                "volIntegrate("+ssf.name()+')',
                ssf.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            inv(fvc::surfaceSum(sqr(mesh_.Sf())/mesh_.magSf()))
          & fvc::surfaceSum((mesh_.Sf()/mesh_.magSf())*ssf),
            zeroGradientFvPatchField<GradType>::typeName
        )
    );

    treconField().correctBoundaryConditions();

    return treconField;
}

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >
>
PLICInterface::reconstruct
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >& tssf
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;
    tmp<GeometricField<GradType, fvPatchField, volMesh> > tvf
    (
        reconstruct(tssf())
    );
    tssf.clear();
    return tvf;
}

#endif
