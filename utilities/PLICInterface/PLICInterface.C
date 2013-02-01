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

\*---------------------------------------------------------------------------*/

#include "PLICInterface.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PLICInterface::PLICInterface
(
    volScalarField& alpha
)
:
    mesh_(alpha.mesh()),
    plicDict_
    (
        IOobject
        (
            "PLICInterfaceDict",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    alpha_(alpha),

    alphaf_
    (
        IOobject
        (
            "alphaf",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphaf", dimless, 0.0)
    ),
    
    phiAlpha_
    (
        IOobject
        (
            "phiAlpha",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("phiAlpha", dimVolume/dimTime, 0.0)
    ),
    
    sumalphaf_
    (
        IOobject
        (
            "sumalphaf",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("sumalphaf",dimless,0.0)
    ),

    iArea_
    (
        IOobject
        (
            "iArea",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("iArea",dimArea, 0.0)
    ),

    iNormal_
    (
        IOobject
        (
            "iNormal",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("iNormal", dimless, vector::zero)
    ),
    
    iPoint_
    (
        IOobject
        (
            "iPoint",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("iPoint", dimless, vector::zero)
    ),
    
    gasC_
    (
        IOobject
        (
            "gasC",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("gasC", dimless, vector::zero)
    ),
    
    liquidC_
    (
        IOobject
        (
            "liquidC",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("liquidC", dimless, vector::zero)
    ),

    alphaMin_(plicDict_.lookup("alphaMin")),
    reconstructTol_(plicDict_.lookup("reconstructTol"))

{
    Foam::Info << "Created immersed boundary" << Foam::endl;

    // Calculate new interface position and fields
    correct();
    
    alpha_.oldTime();
    alphaf_.oldTime();
    iNormal_.oldTime();
    iPoint_.oldTime();
    iArea_.oldTime();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::PLICInterface> Foam::PLICInterface::clone() const
{
    notImplemented("PLICInterface::clone() const");
    return autoPtr<PLICInterface>(NULL);
}

// Perform a correct update of the interface after the mesh is adapted
void Foam::PLICInterface::update()
{
    //Mesh update will simply apply alpha of the parent cell to its child cells
    // which is incorrect for the sharp interface. Fortunately, it will also
    // write iNormal and iPoint to each child cell so that alpha can be 
    // recalculated easily.
    
    // When coarsening, this gets more complicated. To alleviate this, we will
    // require the boundary to always be refined to the maximum level. It can
    // only unrefine when the boundary has moved away.

    forAll(alpha_, cellI)
    {
        if (mag(iNormal_[cellI]) > SMALL && mag(iPoint_[cellI]) > SMALL)
        {
            cuttableCell cc(mesh_, cellI);
            plane p(iPoint_[cellI],iNormal_[cellI]);
            alpha_[cellI] = 1.0 - cc.cut(p);
        }
    }
    alpha_.correctBoundaryConditions();
    
    correct();
}


void Foam::PLICInterface::calculateInterfaceNormal
(
    const volScalarField& intermeds
)
{    
    //Use a smoothing procedure to capture the interface better
        
    // Get gradient
    iNormal_ = fvc::grad(alpha_)*dimensionedScalar("one",dimLength,1.0);
        
    // Do spatial smoothing operation
    // TODO: Make the weights inputtable in dictionary
    surfaceVectorField iNormalf = interpolate(iNormal_);
    iNormal_ = 0.7*fvc::average(iNormalf) + 0.3*iNormal_;
    
    for (label i = 0; i < 1; ++i)
    {
        iNormalf = interpolate(iNormal_);
        iNormal_ = 0.7*fvc::average(iNormalf) + 0.3*iNormal_;
    }
        
    // Normalize and limit iNormal only to intermediate cells
    iNormal_ *= intermeds / (mag(iNormal_) + SMALL);
    iNormal_.correctBoundaryConditions();
}


// Get the outward facing normal vector on faceI relative to cellI
Foam::vector Foam::PLICInterface::outwardNormal
(
    label faceI,
    label cellI
) const
{
    const face& f = mesh_.faces()[faceI];
    const pointField& points = mesh_.points();

    // Calculate face normal
    vector norm = f.normal(points);
    norm /= mag(norm);
    
    // Flip norm if pointed wrong way
    vector vSF = f.centre(points) - mesh_.C()[cellI];
    
    if ((vSF & norm) < 0.0)
    {
        norm *= -1.0;
    }
    
    return norm;
}


// Given alpha, calculate derived interface fields
void Foam::PLICInterface::correct()
{
    // Step 1: Identify intermediate cells based on alpha_ and reconstructTol_
    volScalarField intermeds = pos(alpha_ - reconstructTol_)
                               *pos(1.0 - reconstructTol_ - alpha_);

    // Step 2: Calculate interface normal in intermediate cells
    calculateInterfaceNormal(intermeds);

    // Step 3: Calculate the cut plane and cut area in intermediate cells
    //         setting a_burn to zero in all non-intermediate cells
    forAll(iNormal_, cellI)
    {
        if (intermeds[cellI] > SMALL)
        {
            cuttableCell cc(mesh_, cellI);
            plane p = cc.constructInterface(iNormal_[cellI],1.0-alpha_[cellI]);
            iPoint_[cellI] = p.refPoint();
            iArea_[cellI] = cc.cutArea();
            
            // Save gas and solid portion centroids
            gasC_[cellI] = cc.lostCentroid();
            liquidC_[cellI] = cc.cutCentroid();
            
            // DEBUGGING PURPOSES
            if (iArea_[cellI] < SMALL)
            {
                Info<< "WARNING: Cut area is " << iArea_[cellI]
                    << " with alphaSolid = " << 1.0-alpha_[cellI] << endl;
            }
        }
        else
        {
            iArea_[cellI] = 0.0;
            iPoint_[cellI] = vector::zero;
            if (alpha_[cellI] > 0.5)
            { //gas cell
                gasC_[cellI] = mesh_.C()[cellI];
                liquidC_[cellI] = vector::zero;
            }
            else
            { //liquid cell
                gasC_[cellI] = vector::zero;
                liquidC_[cellI] = mesh_.C()[cellI];
            }
        }
    }
    
    iPoint_.correctBoundaryConditions();
    iNormal_.correctBoundaryConditions();


    // Step 4: Calculate alphaf on all faces, valid only in the homogeneous
    //         regions away from the interface
    alphaf_ = fvc::interpolate(alpha_);

    // Step 5: Use the planes in intermediate cells to correct alphaf 
    //         near the interface. Also catch solid cells that have a partial
    //         burning face as calculated from an intermediate cell cut plane.
    surfaceScalarField hasIntermeds = fvc::interpolate(intermeds);
    
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbor = mesh_.neighbour();

    forAll(alphaf_, faceI)
    {
        label own = owner[faceI];
        label nei = neighbor[faceI];

        // Faces where there is a cut plane in one or both bounding cells
        if (hasIntermeds[faceI] > SMALL)
        {
            cuttableFace cf(mesh_, faceI);

            scalar alphafNei = 1.0;
            scalar alphafOwn = 1.0;

            if (mag(iNormal_[nei]) > SMALL)
            {
                Foam::plane p(iPoint_[nei], iNormal_[nei]);
                alphafNei = cf.cut(p);
            }

            if (mag(iNormal_[own]) > SMALL)
            {
                Foam::plane p(iPoint_[own], iNormal_[own]);
                alphafOwn = cf.cut(p);
            }
            
            alphaf_[faceI] = Foam::min(alphafNei, alphafOwn);
            
            // Catch the cases where:
            //
            //   +---------+
            //   |\        |
            //   | \   g   |
            //   |l \      |
            //   +---======+  <- Face in question
            //   |         |
            //   |    l    |
            //   |         |
            //   +---------+
            //
            //  The value of alphaf (===) for the face is calculated correctly,
            //  but needs to be added to the interface area of the liquid cell
            
            if (alpha_[own] < reconstructTol_.value()
                && alphaf_[faceI] > SMALL)
            {
                //own is solid
                iArea_[own] += mesh_.magSf()[faceI] * alphaf_[faceI];
                iNormal_[own] += outwardNormal(faceI, own)*alphaf_[faceI];
            }
            else if (alpha_[nei] < reconstructTol_.value()
                     && alphaf_[faceI] > SMALL)
            {
                //nei is solid
                iArea_[nei] += mesh_.magSf()[faceI] * alphaf_[faceI];
                iNormal_[nei] += outwardNormal(faceI, nei)*alphaf_[faceI];
            }
        }

        // Catch sharp face-coincident interfaces. In this case, the liquid cell
        // is the one that will be given the interface area
        else if (mag(alpha_[own] - alpha_[nei]) > 0.1)
        {
            alphaf_[faceI] = 0.0;

            //set iNormal and a_burn for this case
            label liquidcell = (alpha_[own] < 0.1) ? own : nei;

            // Increment iArea since a liquid cell could technically have more
            // than one sharp boundary
            iArea_[liquidcell] += mesh_.magSf()[faceI];
                        
            // Increment iNormal_
            iNormal_[liquidcell] += outwardNormal(faceI, liquidcell);
        }
    }
    
    // Now set alphaf on parallel patches
    const volScalarField::GeometricBoundaryField& alphaBf = 
        alpha_.boundaryField();
    const volVectorField::GeometricBoundaryField& iPointBf = 
        iPoint_.boundaryField();
    volVectorField::GeometricBoundaryField& iNormalBf = 
        iNormal_.boundaryField();
    volScalarField::GeometricBoundaryField& iAreaBf = 
        iArea_.boundaryField();
    const surfaceScalarField::GeometricBoundaryField& hasIntermedsBf =
        hasIntermeds.boundaryField();
        
    surfaceScalarField::GeometricBoundaryField& alphafBf =
        alphaf_.boundaryField();
    
    forAll(alphafBf, patchI)
    {
        const fvPatchScalarField& alphaPf = alphaBf[patchI];
        const fvPatchVectorField& iPointPf = iPointBf[patchI];
        
        fvPatchVectorField& iNormalPf = iNormalBf[patchI];
        fvPatchScalarField& iAreaPf = iAreaBf[patchI];
        
        const scalarField& hasIntermedsPf = hasIntermedsBf[patchI];
        scalarField& alphafPf = alphafBf[patchI];

        const labelList& pFaceCells = mesh_.boundary()[patchI].faceCells();
    
        if (alphaPf.coupled()) //returns true for parallel and cyclic patches
        {
            // Get values across parallel patch
            const vectorField iPointPNf(iPointPf.patchNeighbourField());
            const vectorField iNormalPNf(iNormalPf.patchNeighbourField());
            const scalarField alphaPNf(alphaPf.patchNeighbourField());
            const scalarField iAreaPNf(iAreaPf.patchNeighbourField());
            
            //patch face starting IDs in mesh.faces()
            label patchFs = alphaPf.patch().start();
            const fvPatch& meshPf = mesh_.boundary()[patchI];
            
            forAll(alphafPf, pFaceI)
            {
                label pfCellI = pFaceCells[pFaceI];
                
                if (hasIntermedsPf[pFaceI] > SMALL)
                {
                    cuttableFace cf(mesh_, patchFs+pFaceI);
                    
                    scalar alphafNei = 1.0;
                    scalar alphafOwn = 1.0;

                    if (mag(iNormalPNf[pFaceI]) > SMALL)
                    {
                        Foam::plane p(iPointPNf[pFaceI], iNormalPNf[pFaceI]);
                        alphafNei = cf.cut(p);
                    }

                    if (mag(iNormal_[pfCellI]) > SMALL)
                    {
                        Foam::plane p(iPoint_[pfCellI], iNormal_[pfCellI]);
                        alphafOwn = cf.cut(p);
                    }
                    
                    alphafPf[pFaceI] = Foam::min(alphafNei, alphafOwn);
                    
                    // now catch sharp edges
                    if (alpha_[pfCellI] < reconstructTol_.value()
                        && alphafPf[pFaceI] > SMALL)
                    {
                        //own is liquid
                        iArea_[pfCellI] += meshPf.magSf()[pFaceI]
                                           *alphafPf[pFaceI];

                        iNormal_[pfCellI] += 
                                outwardNormal(patchFs+pFaceI, pfCellI)
                                *alphafPf[pFaceI];
                    }
                    
                }
                else if (mag(alpha_[pfCellI] - alphaPNf[pFaceI]) > 0.1)
                {
                    alphafPf[pFaceI] = 0.0;

                    //set iNormal and iArea for this case
                    if (alpha_[pfCellI] < 0.1)
                    {
                        // than one sharp boundary
                        iArea_[pfCellI] += meshPf.magSf()[pFaceI];
                                    
                        // Increment iNormal_
                        iNormal_[pfCellI] += 
                            outwardNormal(patchFs+pFaceI, pfCellI);
                    }
                }
            }
        }
    }
    
    //Re-normalize iNormal (only needed for the cases when it is incremented)
    iNormal_ /= (mag(iNormal_) + VSMALL);
    
    sumalphaf_ = fvc::surfaceSum(alphaf_);
    
    iArea_.correctBoundaryConditions();
    iNormal_.correctBoundaryConditions();
    gasC_.correctBoundaryConditions();
    liquidC_.correctBoundaryConditions();
}



// Calculate alpha that excludes small gas cells
Foam::tmp<Foam::volScalarField> Foam::PLICInterface::alphaCorr() const
{
    return alpha_ * pos(alpha_ - alphaMin_);
}


Foam::tmp<Foam::surfaceScalarField> Foam::PLICInterface::scTransferWeights()
{
    tmp<surfaceScalarField> tw
    (
        new surfaceScalarField
        (
            IOobject
            (
                "tw",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("tw", dimless, 0.0)
        )
    );
    surfaceScalarField& w = tw();
    
    // If negative, alpha is a small cell
    volScalarField alphaShift = alpha_ - alphaMin_;
    
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbor = mesh_.neighbour();

    // Calculate transfer weights
    forAll(w, faceI)
    {
        label own = owner[faceI];
        label nei = neighbor[faceI];

        if (alphaShift[own] * alphaShift[nei] * alphaf_[faceI] < 0.0)
        { //one is small, one is not, and they share a gas boundary

            label sc = (alphaShift[own] < 0.0) ? own : nei;

            w[faceI] = mag
            (
                iNormal_[sc] & mesh_.Sf()[faceI]
            ) * alphaf_[faceI];
            
            alphaf_[faceI] = 0.0;
        }
        
        // if both cells are small, zero out alphaf between then
        if (alphaShift[own] < 0.0 && alphaShift[nei] < 0.0)
        {
            alphaf_[faceI] = 0.0;
        }
    }
    
    
    
    // Now set alphaf on parallel patches
    const volScalarField::GeometricBoundaryField& alphaShiftBf = 
        alphaShift.boundaryField();
    const volVectorField::GeometricBoundaryField& iNormalBf = 
        iNormal_.boundaryField();
        
    surfaceScalarField::GeometricBoundaryField& wBf =
        w.boundaryField();
    surfaceScalarField::GeometricBoundaryField& alphafBf =
        alphaf_.boundaryField();
    
    forAll(alphafBf, patchI)
    {
        const fvPatchScalarField& alphaShiftPf = alphaShiftBf[patchI];
        const fvPatchVectorField& iNormalPf = iNormalBf[patchI];
        
        scalarField& alphafPf = alphafBf[patchI];
        scalarField& wPf = wBf[patchI];
        
        const labelList& pFaceCells = mesh_.boundary()[patchI].faceCells();
    
        if (alphaShiftPf.coupled())
        {
            // Get values across parallel patch
            const scalarField alphaShiftPNf(alphaShiftPf.patchNeighbourField());
            const vectorField iNormalPNf(iNormalPf.patchNeighbourField());
            
            const fvPatch& meshPf = mesh_.boundary()[patchI];
            
            forAll(alphafPf, pFaceI)
            {
                label pfCellI = pFaceCells[pFaceI];
                
                //Calculate weight on small-large cell faces and close
                // small cell faces
                if ( alphaShift[pfCellI] * alphaShiftPNf[pFaceI]
                     * alphafPf[pFaceI] < 0.0)
                {
                
                    vector scNorm = (alphaShift[pfCellI] < 0.0)
                                    ? iNormal_[pfCellI]
                                    : iNormalPNf[pFaceI];

                    wPf[pFaceI] = mag
                    (
                        scNorm & meshPf.Sf()[pFaceI]
                    ) * alphafPf[pFaceI];

                    alphafPf[pFaceI] = 0.0;
                }
                
                //Close faces between small cells
                if (alphaShift[pfCellI] < 0.0 && alphaShiftPNf[pFaceI] < 0.0)
                {
                    alphafPf[pFaceI] = 0.0;
                }
            }
        }
    }
    
    sumalphaf_ = fvc::surfaceSum(alphaf_);
    
    return tw;
}




Foam::tmp<Foam::volScalarField> 
Foam::PLICInterface::getRefinementField
(
    const volVectorField& U
) const
{
    // Force all cells on OR near the interface to refine to the maximum level
    dimensionedScalar C("C",dimLength,1e6);
    tmp<volScalarField> tRefinementField = C*mag(fvc::grad(alpha_));

    //Include curl criteria from Popinet (Gerris), scaled by 0.5, to also
    // refine key fluid flow regions
    tRefinementField().internalField() = max
    (
        tRefinementField().internalField(), 
        Foam::mag(fvc::curl(U)) * Foam::pow(mesh_.V(),1.0/3.0) * 0.5
    );

    return tRefinementField;
}

Foam::tmp<surfaceScalarField> Foam::PLICInterface::phiAlphaFrac
(
    const volVectorField& U
) const
{
    Info<< "Calculating phiAlphaFrac using reconstructed interfaces" << endl;

    //Normal definition applied first, applicable away from interfaces
    tmp<surfaceScalarField> tphiAlphaFrac = fvc::interpolate(alpha_);
    surfaceScalarField& phiAlphaFrac = tphiAlphaFrac();

    //Then adjust at the interface regions
    //  the scheme for this interpolation is set in fvSchemes
    surfaceVectorField Uf = fvc::interpolate(U);
    
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbor = mesh_.neighbour();
    const vectorField& cellCenters = mesh_.C();
    
    scalar dT = mesh_.time().deltaT().value();
    
    //Loop through internal faces (internal to this processor)
    forAll(owner, faceI)
    {
        label own = owner[faceI];
        label nei = neighbor[faceI];

        //figure out which cell is upwind of this face
        vector vON = cellCenters[nei] - cellCenters[own];
        label uwCell;

        if ((vON & Uf[faceI]) > 0.0)
        {
            uwCell = own;
        }
        else
        {
            uwCell = nei;
        }


        if (   iNormal_[uwCell] != vector::zero 
            && mag(Uf[faceI] & mesh_.Sf()[faceI]) > SMALL)
        {
            // Make a cuttable cell by extruding faceI in the direction of the
            // upwind cell's inverse velocity (-Uf*dT)
            cuttableCell pf(mesh_, faceI, -Uf[faceI]*dT);

            // Then cut the extruded face by the interface plane in the upwind
            // cell
            phiAlphaFrac[faceI] = pf.cut( plane(iPoint_[uwCell], iNormal_[uwCell]) );    
        }

        else if (mag(Uf[faceI] & mesh_.Sf()[faceI]) <= SMALL)
        {
            phiAlphaFrac[faceI] = 0.0;
        }

        else if (mag(alpha_[own] - alpha_[nei]) > 0.1)
        {
            phiAlphaFrac[faceI] = alpha_[uwCell];
        }
    }

    //Loop through boundary patches to set phiAlpha at processor patches
    // Assume that the nominal treatment of phiAlpha at boundary patches is ok
    // for now.

    const volScalarField::GeometricBoundaryField& alphaBf = alpha_.boundaryField();
    const volVectorField::GeometricBoundaryField& iPointBf = iPoint_.boundaryField();
    const volVectorField::GeometricBoundaryField& iNormalBf = iNormal_.boundaryField();

    surfaceScalarField::GeometricBoundaryField& phiAlphaFracBf = phiAlphaFrac.boundaryField();
    surfaceVectorField::GeometricBoundaryField& UfBf = Uf.boundaryField();

    forAll(phiAlphaFracBf, patchI)
    {
        const fvPatchScalarField& alphaPf = alphaBf[patchI];
        const fvPatchVectorField& iPointPf = iPointBf[patchI];
        const fvPatchVectorField& iNormalPf = iNormalBf[patchI];

        const vectorField& UfPf = UfBf[patchI];
        scalarField& phiAlphaFracPf = phiAlphaFracBf[patchI];

        const labelList& pFaceCells = mesh_.boundary()[patchI].faceCells();

        if (alphaPf.coupled()) //returns true if a parallel or cyclic patch
        {
            // Get point and normal field values across parallel patch
            const vectorField iPointPNf(iPointPf.patchNeighbourField());
            const vectorField iNormalPNf(iNormalPf.patchNeighbourField());

            forAll(phiAlphaFracPf, pFaceI)
            {
                label pfCellI = pFaceCells[pFaceI];

                //TODO Switch to upwind later

                if ((iNormal_[pfCellI] != vector::zero || iNormalPNf[pFaceI] != vector::zero)
                    && mag(UfPf[pFaceI] & alphaPf.patch().Sf()[pFaceI]) > SMALL)
                {
                    //patch face starting IDs in mesh.faces()
                    label patchFs = alphaPf.patch().start();

                    // Make a cuttable cell by extruding face in direction -Uf*dT
                    cuttableCell pf(mesh_, patchFs+pFaceI, -UfPf[pFaceI]*dT);

                    scalar cutSum = 0.0;
                    label nCuts = 0;

                    if (iNormal_[pfCellI] != vector::zero)
                    {
                        // Cut the polyhedron by the reconstructed interface plane in cell on this proc
                        cutSum += pf.cut( plane(iPoint_[pfCellI], iNormal_[pfCellI]) );
                        nCuts++;
                    }

                    if (iNormalPNf[pFaceI] != vector::zero)
                    {
                        // Cut the polyhedron by the reconstructed interface plane in cell on neighbour proc
                        cutSum += pf.cut( plane(iPointPNf[pFaceI], iNormalPNf[pFaceI]) );
                        nCuts++;
                    }

                    phiAlphaFracPf[pFaceI] = cutSum / nCuts;            
                }
                else if (mag(UfPf[pFaceI] & alphaPf.patch().Sf()[pFaceI]) <= SMALL)
                {
                    phiAlphaFracPf[pFaceI] = 0.0;
                }

                //TODO do upwindy here

            }

        }
        else
        {
            //For other patches, the faces only have one neighbor. We will let
            //the normal boundary conditions on phiAlpha apply and not mess
            //with these.

            /*forAll(phiAlphaPf, pFacei)
            {
                label pfCelli = pFaceCells[pFacei];
                if (iNormal[pfCelli] != vector::zero)
                {
                    // Make an extruded polyhedron from faceI in direction -UdT
                    //Polyhedron pf(mesh, faceI, -Uf[faceI]*runTime.deltaT());

                }
            }*/
        }
    }


    return tphiAlphaFrac;
}


void Foam::PLICInterface::moveInterface
(
    const surfaceScalarField& phi,
    const volVectorField& U   
)
{
    // Get the volume flux of gas phase (alpha) on faces
    phiAlpha_ = phiAlphaFrac(U) * phi;

    MULES::explicitSolve
    (
        geometricOneField(),
        alpha_, 
        phi, 
        phiAlpha_, 
        zeroField(), //Sp
        zeroField(), //Su (include alpha*divU?)
        1,           //alphaMax
        0            //alphaMin
    );
    alpha_.max(0.0);

    // do this outside? just have this return phiAlpha?
    //rhoPhi = phiAlpha*(rho1 - rho2) + phi*rho2;


    Info<< "Gas phase volume fraction = "
        << alpha_.weightedAverage(mesh_.Vsc()).value()
        << "  Min(alpha) = " << min(alpha_).value()
        << "  Max(alpha) = " << max(alpha_).value()
        << endl;

    
    
    //Correct the interface parameters using the new alpha
    //  This updates iNormal, iPoint, iArea, centroids, ...
    correct();
}


        
Foam::tmp<Foam::volScalarField>
Foam::PLICInterface::smallAndLiquidCells() const
{
    return neg(alpha_ - alphaMin_);
}
        
Foam::tmp<Foam::volScalarField> Foam::PLICInterface::smallCells() const
{
    return neg(alpha_ - alphaMin_)*pos(alpha_ - SMALL);
}


Foam::tmp<Foam::volScalarField> Foam::PLICInterface::liquidCells() const
{
    return neg(alpha_ - SMALL);
}

Foam::tmp<Foam::volScalarField> Foam::PLICInterface::gasCells() const
{
    return pos(alpha_ - alphaMin_);
}

Foam::tmp<Foam::volScalarField> Foam::PLICInterface::noCells() const
{
    return neg(alpha_ + 100.0);
}


// ************************************************************************* //
