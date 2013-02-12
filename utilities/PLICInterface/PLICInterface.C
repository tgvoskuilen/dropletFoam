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
    volScalarField& alphaLiquid,
    volScalarField& alphaVapor,
    const volVectorField& U
)
:
    mesh_(alphaLiquid.mesh()),
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
    
    alphaVapor_(alphaVapor),
    alphaLiquid_(alphaLiquid),
    U_(U),
    
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
    
    phiAlphaLiquid_
    (
        IOobject
        (
            "phiAlphaLiquid",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("phiAlphaLiquid", dimVolume/dimTime, 0.0)
    ),
    
    phiAlphaVapor_
    (
        IOobject
        (
            "phiAlphaVapor",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("phiAlphaVapor", dimVolume/dimTime, 0.0)
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
    
    deltaV_
    (
        IOobject
        (
            "deltaV",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("deltaV", dimLength, 0.0)
    ),
    
    deltaL_
    (
        IOobject
        (
            "deltaL",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("deltaL", dimLength, 0.0)
    ),
    
    wL_
    (
        IOobject
        (
            "wL",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("wL", dimless, 0.0)
    ),
    
    wV_
    (
        IOobject
        (
            "wV",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("wV", dimless, 0.0)
    ),
    
    intermeds_
    (
        IOobject
        (
            "intermeds",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("intermeds", dimless, 0.0)
    ),
    
    alphaLsmooth_
    (
        IOobject
        (
            "alphaLsmooth",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphaLsmooth", dimless, 0.0)
    ),
    
    smoothN_
    (
        IOobject
        (
            "smoothN",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("smoothN", dimless, vector::zero)
    ),
    
    
    liquidCells_
    (
        IOobject
        (
            "liquidCells",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("liquidCells", dimless, 0.0)
    ),
    
    smallLiquidCells_
    (
        IOobject
        (
            "smallLiquidCells",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("smallLiquidCells", dimless, 0.0)
    ),
    
    noLiquidCells_
    (
        IOobject
        (
            "noLiquidCells",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("noLiquidCells", dimless, 0.0)
    ),
    
    gasCells_
    (
        IOobject
        (
            "gasCells",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("gasCells", dimless, 0.0)
    ),
    
    smallGasCells_
    (
        IOobject
        (
            "smallGasCells",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("smallGasCells", dimless, 0.0)
    ),
    
    noGasCells_
    (
        IOobject
        (
            "noGasCells",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("noGasCells", dimless, 0.0)
    ),

    kappaI_
    (
        IOobject
        (
            "kappaI",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("kappaI", dimless/dimLength, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    kappaTest_
    (
        IOobject
        (
            "kappaTest",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("kappaTest", dimless/dimLength, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    
    alphaMin_(readScalar(plicDict_.lookup("alphaMin"))),
    reconstructTol_(readScalar(plicDict_.lookup("reconstructTol")))

{
    Foam::Info << "Created PLIC interface" << Foam::endl;

    // Calculate new interface position and fields
    correct();
    
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
    
    forAll(alphaLiquid_, cellI)
    {
        if (mag(iNormal_[cellI]) > SMALL && mag(iPoint_[cellI]) > SMALL)
        {
            cuttableCell cc(mesh_, cellI);
            plane p(iPoint_[cellI],iNormal_[cellI]);
            alphaLiquid_[cellI] = cc.cut(p);
        }
    }
    alphaLiquid_.correctBoundaryConditions();
    
    correct();
}

tmp<surfaceScalarField> Foam::PLICInterface::stf() const
{
    tmp<surfaceScalarField> tstf
    (
        new surfaceScalarField
        (
            IOobject
            (
                "stf",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("stf",dimPressure/dimLength,0.0)
        )
    );
    
    dimensionedScalar deltaN("deltaN", 1e-8/pow(average(mesh_.V()), 1.0/3.0));

    // Cell gradient of alpha
    const volVectorField gradAlpha(fvc::grad(alphaLiquid_));

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    // Face unit interface normal
    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN));

    // Simple expression for curvature
    volScalarField kappaI = -fvc::div(nHatfv & mesh_.Sf());

    tstf() = dimensionedScalar("sigma", dimPressure*dimLength, 0.07)
     * fvc::interpolate(kappaI) * fvc::snGrad(alphaLiquid_);

    return tstf;
}

void Foam::PLICInterface::calculateInterfaceNormal()
{    
    alphaLsmooth_ = alphaLiquid_;
    dimensionedScalar dx = pow(min(mesh_.V()), 1.0/3.0);
    dimensionedScalar dA = dx*dx/4.0; //dTau("dTau",dimArea,dx*dx/4.0);
    
    for(label i = 0; i < 2; i++)
    {
        alphaLsmooth_ += dA * fvc::laplacian(alphaLsmooth_);
    }
    
    dimensionedScalar s("s",dimless/dimLength,SMALL);
    iNormal_ = -fvc::grad(alphaLsmooth_) / (mag(fvc::grad(alphaLsmooth_)) + s) * intermeds_;

    iNormal_.correctBoundaryConditions();
    

    //Use a smoothing procedure to capture the interface better
        
    // Get gradient
    /*iNormal_ = -fvc::grad(alphaLiquid_)*dimensionedScalar("one",dimLength,1.0);
        
    // Do spatial smoothing operation
    // TODO: Make the weights inputtable in dictionary
    surfaceVectorField iNormalf = interpolate(iNormal_);
    iNormal_ = 0.7*fvc::average(iNormalf) + 0.3*iNormal_;
    
    for (label i = 0; i < 3; ++i)
    {
        iNormalf = interpolate(iNormal_);
        iNormal_ = 0.7*fvc::average(iNormalf) + 0.3*iNormal_;
    }
        
    // Normalize and limit iNormal only to intermediate cells
    iNormal_ *= intermeds_ / (mag(iNormal_) + SMALL);
    iNormal_.correctBoundaryConditions();
    */
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

Foam::tmp<Foam::volScalarField> Foam::PLICInterface::AbydeltaLV() const
{
    dimensionedScalar rV("rV",dimless/dimVolume,1.0);
    dimensionedScalar s("s",dimLength,SMALL);
    
    tmp<volScalarField> val = iArea_ / (deltaL_+s)*rV;
    val().internalField() /= mesh_.V();

    return val;
}

Foam::tmp<Foam::volVectorField> Foam::PLICInterface::shearVec
(
    word region, 
    const volVectorField& uL, 
    const volVectorField& uV,
    const volScalarField& muL,
    const volScalarField& muV
) const
{
    // shear vec is (ut - uit)*Ai/delta
    
    tmp<volVectorField> utL = uL - (uL & iNormal_) * iNormal_;
    tmp<volVectorField> utV = uV - (uV & iNormal_) * iNormal_;
    
    dimensionedScalar s("s",dimLength,SMALL);
    tmp<volScalarField> wV = muV/(deltaV_+s);
    tmp<volScalarField> wL = muL/(deltaL_+s);
    
    dimensionedScalar s2("s2",dimMass/dimArea/dimTime,SMALL);
    tmp<volVectorField> uti = (wV()*utV() + wL()*utL())/(wV()+wL()+s2);
    
    dimensionedScalar rV("rV",dimless/dimVolume,0.0); //TEMPORARY
    if( region == "Vapor" )
    {
        tmp<volVectorField> tau = -muV*iArea_/(deltaV_+s)*(utV - uti)*rV;
        tau().internalField() /= mesh_.V();
        return tau;
    }
    else
    {
        tmp<volVectorField> tau = -muL*iArea_/(deltaL_+s)*(utL - uti)*rV;
        tau().internalField() /= mesh_.V();
        return tau;
    }
}



void Foam::PLICInterface::updateKappaTest(bool changeNormals)
{
/*
    alphaLsmooth_ = alphaLiquid_;
    scalar dx = 6e-5;
    dimensionedScalar dTau("dTau",dimArea,dx*dx/4.0);
    
    for(label i = 0; i < 3; i++)
    {
        alphaLsmooth_ += dTau * fvc::laplacian(alphaLsmooth_);
    }
    dimensionedScalar s("s",dimless/dimLength,SMALL);
    
    smoothN_ = fvc::grad(alphaLsmooth_) / (mag(fvc::grad(alphaLsmooth_)) + s);
    
    kappaTest_ = -fvc::div(fvc::interpolate(smoothN_) & mesh_.Sf()) * pos(mag(iNormal_) - SMALL);

    smoothN_ *= pos(mag(iNormal_) - SMALL);
    */

/*
    tmp<surfaceVectorField> nf
    (
        new surfaceVectorField
        (
            IOobject
            (
                "nf",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("nf", dimless, vector::zero)
        )
    );
    

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbor = mesh_.neighbour();
    
    forAll(nf(), faceI)
    {
        label own = owner[faceI];
        label nei = neighbor[faceI];
        

        if( mag(iNormal_[own]) > SMALL && mag(iNormal_[nei]) > SMALL )
        {
            point po = liquidC_[own] + deltaL_[own]*iNormal_[own];
            point pn = liquidC_[nei] + deltaL_[nei]*iNormal_[nei];
            point pf = mesh_.Cf()[faceI];
            vector Af = mesh_.Sf()[faceI];
            
            scalar l1 = ((pf-po) & Af) / mag(Af);
            scalar ltot = ((pn-po) & Af) / mag(Af);
            
            scalar f1 = l1 / ltot;
            
            nf()[faceI] = f1*iNormal_[nei] + (1.0-f1)*iNormal_[own];
        }
    }

    //TODO: Do Parallel patches Here

    kappaTest_ = fvc::div(nf() & mesh_.Sf());
*/

/*
    tmp<surfaceScalarField> kap
    (
        new surfaceScalarField
        (
            IOobject
            (
                "kap",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("kap", dimless/dimLength, 0.0)
        )
    );
    
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbor = mesh_.neighbour();
    
    forAll(kap(), faceI)
    {
        label own = owner[faceI];
        label nei = neighbor[faceI];
        
        if( mag(iNormal_[own]) > SMALL && mag(iNormal_[nei]) > SMALL )
        {
            kap()[faceI] = calcCurvature
            (
                iNormal_[own],
                liquidC_[own] + deltaL_[own]*iNormal_[own],
                iNormal_[nei],
                liquidC_[nei] + deltaL_[nei]*iNormal_[nei]
            );
        }
    }

    //TODO: Do Parallel patches Here


    kappaTest_ = fvc::surfaceSum(kap() * mesh_.magSf()) / fvc::average(mesh_.magSf());
*/

/*
    surfaceScalarField nfDotAf(fvc::interpolate(iNormal_) & mesh_.Sf());
    
    kappaTest_ = fvc::div(nfDotAf) * pos(mag(iNormal_)-SMALL);
*/


    //3D or 2D selection required
    label dims = 2;

    label nb = (dims==3) ? 6 : 3;
    label nA = (nb*(nb+1))/2;

    scalarList A(nA,0.0);
    scalarList b(nb,0.0);
    scalarList x(nb,0.0);
    
    CPCCellToCellStencil wideStencil(mesh_);

    forAll(mesh_.C(), lcellI)
    {
        if( mag(iNormal_[lcellI]) > SMALL )
        {
            const labelList& cell27Stencil = wideStencil[lcellI];
            
            label npts = 0;
            
            scalar dx = Foam::pow(mesh_.V()[lcellI],1.0/3.0);
            
            vector p0 = liquidC_[lcellI] - 2.0*dx*iNormal_[lcellI];
            
            A = 0.0;
            b = 0.0;
            
            vector xdir = vector::zero;
            vector ydir = vector::zero;
            const vector& zdir = iNormal_[lcellI];
            
            //get coordinate system
            forAll(cell27Stencil, cellI)
            {
                if( cell27Stencil[cellI] < mesh_.nCells() && cellI > 0 )
                {
                    label cellJ = cell27Stencil[cellI];
                    if( mag(iNormal_[cellJ]) > SMALL )
                    {   
                        vector r = liquidC_[cellJ]+deltaL_[cellJ]*iNormal_[cellJ]-p0;
                        if( mag( (r/mag(r)) ^ zdir ) > SMALL )
                        {
                            ydir = zdir ^ (r/mag(r));
                            ydir /= mag(ydir);
                            xdir = ydir ^ zdir;
                            break;
                        }
                    }
                }
            }
            
            if( mag(xdir) < SMALL )
            {
                Info<< "ERROR: no directions set" << endl;
            }
            
            //make coordinate system transformation tensor
            Tensor<scalar> T(xdir, ydir, zdir);
            Tensor<scalar> Ti = inv(T);
            
            forAll(cell27Stencil, cellI)
            {
                if( cell27Stencil[cellI] < mesh_.nCells()  )
                {
                    label cellJ = cell27Stencil[cellI];
                    if( mag(iNormal_[cellJ]) > SMALL )
                    {
                        addPoint
                        (
                            liquidC_[cellJ]+deltaL_[cellJ]*iNormal_[cellJ]-p0,
                            T,
                            A,
                            b,
                            1.0, //alphaLiquid_[cellJ]*(1.0-alphaLiquid_[cellJ]),
                            dims
                        );
                        ++npts;
                    }
                }
            }
            
            //solve Ax = b and get kappaTest_[lcellI] if npts is high enough
            // right now just warn/crash if it's not high enough
            if( npts >= nb )
            {
                if( dims == 3 )
                {
                    solveCholesky(A,b,x);
                    kappaTest_[lcellI] = -(x[0]+x[1]); //wrong!
                }
                else
                {
                    solveCholesky(A,b,x);
                    kappaTest_[lcellI] = -2.0*x[0]/Foam::pow(1+x[1]*x[1],1.5);

                    
                    /*if( kappaTest_[lcellI] < 500 || kappaTest_[lcellI] > 1500)
                    {
                        Info<< "norm error at "<<lcellI<< " = " << x[1] << " with k = " << kappaTest_[lcellI] << endl;
                        forAll(cell27Stencil, cellI)
                        {
                            if( cell27Stencil[cellI] < mesh_.nCells()  )
                            {
                                label cellJ = cell27Stencil[cellI];
                                if( mag(iNormal_[cellJ]) > SMALL )
                                {
                                    Info<< "   p = " << liquidC_[cellJ]+deltaL_[cellJ]*iNormal_[cellJ]
                                        << ", n = " << iNormal_[cellJ] 
                                        << ", c = " << mesh_.C()[cellJ]
                                        << "alpha = " << alphaLiquid_[cellJ]
                                        << endl;
                                }
                            }
                        } 
                    }*/
                    
                    if( mag(x[1]) > 0.01 && changeNormals )
                    {
                        vector nn(x[1],0.0,1.0);
                        nn /= mag(nn);
                        vector newn = Ti & nn;
                        iNormal_[lcellI] = newn;
                    }
                }
            }
            else
            {
                Info<< "WARNING: not enough points (found " << npts << ")" << endl;
                kappaTest_[lcellI] = 985.0;
            }
        }
        else
        { //cell is not an interface cell
            kappaTest_[lcellI] = 0.0;
        }
    }
    
    //TODO: Get cells through parallel patches
    Info<< "max kappaTest before smoothing = " << Foam::max(kappaTest_) << endl;
    
    // kappa averaging iterations
    tmp<volScalarField> tkappaIter
    (
        new volScalarField
        (
            IOobject
            (
                "kappaIter",
                mesh_.time().timeName(),
                mesh_
            ),
            kappaTest_
        )
    );
        
    for(label k = 0; k < 5; ++k)
    {
        forAll(mesh_.C(), lcellI)
        {
            if( mag(kappaTest_[lcellI]) > SMALL )
            {
                const labelList& cell27Stencil = wideStencil[lcellI];
                    
                scalar avgKappa = 0.0;
                label n = 0;
                
                //get coordinate system
                forAll(cell27Stencil, cellI)
                {
                    if( cell27Stencil[cellI] < mesh_.nCells() )
                    {
                        label cellJ = cell27Stencil[cellI];
                        if( mag(kappaTest_[cellJ]) > SMALL )
                        {
                            avgKappa += kappaTest_[cellJ];
                            ++n;
                        }
                    }
                }
                
                tkappaIter()[lcellI] = avgKappa/n;
            }
        }
        kappaTest_ = tkappaIter();
        
        Info<< "max kappaTest after iter "<< k 
            <<" = " << Foam::max(kappaTest_) << endl;
    }
    
    /*
    //smooth kappa along interface
    tmp<surfaceScalarField> kapf
    (
        new surfaceScalarField
        (
            IOobject
            (
                "kapf",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("kapf", dimless/dimLength, 0.0)
        )
    );
    tmp<volScalarField> nfaces
    (
        new volScalarField
        (
            IOobject
            (
                "nfaces",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("nfaces", dimless, 0.0)
        )
    );
    

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbor = mesh_.neighbour();
    
    
    for(label i = 0; i < 5; i++)
    {
        kapf() = dimensionedScalar("kapf", dimless/dimLength, 0.0);
        nfaces() = 0;
        
        forAll(kapf(), faceI)
        {
            label own = owner[faceI];
            label nei = neighbor[faceI];
            
            if( mag(kappaTest_[own]) > SMALL && mag(kappaTest_[nei]) > SMALL )
            {
                kapf()[faceI] = 0.5*(kappaTest_[own]+kappaTest_[nei]);
                nfaces()[nei] += 1.0;
                nfaces()[own] += 1.0;
            }
        }
        
        kappaTest_ = fvc::surfaceSum(kapf()) / (nfaces()+SMALL);
        
        Info<< "max kappaTest after iter "<<i<<" = " << Foam::max(kappaTest_) << endl;
    }
    */
}


void Foam::PLICInterface::solveCholesky
(
    scalarList& A,
    scalarList& b,
    scalarList& x
) const
{
    label n = x.size();
    label nA = A.size();
    scalarList L(nA,0.0);
    scalarList y(n,0.0);
    
    //Create L matrix
    label p = 0;
    for(label c = 0; c < n-1; ++c) //iterate matrix columns
    {
        //set diagonal value of L(k,k) = sqrt(A(k,k))
        label pd = p;
        L[p] = sqrt(A[p]); ++p;

        //Fill remainder of column
        for(label r = c+1; r<n; ++r)
        {
            L[p] = A[p]/L[pd]; ++p;
        }

        //Adjust A values accordingly, this destroys A
        for(label cc = c+1; cc<n; ++cc)
        {
            label dicc = cc*n-((cc-1)*cc)/2;

            for(label r = dicc; r < dicc + n - cc; ++r)
            {
                A[r] -= L[pd+cc-c] * L[pd+cc-c+r-dicc];
            }
        }
    }
    L[p] = Foam::sqrt(A[p]);

    //Forward substitute to get y from Ly = b
    y[0] = b[0] / L[0];
    for(label r=1; r<n; ++r)
    {
        label di = r*n - (r*(r-1))/2;
        
        y[r] = b[r] / L[di];
        
        for(label c=0; c<r; ++c)
        {
            label dic = c*n - (c*(c-1))/2;
            y[r] -= L[dic+r-c] * y[c] / L[di];
        }
    }

    //Backward substitute to get x from L^T x = y
    x[n-1] = y[n-1] / L[nA-1];
    
    for(label c=n-2; c>=0; --c)
    {
        label di = c*n - (c*(c-1))/2;
        
        x[c] = y[c] / L[di];
        
        for(label r=n-1; r>c; --r)
        {
            x[c] -= L[di+r-c] * x[r] / L[di];
        }
    }
}

void Foam::PLICInterface::addPoint
(
    const vector& r,
    const Tensor<scalar>& T,
    scalarList& A,
    scalarList& b,
    scalar w,
    label dims
) const
{
   
    //coordinate transformation
    vector d = (T & r);
    scalar wSqr = w*w;
    //Info<< " adding point " << d << " w = " << w << endl;
    
    if( dims == 3 )
    {
        //fit z = ax^2 + by^2 + cxy + dx + ey + f
        
        //TODO
        /*
        A[0] +=
        A[1] +=
        A[2] +=
        A[3] +=
        A[4] +=
        A[5] +=
        A[6] +=
        A[7] +=
        A[8] +=
        A[9] +=
        A[10] +=
        A[11] +=
        A[12] +=
        A[13] +=
        A[14] +=
        A[15] +=
        A[16] +=
        A[17] +=
        A[18] +=
        A[19] +=
        A[20] +=
        
        b[0] +=
        b[1] +=
        b[2] +=
        b[3] +=
        b[4] +=
        b[5] +=
        */
        
    }
    else
    {
        //fit z = ax^2 + bx + c
        
        A[0] += d.x()*d.x()*d.x()*d.x()*wSqr;
        A[1] += d.x()*d.x()*d.x()*wSqr;
        A[2] += d.x()*d.x()*wSqr;
        A[3] += d.x()*d.x()*wSqr;
        A[4] += d.x()*wSqr;
        A[5] += wSqr;
        
        b[0] += d.z()*d.x()*d.x()*wSqr;
        b[1] += d.z()*d.x()*wSqr;
        b[2] += d.z()*wSqr;
    }
}


scalar Foam::PLICInterface::calcCurvature
(
    const vector& n1,
    const point& p1,
    const vector& n2,
    const point& p2
) const
{
    scalar theta = acos(n1 & n2);
    scalar s = mag(p1 - p2);
    return theta / s;
}


void Foam::PLICInterface::updateKappa(bool changeNormals)
{
    /*dimensionedScalar deltaN("deltaN", 1e-8/pow(average(mesh_.V()), 1.0/3.0));
    surfaceVectorField gradAlphaf(fvc::interpolate(-fvc::grad(alphaLiquid_)));
    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN));
    kappaI_ = -fvc::div(nHatfv & mesh_.Sf());*/
    
    //updateKappaTest(changeNormals);
    
    
    
    //kappaI_ = dimensionedScalar("k",dimless/dimLength,1.0/0.00101) * pos(mag(iNormal_)-SMALL);
    //kappaI_ = kappaTest_;
}


void Foam::PLICInterface::calcAlphaf()
{
    // Step 4: Calculate alphaf on all faces, valid only in the homogeneous
    //         regions away from the interface
    alphaf_ = fvc::interpolate(alphaLiquid_);

    // Step 5: Use the planes in intermediate cells to correct alphaf 
    //         near the interface. Also catch solid cells that have a partial
    //         open face as calculated from an intermediate cell cut plane.
    surfaceScalarField hasIntermeds = fvc::interpolate(intermeds_);
    
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

            // If we just use iNormal as a criteria, it can break, since
            //  iNormal can be incremented for cells with coincident boundaries
            //  while iPoint is left at its default (0,0,0)
            if (mag(iNormal_[nei]) > SMALL && mag(iPoint_[nei]) > SMALL)
            {
                Foam::plane p(iPoint_[nei], iNormal_[nei]);
                alphafNei = cf.cut(p);
            }

            if (mag(iNormal_[own]) > SMALL && mag(iPoint_[own]) > SMALL)
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
            //  but needs to be added to the interface area of the liquid cell.
            //  Because this can happen to more than one face of a cell, we
            //  must increment iNormal each time we find one.
            
            if (alphaVapor_[own] < reconstructTol_ 
                && (1.0 - alphaf_[faceI]) > SMALL)
            {                        
                //own is liquid
                iArea_[own] += mesh_.magSf()[faceI] * (1.0-alphaf_[faceI]);
                iNormal_[own] += outwardNormal(faceI,own)*(1.0-alphaf_[faceI]);
            }
            else if (alphaVapor_[nei] < reconstructTol_
                     && (1.0 - alphaf_[faceI]) > SMALL)
            {
                //nei is liquid                
                iArea_[nei] += mesh_.magSf()[faceI] * (1.0-alphaf_[faceI]);
                iNormal_[nei] += outwardNormal(faceI,nei)*(1.0-alphaf_[faceI]);
            }
        }

        // Catch sharp face-coincident interfaces. In this case, the liquid cell
        // is the one that will be given the interface area
        else if (mag(alphaLiquid_[own] - alphaLiquid_[nei]) > 0.1)
        {

            alphaf_[faceI] = 0.0;
            
            //set iNormal and area
            label liquidcell = (alphaLiquid_[own] > 0.5) ? own : nei;
            
            // Increment iArea since a liquid cell could technically have more
            // than one sharp boundary
            iArea_[liquidcell] += mesh_.magSf()[faceI];
                        
            // Increment iNormal_
            iNormal_[liquidcell] += outwardNormal(faceI, liquidcell);
            
            //TODO: Set iPoint_ to face centroid?
        }
    }
    
    // Now set alphaf on parallel patches
    const volScalarField::GeometricBoundaryField& alphaLBf = 
        alphaLiquid_.boundaryField();
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
        const fvPatchScalarField& alphaLPf = alphaLBf[patchI];
        const fvPatchVectorField& iPointPf = iPointBf[patchI];
        
        fvPatchVectorField& iNormalPf = iNormalBf[patchI];
        fvPatchScalarField& iAreaPf = iAreaBf[patchI];
        
        const scalarField& hasIntermedsPf = hasIntermedsBf[patchI];
        scalarField& alphafPf = alphafBf[patchI];

        const labelList& pFaceCells = mesh_.boundary()[patchI].faceCells();
    
        if (alphaLPf.coupled()) //returns true for parallel and cyclic patches
        {
            // Get values across parallel patch
            const vectorField iPointPNf(iPointPf.patchNeighbourField());
            const vectorField iNormalPNf(iNormalPf.patchNeighbourField());
            const scalarField alphaLPNf(alphaLPf.patchNeighbourField());
            const scalarField iAreaPNf(iAreaPf.patchNeighbourField());
            
            //patch face starting IDs in mesh.faces()
            label patchFs = alphaLPf.patch().start();
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
                    if (alphaVapor_[pfCellI] < reconstructTol_
                        && (1.0-alphafPf[pFaceI]) > SMALL)
                    {
                        //own is liquid

                        iArea_[pfCellI] += meshPf.magSf()[pFaceI]
                                           *(1.0-alphafPf[pFaceI]);

                        iNormal_[pfCellI] += 
                                outwardNormal(patchFs+pFaceI, pfCellI)
                                *(1.0-alphafPf[pFaceI]);
                    }
                    
                }
                else if (mag(alphaLiquid_[pfCellI] - alphaLPNf[pFaceI]) > 0.1)
                {
                    alphafPf[pFaceI] = 0.0;

                    //set iNormal and iArea for this case
                    if (alphaLiquid_[pfCellI] > 0.1)
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
    
    sumalphaf_ = fvc::surfaceSum(alphaf_);
    iNormal_ /= (mag(iNormal_)+SMALL);
    iNormal_.correctBoundaryConditions();
    iArea_.correctBoundaryConditions();
}


void Foam::PLICInterface::calcInterfacePlanes()
{
    // Calculate the cut plane and cut area in any cell that has a
    // normal vector defined.
    forAll(iNormal_, cellI)
    {
        if (mag(iNormal_[cellI]) > SMALL)
        {
            cuttableCell cc(mesh_, cellI);
            plane p = cc.constructInterface(iNormal_[cellI], alphaLiquid_[cellI]);
            iPoint_[cellI] = p.refPoint();
            iArea_[cellI] = cc.cutArea();
            
            // Save gas and solid portion centroids
            gasC_[cellI] = cc.lostCentroid();
            liquidC_[cellI] = cc.cutCentroid();
            
            deltaL_[cellI] = mag((liquidC_[cellI]-iPoint_[cellI]) & iNormal_[cellI]);
            deltaV_[cellI] = mag((gasC_[cellI]-iPoint_[cellI]) & iNormal_[cellI]);
            
            // DEBUGGING PURPOSES
            if (iArea_[cellI] < SMALL)
            {
                Info<< "WARNING: Cut area is " << iArea_[cellI]
                    << " with alphaLiquid = " << alphaLiquid_[cellI] << endl;
            }
        }
        else
        {
            iArea_[cellI] = 0.0;
            iPoint_[cellI] = vector::zero;
            deltaL_[cellI] = 0.5*Foam::pow(mesh_.V()[cellI], 1.0/3.0);
            deltaV_[cellI] = 0.5*Foam::pow(mesh_.V()[cellI], 1.0/3.0);
            if (alphaVapor_[cellI] > 0.5)
            { //vapor cell
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
    iArea_.correctBoundaryConditions();
    gasC_.correctBoundaryConditions();
    liquidC_.correctBoundaryConditions();
}

// Given alpha, calculate derived interface fields
void Foam::PLICInterface::correct()
{
    // Clip un-reconstructed alpha portions
    alphaLiquid_ *= pos(alphaLiquid_ - reconstructTol_);
    forAll(alphaLiquid_, cellI)
    {
        if( alphaLiquid_[cellI] > 1.0 - reconstructTol_ )
        {
            alphaLiquid_[cellI] = 1.0;
        }
    }
    alphaLiquid_.correctBoundaryConditions();
    
    // Set alphaVapor
    alphaVapor_ = scalar(1.0) - alphaLiquid_;
    

    
    // Step 1: Identify intermediate cells based on alphaLiquid_ and reconstructTol_
    intermeds_ = pos(alphaLiquid_ - reconstructTol_)
                               *pos(1.0 - reconstructTol_ - alphaLiquid_);

    // Step 2: Calculate interface normal in intermediate cells
    calculateInterfaceNormal();

    calcInterfacePlanes();
    
    // Step 0: Calculate new curvature field (kappaI_) from alphaLiquid_
    updateKappa(false);
        
    //calculate alphaf
    calcAlphaf();

    // Calculate transfer weights on small cells
    calcTransferWeights();
    
    liquidCells_ = liquidCells();
    smallLiquidCells_ = smallLiquidCells();
    noLiquidCells_ = noLiquidCells();
    gasCells_ = gasCells();
    smallGasCells_ = smallGasCells();
    noGasCells_ = noGasCells();
}



// Calculate alpha that excludes small gas cells
Foam::tmp<Foam::volScalarField> Foam::PLICInterface::alphaVaporCorr() const
{
    return alphaVapor_ * pos(alphaVapor_ - alphaMin_);
}

void Foam::PLICInterface::calcTransferWeights()
{
    Info<< "Calculating small cell face linking weights" << endl;
    
    // If negative, alpha is a small cell
    volScalarField alphaShiftV = alphaVapor_ - alphaMin_;
    volScalarField alphaShiftL = alphaLiquid_ - alphaMin_;
    
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbor = mesh_.neighbour();

    // Calculate transfer weights
    forAll(wL_, faceI)
    {
        label own = owner[faceI];
        label nei = neighbor[faceI];
        

        if (alphaShiftL[own] * alphaShiftL[nei] * alphaf_[faceI] < 0.0)
        { //one is small, one is not, and they share a gas boundary

            label sc = (alphaShiftL[own] < 0.0) ? own : nei;

            wL_[faceI] = mag
            (
                iNormal_[sc] & mesh_.Sf()[faceI]
            ) * alphaf_[faceI];
            
            alphaf_[faceI] = 0.0;
        }
        
        if (alphaShiftV[own] * alphaShiftV[nei] * (1.0 - alphaf_[faceI]) < 0.0)
        { //one is small, one is not, and they share a liquid boundary

            label sc = (alphaShiftV[own] < 0.0) ? own : nei;

            wV_[faceI] = mag
            (
                iNormal_[sc] & mesh_.Sf()[faceI]
            ) * (1.0 - alphaf_[faceI]);

            alphaf_[faceI] = 1.0;
        }
        
        // if both cells are small, zero out alphaf between then
        if (alphaShiftL[own] < 0.0 && alphaShiftL[nei] < 0.0)
        {
            alphaf_[faceI] = 0.0;
        }
        
        // if both cells are small, zero out alphaf between then
        if (alphaShiftV[own] < 0.0 && alphaShiftV[nei] < 0.0)
        {
            alphaf_[faceI] = 1.0;
        }
    }
    
    // Now set weights on parallel patches
    const volScalarField::GeometricBoundaryField& alphaShiftLBf = 
        alphaShiftL.boundaryField();
    const volScalarField::GeometricBoundaryField& alphaShiftVBf = 
        alphaShiftV.boundaryField();
    const volVectorField::GeometricBoundaryField& iNormalBf = 
        iNormal_.boundaryField();
        
    surfaceScalarField::GeometricBoundaryField& wLBf =
        wL_.boundaryField();
    surfaceScalarField::GeometricBoundaryField& wVBf =
        wV_.boundaryField();
    surfaceScalarField::GeometricBoundaryField& alphafBf =
        alphaf_.boundaryField();
    
    forAll(alphafBf, patchI)
    {
        const fvPatchScalarField& alphaShiftLPf = alphaShiftLBf[patchI];
        const fvPatchScalarField& alphaShiftVPf = alphaShiftVBf[patchI];
        const fvPatchVectorField& iNormalPf = iNormalBf[patchI];
        
        scalarField& alphafPf = alphafBf[patchI];
        scalarField& wLPf = wLBf[patchI];
        scalarField& wVPf = wVBf[patchI];
        
        const labelList& pFaceCells = mesh_.boundary()[patchI].faceCells();
    
        if (alphaShiftLPf.coupled())
        {
            // Get values across parallel patch
            const scalarField alphaShiftLPNf(alphaShiftLPf.patchNeighbourField());
            const scalarField alphaShiftVPNf(alphaShiftVPf.patchNeighbourField());
            const vectorField iNormalPNf(iNormalPf.patchNeighbourField());
            
            const fvPatch& meshPf = mesh_.boundary()[patchI];
            
            forAll(alphafPf, pFaceI)
            {
                label pfCellI = pFaceCells[pFaceI];
                
                //Calculate weight on small-large cell faces and close
                // small cell faces
                if ( alphaShiftL[pfCellI] * alphaShiftLPNf[pFaceI]
                     * alphafPf[pFaceI] < 0.0)
                {
                
                    vector scNorm = (alphaShiftL[pfCellI] < 0.0)
                                    ? iNormal_[pfCellI]
                                    : iNormalPNf[pFaceI];

                    wLPf[pFaceI] = mag
                    (
                        scNorm & meshPf.Sf()[pFaceI]
                    ) * alphafPf[pFaceI];
                    
                    alphafPf[pFaceI] = 0.0;
                }
                
                if ( alphaShiftV[pfCellI] * alphaShiftVPNf[pFaceI]
                     * (1.0 - alphafPf[pFaceI]) < 0.0)
                {
                
                    vector scNorm = (alphaShiftV[pfCellI] < 0.0)
                                    ? iNormal_[pfCellI]
                                    : iNormalPNf[pFaceI];

                    wVPf[pFaceI] = mag
                    (
                        scNorm & meshPf.Sf()[pFaceI]
                    ) * (1.0 - alphafPf[pFaceI]);
                    
                    alphafPf[pFaceI] = 1.0;
                }
                
                //Close faces between small liquid cells
                if (alphaShiftL[pfCellI] < 0.0 && alphaShiftLPNf[pFaceI] < 0.0)
                {
                    alphafPf[pFaceI] = 0.0;
                }
                
                //Close faces between small gas cells
                if (alphaShiftV[pfCellI] < 0.0 && alphaShiftVPNf[pFaceI] < 0.0)
                {
                    alphafPf[pFaceI] = 1.0;
                }
            }
        }
    }
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
    volScalarField alphaShift = alphaVapor_ - alphaMin_;
    
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
Foam::PLICInterface::getRefinementField() const
{
    // Force all cells on OR near the interface to refine to the maximum level
    dimensionedScalar C("C",dimLength,1e8);
    tmp<volScalarField> tRefinementField = C*mag(fvc::grad(alphaVapor_));

    return tRefinementField;
}

// Find the volume fraction of LIQUID in the flux on each face
Foam::tmp<surfaceScalarField> Foam::PLICInterface::phiAlphaFrac() const
{
    Info<< "Calculating phiAlphaFrac using reconstructed interfaces" << endl;

    //Normal definition applied first, applicable away from interfaces
    tmp<surfaceScalarField> tphiAlphaFrac = fvc::interpolate(alphaLiquid_);
    surfaceScalarField& phiAlphaFrac = tphiAlphaFrac();

    //Then adjust at the interface regions
    //  the scheme for this interpolation is set in fvSchemes
    surfaceVectorField Uf = fvc::interpolate(U_);
    
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


        if (   iNormal_[uwCell] != vector::zero && iPoint_[uwCell] != vector::zero
            && mag(Uf[faceI] & mesh_.Sf()[faceI]) > SMALL)
        {
            // Make a cuttable cell by extruding faceI in the direction of the
            // face's inverse velocity (-Uf*dT)
            cuttableCell pf(mesh_, faceI, -Uf[faceI]*dT);

            // Then cut the extruded face by the interface plane in the upwind
            // cell
            phiAlphaFrac[faceI] = pf.cut( plane(iPoint_[uwCell], iNormal_[uwCell]) );    
        }

        else if (mag(Uf[faceI] & mesh_.Sf()[faceI]) <= SMALL)
        { //phi is SMALL too, no flow through face
            phiAlphaFrac[faceI] = 0.0;
        }

        else if (mag(alphaLiquid_[own] - alphaLiquid_[nei]) > 0.1)
        { //interface coincident with cell boundary
            phiAlphaFrac[faceI] = alphaLiquid_[uwCell];
        }
    }

    //Loop through boundary patches to set phiAlpha at processor patches
    // Assume that the nominal treatment of phiAlpha at boundary patches is ok
    // for now.

    const volScalarField::GeometricBoundaryField& alphaBf = alphaLiquid_.boundaryField();
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

                    if (iNormal_[pfCellI] != vector::zero && iPoint_[pfCellI] != vector::zero)
                    {
                        // Cut the polyhedron by the reconstructed interface plane in cell on this proc
                        cutSum += pf.cut( plane(iPoint_[pfCellI], iNormal_[pfCellI]) );
                        nCuts++;
                    }

                    if (iNormalPNf[pFaceI] != vector::zero && iPointPNf[pFaceI] != vector::zero)
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

tmp<surfaceVectorField> Foam::PLICInterface::liquidUf() const
{

    tmp<surfaceVectorField> tUf = fvc::interpolate(U_);
    surfaceVectorField& Uf = tUf();
    
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbor = mesh_.neighbour();

    
    //Loop through internal faces (internal to this processor)
    forAll(Uf, faceI)
    {
        label own = owner[faceI];
        label nei = neighbor[faceI];
        
        if( alphaLiquid_[own] < reconstructTol_ && alphaLiquid_[nei] > reconstructTol_ )
        {
            Uf[faceI] = U_[nei];
        }
        else if( alphaLiquid_[nei] < reconstructTol_ && alphaLiquid_[own] > reconstructTol_ )
        {
            Uf[faceI] = U_[own];
        }
    }
    
    return tUf;
}

void Foam::PLICInterface::advect
(
    const volScalarField& liquidVolSource,
    const surfaceScalarField& phi
)
{
    // Get the volume flux of liquid phase on faces
    //tmp<surfaceScalarField> tphi = fvc::interpolate(U_) & mesh_.Sf();
    //const surfaceScalarField& phi = tphi();
    
    phiAlphaLiquid_ = phiAlphaFrac() * phi;
    
    volScalarField divU(fvc::div(phi));
    
    Info<< "Max divergence = " << Foam::max(divU).value() << endl;

    volScalarField::DimensionedInternalField Su
    (
        IOobject
        (
            "Su",
            mesh_.time().timeName(),
            mesh_
        ),
        // Divergence term is handled explicitly to be
        // consistent with the explicit transport solution
        //divU*min(alphaLiquid_, scalar(1)) + liquidVolSource
        //divU*pos(alphaLiquid_-0.5) + liquidVolSource // Weymouth 2010
        liquidVolSource
    );  
        
    MULES::explicitSolve
    (
        geometricOneField(),
        alphaLiquid_, 
        phi, 
        phiAlphaLiquid_, 
        zeroField(), //Sp
        Su, //Su
        1,           //alphaMax
        0            //alphaMin
    );
    
    phiAlphaVapor_ = phi - phiAlphaLiquid_;
        
    Info<< "Liquid phase volume fraction pre-correct = "
        << "  Min(alpha) = " << min(alphaLiquid_).value()
        << "  Max(alpha) = " << max(alphaLiquid_).value()
        << endl;
        
    //Correct the interface parameters using the new alphaLiquid
    //  This updates iNormal, iPoint, iArea, centroids, alphaVapor, ...
    correct();

    Info<< "Liquid phase volume fraction post-correct = "
        << alphaLiquid_.weightedAverage(mesh_.Vsc()).value()
        << "  Min(alpha) = " << min(alphaLiquid_).value()
        << "  Max(alpha) = " << max(alphaLiquid_).value()
        << endl;
}


Foam::tmp<Foam::volScalarField> Foam::PLICInterface::noLiquidCells() const
{
    // alphaL = 0
    return neg(alphaLiquid_ - SMALL);
}

Foam::tmp<Foam::volScalarField> Foam::PLICInterface::noGasCells() const
{
    // alphaV = 0
    return neg(alphaVapor_ - SMALL);
}

Foam::tmp<Foam::volScalarField> Foam::PLICInterface::smallLiquidCells() const
{
    // 0 < alphaL < min
    return neg(alphaLiquid_ - alphaMin_)*pos(alphaLiquid_ - SMALL);
}

Foam::tmp<Foam::volScalarField> Foam::PLICInterface::smallGasCells() const
{
    // 0 < alphaV < min
    return neg(alphaVapor_ - alphaMin_)*pos(alphaVapor_ - SMALL);
}
        
Foam::tmp<Foam::volScalarField> Foam::PLICInterface::liquidCells() const
{
    // alphaL > min
    return pos(alphaLiquid_ - alphaMin_);
}

Foam::tmp<Foam::volScalarField> Foam::PLICInterface::gasCells() const
{
    // alphaV > min
    return pos(alphaVapor_ - alphaMin_);
}

Foam::tmp<Foam::volScalarField> Foam::PLICInterface::noCells() const
{
    return neg(alphaVapor_ + 100.0);
}


// ************************************************************************* //
