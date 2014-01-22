{

    //Net mass accumulation
    scalar dM = Foam::sum
    (
        (rho.internalField() - rho.oldTime().internalField()) * mesh.V()
    );

    //Boundary mass net outflux on all non-coupled boundaries
    const surfaceScalarField& rhoPhi = mixture.rhoPhi();
    scalar netOutFlux = 0.0;
    
    forAll(mesh.boundaryMesh(), i)
    {
        if ( !mesh.boundaryMesh()[i].coupled() )
        {
            netOutFlux += Foam::sum
            (
                rhoPhi.boundaryField()[i]
            )*mesh.time().deltaT().value();
        }
    }

    // Add up total accumulation and net outflux over all processors
    Foam::reduce(netOutFlux, sumOp<scalar>());
    Foam::reduce(dM, sumOp<scalar>());

    totalMass = Foam::sum(rho.internalField()*mesh.V());
    massError = dM + netOutFlux;
    
    Foam::reduce(totalMass, sumOp<scalar>());
    
    Info<< "Mass continuity error = "
        << Foam::mag(massError)
        << " of total mass "
        << totalMass
        << endl;
}

