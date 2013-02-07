{

    //Net mass accumulation
    scalar dM = Foam::sum
    (
        (rho.internalField() - rho.oldTime().internalField()) * mesh.V()
    );

    //Boundary mass net outflux
    word patchName = "atmosphere";
    label patchID = mesh.boundaryMesh().findPatchID(patchName);
    scalar netOutFlux = 0.0;

    //patchID = -1 if patch cannot be found
    if (patchID >= 0)
    {
        netOutFlux = Foam::sum
        (
            phi.boundaryField()[patchID]
        )*mesh.time().deltaT().value();
    }

    // Add up total accumulation and net outflux over all processors
    Foam::reduce(netOutFlux, sumOp<scalar>());
    Foam::reduce(dM, sumOp<scalar>());

    Info<< "Mass continuity error = "
        << Foam::mag(dM + netOutFlux)
        << " of total mass "
        << Foam::sum(rho.internalField()*mesh.V())
        << endl;
}

