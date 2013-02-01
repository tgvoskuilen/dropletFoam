// * * * * * ALPHA FIELDS * * * * * * * * * * * * * * * * * * * * * * * * //
    volScalarField alphaVapor
    (
        IOobject
        (
            "alphaVapor",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

