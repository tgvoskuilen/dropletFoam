
    // First create U and phi (mass flux) since thermo looks them up from mesh
    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    U.oldTime();
    
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rho",dimDensity,1.0)
    );
    
    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(rho*U) & mesh.Sf()
    );

    
    Info<< "Creating combustion model\n" << endl;
    autoPtr<combustionModels::rhoChemistryCombustionModel> combustion
    (
        combustionModels::rhoChemistryCombustionModel::New
        (
            mesh
        )
    );

    rhoChemistryModel& chemistry = combustion->pChemistry();
    hsReactionThermo& thermo = chemistry.thermo();
	
	Info<< "Extracting two phase mixture thermo:" << endl;

    hsTwophaseMixtureThermo<reactingMixture<gasThermoPhysics> >& mixture =
        dynamic_cast<hsTwophaseMixtureThermo<reactingMixture<gasThermoPhysics> >& >(thermo.composition());


    Info<< "Setting mixture pointers" << endl;
	mixture.setPtrs( combustion.operator->() );

    volScalarField refinementField
    (
        IOobject
        (
            "refinementField",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("refinementField",dimless,0)
    );
    
    
    rho = mixture.rho();
    
    
    
    // Construct turbulence model
    //   The liquid phase is assumed laminar
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            mixture.vapor().rho(),
            mixture.vapor().U(), 
            mixture.vapor().phi(),
            mixture
        )
    );

    // Set the turbulence into the combustion model
    combustion->setTurbulence(turbulence());


    #include "readGravitationalAcceleration.H"

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    //volScalarField& hs = thermo.hs();
    const volScalarField& T = mixture.T();
    volScalarField& p = mixture.p();
    p.oldTime();
 
    volScalarField hs
    (
        IOobject
        (
            "hs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.hs()
    );
    
    volScalarField Te
    (
        IOobject
        (
            "Te",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo.T()
    );
    Te.oldTime();


    
    
    volScalarField dQ
    (
        IOobject
        (
            "dQ",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("dQ", dimPower, 0.0)
    );




