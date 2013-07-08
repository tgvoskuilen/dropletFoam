
    // First create U and phi (volume flux) since thermo looks them up from mesh
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

    #include "createPhi.H"
    
    
    Info<< "Creating combustion model\n" << endl;
    autoPtr<combustionModels::rhoCombustionModel> reaction
    (
        combustionModels::rhoCombustionModel::New(mesh)
    );

    rhoReactionThermo& thermo = reaction->thermo();
    thermo.validate(args.executable(), "h", "e");
	
	Info<< "Extracting two phase mixture thermo:" << endl;
	
    /*if( !isA<heTwophaseMixtureThermo<reactingMixture<gasThermoPhysics> > >(thermo.composition()) )
    {
        FatalErrorIn
        (
            "reactingInterFoam::createFields.H"
        )   << "Specified thermo must be an hsTwophaseMixtureThermo"
            << exit(FatalError);
    }*/
    
    typedef heTwophaseMixtureThermo
    <
        Foam::rhoReactionThermo,
        Foam::SpecieMixture
        <
            Foam::reactingMixture
            <
                Foam::sutherlandTransport
                <
                    Foam::species::thermo
                    <
                        Foam::janafThermo<Foam::perfectGas<Foam::specie> >, 
                        Foam::sensibleEnthalpy
                    > 
                > 
            > 
        > 
    > reactingTwoPhaseMixtureThermo;
    
    
    reactingTwoPhaseMixtureThermo& mixture = 
        dynamic_cast<reactingTwoPhaseMixtureThermo&>(thermo.composition());
    

    //heTwophaseMixtureThermo<rhoReactionThermo, reactingMixture<gasThermoPhysics> >& mixture =
    //    dynamic_cast<hsTwophaseMixtureThermo<reactingMixture<gasThermoPhysics> >& >(thermo.composition());


    Info<< "Setting mixture pointers" << endl;
	mixture.setPtrs( reaction.operator->() );
    
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
        thermo.rho()
    );
    
    // Construct turbulence model
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U, 
            mixture.rhoPhi(), 
            mixture
        )
    );

    // Set the turbulence into the combustion model
    reaction->setTurbulence(turbulence());


    #include "readGravitationalAcceleration.H"

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    //volScalarField& hs = thermo.hs();
    //const volScalarField& T = mixture.T();
    volScalarField& T = mixture.T();
    T.oldTime(); //needed?
    volScalarField& p = mixture.p();
    p.oldTime(); //needed?
    
    //TODO: Read from dictionary
    dimensionedScalar p0("p0",dimPressure,1e5);
    
    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    //Force p to be consistent with p_rgh (gives better initial conditions
    // than forcing p_rgh to be consistent with p)
    p = p_rgh + p0 + rho*gh;
    
    /*volScalarField hs
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
    );*/
    
    /*volScalarField Te
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
    Te.oldTime();*/

    
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
    
    Info<< "Creating field dpdt\n" << endl;
    volScalarField dpdt("dpdt", fvc::ddt(p));

    Info<< "Creating field kinetic energy K\n" << endl;
    volScalarField K("K", 0.5*magSqr(U));




