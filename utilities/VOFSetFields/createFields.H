
    // Load VOFSetFields dictionary
    IOdictionary VOFDict
    (
        IOobject
        (
            "VOFSetFieldsDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    volScalarField alphaLiquid
    (
        IOobject
        (
            "alphaLiquid",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
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
    
    volScalarField Yremaining
    (
        IOobject
        (
            "Yremaining",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("Yr",dimless,1.0)
    );
    
    
    // Read entries from dictionary
    
    PtrDictionary<droplet> droplets
    (
        VOFDict.lookup("droplets"),
        droplet::iNew(mesh)
    );
    
    HashTable<scalar> defaultSpecieValues(VOFDict.lookup("defaultSpecies"));
    
    
    List<word> specieNames = defaultSpecieValues.toc();
    
    forAllConstIter(PtrDictionary<droplet>, droplets, dropI)
    {
        specieNames.append( dropI().species() );
    }
    
    
    PtrList<volScalarField> species(specieNames.size());
    
    Info<< "Loading species " << specieNames << endl;
    
    for(label i = 0; i < specieNames.size(); ++i)
    {
        //Follow basicMultiComponentMixture
        Info<< "Loading specie " << specieNames[i] << endl;
        
        IOobject header
        (
            specieNames[i],
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        );
        
        if (header.headerOk())
        {
            species.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        specieNames[i],
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )  
            );
        }
        else
        {
            volScalarField Ydefault
            (
                IOobject
                (
                    "Ydefault",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            species.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        specieNames[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    Ydefault
                )
            );
        }
    }
    
    
    alphaLiquid = dimensionedScalar("alphaLiquid",dimless,0.0);
    U = dimensionedVector("U",dimVelocity,vector::zero);
    forAll(species, i)
    {
        species[i] = dimensionedScalar(species[i].name(),dimless,0.0);
    }
    
    

    
