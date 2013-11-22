{
    // Compute and write domain-averaged scalar quantities to a text file

    scalar maxImbalance = 0.0;

    if ( Pstream::parRun() )
    {
        //determine current level of imbalance for parallel runs
        label nGlobalCells = mesh.globalData().nTotalCells();
        scalar idealNCells = scalar(nGlobalCells)/scalar(Pstream::nProcs());
        scalar localImbalance = mag(scalar(mesh.nCells()) - idealNCells);
        Foam::reduce(localImbalance, maxOp<scalar>());
        maxImbalance = localImbalance/idealNCells;
    }

    myFile<< runTime.value() << token::TAB                 // Simulation time
          << runTime.deltaTValue() << token::TAB           // Current time step
          << runTime.elapsedClockTime() << token::TAB      // Wall clock time
          << mesh.globalData().nTotalCells() << token::TAB // Total cells
          << maxImbalance << token::TAB                    // Maximum imbalance
          << totalMass << token::TAB                       // Total mass
          << massError << endl;                            // Mass error

}
