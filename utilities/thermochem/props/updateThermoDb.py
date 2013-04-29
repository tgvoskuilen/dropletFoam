
import sys
sys.path.append('../scripts')

import thermo


thermoData = thermo.read_thermo('thermo.dat')
transportData = thermo.read_transport('transport.dat')

thermo.write_openfoam_database(thermoData, transportData, 'thermo.pydb')
