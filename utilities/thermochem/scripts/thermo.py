
import os
import csv
import time

import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)
    import numpy as np
    from scipy.optimize import curve_fit
    import matplotlib.pyplot as plt

#import ast
# d = ast.literal_eval(str)

def read_Sutherland_coeffs():
    LJcoeffFile = 'SutherlandCoeffs.csv'

    with open(LJcoeffFile,'r') as LJfile:
        LJreader = csv.reader(LJfile)
        coeffs = {}
        for row in LJreader:
            try:
                name = row[0]
                As = float(row[1])
                Ts = float(row[2])
                coeffs[name] = (As, Ts)
            except ValueError:
                pass
            
    return coeffs
    
def molwt(atoms):
    Ws = {'C':  12.011, 
          'H':  1.008, 
          'N':  14.007, 
          'O':  15.999,
          'AR': 39.948,
          'HE': 4.0026,
          'SI': 28.085,
          'E':  0.0} #electrons have essentially 0 mass
          
    W = 0.
    
    for a,n in atoms.iteritems():
        W = W + Ws[a]*float(n)  

    return W
    

def read_thermo_entry(name,lines):
    """
    Read an entry from a CHEMKIN style thermodynamic entry (4 lines) and
    put the data into a dictionary
    """
    # Read atomic composition, these could also be in cols 81-100, but rarely
    atoms = lines[0][24:44]
    a = {}
    while True:
        n = atoms[0:2]    # atom name
        n = n.strip()     # remove spaces
        istr = atoms[3:5] # number of atoms
        try:
            i = int(istr) # istr is a number
            if n and abs(i) > 0:
                try:
                    int(n) # n is not a number
                except ValueError:
                    a[n] = i
                    
        except ValueError:
            break
        atoms = atoms[5:]
            
    # Read Tlow, Tmid, and Thigh
    Tlow = float(lines[0][46:55])
    Thigh = float(lines[0][55:65])
    
    try:
        Tmid = float(lines[0][65:73]) #Tmid is optional
        if Tmid != 1000:
            #print "WARNING: Tmid for specie %s is %5.1f" % (name,Tmid)
            Tmid = 1000.
    except ValueError:
        Tmid = 1000.
            
    
    # Read coefficients
    high=[float(lines[1][0:15]), float(lines[1][15:30]),float(lines[1][30:45]),
          float(lines[1][45:60]),float(lines[1][60:75]),float(lines[2][0:15]),
          float(lines[2][15:30]) ]
    low =[float(lines[2][30:45]), float(lines[2][45:60]),float(lines[2][60:75]),
          float(lines[3][0:15]),float(lines[3][15:30]),float(lines[3][30:45]),
          float(lines[3][45:60]) ]

    data = {'Name':name, 'Atoms':a, 'Tlow':min([200.0, Tlow]), 'Thigh':Thigh,
            'W':molwt(a), 'Tmid':Tmid, 'lowCpCoeffs':low, 'highCpCoeffs':high}

    return data


def atomic_match(a1,a2):
    """
    Return a score indicating how close two atomic sets are to one another
    """
    Ws = {'C':  12.011, 
          'H':  1.008, 
          'N':  14.007, 
          'O':  15.999,
          'AR': 39.948,
          'HE': 4.0026}
          
    scores = []
    
    for atom,num in a1.iteritems():
        if atom in a2:
            scores.append(abs(num - a2[atom])*Ws[atom])
        else:
            scores.append(abs(num)*Ws[atom])

    for atom,num in a2.iteritems():
        if atom not in a1:
            scores.append(abs(num)*Ws[atom])

    return sum(scores)
    
    
    
def read_thermo(thermoFile):

    with open(thermoFile,'r') as f:
        lines = f.readlines()

    for i,line in enumerate(lines):
        lines[i] = line.rstrip() #remove line end char
        
    data = {}
    currentEntry = []

    for line in lines:
        try:
            ecstr = line[79]
            try:
                ec = int(ecstr)
            except ValueError:
                print "Could not read end char from line: '%s'" % line
                raise
                
            currentEntry.append(line)
            if ec == 4:
                name = currentEntry[0][0:18]
                name = name.strip().split(' ')[0]
                try:
                    data[name] = read_thermo_entry(name,currentEntry)
                except KeyError:
                    print "\n *** Error reading entry '%s' ***\n" % name
                finally:
                    currentEntry = []
            
        except IndexError:
            pass

    return data
    
    
def Sutherland(T,Ts,As):
    """
    Sutherland transport model for viscosity
    """
    return As*T**1.5/(T+Ts)
    
    
def calc_NASA_mu(Tmax,muCoeffs):
    """
    muCoeffs is an array of dictionaries, each of which has T1 and T2 and
    A,B,C,D for the fit ln(mu) = A*ln(T) + B/T + C/T**2 + D
    """

    Tmin = min([x['T1'] for x in muCoeffs])
    
    T = np.linspace(Tmin,Tmax,500)
    mu = np.zeros_like(T)
    
    for c in muCoeffs:
        muC = 1e-7*np.exp(c['A']*np.log(T) + c['B']/T + c['C']/T/T + c['D'])
        for i in np.ndindex(T.shape):
            if T[i] >= c['T1'] and T[i] < c['T2']:
                mu[i] = muC[i]
    
    return T,mu


    
def read_transport_entry(lines):
    name1 = lines[0][0:15].strip()
    name2 = lines[0][15:30].strip()
    
    name = name1 if not name2 else name1+'-'+name2
    species = [name1]
    if name2:
        species.append(name2)

    data = {'Name': name, 'Species': species, 
            'Viscosity': [], 'Conductivity': [] }

    lines = lines[1:]
    
    for line in lines:
        l = line.strip()
        
        entry = {'T1': float(l[1:8]),
                 'T2': float(l[8:17]),
                 'A' : float(l[19:34]),
                 'B' : float(l[34:49]),
                 'C' : float(l[49:64]),
                 'D' : float(l[64:79])}
        
        if l[0] == 'V':
            data['Viscosity'].append(entry)
        elif l[0] == 'C':
            data['Conductivity'].append(entry)
        else:
            print "Unknown line format:", l
            raise ValueError


    if len(species) == 1:
        #fit viscosity data to sutherland model
        Tmax = 3000
        T,muNASA = calc_NASA_mu(Tmax,data['Viscosity'])
        
        opt,cov = curve_fit(Sutherland, T, muNASA, p0=[200,1e-5])
        
        err = 0.5*cov[0][0]**0.5/opt[0] + 0.5*cov[1][1]**0.5/opt[1]
        
        if err > 0.05:
            print "Fit error for %s is %5.4f (%5.2f, %4.3e)" % (name, err, opt[0], opt[1])
            muFit = Sutherland(T,opt[0],opt[1])
            f = plt.figure()
            ax = f.add_subplot(111)
            ax.plot(T,muNASA,'ro')
            ax.plot(T,muFit,'-k')
            f.suptitle(name)
            f.show()
            raw_input("This fit is off by more than 5%...")
            plt.close(f)

        data['Ts'] = opt[0]
        data['As'] = opt[1]

    return data


def read_transport(transportFile):
    
    with open(transportFile,'r') as f:
        lines = f.readlines()

    for i,line in enumerate(lines):
        lines[i] = line.rstrip() #remove line end char
        
    data = {}
    currentEntry = []

    for line in lines:
        if line:
            if line[0] != '!':
                if line[0] != ' ' and currentEntry:
                    entry = read_transport_entry(currentEntry)
                    data[entry['Name']] = entry
                    currentEntry = []
                    
                currentEntry.append(line)
                
    if currentEntry:
        entry = read_transport_entry(currentEntry)
        data[entry['Name']] = entry

    return data
    
    
    
def write_openfoam_database(thermo, transport, dbfilename):
    """
    Write python-readable database of openfoam-required inputs for items
    found in thermo and transport.
    """
    dataSet = {}
    
    for name,data in thermo.iteritems():
        dataSet[name] = data
        if name in transport:
            dataSet[name]['As'] = transport[name]['As']
            dataSet[name]['Ts'] = transport[name]['Ts']
        else:
            #TODO
            dataSet[name]['As'] = 1e-5
            dataSet[name]['Ts'] = 200
    
    with open(dbfilename,'w') as dbfile:
        for name,data in dataSet.iteritems():
            dbfile.write(str(data)+"\n")
    


    
