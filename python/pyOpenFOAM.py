
import os
import re
import string

# Database of sutherland transport coefficients
# http://www.osti.gov/bridge/servlets/purl/10181503-tLKFUg/10181503.pdf
sutherlandData = [{'Specie':'N2',      'As':3.61e-7,   'Ts':-9.549},
                  {'Specie':'O2',      'As':2.29e-6,   'Ts':164.4},
                  {'Specie':'H2',      'As':1.96e-7,   'Ts':2.187},
                  {'Specie':'H',       'As':3.95e-7,   'Ts':0},
                  {'Specie':'O',       'As':1.15e-6,   'Ts':0},
                  {'Specie':'OH',      'As':1.10e-6,   'Ts':0},
                  {'Specie':'H2O',     'As':1.60e-6,   'Ts':0} ]

#Alternately, use Lennard-Jones parameters



class transportProperty(object):
    """
    Class for each transport property. Defaults to gas properties if entries
    are omitted.
    """
    def __init__(self, name, lines):
        self.name = name

        # The viscosity model is only used for liquids
        self.transportModel = 'Newtonian'
        self.nu = 8.91e-7

        # Sutherland transport properties
        self.As = 1.67212e-06
        self.Ts = 170.672

        # Set defaults
        self.rho0 = 0
        self.Z0 = 1
        self.kth = 0.1
        self.eps = 0.0
        self.radius = [0, 0, 0]
        self.center = [0, 0, 0]
        self.U0 = [0, 0, 0]
        self.delVapor = [2e-5, 2e-5, 2e-5]
        self.dimensions = 2
        
        self.sigma0 = 0
        self.sigma_a = 1
        
        self.vapor_name = None
        self.vapor_layer = []
        
        self.Lb = 0
        self.Tb = 0
        self.L_a = 0.38
        self.betaM = 0.5
        
        self.Tc = 0
        self.Pc = 0
        
        self.PvCoeffs = [0, 0, 0]
        
        self.D = 2e-5

        # Read values to overwrite defaults
        for i,line in enumerate(lines):
            if 'transportModel' in line:
                self.transportModel = string.split(line)[1]

            if 'vapor' in line:
                self.vapor_name = string.split(line)[1]
                
            if 'surrounding' in line:
                value = float(string.split(line)[1])
                name = string.split(line)[2]
                self.vapor_layer.append( (value, name) )

            if 'Lb' in line:
                self.Lb = float(string.split(line)[1])
                
            if 'Tb' in line:
                self.Tb = float(string.split(line)[1])
                
            if 'powaL' in line:
                self.L_a = float(string.split(line)[1])
                
            if 'betaM' in line:
                self.betaM = float(string.split(line)[1])
                  
            if 'Tc' in line:
                self.Tc = float(string.split(line)[1])
                
            if 'Pc' in line:
                self.Pc = float(string.split(line)[1])
                                     
            if 'PvCoeffs' in line:
                self.PvCoeffs = [float(s) for s in string.split(line)[1:4]]
                
            if 'nu' in line:
                self.nu = float(string.split(line)[1])

            if 'D' in line:
                self.D = float(string.split(line)[1])

            if 'rho0' in line:
                self.rho0 = float(string.split(line)[1])

            if 'Z0' in line:
                self.Z0 = float(string.split(line)[1])

            if 'kth' in line:
                self.kth = float(string.split(line)[1])

            if 'eps' in line:
                self.eps = float(string.split(line)[1])

            if 'radius' in line:
                self.radius = [float(s) for s in string.split(line)[1:4]]

            if 'delVapor' in line:
                self.delVapor = [float(s) for s in string.split(line)[1:4]]
                
            if 'dimensions' in line:
                self.dimensions = int(string.split(line)[1])

            if 'U0' in line:
                self.U0 = [float(s) for s in string.split(line)[1:4]]

            if 'center' in line:
                self.center = [float(s) for s in string.split(line)[1:4]]

            if 'sigma0' in line:
                self.sigma0 = float(string.split(line)[1])
                
            if 'powasigma' in line:
                self.sigma_a = float(string.split(line)[1])
                
                
                
    def write(self, f, name=None):
        """
        Write the transportProperty dictionary entry for this object
        """
        if name is None:
            name = self.name

        f.write('    '+name+'\n    {\n')
        f.write('       transportModel  '+self.transportModel+';\n')
        f.write('       nu              nu [ 0 2 -1 0 0 0 0 ] '+str(self.nu)+';\n\n')
        f.write('       rho0            rho0 [ 1 -3 0 0 0 0 0 ] '+str(self.rho0)+';\n')
        f.write('       Z0              Z0 [ 0 0 0 0 0 0 0 ] '+str(self.Z0)+';\n\n')
        f.write('       kth             kth [ 1 1 -3 -1 0 0 0] '+str(self.kth)+';\n\n')
        f.write('       eps             eps [ 0 0 0 0 0 0 0] '+str(self.eps)+';\n\n')
        f.write('       D               D [ 0 2 -1 0 0 0 0] '+str(self.D)+';\n\n')
        
        f.write('       surrounding\n       (\n');
        if self.vapor_layer:
            for vap in self.vapor_layer:
                f.write('           '+vap[1]+' '+str(vap[0])+'\n')
        elif self.vapor_name is not None:
            f.write('           '+self.vapor_name+' 1.0\n')
        f.write('       );\n')
            
        if self.vapor_name is not None:
            if (self.Lb == 0 or self.Tc == 0 or self.Pc == 0 or 
                self.Tb == 0 or self.PvCoeffs[0] == 0):
                raise ValueError("Critical evaporation components not set")
            

                
            f.write('       vapor           '+self.vapor_name+';\n')
            f.write('       Lb              Lb [ 1 2 -2 0 -1 0 0] '+str(self.Lb)+';\n')
            f.write('       Tb              Tb [ 0 0 0 1 0 0 0] '+str(self.Tb)+';\n')
            f.write('       La              '+str(self.L_a)+';\n')
            f.write('       betaM           '+str(self.betaM)+';\n\n')
            
            f.write('       Tc              Tc [ 0 0 0 1 0 0 0] '+str(self.Tc)+';\n')
            f.write('       Pc              Pc [ 1 -1 -2 0 0 0 0] '+str(self.Pc)+';\n')
            f.write('       PvCoeffs        ('+str(self.PvCoeffs[0])+' '+str(self.PvCoeffs[1])+' '+str(self.PvCoeffs[2])+');\n\n')
            
            f.write('       sigma0          sigma0 [ 1 0 -2 0 0 0 0] '+str(self.sigma0)+';\n')
            f.write('       sigmaa          '+str(self.sigma_a)+';\n\n')
            
            

        if sum(self.radius) > 0.0:
            f.write('       radius          ('+str(self.radius[0])+' '+str(self.radius[1])+' '+str(self.radius[2])+');\n')
            f.write('       center          ('+str(self.center[0])+' '+str(self.center[1])+' '+str(self.center[2])+');\n')
            f.write('       U0              ('+str(self.U0[0])+' '+str(self.U0[1])+' '+str(self.U0[2])+');\n')
            f.write('       delVapor        ('+str(self.delVapor[0])+' '+str(self.delVapor[1])+' '+str(self.delVapor[2])+');\n')
            f.write('       dimensions      '+str(self.dimensions)+';\n')

        f.write('    }\n\n')
        

def read_transport_db():
    """
    Reads the properties.db text file into a Python dictionary
    """
    db_f = open('chemkin/properties.db','r').readlines()

    for i, s in enumerate(db_f):
        db_f[i] = s.rstrip()

    start = 0
    end = 0
    db = {}
    db['Default'] = transportProperty('Default',[])

    for i,s in enumerate(db_f):
        if s == '{':
            start = i-1
        if s == '}':
            end = i
            db[db_f[start]] = transportProperty(db_f[start],db_f[start+1:end])

    return db


def add_species_to_thermo_dict(species, default):
    """
    Takes the list of species and adds all non liquid ones (no "L") to the
    "Vapor" subspecies list
    """
    f_thermoDict = open('constant/thermophysicalProperties','r')
    thermo_lines = f_thermoDict.readlines()
    f_thermoDict.close()

    in_vapor = False
    in_subspecies = False
    ss_line_start = 0
    ss_line_end = 0

    for i,line in enumerate(thermo_lines):
        if "Vapor" in line:
            in_vapor = True

        if "subspecies" in line and in_vapor:
            in_subspecies = True
            ss_line_start = i+2

        if ");" in line and in_vapor and in_subspecies:
            ss_line_end = i
            break

    new_lines = []
    for specie in species:
        if "L" not in specie:
            value = 1.0 if specie == default else 0.0
            new_lines.append("           %s   %f\n" % (specie, value))

    new_thermo = thermo_lines[:ss_line_start] + new_lines  \
                + thermo_lines[ss_line_end:]
    
    for line in new_thermo:
        print line


def fix_thermo_transport_props():
    """
    Fix the sutherland transport properties in "thermo" which openfoam always
    puts to the same constant values
    """
    f_thermo = open('constant/thermo','r')
    thermo_lines = f_thermo.readlines()
    f_thermo.close()

    new_thermo = []
    currentSpecie = ''
    hasNewData = False
    
    for line in thermo_lines:
        if line[0] != '{' and line[0] != '}' and line[0] != ' ':
            currentSpecie = line.strip()
            matches = [x for x in sutherlandData if 
                x['Specie'] == currentSpecie]

            hasNewData = len(matches)==1

            if hasNewData:
                As = matches[0]['As']
                Ts = matches[0]['Ts']
                print "Updating Sutherland coefficients for ", currentSpecie


        if 'As' in line and hasNewData:
            new_thermo.append('        As              %e;\n' % As)

        elif 'Ts' in line and hasNewData:
            new_thermo.append('        Ts              %f;\n' % Ts)
            hasNewData = False

        else:
            new_thermo.append(line)

    f_thermo = open('constant/thermo','w')
    f_thermo.writelines(new_thermo)
    f_thermo.close()


def make_thermo_properties(species, Ti=280):
    """
    Make the thermophysicalProperties dict using the properties.db file
    """
    db = read_transport_db()

    f_tp = open('constant/thermophysicalProperties','w')
    f_tp.write(\
r"""/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.1.1                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      thermophysicalProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

thermoType hsMultiphaseMixtureThermo<reactingMixture<gasThermoPhysics>>;

chemistryReader foamChemistryReader;

foamChemistryFile "$FOAM_CASE/constant/reactions";

foamChemistryThermoFile "$FOAM_CASE/constant/thermo";

phases
(
""")

    # Write entries for each specie
    for specie in species:
        if specie in db:
            db[specie].write(f_tp)
        else:
            db['Default'].write(f_tp, name=specie)

    f_tp.write(');\n\n')

    # Write Ti
    f_tp.write('Ti    Ti [0 0 0 1 0 0 0] %f;\n\n' % Ti)

    # Write sigma table for all non-zero pairs
    """f_tp.write('sigmas\n')
    f_tp.write('(\n')

    for i in range(0,len(species)):
        for j in range(i+1,len(species)):

            if species[i] in db:
                tp1 = db[species[i]]
            else:
                tp1 = db['Default']

            if species[j] in db:
                tp2 = db[species[j]]
            else:
                tp2 = db['Default']

            # this is a rather case-specific set of criteria
            if (tp1.sigma != 0.0 or tp2.sigma != 0.0) and tp1.sigma != tp2.sigma:
                f_tp.write('  ('+species[i]+' '+species[j]+') '+str(max([tp1.sigma,tp2.sigma]))+'\n')

    f_tp.write(');\n\n')"""
    f_tp.write('// ************************************************************************* //')
    f_tp.close()


def make_zero_dir(species, default):
    """
    Make the "0" directory for this list of species
    """
    os.system('rm -rf 0')
    os.system('cp -r 0.org 0')
    os.system('rm -f 0/specieX')
    os.system('rm -f 0/alphaX')

    for specie in species:
        if specie == default:
            value = 1.0
        else:
            value = 0.0
        os.system("sed"+
                  " -e s/specieX/"+specie+"/g"+
                  " -e s/DEFAULT_VALUE/"+str(value)+"/g"+
                  " 0.org/specieX > 0/"+specie)
        os.system("sed"+
                  " -e s/alphaX/alpha"+specie+"/g"+
                  " -e s/DEFAULT_VALUE/"+str(value)+"/g"+
                  " 0.org/alphaX > 0/alpha"+specie)


def read_species():
    """
    Reads the list of species from the "reactions" file
    """
    start = 0
    ns = 0
    lines = open('constant/reactions','r').readlines()
    for i,s in enumerate(lines):
        if "species" in s:
            ns = int(lines[i+1].rstrip())
            start = i+2
        lines[i] = s.rstrip()

    return lines[start+1:start+ns+1]
    


def get_proc_dirs():
    """
    Returns a list of all processor directories. For example:
     ['processor0','processor1','processor2']
    """
    proc_dirs = []

    dirs = [ d for d in os.listdir(os.getcwd()) 
             if os.path.isdir(os.path.join(os.getcwd(), d)) ]

    for dirname in dirs:
        if re.match('processor[0-9]+',dirname):
            proc_dirs.append(dirname)
    
    return proc_dirs



def touch_foam_files(name):
    """ Creates empty .foam files for paraview """
    foamname = name+'.foam'
    os.system('touch '+foamname)

    for dirname in get_proc_dirs():
        os.system('touch '+dirname+'/'+foamname)


def read_inputs(argv):
    """ Reads the command line input for np """
    if len(argv) > 1:
        try:
            return int(argv[1])
        except ValueError:
            return 1
    else:
        return 1


def get_imbalance():
    """
    Gets the current parallel imbalance at the latest time.
    The imbalance is calculated from the filesize of each processor's
    most recent time step.
    """
    np = num_procs() 
    time_dirs = get_proc_times(0)

    if time_dirs is not None:
        time_values = [(float(td),td) for td in time_dirs]
        max_time_dir = max(time_values)[1]

        proc_sizes = []

        for p in range(0,np):
            dir_path = os.path.join(os.getcwd(),"processor"+str(p),max_time_dir)
            dir_size = get_dir_size(dir_path)
            proc_sizes.append(dir_size)

        proc_load = [abs(100.*np*s/float(sum(proc_sizes))-100.) for s in proc_sizes]

        return proc_load

    else:
        return None


def num_procs():
    """ Returns the number of processors used in this case """
    return length(get_proc_dirs())



def get_dir_size(start_path = '.'):
    """
    General function to recursively get the size of files in a directory
    """
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(start_path):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            total_size += os.path.getsize(fp)
    return total_size


def get_proc_times(proc=0, isParallel=True):
    """
    Returns a list of the time folders in a given processor
        (default folder: processor0). You can also get serial runs by setting
        isParallel to False (then proc argument is ignored)
    """
    if isParallel:
        proc_path = os.path.join(os.getcwd(), "processor"+str(proc))
    else:
        proc_path = os.getcwd()

    if os.path.isdir(proc_path):

        proc_dirs = [ d for d in os.listdir(proc_path) 
                      if os.path.isdir(os.path.join(proc_path, d)) ]

        time_dirs = []

        for dirname in proc_dirs:
            try:
                t = float(dirname)
                time_dirs.append(dirname)
            except ValueError:
                pass

        time_dirs.sort(key=float)

        return time_dirs

    else:
        return None


def get_sorted_time_folders(isParallel=True):
    """
    Returns a sorted list of time folders with the lowest time folder
    omitted. This can be passed to the parallel reconstruction routines.
    """
    time_dirs = get_proc_times(proc=0,isParallel=isParallel)
    if time_dirs is not None:
        time_dirs.sort(key=float)
        time_dirs.pop(0)
        return time_dirs
    else:
        return None


def get_proc_pairing(np, ratio=1.0):
    """
    Get the pair of factors of np that most closely gives the desired aspect
    ratio (ratio = num1 / num2)
    """
    f = factors(np)
    fr = [abs(float(x[0])/float(x[1])-ratio) for x in f]
    return f[fr.index(min(fr))]



def set_decompose_par(np, method="simple", ratio=1.0):
    """
    Sets the decomposeParDict file using the template .org file and two
    inputs (np and method).
    """
    pair = get_proc_pairing(np, ratio)
    os.system('rm -f system/decomposeParDict')
    os.system("sed"+
              " -e s/NUMPROCS/"+str(np)+"/g"+
              " -e s/DECOMPOSEPAR_METHOD/"+method+"/g"+
              " -e s/VALUE1/"+str(pair[0])+"/g"+
              " -e s/VALUE2/"+str(pair[1])+"/g"+
              " system/decomposeParDict.org > system/decomposeParDict")


def run(program, np=1, args="", hide=True):
    """
    General routine for running OpenFOAM programs, either in serial or parallel
    If np is set as 1, a serial run is used. Otherwise, mpirun is called.
    """ 
    if np > 1:
        cmd = "mpirun -np " + str(np) + " " + program + " -parallel " + args
    else:
        cmd = program + " " + args
        
    if hide:
        cmd = cmd + " > log." + program

    print " - Running Command: %s" % cmd
    os.system(cmd)


def factors(n):
    """ Return a list of factor pairs of the integer n"""
    fs = filter(lambda i: n % i == 0, range(1, n + 1))
    return [(x, n/x) for x in fs]



if __name__ == "__main__":
    print "This is pyOpenFOAM"




