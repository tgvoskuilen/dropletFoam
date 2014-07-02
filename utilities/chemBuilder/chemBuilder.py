
import inspect
import os

class Specie(object):
    def __init__(self,name):
        self.name = name
        self.isLiquid = (name[-1] == 'L')
        
    def __str__(self):
        return self.name
        
        

class Reaction(object):
    def __init__(self, lineStart, lineEnd, allLines):
        self.lines = allLines[lineStart:lineEnd]

    def set_zero_thirdbody(self,specie_list):
        if any(["THIRDBODY" in x for x in self.lines]):
            tbi = [i for i,x in enumerate(self.lines) if 'THIRDBODY' in x][0]
            
            self.lines.pop(tbi)

            for s in specie_list:
                self.lines.append('    %s/ 0.0/' % s)

    def __str__(self):
        return '\n'.join(self.lines)
    
class ChemkinInput(object):
    def __init__(self, filename=None, extra_species=None):
        self.elements = []
        if extra_species is not None:
            self.species = extra_species
        else:
            self.species = []
        self.reactions = []
            
        if filename is not None:
            basepath = os.path.dirname(inspect.getfile(ChemkinInput))
            
            # Read the input file
            with open(os.path.join(basepath,filename+'.inp'),'r') as f:
                lines = f.readlines()
                
            # Remove the carriage return characters from the lines
            lines = [x.rstrip() for x in lines]
            
            # Get section start lines
            es = lines.index('ELEMENTS')+1
            ss = lines.index('SPECIES')+1
            rs = lines.index('REACTIONS')+1
            
            # Get section end lines
            ends = [i for i,x in enumerate(lines) if 'END' in x]
            
            ee = min([x for x in ends if x >= es])
            se = min([x for x in ends if x >= ss])
            re = min([x for x in ends if x >= rs])
            
            # Read elements (can be on one or multiple lines)
            for i in range(es,ee):
                l = lines[i].split('!')[0]  # get only parts before the '!'
                entries = l.split(' ')
                for e in entries:
                    if e not in self.elements:
                        self.elements.append(e)
                        
            # Read species (can be on one or multiple lines)
            for i in range(ss,se):
                l = lines[i].split('!')[0]  # get only parts before the '!'
                entries = l.split(' ')
                for e in entries:
                    if e not in self.species:
                        self.species.append(e)
             
            # Read reactions - this is trickier, since they can be multi-line
            # I will assume that a line is the start of a reaction if it has an
            # '=' on it (and not in a comment)
            
            crs = None
            
            for i in range(rs,re):
                l = lines[i].split('!')[0]  # get only parts before the '!'
                if '=' in l:
                    if crs is not None:
                        self.reactions.append( Reaction(crs,i,lines) )
                    crs = i

            if crs is not None:
                self.reactions.append( Reaction(crs,re,lines) )

    def get_liquid_species(self):
        return [x for x in self.species if x[-1] == 'L']


    def set_zero_thirdbody(self,specie_list):
        specie_list += self.get_liquid_species()
        for r in self.reactions:
            r.set_zero_thirdbody(specie_list)

    def __add__(self,other):
        c = ChemkinInput() # empty
        
        # make element list (no duplicates)
        c.elements = list(set(self.elements+other.elements))
        c.elements.sort()
        
        # make specie list (no duplicates)
        c.species = list(set(self.species+other.species))
        c.species.sort()
        
        # make reactions list (duplicates not checked for)
        c.reactions = self.reactions + other.reactions
        
        return c
        
    def __str__(self):
        s = 'ELEMENTS\n'+' '.join(self.elements)+'\nEND\n\nSPECIES\n' + \
            ('! %d species\n'%len(self.species)) + \
            '\n'.join(self.species) + '\nEND\n\nREACTIONS\n' + \
            ('! %d reactions\n'%len(self.reactions))
        
        for r in self.reactions:
            s = s + str(r) + '\n'
        
        s = s + 'END'
        
        return s

    def write_file(self, filename):
        with open(filename,'w') as f:
            f.write(str(self))

def build_set(gas=None, liquid=None, extra_species=None):

    # first load both gas and liquid sets
    if gas is not None:
        gasChem = ChemkinInput(gas)
        
    if liquid is not None or extra_species is not None:
        liquidChem = ChemkinInput(liquid, extra_species)
            
        if gas is not None:
            # then insert liquid species (ending with 'L') into gas 
            # thirdbody reactions
            gasChem.set_zero_thirdbody(liquidChem.get_liquid_species())

        return liquidChem + gasChem
        
    elif gas is not None:
        return gasChem
        
    else:
        raise ValueError('Invalid chemistry building inputs')
    
    
if __name__ == "__main__":

    rxnset = build_set(gas='GasChem2',extra_species=["CH3NHNH2L","HNO3L"])
    
    rxnset.write_file('chem.inp')
    

