
import os
import gzip
import string
import csv

def get_times():
    """
    Returns a sorted list of the time folders
    """
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
        time_dirs.pop(0) #remove zero time folder

        return time_dirs

    else:
        return None
        
def get_fields(folder):
    """
    Get a list of the field names from the first populated folder
    """
    time_path = os.path.join(os.getcwd(), folder)
    
    if os.path.isdir(time_path):
        gzFiles = [f for f in os.listdir(time_path)
                   if f.endswith('.gz')]
        
        fields = ['Time']
        for f in gzFiles:
            n,e = os.path.splitext(f)
            fields.append(n)
            
        return fields
    
    else:
        return None


def get_field_data(folder,fieldList):
    """
    Unzip and read the data from all the specified fields in a given folder
    """
    time_path = os.path.join(os.getcwd(), folder)
    
    if os.path.isdir(time_path):
        values = [float(folder)]
        for i,f in enumerate(fieldList):
            if i > 0:
                filePath = os.path.join(time_path,f+'.gz')
                fz = gzip.open(filePath,'rb')
                content = fz.read()
                fz.close()
                
                loc1 = string.find(content,'internalField')
                chop1 = content[loc1:]
                loc2 = string.find(chop1,';')
                chop2 = chop1[13:loc2]
                if "nonuniform" not in chop2:
                    values.append(float(string.split(chop2)[1]))
                else:
                    values.append(0.)
            
        return values
    
    else:
        return None



if __name__ == "__main__":
    times = get_times()
    fields = get_fields(times[0])
    data = []
    
    for i,t in enumerate(times):
        print "Reading time %d of %d" % (i+1,len(times))
        data.append(get_field_data(t,fields))
    
    with open('RChem2S_Data.csv','w') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(fields)
        for row in data:
            writer.writerow(row)
            
            
    
