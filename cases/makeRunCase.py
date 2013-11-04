
import argparse
import shutil
import os

parser = argparse.ArgumentParser(description='Move case to run folder')

# Required input
parser.add_argument('caseName', help='Name of the case to copy')

# Optional inputs
parser.add_argument('-f', dest='force', action='store_true', 
       help='Force case creation - WARNING - ' \
        +'This will overwrite existing case and results!!!')
        
parser.add_argument('runName', nargs='?', action='store',
                    default=None, help='Name of run folder')

# Get inputs
args = parser.parse_args()
force = args.force
caseName = args.caseName

if args.runName is not None:
    runName = args.runName
else:
    runName = caseName
    
# Check if destination folder exists
runPath = os.path.join(os.getcwd(),'../run')
newPath = os.path.join(runPath,runName)
basePath = os.path.join(os.getcwd(),caseName)

if not os.path.isdir(basePath):
    raise IOError('Specified case not found')

if not os.path.isdir(runPath):
    os.mkdir(runPath)
    
if os.path.isdir(newPath):
    if not force:
        raise IOError('Destination case folder already exists')
    else:
        shutil.rmtree(newPath)
        
# Copy case to run folder
shutil.copytree(basePath, newPath)




