


# Get the meshBuilder libraryA
import sys
toolPath = '../../../foamTools/python'
sys.path.append(toolPath)
from meshBuilder import *

# Set the mesh base unit to millimeters (conversion factor to meters)
baseUnit = 0.001

# target cell size (in base units, mm)
cellSize = 1

# Definitions:
#  A is the fuel jet side
#  B is the ox jet side

# Dimensions (in base units, mm)
wi = 2.375 # this is half the total width between the jet center points
wo = 15. # this is half the total domain width
DiA = 1.6 # actual size from jacob's paper (1.32 and 1.52)
hiA = 0.35 # half-height of square in inlet jet center
DiB = 1.2 # actual size from jacob's paper (1.32 and 1.52)
hiB = 0.25 # half-height of square in inlet jet center
Do = 3.05 # outer jet diameter
t = 20 # thickness in z dir
h = 40 # height of square block
hCL = 25 # distance up center for skewed bottom block
Li = 10 # inlet jet channel length

# this is the peak angle, the jet angle is 180 - theta
theta = math.pi * 120. / 180.  # radians (60 degrees)

# mesh sizes
nzOut = 19 # number in z direction on outside
ndia = 21 # x and y size of diamand at the top
nDo = 20  # number across Do
nzIn = nDo  # number in z direction in inner block
nSpr = 45 # number in spreading triangle part
nDS = 50  # number in downstream length portion
nOD = 10   # number of radial elements between Di and Do
nID = 6 # number of radial elements in curved portions of smaller cylinder
nIL = 50 # number of elements in inlet jets length

# mesh grading factors
gS = 1./3. # side part grading in z direction, refined towards jets
rgS = 1./gS





rt2 = math.sqrt(2.0)

# Calculate some intermediate position points on "roof" plane
cp = wi/math.sin(theta/2.) # center point of inlet

# length of points along "roof plane"
rLA = [cp - Do/(2.*rt2), \
       cp - DiA/(2.*rt2), \
       cp - hiA, \
       cp + hiA, \
       cp + DiA/(2.*rt2), \
       cp + Do/(2.*rt2), \
       wo/math.sin(theta/2.) ]
       
rLB = [cp - Do/(2.*rt2), \
       cp - DiB/(2.*rt2), \
       cp - hiB, \
       cp + hiB, \
       cp + DiB/(2.*rt2), \
       cp + Do/(2.*rt2), \
       wo/math.sin(theta/2.) ]
       
# global coordinates of points on "roof plane"
rxA = [L * math.sin(theta/2.) for L in rLA]
ryA = [L * math.cos(theta/2.) for L in rLA]

rxB = [-L * math.sin(theta/2.) for L in rLB]
ryB = [L * math.cos(theta/2.) for L in rLB]

# global coordinates of points on "roof plane" offset by Li
rxAT = [L * math.sin(theta/2.) + Li * math.cos(theta/2.) for L in rLA]
ryAT = [L * math.cos(theta/2.) - Li * math.sin(theta/2.) for L in rLA]

rxBT = [-L * math.sin(theta/2.) - Li * math.cos(theta/2.) for L in rLB]
ryBT = [L * math.cos(theta/2.) - Li * math.sin(theta/2.) for L in rLB]

# calculate some heights down center line (average of two sides)
cL = [0.5*(La+Lb) / math.cos(theta/2.) for La,Lb in zip(rLA,rLB)]
cL.pop()
# TODO: Check that hCL here does not violate contiguity of the "hat" portion
cL.append( ryA[-1] + h - hCL )
cL.append( ryA[-1] + h )


# calculate z direction distances
zsA = [-Do/(2.*rt2)-t, -Do/(2.*rt2), -DiA/(2.*rt2), -hiA, hiA, DiA/(2.*rt2), Do/(2.*rt2), Do/(2.*rt2)+t]
zsB = [-Do/(2.*rt2)-t, -Do/(2.*rt2), -DiB/(2.*rt2), -hiB, hiB, DiB/(2.*rt2), Do/(2.*rt2), Do/(2.*rt2)+t]
zsC = [0.5*(zA+zB) for zA,zB in zip(zsA,zsB)]

# Create the individual blocks in the mesh
blocks = []

# ----------------------------------------------------------------  
# mesh base
blocks.append( NonIsoBlock(Point(0,0,zsC[0]),Point(rxA[0],ryA[0],zsA[0]),Point(0,cL[0],zsC[0]),Point(rxB[0],ryB[0],zsB[0]), \
                           Point(0,0,zsC[1]),Point(rxA[0],ryA[0],zsA[1]),Point(0,cL[0],zsC[1]),Point(rxB[0],ryB[0],zsB[1]), \
                           (1,1,1,1,1,1,1,1,gS,gS,gS,gS), size=(ndia,ndia,nzOut) ) )
                           
blocks.append( NonIsoBlock(Point(rxB[0],ryB[0],zsB[0]),Point(0,cL[0],zsC[0]),Point(0,cL[5],zsC[0]),Point(rxB[5],ryB[5],zsB[0]), \
                           Point(rxB[0],ryB[0],zsB[1]),Point(0,cL[0],zsC[1]),Point(0,cL[5],zsC[1]),Point(rxB[5],ryB[5],zsB[1]), \
                           (1,1,1,1,1,1,1,1,gS,gS,gS,gS), size=(ndia,nDo,nzOut) ) )
                           
blocks.append( NonIsoBlock(Point(rxA[0],ryA[0],zsA[0]),Point(rxA[5],ryA[5],zsA[0]),Point(0,cL[5],zsC[0]),Point(0,cL[0],zsC[0]), \
                           Point(rxA[0],ryA[0],zsA[1]),Point(rxA[5],ryA[5],zsA[1]),Point(0,cL[5],zsC[1]),Point(0,cL[0],zsC[1]), \
                           (1,1,1,1,1,1,1,1,gS,gS,gS,gS), size=(nDo,ndia,nzOut) ) )
                           
blocks.append( NonIsoBlock(Point(rxB[5],ryB[5],zsB[0]),Point(0,cL[5],zsC[0]),Point(0,cL[6],zsC[0]),Point(rxB[6],ryB[6],zsB[0]), \
                           Point(rxB[5],ryB[5],zsB[1]),Point(0,cL[5],zsC[1]),Point(0,cL[6],zsC[1]),Point(rxB[6],ryB[6],zsB[1]), \
                           (1,0.5,0.5,1,1,4,4,1,gS,gS,gS,gS), size=(ndia,nSpr,nzOut) ) )
                           
blocks.append( NonIsoBlock(Point(rxA[5],ryA[5],zsA[0]),Point(rxA[6],ryA[6],zsA[0]),Point(0,cL[6],zsC[0]),Point(0,cL[5],zsC[0]), \
                           Point(rxA[5],ryA[5],zsA[1]),Point(rxA[6],ryA[6],zsA[1]),Point(0,cL[6],zsC[1]),Point(0,cL[5],zsC[1]), \
                           (1,4,4,1,1,0.5,0.5,1,gS,gS,gS,gS), size=(nSpr,ndia,nzOut) ) )
           
blocks.append( NonIsoBlock(Point(rxB[6],ryB[6],zsB[0]),Point(0,cL[6],zsC[0]),Point(0,cL[7],zsA[0]),Point(rxB[6],cL[7],zsB[0]), \
                           Point(rxB[6],ryB[6],zsB[1]),Point(0,cL[6],zsC[1]),Point(0,cL[7],zsA[1]),Point(rxB[6],cL[7],zsB[1]), \
                           (0.5,0.5,0.5,0.5,1,1,1,1,gS,gS,gS,gS), size=(ndia,nDS,nzOut) ) )
                           
blocks.append( NonIsoBlock(Point(rxA[6],ryA[6],zsA[0]),Point(rxA[6],cL[7],zsA[0]),Point(0,cL[7],zsC[0]),Point(0,cL[6],zsC[0]), \
                           Point(rxA[6],ryA[6],zsA[1]),Point(rxA[6],cL[7],zsA[1]),Point(0,cL[7],zsC[1]),Point(0,cL[6],zsC[1]), \
                           (1,1,1,1,0.5,0.5,0.5,0.5,gS,gS,gS,gS), size=(nDS,ndia,nzOut) ) )


# ----------------------------------------------------------------                                  
# mesh center
blocks.append( NonIsoBlock(Point(0,0,zsC[1]),Point(rxA[0],ryA[0],zsA[1]),Point(0,cL[0],zsC[1]),Point(rxB[0],ryB[0],zsB[1]), \
                           Point(0,0,zsC[6]),Point(rxA[0],ryA[0],zsA[6]),Point(0,cL[0],zsC[6]),Point(rxB[0],ryB[0],zsB[6]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(ndia,ndia,nzIn) ) )

# cylinder components - outer cylinder
blocks.append( NonIsoBlock(Point(rxB[0],ryB[0],zsB[1]),Point(0,cL[0],zsC[1]),Point(0,cL[5],zsC[1]),Point(rxB[5],ryB[5],zsB[1]), \
                           Point(rxB[1],ryB[1],zsB[2]),Point(0,cL[1],zsC[2]),Point(0,cL[4],zsC[2]),Point(rxB[4],ryB[4],zsB[2]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(ndia,nDo,nOD) ) )
                           
blocks.append( NonIsoBlock(Point(rxA[0],ryA[0],zsA[1]),Point(rxA[5],ryA[5],zsA[1]),Point(0,cL[5],zsC[1]),Point(0,cL[0],zsC[1]), \
                           Point(rxA[1],ryA[1],zsA[2]),Point(rxA[4],ryA[4],zsA[2]),Point(0,cL[4],zsC[2]),Point(0,cL[1],zsC[2]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(nDo,ndia,nOD) ) )

blocks.append( NonIsoBlock(Point(rxB[1],ryB[1],zsB[5]),Point(0,cL[1],zsC[5]),Point(0,cL[4],zsC[5]),Point(rxB[4],ryB[4],zsB[5]), \
                           Point(rxB[0],ryB[0],zsB[6]),Point(0,cL[0],zsC[6]),Point(0,cL[5],zsC[6]),Point(rxB[5],ryB[5],zsB[6]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(ndia,nDo,nOD) ) )
                           
blocks.append( NonIsoBlock(Point(rxA[1],ryA[1],zsA[5]),Point(rxA[4],ryA[4],zsA[5]),Point(0,cL[4],zsC[5]),Point(0,cL[1],zsC[5]), \
                           Point(rxA[0],ryA[0],zsA[6]),Point(rxA[5],ryA[5],zsA[6]),Point(0,cL[5],zsC[6]),Point(0,cL[0],zsC[6]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(nDo,ndia,nOD) ) )

blocks.append( NonIsoBlock(Point(rxB[0],ryB[0],zsB[1]),Point(0,cL[0],zsC[1]),Point(0,cL[1],zsC[2]),Point(rxB[1],ryB[1],zsB[2]), \
                           Point(rxB[0],ryB[0],zsB[6]),Point(0,cL[0],zsC[6]),Point(0,cL[1],zsC[5]),Point(rxB[1],ryB[1],zsB[5]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(ndia,nOD,nDo) ) )
                           
blocks.append( NonIsoBlock(Point(rxB[4],ryB[4],zsB[2]),Point(0,cL[4],zsC[2]),Point(0,cL[5],zsC[1]),Point(rxB[5],ryB[5],zsB[1]), \
                           Point(rxB[4],ryB[4],zsB[5]),Point(0,cL[4],zsC[5]),Point(0,cL[5],zsC[6]),Point(rxB[5],ryB[5],zsB[6]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(ndia,nOD,nDo) ) )
                           
blocks.append( NonIsoBlock(Point(rxA[0],ryA[0],zsA[6]),Point(0,cL[0],zsC[6]),Point(0,cL[1],zsC[5]),Point(rxA[1],ryA[1],zsA[5]), \
                           Point(rxA[0],ryA[0],zsA[1]),Point(0,cL[0],zsC[1]),Point(0,cL[1],zsC[2]),Point(rxA[1],ryA[1],zsA[2]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(ndia,nOD,nDo) ) )
                           
blocks.append( NonIsoBlock(Point(rxA[4],ryA[4],zsA[5]),Point(0,cL[4],zsC[5]),Point(0,cL[5],zsC[6]),Point(rxA[5],ryA[5],zsA[6]), \
                           Point(rxA[4],ryA[4],zsA[2]),Point(0,cL[4],zsC[2]),Point(0,cL[5],zsC[1]),Point(rxA[5],ryA[5],zsA[1]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(ndia,nOD,nDo) ) )

# inner cylinder
blocks.append( NonIsoBlock(Point(rxB[1],ryB[1],zsB[2]),Point(0,cL[1],zsC[2]),Point(0,cL[4],zsC[2]),Point(rxB[4],ryB[4],zsB[2]), \
                           Point(rxB[2],ryB[2],zsB[3]),Point(0,cL[2],zsC[3]),Point(0,cL[3],zsC[3]),Point(rxB[3],ryB[3],zsB[3]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(ndia,nDo,nID) ) )
                           
blocks.append( NonIsoBlock(Point(rxB[2],ryB[2],zsB[4]),Point(0,cL[2],zsC[4]),Point(0,cL[3],zsC[4]),Point(rxB[3],ryB[3],zsB[4]), \
                           Point(rxB[1],ryB[1],zsB[5]),Point(0,cL[1],zsC[5]),Point(0,cL[4],zsC[5]),Point(rxB[4],ryB[4],zsB[5]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(ndia,nDo,nID) ) )
                           
blocks.append( NonIsoBlock(Point(rxB[1],ryB[1],zsB[2]),Point(0,cL[1],zsC[2]),Point(0,cL[2],zsC[3]),Point(rxB[2],ryB[2],zsB[3]), \
                           Point(rxB[1],ryB[1],zsB[5]),Point(0,cL[1],zsC[5]),Point(0,cL[2],zsC[4]),Point(rxB[2],ryB[2],zsB[4]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(ndia,nID,nDo) ) )
                           
blocks.append( NonIsoBlock(Point(rxB[3],ryB[3],zsB[3]),Point(0,cL[3],zsC[3]),Point(0,cL[4],zsC[2]),Point(rxB[4],ryB[4],zsB[2]), \
                           Point(rxB[3],ryB[3],zsB[4]),Point(0,cL[3],zsC[4]),Point(0,cL[4],zsC[5]),Point(rxB[4],ryB[4],zsB[5]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(ndia,nID,nDo) ) )

blocks.append( NonIsoBlock(Point(rxA[1],ryA[1],zsA[2]),Point(rxA[4],ryA[4],zsA[2]),Point(0,cL[4],zsC[2]),Point(0,cL[1],zsC[2]), \
                           Point(rxA[2],ryA[2],zsA[3]),Point(rxA[3],ryA[3],zsA[3]),Point(0,cL[3],zsC[3]),Point(0,cL[2],zsC[3]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(nDo,ndia,nID) ) )
                                             
blocks.append( NonIsoBlock(Point(rxA[2],ryA[2],zsA[4]),Point(rxA[3],ryA[3],zsA[4]),Point(0,cL[3],zsC[4]),Point(0,cL[2],zsC[4]), \
                           Point(rxA[1],ryA[1],zsA[5]),Point(rxA[4],ryA[4],zsA[5]),Point(0,cL[4],zsC[5]),Point(0,cL[1],zsC[5]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(nDo,ndia,nID) ) )
                           
blocks.append( NonIsoBlock(Point(rxA[1],ryA[1],zsA[2]),Point(rxA[2],ryA[2],zsA[3]),Point(0,cL[2],zsC[3]),Point(0,cL[1],zsC[2]), \
                           Point(rxA[1],ryA[1],zsA[5]),Point(rxA[2],ryA[2],zsA[4]),Point(0,cL[2],zsC[4]),Point(0,cL[1],zsC[5]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(nID,ndia,nDo) ) )

blocks.append( NonIsoBlock(Point(rxA[3],ryA[3],zsA[3]),Point(rxA[4],ryA[4],zsA[2]),Point(0,cL[4],zsC[2]),Point(0,cL[3],zsC[3]), \
                           Point(rxA[3],ryA[3],zsA[4]),Point(rxA[4],ryA[4],zsA[5]),Point(0,cL[4],zsC[5]),Point(0,cL[3],zsC[4]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(nID,ndia,nDo) ) )
                
# cores
blocks.append( NonIsoBlock(Point(rxB[2],ryB[2],zsB[3]),Point(0,cL[2],zsC[3]),Point(0,cL[3],zsC[3]),Point(rxB[3],ryB[3],zsB[3]), \
                           Point(rxB[2],ryB[2],zsB[4]),Point(0,cL[2],zsC[4]),Point(0,cL[3],zsC[4]),Point(rxB[3],ryB[3],zsB[4]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(ndia,nDo,nDo) ) )
                           
blocks.append( NonIsoBlock(Point(rxA[2],ryA[2],zsA[3]),Point(rxA[3],ryA[3],zsA[3]),Point(0,cL[3],zsC[3]),Point(0,cL[2],zsC[3]), \
                           Point(rxA[2],ryA[2],zsA[4]),Point(rxA[3],ryA[3],zsA[4]),Point(0,cL[3],zsC[4]),Point(0,cL[2],zsC[4]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(nDo,ndia,nDo) ) )

# inlet channel cylinders
blocks.append( NonIsoBlock(Point(rxBT[1],ryBT[1],zsB[2]),Point(rxB[1],ryB[1],zsB[2]),Point(rxB[4],ryB[4],zsB[2]),Point(rxBT[4],ryBT[4],zsB[2]), \
                           Point(rxBT[2],ryBT[2],zsB[3]),Point(rxB[2],ryB[2],zsB[3]),Point(rxB[3],ryB[3],zsB[3]),Point(rxBT[3],ryBT[3],zsB[3]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(nIL,nDo,nID) ) )
                           
blocks.append( NonIsoBlock(Point(rxBT[2],ryBT[2],zsB[4]),Point(rxB[2],ryB[2],zsB[4]),Point(rxB[3],ryB[3],zsB[4]),Point(rxBT[3],ryBT[3],zsB[4]), \
                           Point(rxBT[1],ryBT[1],zsB[5]),Point(rxB[1],ryB[1],zsB[5]),Point(rxB[4],ryB[4],zsB[5]),Point(rxBT[4],ryBT[4],zsB[5]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(nIL,nDo,nID) ) )
                           
blocks.append( NonIsoBlock(Point(rxBT[1],ryBT[1],zsB[2]),Point(rxB[1],ryB[1],zsB[2]),Point(rxB[2],ryB[2],zsB[3]),Point(rxBT[2],ryBT[2],zsB[3]), \
                           Point(rxBT[1],ryBT[1],zsB[5]),Point(rxB[1],ryB[1],zsB[5]),Point(rxB[2],ryB[2],zsB[4]),Point(rxBT[2],ryBT[2],zsB[4]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(nIL,nID,nDo) ) )
                           
blocks.append( NonIsoBlock(Point(rxBT[3],ryBT[3],zsB[3]),Point(rxB[3],ryB[3],zsB[3]),Point(rxB[4],ryB[4],zsB[2]),Point(rxBT[4],ryBT[4],zsB[2]), \
                           Point(rxBT[3],ryBT[3],zsB[4]),Point(rxB[3],ryB[3],zsB[4]),Point(rxB[4],ryB[4],zsB[5]),Point(rxBT[4],ryBT[4],zsB[5]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(nIL,nID,nDo) ) )
                           
blocks.append( NonIsoBlock(Point(rxBT[2],ryBT[2],zsB[3]),Point(rxB[2],ryB[2],zsB[3]),Point(rxB[3],ryB[3],zsB[3]),Point(rxBT[3],ryBT[3],zsB[3]), \
                           Point(rxBT[2],ryBT[2],zsB[4]),Point(rxB[2],ryB[2],zsB[4]),Point(rxB[3],ryB[3],zsB[4]),Point(rxBT[3],ryBT[3],zsB[4]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(nIL,nDo,nDo) ) )
                           


blocks.append( NonIsoBlock(Point(rxAT[1],ryAT[1],zsA[2]),Point(rxAT[4],ryAT[4],zsA[2]),Point(rxA[4],ryA[4],zsA[2]),Point(rxA[1],ryA[1],zsA[2]), \
                           Point(rxAT[2],ryAT[2],zsA[3]),Point(rxAT[3],ryAT[3],zsA[3]),Point(rxA[3],ryA[3],zsA[3]),Point(rxA[2],ryA[2],zsA[3]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(nDo,nIL,nID) ) )
                                             
blocks.append( NonIsoBlock(Point(rxAT[2],ryAT[2],zsA[4]),Point(rxAT[3],ryAT[3],zsA[4]),Point(rxA[3],ryA[3],zsA[4]),Point(rxA[2],ryA[2],zsA[4]), \
                           Point(rxAT[1],ryAT[1],zsA[5]),Point(rxAT[4],ryAT[4],zsA[5]),Point(rxA[4],ryA[4],zsA[5]),Point(rxA[1],ryA[1],zsA[5]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(nDo,nIL,nID) ) )
                           
blocks.append( NonIsoBlock(Point(rxAT[1],ryAT[1],zsA[2]),Point(rxAT[2],ryAT[2],zsA[3]),Point(rxA[2],ryA[2],zsA[3]),Point(rxA[1],ryA[1],zsA[2]), \
                           Point(rxAT[1],ryAT[1],zsA[5]),Point(rxAT[2],ryAT[2],zsA[4]),Point(rxA[2],ryA[2],zsA[4]),Point(rxA[1],ryA[1],zsA[5]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(nID,nIL,nDo) ) )

blocks.append( NonIsoBlock(Point(rxAT[3],ryAT[3],zsA[3]),Point(rxAT[4],ryAT[4],zsA[2]),Point(rxA[4],ryA[4],zsA[2]),Point(rxA[3],ryA[3],zsA[3]), \
                           Point(rxAT[3],ryAT[3],zsA[4]),Point(rxAT[4],ryAT[4],zsA[5]),Point(rxA[4],ryA[4],zsA[5]),Point(rxA[3],ryA[3],zsA[4]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(nID,nIL,nDo) ) )
                           
blocks.append( NonIsoBlock(Point(rxAT[2],ryAT[2],zsA[3]),Point(rxAT[3],ryAT[3],zsA[3]),Point(rxA[3],ryA[3],zsA[3]),Point(rxA[2],ryA[2],zsA[3]), \
                           Point(rxAT[2],ryAT[2],zsA[4]),Point(rxAT[3],ryAT[3],zsA[4]),Point(rxA[3],ryA[3],zsA[4]),Point(rxA[2],ryA[2],zsA[4]), \
                           (1,1,1,1,1,1,1,1,1,1,1,1), size=(nDo,nIL,nDo) ) )
                           

# body blocks
blocks.append( NonIsoBlock(Point(rxB[5],ryB[5],zsB[1]),Point(0,cL[5],zsC[1]),Point(0,cL[6],zsC[1]),Point(rxB[6],ryB[6],zsB[1]), \
                           Point(rxB[5],ryB[5],zsB[6]),Point(0,cL[5],zsC[6]),Point(0,cL[6],zsC[6]),Point(rxB[6],ryB[6],zsB[6]), \
                           (1,0.5,0.5,1,1,4,4,1,1,1,1,1), size=(ndia,nSpr,nzIn) ) )
                           
blocks.append( NonIsoBlock(Point(rxA[5],ryA[5],zsA[1]),Point(rxA[6],ryA[6],zsA[1]),Point(0,cL[6],zsC[1]),Point(0,cL[5],zsC[1]), \
                           Point(rxA[5],ryA[5],zsA[6]),Point(rxA[6],ryA[6],zsA[6]),Point(0,cL[6],zsC[6]),Point(0,cL[5],zsC[6]), \
                           (1,4,4,1,1,0.5,0.5,1,1,1,1,1), size=(nSpr,ndia,nzIn) ) )
           
           
blocks.append( NonIsoBlock(Point(rxB[6],ryB[6],zsB[1]),Point(0,cL[6],zsC[1]),Point(0,cL[7],zsC[1]),Point(rxB[6],cL[7],zsB[1]), \
                           Point(rxB[6],ryB[6],zsB[6]),Point(0,cL[6],zsC[6]),Point(0,cL[7],zsC[6]),Point(rxB[6],cL[7],zsB[6]), \
                           (0.5,0.5,0.5,0.5,1,1,1,1,1,1,1,1), size=(ndia,nDS,nzIn) ) )
                           
blocks.append( NonIsoBlock(Point(rxA[6],ryA[6],zsA[1]),Point(rxA[6],cL[7],zsA[1]),Point(0,cL[7],zsC[1]),Point(0,cL[6],zsC[1]), \
                           Point(rxA[6],ryA[6],zsA[6]),Point(rxA[6],cL[7],zsA[6]),Point(0,cL[7],zsC[6]),Point(0,cL[6],zsC[6]), \
                           (1,1,1,1,0.5,0.5,0.5,0.5,1,1,1,1), size=(nDS,ndia,nzIn) ) )


# ----------------------------------------------------------------  
# mesh top
blocks.append( NonIsoBlock(Point(0,0,zsC[6]),Point(rxA[0],ryA[0],zsA[6]),Point(0,cL[0],zsC[6]),Point(rxB[0],ryB[0],zsB[6]), \
                           Point(0,0,zsC[7]),Point(rxA[0],ryA[0],zsA[7]),Point(0,cL[0],zsC[7]),Point(rxB[0],ryB[0],zsB[7]), \
                           (1,1,1,1,1,1,1,1,rgS,rgS,rgS,rgS), size=(ndia,ndia,nzOut) ) )
                           
blocks.append( NonIsoBlock(Point(rxB[0],ryB[0],zsB[6]),Point(0,cL[0],zsC[6]),Point(0,cL[5],zsC[6]),Point(rxB[5],ryB[5],zsB[6]), \
                           Point(rxB[0],ryB[0],zsB[7]),Point(0,cL[0],zsC[7]),Point(0,cL[5],zsC[7]),Point(rxB[5],ryB[5],zsB[7]), \
                           (1,1,1,1,1,1,1,1,rgS,rgS,rgS,rgS), size=(ndia,nDo,nzOut) ) )
                           
blocks.append( NonIsoBlock(Point(rxA[0],ryA[0],zsA[6]),Point(rxA[5],ryA[5],zsA[6]),Point(0,cL[5],zsC[6]),Point(0,cL[0],zsC[6]), \
                           Point(rxA[0],ryA[0],zsA[7]),Point(rxA[5],ryA[5],zsA[7]),Point(0,cL[5],zsC[7]),Point(0,cL[0],zsC[7]), \
                           (1,1,1,1,1,1,1,1,rgS,rgS,rgS,rgS), size=(nDo,ndia,nzOut) ) )
                           
                           
blocks.append( NonIsoBlock(Point(rxB[5],ryB[5],zsB[6]),Point(0,cL[5],zsC[6]),Point(0,cL[6],zsC[6]),Point(rxB[6],ryB[6],zsB[6]), \
                           Point(rxB[5],ryB[5],zsB[7]),Point(0,cL[5],zsC[7]),Point(0,cL[6],zsC[7]),Point(rxB[6],ryB[6],zsB[7]), \
                           (1,0.5,0.5,1,1,4,4,1,rgS,rgS,rgS,rgS), size=(ndia,nSpr,nzOut) ) )
                           
blocks.append( NonIsoBlock(Point(rxA[5],ryA[5],zsA[6]),Point(rxA[6],ryA[6],zsA[6]),Point(0,cL[6],zsC[6]),Point(0,cL[5],zsC[6]), \
                           Point(rxA[5],ryA[5],zsA[7]),Point(rxA[6],ryA[6],zsA[7]),Point(0,cL[6],zsC[7]),Point(0,cL[5],zsC[7]), \
                           (1,4,4,1,1,0.5,0.5,1,rgS,rgS,rgS,rgS), size=(nSpr,ndia,nzOut) ) )
           
           
blocks.append( NonIsoBlock(Point(rxB[6],ryB[6],zsB[6]),Point(0,cL[6],zsC[6]),Point(0,cL[7],zsC[6]),Point(rxB[6],cL[7],zsB[6]), \
                           Point(rxB[6],ryB[6],zsB[7]),Point(0,cL[6],zsC[7]),Point(0,cL[7],zsC[7]),Point(rxB[6],cL[7],zsB[7]), \
                           (0.5,0.5,0.5,0.5,1,1,1,1,rgS,rgS,rgS,rgS), size=(ndia,nDS,nzOut) ) )
                           
blocks.append( NonIsoBlock(Point(rxA[6],ryA[6],zsA[6]),Point(rxA[6],cL[7],zsA[6]),Point(0,cL[7],zsC[6]),Point(0,cL[6],zsC[6]), \
                           Point(rxA[6],ryA[6],zsA[7]),Point(rxA[6],cL[7],zsA[7]),Point(0,cL[7],zsC[7]),Point(0,cL[6],zsC[7]), \
                           (1,1,1,1,0.5,0.5,0.5,0.5,rgS,rgS,rgS,rgS), size=(nDS,ndia,nzOut) ) )
                           


mesh = BlockMesh(blocks)

# tag the inlet and outlet faces based on a plane point and normal (all untagged faces are 'walls')
#mesh.tagFaces((0,0,0),(0,1,0),"outlet")
mesh.tagFaces((0,cL[7],0),(0,1,0),"atmosphere")
mesh.tagFaces((rxA[6],0,0),(1,0,0),"atmosphere")
mesh.tagFaces((rxB[6],0,0),(1,0,0),"atmosphere")
mesh.tagFaces((0,0,zsA[0]),(0,0,1),"atmosphere")
mesh.tagFaces((0,0,zsA[7]),(0,0,1),"atmosphere")
mesh.tagFaces((rxA[3],ryA[3],0),(math.cos(theta/2.),-math.sin(theta/2.),0),"atmosphere")
mesh.tagFaces((rxB[3],ryB[3],0),(-math.cos(theta/2.),-math.sin(theta/2.),0),"atmosphere")

# get injector faces as "walls"
mesh.tagFaces([Point(rxA[0],ryA[0],zsA[1]), Point(rxA[1],ryA[1],zsA[2]), Point(rxA[4],ryA[4],zsA[2])], None, "walls")
mesh.tagFaces([Point(rxA[0],ryA[0],zsA[1]), Point(rxA[0],ryA[0],zsA[6]), Point(rxA[1],ryA[1],zsA[2])], None, "walls")
mesh.tagFaces([Point(rxA[5],ryA[5],zsA[6]), Point(rxA[4],ryA[4],zsA[5]), Point(rxA[4],ryA[4],zsA[2])], None, "walls")
mesh.tagFaces([Point(rxA[5],ryA[5],zsA[6]), Point(rxA[0],ryA[0],zsA[6]), Point(rxA[1],ryA[1],zsA[5])], None, "walls")

mesh.tagFaces([Point(rxB[0],ryB[0],zsB[1]), Point(rxB[1],ryB[1],zsB[2]), Point(rxB[4],ryB[4],zsB[2])], None, "walls")
mesh.tagFaces([Point(rxB[0],ryB[0],zsB[1]), Point(rxB[0],ryB[0],zsB[6]), Point(rxB[1],ryB[1],zsB[2])], None, "walls")
mesh.tagFaces([Point(rxB[5],ryB[5],zsB[6]), Point(rxB[4],ryB[4],zsB[5]), Point(rxB[4],ryB[4],zsB[2])], None, "walls")
mesh.tagFaces([Point(rxB[5],ryB[5],zsB[6]), Point(rxB[0],ryB[0],zsB[6]), Point(rxB[1],ryB[1],zsB[5])], None, "walls")


mesh.tagFaces((rxAT[3],ryAT[3],0),(math.cos(theta/2.),-math.sin(theta/2.),0),"FuelInlet")
mesh.tagFaces((rxBT[3],ryBT[3],0),(-math.cos(theta/2.),-math.sin(theta/2.),0),"OxInlet")


# Make point pairs on a given cylinder have arcs
mesh.makeArcs(base=(cp*math.sin(theta/2.),cp*math.cos(theta/2.),0), \
              radius=Do/2., \
              direction=(math.cos(theta/2.),-math.sin(theta/2.),0))
              
mesh.makeArcs(base=(-cp*math.sin(theta/2.),cp*math.cos(theta/2.),0), \
              radius=Do/2., \
              direction=(-math.cos(theta/2.),-math.sin(theta/2.),0))
              
mesh.makeArcs(base=(cp*math.sin(theta/2.),cp*math.cos(theta/2.),0), \
              radius=DiA/2., \
              direction=(math.cos(theta/2.),-math.sin(theta/2.),0))
              
mesh.makeArcs(base=(-cp*math.sin(theta/2.),cp*math.cos(theta/2.),0), \
              radius=DiB/2., \
              direction=(-math.cos(theta/2.),-math.sin(theta/2.),0))
              
              
mesh.removeArcsOnPlane(point=(0,0,0), normal=(1,0,0))

# set the mesh scale in cells/arb unit
mesh.scale = 1. / cellSize

# Set the mesh base unit
mesh.baseUnit = baseUnit

f = open('constant/polyMesh/blockMeshDict','w')
f.write(str(mesh))
f.close()


