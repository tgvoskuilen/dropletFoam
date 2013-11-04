from paraview.simple import *

#import os
#import csv

source = FindSource('hypergol.foam')

times = source.TimestepValues

PlotOverLine1 = PlotOverLine( Source="High Resolution Line Source" )

PlotOverLine1.Source.Point1 = [0, 0.0, 0.003]
PlotOverLine1.Source.Point2 = [0, 0.01, 0.003]

writer = CreateWriter('/scratch/chaf124lnx01/a/tvoskuil/OpenFOAM/tvoskuil-2.1.x/repos/dropletFoam/cases/vialTest2/LineOutput.csv', PlotOverLine1)
writer.FieldAssociation = "Points"
writer.WriteAllTimeSteps = True
writer.UpdatePipeline()
del writer
        
        
"""
view = GetActiveView()

for i,time in enumerate(times):
    if (i > 3 and i < 8) or i == 300:
        view.ViewTime = time
        render=Render()
        filename = "time_%s_ms.csv" % str(time)
        
        print "Exporting to %s" % filename

        writer = CreateWriter(filename, PlotOverLine1)
        writer.FieldAssociation = "Points"
        writer.WriteAllTimeSteps = True
        writer.UpdatePipeline()
        del writer
"""
