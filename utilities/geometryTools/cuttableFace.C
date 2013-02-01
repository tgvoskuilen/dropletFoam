/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::cuttableFace

Description


SourceFiles
    cuttableFace.C

\*---------------------------------------------------------------------------*/

#include "cuttableFace.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cuttableFace::cuttableFace(const fvMesh& mesh, label faceI)
:
    mesh_(mesh),
    faceID_(faceI)
{}


// * * * * * * * * * * * * * * * * Public Methods  * * * * * * * * * * * * * //
Foam::scalar Foam::cuttableFace::cut(const plane& p) const
{
    const face& f = mesh_.faces()[faceID_];
    const pointField& points = mesh_.points();
    
    scalar faceArea = Foam::mag(f.normal(points));

    //scalar tol = Foam::sqrt(faceArea) * 1e-3;

    // Classify each point of the face as either
    // -1 (kept side), 0 (coplanar), or 1 (opposite)
    // Points with a -1 are opposite the plane normal (in the solid)
    Foam::labelList pointState(f.size()); 

    forAll(f, pointI)
    {
        Foam::vector vp = points[f[pointI]] - p.refPoint();
        bool inFront = (vp & p.normal()) > 0.0;
        bool coplanar = p.distance(points[f[pointI]]) <= SMALL;
        
        pointState[pointI] = (coplanar) ? 0 : ((inFront) ? 1 : -1);
    }
    
    //exit here if face is not cut
    if (min(pointState) > -1) //has only 0's and 1's
    { //face is completely on the gas side
        return 1.0; 
    }
    else if (max(pointState) < 1) //has only -1's and 0's
    { //face is completely solid
        return 0.0; 
    }
    
    // If we're still here, it's a cut face. Find the new face's points.
    Foam::DynamicList<Foam::point> newFacePoints(f.size());
    forAll(f, pointI)
    {
        // vs,ve are global vertex indices
        label vs = f[pointI];
        label ve = f.nextLabel(pointI);
        
        // ps,pe are local (to this face) vertex indices
        label ps = pointI;
        label pe = (pointI == f.size()-1) ? 0 : pointI+1;
        
        if (pointState[ps] * pointState[pe] == -1) //edge is cut
        {
            //find intersection point (pI)
            Foam::edge e(vs,ve);
            scalar s = p.normalIntersect( points[vs], e.vec(points) );
            Foam::point pI = points[vs] + s*e.vec(points);

            if (pointState[ps] == 1)
            {
                newFacePoints.append(points[vs]);
                newFacePoints.append(pI);
            }
            else
            {
                newFacePoints.append(pI);
            }
        }
        else if (pointState[ps] >= 0) //edge starts on the kept side
        {
            newFacePoints.append(points[vs]);
        }
    }

    point center = sum(newFacePoints)/newFacePoints.size();
    
    //Find cut face area   
    scalar area = 0.0;
    for(label i=0; i<newFacePoints.size(); ++i)
    {
        label j = (i==newFacePoints.size()-1) ? 0 : i+1;
        vector vCi = newFacePoints[i] - center;
        vector vCj = newFacePoints[j] - center;
        area += 0.5*mag(vCi ^ vCj);
    }
    
    return area/faceArea;
}



