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
    Foam::cuttableCell

Description


SourceFiles
    cuttableCell.C

\*---------------------------------------------------------------------------*/

#include "cuttableCell.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cuttableCell::cuttableCell(const Foam::fvMesh& mesh, label cellI)
 : points_(mesh.cellPoints()[cellI].size()),
   faces_(mesh.cells()[cellI].size()),
   pointEqualTol_(1e-5),
   baseVol_(mesh.V()[cellI]),
   centroid_(mesh.C()[cellI]),
   cutArea_(0.0)
{
    Foam::labelList pointMap(points_.size());
    
    const Foam::pointField& meshPoints = mesh.points();
    const Foam::labelList& cellPoints = mesh.cellPoints()[cellI];
    
    // Gather points
    forAll(cellPoints, pointI)
    {
        points_[pointI] = meshPoints[cellPoints[pointI]];
        pointMap[pointI] = cellPoints[pointI];
    }

    cutCentroid_ = centroid_;
    lostCentroid_ = vector::zero;
    
    // Map faces
    const Foam::cell& faces = mesh.cells()[cellI];
    forAll(faces, faceI)
    {
        
        const Foam::labelList& facePoints = mesh.faces()[faces[faceI]];
        Foam::face mappedFace(facePoints.size());
        
        forAll(facePoints, pointI)
        {
            label pOld = facePoints[pointI];
            label pNew = -1;
            
            forAll(pointMap, p)
            {
                if (pointMap[p] == pOld)
                {
                    pNew = p;
                    break;
                }
            }
            
            if (pNew == -1)
            {
                Info<<"Could not find point " << pOld
                    << " in " << pointMap << endl;
            }
            
            mappedFace[pointI] = pNew;
        }
        faces_[faceI] = mappedFace;
    }
}

// For the version with a vectorField for d, the only things that will change
//  are the calculation of baseVol_ and the addition of d[pointI] rather than
//  d in the point gathering loop. A check to make sure d.size() == f.size()
//  would also be good.
Foam::cuttableCell::cuttableCell(const Foam::fvMesh& mesh, label faceI, Foam::vector d)
 : points_(mesh.faces()[faceI].size()*2),
   faces_(mesh.faces()[faceI].size()+2),
   pointEqualTol_(1e-5),
   baseVol_(Foam::mag(mesh.faces()[faceI].normal(mesh.points()) & d)),
   centroid_(Foam::vector::zero),
   cutArea_(0.0)
{
    Foam::labelList pointMap(points_.size());
    
    const Foam::pointField& meshPoints = mesh.points();
    const Foam::face& f = mesh.faces()[faceI];
        
    // base and top faces are the same size
    Foam::face base(f.size());
    Foam::face top(f.size());

    //Gather all points and make faces base (original) and top (extruded)
    forAll(f, pointI)
    {
        points_[pointI] = meshPoints[f[pointI]];
        points_[pointI + f.size()] = meshPoints[f[pointI]] + d;
        pointMap[pointI] = f[pointI];
        base[pointI] = pointI;
        top[pointI] = pointI + f.size();
    }
    centroid_ = f.centre(meshPoints) + 0.5*d;

    faces_[0] = base;
    faces_[1] = top;

    //Now loop around the base and make the side faces
    forAll(base, pointI)
    {
        Foam::face side(4);
        label ps = base[pointI];
        label pe = base.nextLabel(pointI);
        
        side[0] = ps;
        side[1] = pe;
        side[2] = pe + base.size();
        side[3] = ps + base.size();
        
        faces_[pointI+2] = side;
    }
}


// * * * * * * * * * * * * * * * * Methods * * * * * * * * * * * * * * * * * //

Foam::plane Foam::cuttableCell::constructInterface
(
    const Foam::vector& pn,
    Foam::scalar alpha
)
{
    // plane refPoint = centroid + d*n, want to find scalar 'd'
    Foam::vector n = pn / mag(pn); //make sure plane normal is normalized
    
    //find dmin and dmax
    scalar dmin = GREAT;
    scalar dmax = -GREAT;
    
    forAll(points_, p)
    {
        scalar dp = (n & (points_[p] - centroid_));
        dmin = (dp < dmin) ? dp : dmin;
        dmax = (dp > dmax) ? dp : dmax;
    }
    
    //iterate using the secant method
    scalar dL = dmin;
    scalar dH = dmax;
    scalar dM = dL;
    //scalar fL = -alpha;
    //scalar fH = 1.0-alpha;
    scalar res = 1.0; //mag(f2);
    label iters = 0;
    scalar resmin = 1e-4;
    
    scalar dspan = dmax - dmin;
    scalar dtol = dspan / 1e4;
    

    if (Foam::mag(alpha) <= resmin)
    {
        dM = dL;
        Foam::plane pl(centroid_ + n*(dM+dtol),n);
        cut(pl); //we must still call cut so the cut area is calculated
    }
    else if (Foam::mag(alpha-1.0) <= resmin)
    {
        dM = dH;
        Foam::plane pl(centroid_ + n*(dM-dtol),n);
        cut(pl); //we must still call cut so the cut area is calculated
    }
    else
    {
        DynamicList<scalar> ds(10);
        DynamicList<scalar> fs(10);
        
        ds.append(dL);
        ds.append(dH);
        fs.append(-alpha);
        fs.append(1.0-alpha);
        
        while (res > resmin)
        {
            // Bisection method
            dM = 0.5*(dL + dH);
            Foam::plane pl(centroid_ + n*dM,n);
            scalar fM = cut(pl) - alpha;
            
            res = Foam::mag(fM);
            
            ds.append(dM);
            fs.append(fM);

            dL = (fM < 0.0) ? dM : dL;
            dH = (fM < 0.0) ? dH : dM;
            
            

            if (Foam::mag(dL - dH) < SMALL)
            {
                Info<< "\nWARNING: Bisection method failed" << endl;
                Info<< "  ds = " << ds << endl;
                Info<< "  fs = " << fs << endl;
                Info<< "  n = " << n << endl;
                Info<< "  alpha = " << alpha << endl;
                Info<< "  cell = " << points_ << endl;
                Info<< "  centroid = " << centroid_ << endl;
                FatalError << "bisection method failed " << abort(FatalError);
                break;
            }

    /*
            // Secant method
            scalar d3 = 0.0;

            // Deal with inflection points (highly unlikely..., but possible)
            if (mag(f2 - f1) < SMALL)
            {
                d3 = 0.9*d2 + 0.1*d1;
            }
            else
            {
                d3 = d2 - f2*(d2-d1)/(f2-f1);

                //keep d3 from going out of bounds
                d3 = (d3 > dmax) ? dmax : d3;
                d3 = (d3 < dmin) ? dmin : d3;
            }

            Foam::plane pl(centroid_ + n*d3,n);
            scalar f3 = cut(pl) - alpha;

            f1 = f2;
            f2 = f3;
            d1 = d2;
            d2 = d3;

            res = Foam::mag(f2);
    */

            if (++iters > 1000)
            {
                Info<< "\nWARNING: Interface construction has stalled "
                    << "with alpha = " << alpha
                    << " normal = " << n 
                    << " and res = " << res << endl << endl;
                break;
            }
        }
    }

    return Foam::plane(centroid_ + n*dM, n);
}

void Foam::cuttableCell::reduceCutPoints
(
    Foam::DynamicList<Foam::point>& cutPoints
) const
{
    Foam::DynamicList<Foam::point> reducedPoints(cutPoints.size());
    
    // eliminate duplicates first
    reducedPoints.append( cutPoints[0] );
    
    for(label p=1; p<cutPoints.size(); ++p)
    {
        bool inList = false;
        for(label j=0; j<reducedPoints.size(); ++j)
        {
            if (mag(reducedPoints[j] - cutPoints[p]) < SMALL)
            {
                inList = true;
                break;
            }
        }
        
        if (!inList)
        {
            reducedPoints.append( cutPoints[p] );
        }
    }
    
    // Then sort points around cut face
    if (reducedPoints.size() > 3)
    {
        point origin = sum(reducedPoints) / reducedPoints.size();
                
        Foam::SortableList<scalar> angles( reducedPoints.size() );
        
        angles[0] = 1.0; //measure from this vertex

        Foam::vector vC0 = reducedPoints[0] - origin;
        
        Foam::vector n = (vC0 ^ (reducedPoints[1] - origin));
        
        if (mag(n) < SMALL)
        {
            n = (vC0 ^ (reducedPoints[2] - origin));
        }
        
        for(label i=1; i<reducedPoints.size(); ++i)
        {
            Foam::vector vCN = reducedPoints[i] - origin;
            angles[i] = (vCN & vC0) / (mag(vCN)*mag(vC0)); // = cos(theta)
            
            if ((n & (vC0 ^ vCN)) < 0.0)
            {
                angles[i] = -angles[i]-2.0;
            }
        }
        
        angles.sort();
        Foam::labelList idx = angles.indices();
        cutPoints.resize(reducedPoints.size());
        forAll(idx, i)
        {
            cutPoints[i] = reducedPoints[idx[i]];
        }
    }
    else
    {
        cutPoints = reducedPoints;
    }
}

bool Foam::cuttableCell::makePoints
(
    Foam::DynamicList<Foam::point>& pointList,
    Foam::point& centroid,
    Foam::vector& areaVec
) const
{
    if (pointList.size() == 0)
    {
        Info<< "ERROR: makePoints obtained an empty pointList" << endl;
        return false;
    }

    //Find any old point inside the polygon
    point origin = sum(pointList) / pointList.size();
    
    
    centroid = Foam::vector::zero;
    scalar area = 0.0;
    Foam::vector vCi = pointList[0] - origin;
    Foam::vector vCj = pointList[1] - origin;
    areaVec = (vCi ^ vCj);
    
    //if point 1 happens to be diagonally opposite to point 0, the area vector
    // will be zero. Just pick the next point then.
    if (mag(areaVec) < VSMALL)
    {
        vCj = pointList[2] - origin;
        areaVec = (vCi ^ vCj);
    }
    
    //Find area and centroid
    for (label i=0; i<pointList.size(); i++)
    {
        label j = (i == pointList.size()-1) ? 0 : i+1;
        vCi = pointList[i] - origin;
        vCj = pointList[j] - origin;
        scalar tarea = 0.5*mag(vCi ^ vCj);
        area += tarea;
        centroid += (1.0/3.0)*(3.0*origin + vCi + vCj) * tarea;
    }
    
    if (mag(area) < VSMALL || mag(areaVec) < VSMALL)
    {
        Info<< "\n--> Bad point list:\n" << pointList << endl;
        Info<< "\n--> Area:" << mag(area) << endl;
        //FatalError << "makePoints encountered zero area" << abort(FatalError);
        return false;
    }
    
    centroid /= area;
    areaVec /= mag(areaVec);
    areaVec *= area;
    
    return true;
}

// Cuts a face and returns the centroid and normal of the resulting cut face
// on the back side of the cutting plane (opposite the plane normal).
// Returns true if the face has any parts on the kept side, and false if the
// face is totally removed. Any cut points are added to the cutPoints field
bool Foam::cuttableCell::cutFace
(
    const Foam::plane& cutPlane,
    const label faceI,
    Foam::point& facePoint,
    Foam::vector& faceArea,
    Foam::DynamicList<Foam::point>& cutPoints,
    const labelList& vertStates
) const
{
    const face& f = faces_[faceI];

    facePoint = f.centre(points_);
    faceArea = f.normal(points_);

    //-1 (kept side), 0 (coplanar), 1 (opposite)
    Foam::labelList pointState(f.size()); 

    forAll(f, pointI)
    {
        pointState[pointI] = vertStates[f[pointI]];
    }
    
    //exit here if face is not cut
    if (min(pointState) > -1) //has only 0's and 1's
    { //face is completely lost
        return false; 
    }
    else if (max(pointState) < 1) //has only -1's and 0's
    { //face is kept in tact, any coplanar points are considered "cut"
    
        forAll(pointState, p)
        {
            if (pointState[p] == 0)
            {
                cutPoints.append( points_[ f[p] ] );
            }
        }
        
        if (min(pointState) == 0) //fully co-planar cut
        { //if the face is fully co-planar, don't return both the face, and
          // the cut points or it is counted twice
            return false;
        }
        else
        {
            return true;
        }
    }
    
    // If we're still here, it's a cut face. Find the new face's points.
    Foam::DynamicList<Foam::point> newFacePoints(f.size());
    forAll(f, pointI)
    {
        label vs = f[pointI];
        label ve = f.nextLabel(pointI);
        label ps = pointI;
        label pe = (pointI == f.size()-1) ? 0 : pointI+1;
        
        if (pointState[ps] * pointState[pe] == -1) //edge is cut
        {
            //find intersection point
            Foam::edge e(vs,ve);
            scalar s = cutPlane.normalIntersect( points_[vs], e.vec(points_) );
            Foam::point pI = points_[vs] + s*e.vec(points_);

            if (pointState[ps] == -1)
            {
                newFacePoints.append(points_[vs]);
                newFacePoints.append(pI);
            }
            else
            {
                newFacePoints.append(pI);
            }
            cutPoints.append(pI);
        }
        else if (pointState[ps] == 0) //edge starts on the cutting plane
        {
            cutPoints.append(points_[vs]);
            newFacePoints.append(points_[vs]);
        }
        else if (pointState[ps] == -1) //edge is entirely kept
        {
            newFacePoints.append(points_[vs]);
        }
    }
    //Now figure out face centroid and area vector
    // points should still be in order in newFacePoints, so no need to sort
    if (!makePoints(newFacePoints, facePoint, faceArea))
    {
        Info<<"\n\n--> Points: " << points_ << endl;
        Info<<"\n\n--> Point state: " << vertStates << endl;
        Info<<"\n\n--> Cut points: " << cutPoints << endl;
        Info<<"\n\n--> Face: " << f << endl;
        Info<<"\n\n--> Face Center: " << f.centre(points_) << endl;
        Info<<"\n\n--> Face Area: " << f.normal(points_) << endl;
        Info<<"\n\n--> Plane: " << cutPlane << endl;
        FatalError << "makePoints encountered zero area" << abort(FatalError);
    }

    return true;
}


scalar Foam::cuttableCell::cut(const Foam::plane& cutPlane)
{    
    Foam::DynamicList<Foam::point> cutPoints(faces_.size());
    Foam::DynamicList<Foam::point> facePoints(faces_.size());
    Foam::DynamicList<Foam::vector> faceAreas(faces_.size());
    
    labelList vertStates(points_.size());
    scalar tol = Foam::pow(baseVol_, 1.0/3.0) * pointEqualTol_;
    forAll(points_, p)
    {
        Foam::vector vp = points_[p] - cutPlane.refPoint();
        bool inFront = (vp & cutPlane.normal()) > 0.0;
        bool coplanar = cutPlane.distance(points_[p]) <= tol;
        
        vertStates[p] = (coplanar) ? 0 : ((inFront) ? 1 : -1);
    }
    
    forAll(faces_, faceI)
    {
        Foam::point facePoint = Foam::vector::zero;
        Foam::vector faceArea = Foam::vector::zero;
        
        if (cutFace(cutPlane,faceI,facePoint,faceArea,cutPoints,vertStates))
        {
            facePoints.append( facePoint );
            faceAreas.append( faceArea );
        }
    }

    cutArea_ = 0.0; //assume no cut or degenerate cut

    if (cutPoints.size() > 0)
    {
        reduceCutPoints( cutPoints );
    }
    
    if (cutPoints.size() < 3) //degenerate edge or point cut
    {
        if (facePoints.size() > 0)
        {
            return 1.0;
        }
        else
        {
            return 0.0;
        }
    }
    
    //Get cut face centroid and area vector
    Foam::point cutPoint = Foam::vector::zero;
    Foam::vector cutNormal = Foam::vector::zero;
    makePoints(cutPoints, cutPoint, cutNormal);
    cutArea_ = mag(cutNormal);

    //Assemble cut portion centroid
    cutCentroid_ = vector::zero;
    scalar volume = 0.0;
    
    // Pick a point inside the cut polyhedron
    point origin = (sum(facePoints)+cutPoint)/(facePoints.size()+1);
    
    //Calculate volume by adding pyramid volumes
    scalar pvolume = Foam::mag((1.0/3.0)*(cutNormal & (origin - cutPoint)));
    cutCentroid_ += (0.25*origin + 0.75*cutPoint)*pvolume;
    volume += pvolume;
    
    forAll(facePoints, fp)
    {
        scalar pvolume = Foam::mag
        (
            (1.0/3.0)*(faceAreas[fp] & (origin - facePoints[fp]))
        );
        cutCentroid_ += (0.25*origin + 0.75*facePoints[fp])*pvolume;
        volume += pvolume;
    }
    cutCentroid_ /= volume;
    
    scalar alpha = volume/baseVol_;
    
    lostCentroid_ = centroid_ + 
                    (alpha - 1.0)/(alpha + SMALL)*(cutCentroid_ - centroid_);
    
    return alpha;
}



