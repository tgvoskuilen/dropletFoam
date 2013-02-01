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

#ifndef cuttableCell_H
#define cuttableCell_H

#include "fvCFD.H"
#include "plane.H"
#include "SortableList.H"

namespace Foam
{

//- Create a cell that can be cut by a plane and have the cut volume and area
//  returned.
class cuttableCell 
{
    // Private attributes

        //- Points that make the cell
        pointField points_;

        //- Faces that make the cell, referenced to points_
        faceList faces_;
        
        //- Tolerance for point equality (to catch degenerate cuts)
        scalar pointEqualTol_;

        //- Volume of the uncut cell
        scalar baseVol_;

        //- Centroid of the uncut cell
        point centroid_;
        
        //- Centroid of the cut cell
        point cutCentroid_;
        point lostCentroid_;

        //- Area of the last cut done on the cell
        scalar cutArea_;
        
    // Private methods

        //- Cut a face with a plane and return the cut face plane components.
        //   Also adds any cut points to the cutPoints list
        bool cutFace
        (
            const plane& cutPlane,
            const label faceI,
            point& facePoint,
            vector& faceNormal,
            DynamicList<point>& cutPoints,
            const labelList& vertStates
        ) const;
        
        //- Given a list of co-planar sorted points, 
        //    find the centroid and area vector
        bool makePoints
        (
            DynamicList<point>& pointList,
            point& centroid,
            vector& areaVec
        ) const;

        //- Given the list of cut points, remove duplicates and sort in order
        //   around the cut face
        void reduceCutPoints
        (
            DynamicList<point>& cutPoints
        ) const;
    
    
    
    public:
    
    // Constructors

        //- Construct from a cell in the mesh
        cuttableCell(const fvMesh& mesh, label cellI);

        //- Construct by extruding a face in the mesh in direction d
        cuttableCell(const fvMesh& mesh, label faceI, vector d);
        
    // Methods

        //- Cut the cell by the plane p and return the cut volume. The cut
        //    volume is always the volume opposite the plane normal direction
        scalar cut(const plane& p);

        //- Find the plane with normal n that cuts the cell with a volume
        //    fraction alpha
        plane constructInterface(const vector& n, scalar alpha);
        
        //- Return the area of the last cut done on the cuttableCell
        inline scalar cutArea() const { return cutArea_; }
        
        inline point cutCentroid() const { return cutCentroid_; }
        inline point lostCentroid() const { return lostCentroid_; }
};

};

#endif


