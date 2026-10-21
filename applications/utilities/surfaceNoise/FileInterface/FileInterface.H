/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::FileInterface

Description
    Reader of surface pressure time history and writer of complex pressure
    amplitudes for the surfaceNoise utility.

    Input file is a text table written by surfacePressureSampling, one line
    per time step:
    \verbatim
    t   p_1  p_2  ...  p_N
    \endverbatim

    Blank lines are skipped. A short last line (interrupted write) is
    dropped with a warning, any other line of wrong length is an error.
    Times must be strictly ascending.

    Output is a gmsh $NodeData block, values are written as Re+Imj.

SourceFiles
    FileInterface.C

\*---------------------------------------------------------------------------*/

#ifndef FileInterface_H
#define FileInterface_H

#include "fileName.H"
#include "scalarList.H"
#include "complex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class FileInterface Declaration
\*---------------------------------------------------------------------------*/

class FileInterface
{
    // Private Data

        //- Name of the file to read or write
        fileName name_;


public:

    // Constructors

        //- Construct from file name
        explicit FileInterface(const fileName& name);


    //- Destructor
    ~FileInterface() = default;


    // Member Functions

        //- Read pressure time history.
        //  Pressure is stored per point: pressure[pointI][timeI]
        void read(scalarList& times, List<scalarList>& pressure) const;

        //- Write complex amplitudes at the given frequency
        void write(const scalar frequency, const UList<complex>& data) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
