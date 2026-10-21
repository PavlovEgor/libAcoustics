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
    Foam::Interpolator

Description
    Abstract base class for interpolation of a time signal from ascending
    sampled times to the uniform grid t0 + j*dt, j = 0..n-1.

    The sampled times and the grid are the same for all signals, so
    everything that depends on them only is computed once on construction.
    Interval search is done here, schemes add their own coefficients.

    Points of the grid outside the sampled range take the end values.

SourceFiles
    Interpolator.C

\*---------------------------------------------------------------------------*/

#ifndef Interpolator_H
#define Interpolator_H

#include "scalarList.H"
#include "labelList.H"
#include "autoPtr.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class Interpolator Declaration
\*---------------------------------------------------------------------------*/

class Interpolator
{
protected:

    // Protected Data

        //- Lower sampled index of the interval for each grid point
        labelList lower_;

        //- Position of each grid point in its interval, 0..1
        scalarList u_;


public:

    //- Runtime type information
    TypeName("Interpolator");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            Interpolator,
            times,
            (
                const UList<scalar>& tSampled,
                const scalar t0,
                const scalar dt,
                const label n
            ),
            (tSampled, t0, dt, n)
        );


    // Constructors

        //- Construct from sampled times and the uniform grid
        Interpolator
        (
            const UList<scalar>& tSampled,
            const scalar t0,
            const scalar dt,
            const label n
        );

        //- No copy construct
        Interpolator(const Interpolator&) = delete;

        //- No copy assignment
        void operator=(const Interpolator&) = delete;


    // Selectors

        //- Return a reference to the selected interpolation scheme
        static autoPtr<Interpolator> New
        (
            const word& type,
            const UList<scalar>& tSampled,
            const scalar t0,
            const scalar dt,
            const label n
        );


    //- Destructor
    virtual ~Interpolator() = default;


    // Member Functions

        //- Number of points of the uniform grid
        label size() const noexcept
        {
            return lower_.size();
        }

        //- Interpolate sampled values to the uniform grid,
        //- result has size() elements
        virtual void interpolate
        (
            const UList<scalar>& values,
            UList<scalar>& result
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
