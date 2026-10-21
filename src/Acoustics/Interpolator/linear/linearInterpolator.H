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
    Foam::linearInterpolator

Description
    Piecewise linear interpolation.

    Attenuates amplitude of the frequency f approximately by
    sinc^2(f*dt) of the sampled step dt: -0.3 dB at f*dt = 0.1.

SourceFiles
    linearInterpolator.C

\*---------------------------------------------------------------------------*/

#ifndef linearInterpolator_H
#define linearInterpolator_H

#include "Interpolator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class linearInterpolator Declaration
\*---------------------------------------------------------------------------*/

class linearInterpolator
:
    public Interpolator
{
public:

    //- Runtime type information
    TypeName("linear");


    // Constructors

        //- Construct from sampled times and the uniform grid
        linearInterpolator
        (
            const UList<scalar>& tSampled,
            const scalar t0,
            const scalar dt,
            const label n
        );


    //- Destructor
    virtual ~linearInterpolator() = default;


    // Member Functions

        //- Interpolate sampled values to the uniform grid
        virtual void interpolate
        (
            const UList<scalar>& values,
            UList<scalar>& result
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
