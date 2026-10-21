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
    Foam::cubicSplineInterpolator

Description
    Natural cubic spline interpolation on non-uniform sampled times.

    Second derivatives M at the sampled points satisfy the tridiagonal system
    \verbatim
        h[i-1] M[i-1] + 2 (h[i-1] + h[i]) M[i] + h[i] M[i+1]
      = 6 ((y[i+1] - y[i])/h[i] - (y[i] - y[i-1])/h[i-1]),

        M[0] = M[n-1] = 0
    \endverbatim
    The matrix depends on the sampled times only, so it is factorised once
    and each signal costs one forward and one backward substitution.

    Amplitude error is below 0.1 dB up to f*dt = 0.2 of the sampled step dt.
    Requires strictly ascending sampled times.

    The second derivatives are kept in a workspace member, so one object
    must not be used from several threads at once.

SourceFiles
    cubicSplineInterpolator.C

\*---------------------------------------------------------------------------*/

#ifndef cubicSplineInterpolator_H
#define cubicSplineInterpolator_H

#include "../Interpolator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class cubicSplineInterpolator Declaration
\*---------------------------------------------------------------------------*/

class cubicSplineInterpolator
:
    public Interpolator
{
    // Private Data

        //- Sampled intervals
        scalarList h_;

        //- Inverse of the sampled intervals
        scalarList invH_;

        //- Super-diagonal of the factorised system
        scalarList cPrime_;

        //- Inverse of the diagonal of the factorised system
        scalarList invDiag_;

        //- Weight of the second derivative at the lower sampled point
        scalarList wLower_;

        //- Weight of the second derivative at the upper sampled point
        scalarList wUpper_;

        //- Second derivatives at the sampled points (workspace)
        mutable scalarList secondDeriv_;


public:

    //- Runtime type information
    TypeName("cubicSpline");


    // Constructors

        //- Construct from sampled times and the uniform grid
        cubicSplineInterpolator
        (
            const UList<scalar>& tSampled,
            const scalar t0,
            const scalar dt,
            const label n
        );


    //- Destructor
    virtual ~cubicSplineInterpolator() = default;


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
