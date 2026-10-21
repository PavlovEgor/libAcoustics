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

\*---------------------------------------------------------------------------*/

#include "cubicSplineInterpolator.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cubicSplineInterpolator, 0);
    addToRunTimeSelectionTable(Interpolator, cubicSplineInterpolator, times);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cubicSplineInterpolator::cubicSplineInterpolator
(
    const UList<scalar>& tSampled,
    const scalar t0,
    const scalar dt,
    const label n
)
:
    Interpolator(tSampled, t0, dt, n),
    h_(tSampled.size() - 1),
    invH_(h_.size()),
    cPrime_(tSampled.size(), Zero),
    invDiag_(tSampled.size(), Zero),
    wLower_(n),
    wUpper_(n),
    secondDeriv_(tSampled.size(), Zero)
{
    forAll(h_, i)
    {
        h_[i] = tSampled[i+1] - tSampled[i];

        if (h_[i] <= 0)
        {
            FatalErrorInFunction
                << "Sampled times must be strictly ascending, got "
                << tSampled[i] << " and " << tSampled[i+1]
                << " at indices " << i << " and " << i+1
                << exit(FatalError);
        }

        invH_[i] = 1/h_[i];
    }

    // Forward elimination of the matrix, M[0] = 0 gives cPrime_[0] = 0
    for (label i = 1; i < tSampled.size() - 1; ++i)
    {
        invDiag_[i] = 1/(2*(h_[i-1] + h_[i]) - h_[i-1]*cPrime_[i-1]);
        cPrime_[i] = h_[i]*invDiag_[i];
    }

    forAll(wLower_, j)
    {
        const scalar h2 = sqr(h_[lower_[j]])/6;
        const scalar u = u_[j];

        wLower_[j] = h2*(pow3(1 - u) - (1 - u));
        wUpper_[j] = h2*(pow3(u) - u);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cubicSplineInterpolator::interpolate
(
    const UList<scalar>& values,
    UList<scalar>& result
) const
{
    scalarList& M = secondDeriv_;

    const label lastI = M.size() - 1;

    // Forward substitution, M[0] = 0
    scalar slopeLower = (values[1] - values[0])*invH_[0];

    for (label i = 1; i < lastI; ++i)
    {
        const scalar slopeUpper = (values[i+1] - values[i])*invH_[i];

        M[i] =
            (6*(slopeUpper - slopeLower) - h_[i-1]*M[i-1])*invDiag_[i];

        slopeLower = slopeUpper;
    }

    // Backward substitution, M[lastI] = 0
    for (label i = lastI - 1; i > 0; --i)
    {
        M[i] -= cPrime_[i]*M[i+1];
    }

    forAll(result, j)
    {
        const label i = lower_[j];
        const scalar y0 = values[i];

        result[j] =
            y0 + u_[j]*(values[i+1] - y0)
          + wLower_[j]*M[i] + wUpper_[j]*M[i+1];
    }
}


// ************************************************************************* //
