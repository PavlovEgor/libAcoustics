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

#include "Interpolator.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Interpolator, 0);
    defineRunTimeSelectionTable(Interpolator, times);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Interpolator::Interpolator
(
    const UList<scalar>& tSampled,
    const scalar t0,
    const scalar dt,
    const label n
)
:
    lower_(n),
    u_(n)
{
    if (tSampled.size() < 2)
    {
        FatalErrorInFunction
            << "At least 2 sampled times are required, got "
            << tSampled.size() << exit(FatalError);
    }

    const label lastI = tSampled.size() - 2;

    label i = 0;

    forAll(lower_, j)
    {
        const scalar t = t0 + j*dt;

        while (i < lastI && tSampled[i+1] < t)
        {
            ++i;
        }

        const scalar h = tSampled[i+1] - tSampled[i];

        lower_[j] = i;
        u_[j] =
        (
            h > 0
          ? min(max((t - tSampled[i])/h, scalar(0)), scalar(1))
          : scalar(0)
        );
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::Interpolator> Foam::Interpolator::New
(
    const word& type,
    const UList<scalar>& tSampled,
    const scalar t0,
    const scalar dt,
    const label n
)
{
    Info<< "Selecting interpolation scheme " << type << endl;

    auto* ctorPtr = timesConstructorTable(type);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "Interpolator",
            type,
            *timesConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<Interpolator>(ctorPtr(tSampled, t0, dt, n));
}


// ************************************************************************* //
