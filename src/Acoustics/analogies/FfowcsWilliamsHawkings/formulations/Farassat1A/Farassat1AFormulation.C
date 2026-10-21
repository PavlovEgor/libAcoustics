/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "Farassat1AFormulation.H"
#include "FfowcsWilliamsHawkings.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(Farassat1AFormulation, 0);

    addToRunTimeSelectionTable
    (
        fwhFormulation,
        Farassat1AFormulation,
        analogy
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::Farassat1AFormulation::calculateAcousticPressure
(
    const vectorField& Sf,
    const vectorField& uS,
    const scalarField& rhoS,
    const scalarField& pS,
    label iObs,
    label iSurf,
    scalar ct
)
{
    const scalar c0 = fwh_.c0();
    const scalar rhoRef = fwh_.rhoRef();
    const vectorField& vS = fwh_.surfaceVelocities()[iSurf];

    forAll(Sf, iFace)
    {
        const vector& r = robs_[iObs][iSurf][iFace];
        const scalar magr = magrobs_[iObs][iSurf][iFace];
        const vector rh = r / magr;
        const scalar dS = mag(Sf[iFace]);
        const vector& n = ni_[iSurf].value(iFace);
        const vector M = vS[iFace] / c0;
        const scalar Mr = M & rh;
        const scalar magM = mag(M);

        const vector U =
            (1.0 - rhoS[iFace] / rhoRef) * vS[iFace]
          + rhoS[iFace] * uS[iFace] / rhoRef;

        const tensor Pf =
            pS[iFace]*tensor::I
          + rhoS[iFace]*uS[iFace]*(uS[iFace] - vS[iFace]);

        const vector L = Pf & n;
        const scalar lM = L & M;
        const scalar lr = L & rh;
        const scalar Un = U & n;

        const scalar OneByOneMr = 1.0 / (1.0 - Mr);
        const scalar OneByOneMrSq = OneByOneMr*OneByOneMr;

        Un_[iObs][iSurf].value(iFace) = Un;
        L_[iObs][iSurf].value(iFace) = L;
        M_[iObs][iSurf].value(iFace) = M;

        const scalar dotlr = L_[iObs][iSurf].dot(ct, iFace) & rh;
        const scalar dotMr = M_[iObs][iSurf].dot(ct, iFace) & rh;
        const scalar dotUn = Un_[iObs][iSurf].dot(ct, iFace);

        // thickness source
        qds_[iObs][iSurf][iFace].first().append(tobs_[iObs][iSurf][iFace]);
        qds_[iObs][iSurf][iFace].second().append
        (
            (
                rhoRef * (dotUn) * OneByOneMrSq / magr
                +
                rhoRef * Un * (magr * dotMr + c0 * (Mr - magM*magM)) *
                OneByOneMrSq * OneByOneMr / magr / magr
            )*dS
        );

        // loading source
        const scalar fpart1 = dotlr * (dS / magr / c0) * OneByOneMrSq;
        const scalar fpart2 = (lr - lM) * (dS / magr / magr) * OneByOneMrSq;
        const scalar fpart3 =
            lr * (dS / magr / magr / c0) * OneByOneMrSq * OneByOneMr *
            (magr * dotMr + c0 * Mr - c0 * magM * magM);

        fds_[iObs][iSurf][iFace].first().append(tobs_[iObs][iSurf][iFace]);
        fds_[iObs][iSurf][iFace].second().append(fpart1 + fpart2 + fpart3);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::Farassat1AFormulation::Farassat1AFormulation
(
    const FfowcsWilliamsHawkings& fwh
)
:
    fwhFormulation(fwh),
    Un_(0),
    L_(0),
    M_(0),
    intDotQdS_(fwh.getSoundObservers().size(), 0.0),
    intFdS_(fwh.getSoundObservers().size(), 0.0)
{
    const PtrList<sampledSurface>& surfaces = fwh_.controlSurfaces();
    const label nObs = fwh_.getSoundObservers().size();
    const bool varTimeStep = fwh_.time().isAdjustTimeStep();

    L_.resize(nObs);
    M_.resize(nObs);
    Un_.resize(nObs);

    forAll(L_, iObs)
    {
        L_[iObs].resize(surfaces.size());
        M_[iObs].resize(surfaces.size());
        Un_[iObs].resize(surfaces.size());

        forAll(surfaces, iSurf)
        {
            const label nFaces = surfaces[iSurf].Cf().size();

            L_[iObs][iSurf].setTimeStepVariable(varTimeStep);
            L_[iObs][iSurf].resize(nFaces);
            M_[iObs][iSurf].setTimeStepVariable(varTimeStep);
            M_[iObs][iSurf].resize(nFaces);
            Un_[iObs][iSurf].setTimeStepVariable(varTimeStep);
            Un_[iObs][iSurf].resize(nFaces);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::functionObjects::Farassat1AFormulation::observerAcousticPressure
(
    const vectorField& Sf,
    const vectorField& uS,
    const scalarField& rhoS,
    const scalarField& pS,
    label iObs,
    label iSurf,
    scalar ct
)
{
    calculateAcousticPressure(Sf, uS, rhoS, pS, iObs, iSurf, ct);

    // slightly increase time to get inside of time step
    const scalar ct1 = ct + fwh_.time().deltaTValue()*1.0e-6;

    // acoustic pressure, zero if the source didn't reach the observer
    intDotQdS_[iObs] = 0.0;
    intFdS_[iObs] = 0.0;

    forAll(qds_[iObs][iSurf], iFace)
    {
        scalar q, f;
        sourcesAt(iObs, iSurf, iFace, ct1, q, f);

        intDotQdS_[iObs] += q;
        intFdS_[iObs] += f;
    }

    reduce(intDotQdS_[iObs], sumOp<scalar>());
    reduce(intFdS_[iObs], sumOp<scalar>());

    const scalar coeff1 = 1. / 4. / constant::mathematical::pi;

    return (intDotQdS_[iObs] + intFdS_[iObs])*coeff1;
}


// ************************************************************************* //
