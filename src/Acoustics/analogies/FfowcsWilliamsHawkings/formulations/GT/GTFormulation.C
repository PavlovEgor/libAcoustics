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

#include "GTFormulation.H"
#include "FfowcsWilliamsHawkings.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(GTFormulation, 0);

    addToRunTimeSelectionTable
    (
        fwhFormulation,
        GTFormulation,
        analogy
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::GTFormulation::calculateAcousticPressure
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
    const vector& U0 = fwh_.U0();
    const vectorField& vS = fwh_.surfaceVelocities()[iSurf];

    const scalar M0 = mag(U0) / c0;
    const scalar sqrBetta = 1 - sqr(M0);

    forAll(Sf, iFace)
    {
        const vector& r = robs_[iObs][iSurf][iFace];

        const scalar Rz = sqrt(sqr(r[0]) + sqrBetta * (sqr(r[1]) + sqr(r[2])));
        const scalar R = -M0 * r[0] / sqrBetta + Rz / sqrBetta;
        const vector Rh
        (
            (-M0 * Rz + r[0]) / sqrBetta / R,
            r[1] / R,
            r[2] / R
        );

        const scalar dS = mag(Sf[iFace]);
        const vector& n = ni_[iSurf].value(iFace);
        const vector M = vS[iFace] / c0;
        const scalar magM = mag(M);

        const vector Qi =
            (-rhoRef * U0 + rhoS[iFace] * (uS[iFace] + U0));
        const tensor Lij =
            pS[iFace] * tensor::I + rhoS[iFace] * uS[iFace] * (uS[iFace] + U0);
        const scalar Qn = Qi & n;
        const vector Li = Lij & n;

        const scalar LR = Li & Rh;
        const scalar LM = Li & M;
        const scalar MR = M & Rh;
        const scalar OneByOneMr = 1.0 / (1.0 - MR);
        const scalar OneByOneMrSq = OneByOneMr * OneByOneMr;

        Q_[iObs][iSurf].value(iFace) = Qi;
        L_[iObs][iSurf].value(iFace) = Li;

        const scalar dotQn = Q_[iObs][iSurf].dot(ct, iFace) & n;
        const scalar dotLR = L_[iObs][iSurf].dot(ct, iFace) & Rh;

        // thickness source
        qds_[iObs][iSurf][iFace].first().append(tobs_[iObs][iSurf][iFace]);
        qds_[iObs][iSurf][iFace].second().append
        (
            (
                (dotQn) * OneByOneMrSq / R
                +
                Qn * c0 * (MR - magM*magM) * OneByOneMrSq * OneByOneMr / R / R
            ) * dS
        );

        // loading source
        const scalar fpart1 = (dotLR) / R / c0 * OneByOneMrSq * dS;
        const scalar fpart2 = (LR - LM) / R / R * OneByOneMrSq * dS;
        const scalar fpart3 =
            LR * (MR - magM * magM) / R / R * OneByOneMrSq * OneByOneMr * dS;

        fds_[iObs][iSurf][iFace].first().append(tobs_[iObs][iSurf][iFace]);
        fds_[iObs][iSurf][iFace].second().append(fpart1 + fpart2 + fpart3);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::GTFormulation::GTFormulation
(
    const FfowcsWilliamsHawkings& fwh
)
:
    fwhFormulation(fwh),
    Q_(0),
    L_(0),
    intDotQdS_(fwh.getSoundObservers().size(), 0.0),
    intFdS_(fwh.getSoundObservers().size(), 0.0)
{
    const PtrList<sampledSurface>& surfaces = fwh_.controlSurfaces();
    const label nObs = fwh_.getSoundObservers().size();
    const bool varTimeStep = fwh_.time().isAdjustTimeStep();

    L_.resize(nObs);
    Q_.resize(nObs);

    forAll(L_, iObs)
    {
        L_[iObs].resize(surfaces.size());
        Q_[iObs].resize(surfaces.size());

        forAll(surfaces, iSurf)
        {
            const label nFaces = surfaces[iSurf].Cf().size();

            L_[iObs][iSurf].setTimeStepVariable(varTimeStep);
            L_[iObs][iSurf].resize(nFaces);
            Q_[iObs][iSurf].setTimeStepVariable(varTimeStep);
            Q_[iObs][iSurf].resize(nFaces);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::functionObjects::GTFormulation::observerAcousticPressure
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

    // acoustic pressure of all surfaces,
    // zero if the source didn't reach the observer
    intDotQdS_[iObs] = 0.0;
    intFdS_[iObs] = 0.0;

    forAll(qds_[iObs], surfI)
    {
        forAll(qds_[iObs][surfI], iFace)
        {
            intDotQdS_[iObs] += valueAt(qds_, iObs, surfI, iFace, ct1);
            intFdS_[iObs] += valueAt(fds_, iObs, surfI, iFace, ct1);
        }
    }

    reduce(intDotQdS_[iObs], sumOp<scalar>());
    reduce(intFdS_[iObs], sumOp<scalar>());

    const scalar coeff1 = 1. / 4. / constant::mathematical::pi;

    return (intDotQdS_[iObs] + intFdS_[iObs])*coeff1;
}


// ************************************************************************* //
