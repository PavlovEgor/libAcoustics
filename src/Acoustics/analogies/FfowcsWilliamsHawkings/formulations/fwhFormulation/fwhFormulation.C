/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
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

#include "fwhFormulation.H"
#include "FfowcsWilliamsHawkings.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fwhFormulation, 0);
    defineRunTimeSelectionTable(fwhFormulation, analogy);
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fwhFormulation::initialize()
{
    const List<SoundObserver>& observers = fwh_.getSoundObservers();
    const PtrList<sampledSurface>& surfaces = fwh_.controlSurfaces();

    const label nObs = observers.size();
    const label nSurf = surfaces.size();

    // source histories and observer times
    qds_.resize(nObs);
    fds_.resize(nObs);
    lastIndex_.resize(nObs);
    tobs_.resize(nObs);
    robs_.resize(nObs);
    magrobs_.resize(nObs);

    forAll(observers, iObs)
    {
        qds_[iObs].resize(nSurf);
        fds_[iObs].resize(nSurf);
        lastIndex_[iObs].resize(nSurf);
        tobs_[iObs].resize(nSurf);
        robs_[iObs].resize(nSurf);
        magrobs_[iObs].resize(nSurf);

        forAll(surfaces, iSurf)
        {
            const label nFaces = surfaces[iSurf].Cf().size();

            qds_[iObs][iSurf].resize(nFaces);
            fds_[iObs][iSurf].resize(nFaces);
            lastIndex_[iObs][iSurf].resize(nFaces, 0);
            tobs_[iObs][iSurf].resize(nFaces, 0.0);
            robs_[iObs][iSurf].resize(nFaces);
            magrobs_[iObs][iSurf].resize(nFaces);
        }
    }

    tauMax_.resize(nObs, 0.0);
    rMax_.resize(nObs, 0.0);

    if (fwh_.fixedResponseDelay())
    {
        tauMin_.resize(nObs);
        rMin_.resize(nObs);
    }

    calcDistances();

    // normals
    ni_.resize(nSurf);
    nl_.resize(nSurf);

    forAll(surfaces, iSurf)
    {
        const label nFaces = surfaces[iSurf].Cf().size();

        ni_[iSurf].resize(nFaces);
        nl_[iSurf].resize(nFaces);
    }

    calcNormals();
}


void Foam::functionObjects::fwhFormulation::calcDistances()
{
    const List<SoundObserver>& observers = fwh_.getSoundObservers();
    const PtrList<sampledSurface>& surfaces = fwh_.controlSurfaces();

    // effective distance in the uniform flow U0
    const scalar M0 = mag(fwh_.U0())/fwh_.c0();
    const scalar sqrBeta = 1 - sqr(M0);

    forAll(observers, iObs)
    {
        rMax_[iObs] = 0.0;
        scalar rMin = GREAT;

        forAll(surfaces, iSurf)
        {
            const vectorField& Cf = surfaces[iSurf].Cf();

            forAll(Cf, i)
            {
                const vector r = observers[iObs].position() - Cf[i];

                robs_[iObs][iSurf][i] = r;

                const scalar R =
                    sqrt(sqr(r[0]) + sqrBeta*(sqr(r[1]) + sqr(r[2])));

                const scalar magr = (-M0*r[0] + R)/sqrBeta;

                magrobs_[iObs][iSurf][i] = magr;

                if (magr > rMax_[iObs])
                {
                    rMax_[iObs] = magr;
                }

                if (magr < rMin)
                {
                    rMin = magr;
                }
            }
        }

        reduce(rMax_[iObs], maxOp<scalar>());
        tauMax_[iObs] = rMax_[iObs]/fwh_.c0();

        if (tauMin_.size())
        {
            reduce(rMin, minOp<scalar>());

            // no faces: zero delay
            rMin_[iObs] = (rMin < GREAT) ? rMin : rMax_[iObs];
            tauMin_[iObs] = rMin_[iObs]/fwh_.c0();
        }
    }
}


void Foam::functionObjects::fwhFormulation::calcNormals()
{
    const PtrList<sampledSurface>& surfaces = fwh_.controlSurfaces();

    forAll(ni_, iSurf)
    {
        const vectorField& Sf = surfaces[iSurf].Sf();
        const vectorField& Cf = surfaces[iSurf].Cf();

        // centre of the surface
        vector Cs = gSum(Cf);
        scalar surfSize = scalar(Cf.size());
        reduce(surfSize, sumOp<scalar>());
        Cs /= surfSize;

        // normals point outwards from the centre
        forAll(ni_[iSurf], iFace)
        {
            ni_[iSurf].value(iFace) = Sf[iFace]/mag(Sf[iFace]);

            nl_[iSurf][iFace] =
                ((Cf[iFace] - Cs) & ni_[iSurf].value(iFace)) > 0 ? 1.0 : -1.0;

            ni_[iSurf].value(iFace) *= nl_[iSurf][iFace];
        }
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fwhFormulation::sourcesAt
(
    label iObs,
    label iSurf,
    label iFace,
    scalar tau,
    scalar& q,
    scalar& f
)
{
    // qds and fds have the same times
    const DynamicList<scalar>& t = qds_[iObs][iSurf][iFace].first();
    const DynamicList<scalar>& qv = qds_[iObs][iSurf][iFace].second();
    const DynamicList<scalar>& fv = fds_[iObs][iSurf][iFace].second();

    const label n = t.size();

    q = 0;
    f = 0;

    // zero outside of the history, exact values at its bounds
    if (n < 1 || tau < t[0] || tau > t[n-1])
    {
        return;
    }
    if (tau == t[0])
    {
        q = qv[0];
        f = fv[0];
        return;
    }
    if (tau == t[n-1])
    {
        q = qv[n-1];
        f = fv[n-1];
        return;
    }

    // t[0] < tau < t[n-1]: interval t[li] <= tau < t[li+1] from the last one
    label& li = lastIndex_[iObs][iSurf][iFace];
    li = min(max(li, label(0)), n - 2);

    while (t[li] > tau)
    {
        --li;
    }
    while (t[li+1] <= tau)
    {
        ++li;
    }

    if (tau == t[li])
    {
        q = qv[li];
        f = fv[li];
        return;
    }

    const label ui = li + 1;

    q = qv[li] + ((qv[ui] - qv[li]) / (t[ui] - t[li])) * (tau - t[li]);
    f = fv[li] + ((fv[ui] - fv[li]) / (t[ui] - t[li])) * (tau - t[li]);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fwhFormulation::fwhFormulation
(
    const FfowcsWilliamsHawkings& fwh
)
:
    fwh_(fwh),
    fwhProbeI_(0),
    qds_(0),
    fds_(0),
    tobs_(0),
    robs_(0),
    magrobs_(0),
    ni_(0),
    nl_(0),
    rMax_(0),
    tauMax_(0),
    rMin_(0),
    tauMin_(0),
    lastIndex_(0)
{
    initialize();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::functionObjects::fwhFormulation>
Foam::functionObjects::fwhFormulation::New
(
    const word& type,
    const FfowcsWilliamsHawkings& fwh
)
{
    Info<< "Selecting FW-H formulation " << type << endl;

    auto* ctorPtr = analogyConstructorTable(type);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "fwhFormulation",
            type,
            *analogyConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<fwhFormulation>(ctorPtr(fwh));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::fwhFormulation::update()
{
    const scalar ct = fwh_.time().value();

    // moving surfaces: update the geometry
    if (mag(fwh_.Ufwh()) > SMALL || fwh_.nonUniformSurfaceMotion())
    {
        calcDistances();
        calcNormals();
    }

    // observer times of the current sources
    forAll(tobs_, iObs)
    {
        forAll(tobs_[iObs], iSurf)
        {
            forAll(tobs_[iObs][iSurf], i)
            {
                tobs_[iObs][iSurf][i] =
                    ct + magrobs_[iObs][iSurf][i]/fwh_.c0();

                if (tauMin_.size())
                {
                    tobs_[iObs][iSurf][i] -=
                        (tauMin_[iObs] - fwh_.responseDelay());
                }
            }
        }
    }
}


void Foam::functionObjects::fwhFormulation::clearExpiredData()
{
    scalar ct = fwh_.time().value();
    reduce(ct, minOp<scalar>());

    fwhProbeI_++;

    if (fwhProbeI_ % fwh_.cleanFreq() != 0)
    {
        return;
    }

    fwhProbeI_ = 0;

    forAll(qds_, iObs)
    {
        // sources emitted before expiredTime have reached the observer
        const scalar expiredTime =
            tauMin_.size()
          ? ct - (tauMax_[iObs] - tauMin_[iObs] + fwh_.responseDelay())
          : ct - tauMax_[iObs];

        forAll(qds_[iObs], iSurf)
        {
            forAll(qds_[iObs][iSurf], iFace)
            {
                // -1 if nothing found, 0 .. size-1 for indices to remove
                const label expiredIndex =
                    findExpiredIndex(qds_[iObs][iSurf][iFace], expiredTime);

                if (expiredIndex > -1)
                {
                    qds_[iObs][iSurf][iFace] =
                        getNewPointData(qds_[iObs][iSurf][iFace], expiredIndex);
                    fds_[iObs][iSurf][iFace] =
                        getNewPointData(fds_[iObs][iSurf][iFace], expiredIndex);
                }
            }
        }
    }
}


Foam::label Foam::functionObjects::fwhFormulation::findExpiredIndex
(
    const pointTimeData& timeData,
    scalar expiredTime
) const
{
    label ui = timeData.first().size();
    label li = 0;
    label mi = (ui + li) / 2;

    if (ui < 1)
    {
        return -1;
    }
    else
    {
        if ((timeData.first()[ui-1] < expiredTime) && (mi > 1))
        {
            Info << "WARNING: ExpiredTime is out of bounds." << endl;
            return -1;
        }
    }
    if (mi < 1)
    {
        if (timeData.first()[mi] < expiredTime)
        {
            return 0;
        }
        else
        {
            return -1;
        }
    }
    else
    {
        if (timeData.first()[li] >= expiredTime)
        {
            return -1;
        }
        do
        {
            if (timeData.first()[mi] < expiredTime)
            {
                li = mi;
                mi = (ui + li) / 2;
            }
            else
            {
                ui = mi;
                mi = (ui + li) / 2;
            }
            if (expiredTime == timeData.first()[mi])
            {
                return (mi - 1);
            }
        }
        while
        (
            ((ui - li) != 1)
        );
    }

    if (ui < timeData.size())
    {
        if (timeData.first()[ui] < expiredTime)
        {
            return ui;
        }
    }

    return li;
}


Foam::scalar Foam::functionObjects::fwhFormulation::valueAt
(
    const surfaceTimeData& data,
    label iObs,
    label iSurf,
    label iFace,
    scalar tau
) const
{
    const pointTimeData& timeData = data[iObs][iSurf][iFace];

    label ui = timeData.first().size();
    label li = 0;
    label mi = (ui + li) / 2;
    //exit if size is zero
    if (ui < 1)
    {
        return 0.0;
    }
    //exit if tau out of bounds
    if (
        (tau < timeData.first()[li])
        ||
        (tau > timeData.first()[ui-1])
       )
    {
        return 0.0;
    }
    //exit if tau is at bounds
    if (tau == timeData.first()[li])
    {
        return timeData.second()[li];
    }
    if (tau == timeData.first()[ui-1])
    {
        return timeData.second()[ui-1];
    }
    //check that size is sufficient
    if (mi >= 1)
    {
        do
        {
            if (tau < timeData.first()[mi])
            {
                ui = mi;
                mi = (ui + li) / 2;
            }
            else if (tau > timeData.first()[mi])
            {
                li = mi;
                mi = (ui + li) / 2;
            }
            else //(tau == timeData.first()[mi])
            {
                return timeData.second()[mi];
            }
        }
        while
        (
            ((ui - li) != 1)
        );
    }
    //return value between li and ui
    return timeData.second()[li] +
    (
        (timeData.second()[ui] - timeData.second()[li])
        /
        (timeData.first()[ui] - timeData.first()[li])
    ) * (tau - timeData.first()[li]);
}


Foam::functionObjects::fwhFormulation::pointTimeData
Foam::functionObjects::fwhFormulation::getNewPointData
(
    const pointTimeData& timeData,
    label expiredIndex
) const
{
    if (expiredIndex < 0)
    {
        return timeData;
    }

    const label newsize = timeData.first().size() - (expiredIndex + 1);

    pointTimeData newPointData;
    newPointData.first().resize(newsize);
    newPointData.second().resize(newsize);

    for (label iTime = expiredIndex + 1; iTime < timeData.first().size(); ++iTime)
    {
        newPointData.first()[iTime-(expiredIndex+1)] = timeData.first()[iTime];
        newPointData.second()[iTime-(expiredIndex+1)] = timeData.second()[iTime];
    }

    return newPointData;
}


// ************************************************************************* //
