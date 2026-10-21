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

#include "FfowcsWilliamsHawkings.H"
#include "fwhFormulation.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(FfowcsWilliamsHawkings, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        FfowcsWilliamsHawkings,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::FfowcsWilliamsHawkings::initialize()
{
    if (nonUniformSurfaceMotion_)
    {
        Cf0_.resize(controlSurfaces_.size());

        forAll(controlSurfaces_, iSurf)
        {
            Cf0_[iSurf] = controlSurfaces_[iSurf].Cf();
        }
    }

    vS_.resize(controlSurfaces_.size());

    forAll(controlSurfaces_, iSurf)
    {
        vS_[iSurf].resize(controlSurfaces_[iSurf].Cf().size());
        vS_[iSurf] = vector::zero;
    }

    fwhFormulationPtr_ = fwhFormulation::New(formulationType_, *this);
}


void Foam::functionObjects::FfowcsWilliamsHawkings::updateSurfaceCf()
{
    if (nonUniformSurfaceMotion_)
    {
        expire();
        update();

        forAll(controlSurfaces_, iSurf)
        {
            forAll(Cf0_[iSurf], iF)
            {
                Cf0_[iSurf][iF] = controlSurfaces_[iSurf].Cf()[iF];
            }
        }
    }
}


void Foam::functionObjects::FfowcsWilliamsHawkings::correct()
{
    if (nonUniformSurfaceMotion_)
    {
        expire();
    }
    update();

    const scalar dt = time().deltaTValue();
    const scalar ct = time().value();

    // surface velocities
    forAll(controlSurfaces_, iSurf)
    {
        if (nonUniformSurfaceMotion_)
        {
            vS_[iSurf] = (controlSurfaces_[iSurf].Cf() - Cf0_[iSurf])/dt;

            // positions for the next time step, also at every probe step
            Cf0_[iSurf] = controlSurfaces_[iSurf].Cf();
        }
        else
        {
            vS_[iSurf] = Ufwh_;
        }
    }

    // formulation-specific data
    fwhFormulationPtr_->update();

    // flow on the control surfaces
    List<vectorField> uS(controlSurfaces_.size());
    List<scalarField> rhoS(controlSurfaces_.size());
    List<scalarField> pS(controlSurfaces_.size());

    forAll(controlSurfaces_, iSurf)
    {
        const sampledSurface& surf = controlSurfaces_[iSurf];

        if (surf.interpolate())
        {
            WarningInFunction
                << "Interpolation for surface " << surf.name()
                << " is on, turn it off" << endl;
        }

        uS[iSurf] = surfaceVelocity(surf);
        rhoS[iSurf] = surfaceDensity(surf);
        pS[iSurf] = surfacePressure(surf)() - pInf_;
    }

    // acoustic pressure is the sum over all control surfaces
    forAll(observers_, iObs)
    {
        scalar oap = 0;

        forAll(controlSurfaces_, iSurf)
        {
            oap += fwhFormulationPtr_->observerAcousticPressure
            (
                controlSurfaces_[iSurf].Sf(),
                uS[iSurf],
                rhoS[iSurf],
                pS[iSurf],
                iObs,
                iSurf,
                ct
            );
        }

        if (Pstream::master())
        {
            observers_[iObs].append(ct, oap);
            Log<< "OAP = " << oap << nl;
        }
    }

    // remove old data if needed
    fwhFormulationPtr_->clearExpiredData();
}


Foam::tmp<Foam::scalarField>
Foam::functionObjects::FfowcsWilliamsHawkings::surfaceDensity
(
    const sampledSurface& surface
) const
{
    return sampleOrInterpolate<scalar>(rho()(), surface);
}


Foam::tmp<Foam::vectorField>
Foam::functionObjects::FfowcsWilliamsHawkings::surfaceVelocity
(
    const sampledSurface& surface
) const
{
    const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

    return sampleOrInterpolate<vector>(U, surface);
}


Foam::tmp<Foam::scalarField>
Foam::functionObjects::FfowcsWilliamsHawkings::surfacePressure
(
    const sampledSurface& surface
) const
{
    const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);

    tmp<scalarField> pSampled = sampleOrInterpolate<scalar>(p, surface);

    if (p.dimensions() != dimPressure)
    {
        pSampled.ref() *= rhoRef_;
    }

    return pSampled;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::FfowcsWilliamsHawkings::FfowcsWilliamsHawkings
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    AcousticAnalogy(name, runTime, dict),
    formulationType_(word::null),
    fixedResponseDelay_(true),
    responseDelay_(0.0),
    Ufwh_(Zero),
    U0_(Zero),
    nonUniformSurfaceMotion_(false),
    Cf0_(0),
    vS_(0),
    pInf_(0.0),
    fwhFormulationPtr_(nullptr),
    cleanFreq_(100),
    interpolationScheme_(word::null),
    controlSurfaces_(0)
{
    read(dict);
    update();
    initialize();
}


Foam::functionObjects::FfowcsWilliamsHawkings::FfowcsWilliamsHawkings
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    AcousticAnalogy(name, obr, dict),
    formulationType_(word::null),
    fixedResponseDelay_(true),
    responseDelay_(0.0),
    Ufwh_(Zero),
    U0_(Zero),
    nonUniformSurfaceMotion_(false),
    Cf0_(0),
    vS_(0),
    pInf_(0.0),
    fwhFormulationPtr_(nullptr),
    cleanFreq_(100),
    interpolationScheme_(word::null),
    controlSurfaces_(0)
{
    read(dict);
    update();
    initialize();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::FfowcsWilliamsHawkings::~FfowcsWilliamsHawkings()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::FfowcsWilliamsHawkings::read(const dictionary& dict)
{
    if (!AcousticAnalogy::read(dict))
    {
        return false;
    }

    dict.readEntry("formulationType", formulationType_);

    if (!fwhFormulation::analogyConstructorTable(formulationType_))
    {
        FatalIOErrorInLookup
        (
            dict,
            "formulationType",
            formulationType_,
            *fwhFormulation::analogyConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    Ufwh_ = dict.getOrDefault<vector>("Ufwh", Zero);

    U0_ = dict.getOrDefault<vector>("U0", Zero);

    dict.readEntry("pInf", pInf_);

    dict.readEntry("interpolationScheme", interpolationScheme_);

    dict.readIfPresent("nonUniformSurfaceMotion", nonUniformSurfaceMotion_);

    dict.readIfPresent("fixedResponseDelay", fixedResponseDelay_);

    if (!fixedResponseDelay_ || !dict.readIfPresent("responseDelay", responseDelay_))
    {
        fixedResponseDelay_ = false;
    }

    cleanFreq_ = dict.getOrDefault<label>("cleanFreq", 100);

    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    PtrList<sampledSurface> newList
    (
        dict.lookup("surfaces"),
        sampledSurface::iNew(mesh)
    );

    controlSurfaces_.transfer(newList);

    // Ensure all surfaces are expired
    expire();

    if (controlSurfaces_.size())
    {
        Log << "Function object " << name() << ":" << nl
            << "    Reading Ffowcs Williams - Hawkings analogy control"
            << " surface description:" << nl;

        for (const sampledSurface& surf : controlSurfaces_)
        {
            Log << "        " << surf.name() << nl;
        }
        Log << endl;
    }

    if (Pstream::master() && debug)
    {
        Pout<< "FWH control surfaces additional info:" << nl << "(" << nl;

        for (const sampledSurface& surf : controlSurfaces_)
        {
            Pout<< "  " << surf << endl;
        }
        Pout<< ")" << endl;
    }

    return true;
}


bool Foam::functionObjects::FfowcsWilliamsHawkings::expire()
{
    bool justExpired = false;

    for (sampledSurface& surf : controlSurfaces_)
    {
        if (surf.expire())
        {
            justExpired = true;
        }
    }

    // true if any surfaces just expired
    return justExpired;
}


bool Foam::functionObjects::FfowcsWilliamsHawkings::needsUpdate() const
{
    for (const sampledSurface& surf : controlSurfaces_)
    {
        if (surf.needsUpdate())
        {
            return true;
        }
    }

    return false;
}


bool Foam::functionObjects::FfowcsWilliamsHawkings::update()
{
    if (!needsUpdate())
    {
        return false;
    }

    bool updated = false;

    for (sampledSurface& surf : controlSurfaces_)
    {
        if (surf.update())
        {
            updated = true;
        }
    }

    return updated;
}


// ************************************************************************* //
