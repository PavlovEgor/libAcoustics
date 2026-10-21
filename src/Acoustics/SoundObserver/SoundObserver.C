/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

#include "SoundObserver.H"
#include "Interpolator.H"
#include "amplitudeSpectrum.H"
#include "SubList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SoundObserver::SoundObserver()
:
    SoundObserver(word::null, Zero, 2e-5, 1024)
{}


Foam::SoundObserver::SoundObserver
(
    const word& name,
    const vector& position,
    const scalar pRef,
    const label fftFreq
)
:
    name_(name),
    position_(position),
    pRef_(pRef),
    fftFreq_(fftFreq),
    times_(),
    pressure_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SoundObserver::append(const scalar t, const scalar p)
{
    times_.push_back(t);
    pressure_.push_back(p);
}


bool Foam::SoundObserver::fftDue() const
{
    const label n = pressure_.size();

    return fftFreq_ > 0 && n > 1 && n % fftFreq_ == 0;
}


void Foam::SoundObserver::fft
(
    const word& interpolationScheme,
    scalarField& freq,
    scalarField& amp,
    scalarField& spl
) const
{
    freq.clear();
    amp.clear();
    spl.clear();

    if (pressure_.size() < 2 || times_.back() <= times_.front())
    {
        return;
    }

    // uniform grid over the end of the record with the mean time step
    amplitudeSpectrum spectrum(times_);

    const label start = spectrum.start();
    const label n = spectrum.size();
    const scalar dt = spectrum.dt();

    // samples with the constant time step are used as they are
    bool uniform = true;

    for (label i = start + 1; uniform && i < times_.size(); ++i)
    {
        uniform = mag(times_[i] - times_[i-1] - dt) < 1e-6*dt;
    }

    if (uniform)
    {
        spectrum.transform(SubList<scalar>(pressure_, n, start));
    }
    else
    {
        // otherwise they are interpolated to the grid. The scheme is
        // selected without the New() message, the update is frequent
        auto* ctorPtr =
            Interpolator::timesConstructorTable(interpolationScheme);

        if (!ctorPtr)
        {
            FatalErrorInLookup
            (
                "Interpolator",
                interpolationScheme,
                *Interpolator::timesConstructorTablePtr_
            ) << exit(FatalError);
        }

        autoPtr<Interpolator> interpolator
        (
            ctorPtr(times_, spectrum.t0(), dt, n)
        );

        scalarList pUniform(n);
        interpolator->interpolate(pressure_, pUniform);

        spectrum.transform(pUniform);
    }

    const label nFreq = spectrum.nFrequencies();

    freq.resize(nFreq);
    amp.resize(nFreq);
    spl.resize(nFreq);

    forAll(freq, k)
    {
        freq[k] = spectrum.frequency(k);
        amp[k] = mag(spectrum.amplitude(k));
        spl[k] =
            (amp[k] > SMALL) ? 20*Foam::log10(spectrum.rms(k)/pRef_) : 0;
    }
}


// ************************************************************************* //
