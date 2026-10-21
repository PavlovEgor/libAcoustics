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

#include <fftw3.h>
#include <algorithm>

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

    const label n = pressure_.size();

    if (n < 2 || times_.back() <= times_.front())
    {
        return;
    }

    // uniform FFT grid with the mean sampled time step
    const scalar dt = (times_.back() - times_.front())/(n - 1);

    // samples with the constant time step are used as they are
    bool uniform = true;

    for (label i = 1; uniform && i < n; ++i)
    {
        uniform = mag(times_[i] - times_[i-1] - dt) < 1e-6*dt;
    }

    // otherwise they are interpolated to the uniform grid
    autoPtr<Interpolator> interpolator;

    if (!uniform)
    {
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

        interpolator.reset(ctorPtr(times_, times_.front(), dt, n));
    }

    // real-to-complex FFT: n/2 + 1 frequencies up to the Nyquist one
    const label nFreq = n/2 + 1;

    double* fftIn = fftw_alloc_real(n);
    fftw_complex* fftOut = fftw_alloc_complex(nFreq);

    fftw_plan fftPlan =
        fftw_plan_dft_r2c_1d(n, fftIn, fftOut, FFTW_ESTIMATE);

    if (interpolator)
    {
        UList<scalar> pUniform(fftIn, n);

        interpolator->interpolate(pressure_, pUniform);
    }
    else
    {
        std::copy(pressure_.cbegin(), pressure_.cend(), fftIn);
    }

    fftw_execute(fftPlan);

    freq.resize(nFreq);
    amp.resize(nFreq);
    spl.resize(nFreq);

    // n samples cover the period n*dt
    const scalar period = n*dt;

    forAll(freq, k)
    {
        // a single value at zero and Nyquist frequencies: its amplitude
        // is |X|/n and the rms is equal to it, otherwise 2|X|/n and
        // amplitude/sqrt(2)
        const bool single = (k == 0 || 2*k == n);

        freq[k] = k/period;

        amp[k] =
            (single ? 1.0 : 2.0)/n
           *Foam::sqrt(sqr(fftOut[k][0]) + sqr(fftOut[k][1]));

        const scalar pRms = single ? amp[k] : amp[k]/Foam::sqrt(2.0);

        spl[k] = (amp[k] > SMALL) ? 20*Foam::log10(pRms/pRef_) : 0;
    }

    fftw_destroy_plan(fftPlan);
    fftw_free(fftIn);
    fftw_free(fftOut);
}


// ************************************************************************* //
