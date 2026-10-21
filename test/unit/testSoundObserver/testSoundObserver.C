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

Application
    testSoundObserver

Description
    Checks the pressure history and the spectrum of SoundObserver:
    tone amplitudes, zero and Nyquist frequencies, SPL from the rms pressure
    and the variable time step.

\*---------------------------------------------------------------------------*/

#include "catch.hpp"

#include "SoundObserver.H"
#include "mathematicalConstants.H"

using namespace Foam;
using constant::mathematical::twoPi;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Largest amplitude outside of the given bins
scalar maxOtherAmp(const scalarField& amp, const labelList& bins)
{
    scalar maxAmp = 0;

    forAll(amp, k)
    {
        if (!bins.found(k))
        {
            maxAmp = max(maxAmp, amp[k]);
        }
    }

    return maxAmp;
}


TEST_CASE("Checks SoundObserver history and spectrum", "[testSoundObserver]")
{
    const scalar pRef = 2e-5;
    const label n = 200;
    const scalar dt = 1e-4;
    const scalar t0 = 0.3;

    // frequency step 1/(n*dt) = 50 Hz
    const scalar df = 1/(n*dt);

    SECTION("History and update")
    {
        SoundObserver obs("R-A", vector(0, 0, 1), pRef, n);

        REQUIRE(obs.apressure() == 0);
        REQUIRE(!obs.fftDue());

        for (label i = 0; i < n; ++i)
        {
            REQUIRE(!obs.fftDue());
            obs.append(t0 + i*dt, i);
        }

        REQUIRE(obs.fftDue());
        REQUIRE(obs.apressure() == n - 1);
        REQUIRE(obs.times().size() == n);
        REQUIRE(obs.pressure().size() == n);

        obs.append(t0 + n*dt, 0);
        REQUIRE(!obs.fftDue());

        // copies keep the history
        List<SoundObserver> observers;
        observers.append(obs);
        REQUIRE(observers[0].pressure().size() == n + 1);

        // no update without fftFreq
        SoundObserver noFft("R-B", Zero, pRef, 0);
        noFft.append(0, 1);
        noFft.append(dt, 1);
        REQUIRE(!noFft.fftDue());
    }

    SECTION("Constant time step")
    {
        SoundObserver obs("R-A", Zero, pRef, n);

        // mean, two tones on the bins 7 and 31, Nyquist component
        for (label i = 0; i < n; ++i)
        {
            const scalar t = i*dt;

            obs.append
            (
                t0 + t,
                1.5
              + 2.0*Foam::cos(twoPi*7*df*t + 0.4)
              + 0.5*Foam::sin(twoPi*31*df*t)
              + 0.25*(i % 2 ? -1 : 1)
            );
        }

        scalarField freq, amp, spl;
        obs.fft("linear", freq, amp, spl);

        // n/2 + 1 frequencies up to the Nyquist one
        REQUIRE(freq.size() == n/2 + 1);
        REQUIRE(mag(freq[1] - df) < 1e-9*df);
        REQUIRE(mag(freq.last() - 0.5/dt) < 1e-9*df);

        REQUIRE(mag(amp[0] - 1.5) < 1e-12);
        REQUIRE(mag(amp[7] - 2.0) < 1e-12);
        REQUIRE(mag(amp[31] - 0.5) < 1e-12);
        REQUIRE(mag(amp[n/2] - 0.25) < 1e-12);
        REQUIRE(maxOtherAmp(amp, labelList({0, 7, 31, n/2})) < 1e-12);

        // SPL from the rms pressure: amplitude/sqrt(2) for tones,
        // the value itself at zero and Nyquist frequencies
        REQUIRE(mag(spl[7] - 20*Foam::log10(2.0/Foam::sqrt(2.0)/pRef)) < 1e-9);
        REQUIRE(mag(spl[0] - 20*Foam::log10(1.5/pRef)) < 1e-9);
        REQUIRE(mag(spl[n/2] - 20*Foam::log10(0.25/pRef)) < 1e-9);
    }

    SECTION("Variable time step")
    {
        // time step varies by +-30%
        scalarList times(2*n);
        times[0] = t0;

        for (label i = 1; i < times.size(); ++i)
        {
            times[i] = times[i-1] + dt*(1 + 0.3*Foam::sin(1.7*i));
        }

        // tone on the bin 5 of the uniform grid with the mean time step
        const scalar dtMean = (times.last() - times.first())/(times.size() - 1);
        const scalar f = 5/(times.size()*dtMean);

        SoundObserver obs("R-A", Zero, pRef, times.size());

        for (const scalar t : times)
        {
            obs.append(t, Foam::cos(twoPi*f*(t - t0)));
        }

        scalarField freq, amp, spl;

        obs.fft("linear", freq, amp, spl);

        REQUIRE(mag(freq[5] - f) < 1e-9*f);
        REQUIRE(mag(amp[5] - 1) < 1e-3);
        REQUIRE(maxOtherAmp(amp, labelList({5})) < 1e-3);

        obs.fft("cubicSpline", freq, amp, spl);

        REQUIRE(mag(amp[5] - 1) < 1e-5);
        REQUIRE(maxOtherAmp(amp, labelList({5})) < 1e-5);
    }

    SECTION("Not enough samples")
    {
        SoundObserver obs("R-A", Zero, pRef, 1);
        scalarField freq, amp, spl;

        obs.fft("linear", freq, amp, spl);
        REQUIRE(freq.empty());

        obs.append(t0, 1);
        obs.fft("linear", freq, amp, spl);
        REQUIRE(freq.empty());
        REQUIRE(amp.empty());
        REQUIRE(spl.empty());
    }
}

// ************************************************************************* //
