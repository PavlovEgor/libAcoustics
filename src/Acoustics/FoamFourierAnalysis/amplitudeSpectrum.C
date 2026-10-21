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

#include "amplitudeSpectrum.H"
#include "error.H"

#include <algorithm>

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::label Foam::amplitudeSpectrum::fftSize(label n)
{
    for (; n > 1; --n)
    {
        label m = n;

        for (const label p : {2, 3, 5, 7})
        {
            while (m % p == 0)
            {
                m /= p;
            }
        }

        if (m == 1)
        {
            break;
        }
    }

    return n;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::amplitudeSpectrum::amplitudeSpectrum(const UList<scalar>& times)
:
    start_(0),
    t0_(0),
    period_(0),
    dt_(0),
    n_(0),
    in_(nullptr),
    out_(nullptr),
    plan_(nullptr)
{
    const label nSamples = times.size();

    if (nSamples < 2)
    {
        FatalErrorInFunction
            << "At least 2 time samples are required, got " << nSamples
            << exit(FatalError);
    }

    n_ = fftSize(nSamples - 1);
    start_ = nSamples - 1 - n_;
    t0_ = times[start_];

    period_ = times.back() - t0_;

    if (period_ <= 0)
    {
        FatalErrorInFunction
            << "Sampled times must be ascending, got duration " << period_
            << exit(FatalError);
    }

    dt_ = period_/n_;

    in_ = fftw_alloc_real(n_);
    out_ = fftw_alloc_complex(nFrequencies());

    plan_ = fftw_plan_dft_r2c_1d(n_, in_, out_, FFTW_ESTIMATE);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::amplitudeSpectrum::~amplitudeSpectrum()
{
    fftw_destroy_plan(plan_);
    fftw_free(in_);
    fftw_free(out_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::amplitudeSpectrum::transform(const UList<scalar>& values)
{
    if (values.size() != n_)
    {
        FatalErrorInFunction
            << "Expected " << n_ << " values on the grid, got "
            << values.size() << exit(FatalError);
    }

    std::copy(values.cbegin(), values.cend(), in_);

    fftw_execute(plan_);
}


Foam::complex Foam::amplitudeSpectrum::amplitude(const label k) const
{
    const scalar norm = (k == 0 || 2*k == n_) ? 1.0/n_ : 2.0/n_;

    return complex(out_[k][0], out_[k][1])*norm;
}


Foam::scalar Foam::amplitudeSpectrum::rms(const label k) const
{
    const scalar amp = mag(amplitude(k));

    return (k == 0 || 2*k == n_) ? amp : amp/Foam::sqrt(2.0);
}


// ************************************************************************* //
