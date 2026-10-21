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

Class
    Foam::amplitudeSpectrum

Description
    Amplitude spectrum of a sampled signal by the real-to-complex FFTW
    transform.

    The uniform FFT grid covers the end of the sampled record:
    \verbatim
        t0 + i*dt,  i = 0 .. size()-1,  t0 = times[start()],
        dt = (times.last() - t0)/size()
    \endverbatim
    The FFT size is the largest number not above times.size() - 1 with
    prime factors 2, 3, 5, 7 only: FFTW is up to ten times faster for such
    sizes than for sizes with large prime factors. The first start()
    samples are dropped for that. For the constant time step the grid
    points are the samples start() .. start() + size() - 1.

    One plan is used for any number of transforms of the same size.

SourceFiles
    amplitudeSpectrum.C

SeeAlso
    Foam::Interpolator

\*---------------------------------------------------------------------------*/

#ifndef amplitudeSpectrum_H
#define amplitudeSpectrum_H

#include "scalarList.H"
#include "complex.H"

#include <fftw3.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class amplitudeSpectrum Declaration
\*---------------------------------------------------------------------------*/

class amplitudeSpectrum
{
    // Private Data

        //- Index of the first sample on the grid
        label start_;

        //- Start time of the grid
        scalar t0_;

        //- Period of the analysed signal: times.last() - t0
        scalar period_;

        //- Time step of the grid
        scalar dt_;

        //- FFT size
        label n_;

        //- FFT input
        double* in_;

        //- FFT output, n/2 + 1 frequencies up to the Nyquist one
        fftw_complex* out_;

        //- FFT plan
        fftw_plan plan_;


public:

    // Constructors

        //- Construct the grid and the plan for the sampled times.
        //  The times must be ascending, at least 2 samples
        explicit amplitudeSpectrum(const UList<scalar>& times);

        //- No copy construct
        amplitudeSpectrum(const amplitudeSpectrum&) = delete;

        //- No copy assignment
        void operator=(const amplitudeSpectrum&) = delete;


    //- Destructor
    ~amplitudeSpectrum();


    // Static Member Functions

        //- Largest size not above n with prime factors 2, 3, 5, 7 only
        static label fftSize(label n);


    // Member Functions

        // Grid

            //- Index of the first sample on the grid
            label start() const noexcept
            {
                return start_;
            }

            //- Start time of the grid
            scalar t0() const noexcept
            {
                return t0_;
            }

            //- Time step of the grid
            scalar dt() const noexcept
            {
                return dt_;
            }

            //- FFT size: number of grid points
            label size() const noexcept
            {
                return n_;
            }

            //- Period of the analysed signal
            scalar period() const noexcept
            {
                return period_;
            }


        // Spectrum

            //- Number of frequencies up to the Nyquist one
            label nFrequencies() const noexcept
            {
                return n_/2 + 1;
            }

            //- Frequency k
            scalar frequency(const label k) const
            {
                return k/period_;
            }

            //- Transform values on the grid, size() elements
            void transform(const UList<scalar>& values);

            //- Complex amplitude at frequency k of the last transform:
            //  2X/n, but X/n at zero and Nyquist frequencies
            complex amplitude(const label k) const;

            //- Rms value at frequency k of the last transform:
            //  amplitude/sqrt(2), but the amplitude itself at zero and
            //  Nyquist frequencies
            scalar rms(const label k) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
