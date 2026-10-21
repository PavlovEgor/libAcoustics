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

Class
    Foam::SoundObserver

Description
    Acoustic observer: time history of the acoustic pressure at a point
    and its amplitude spectrum.

    Observer entry of an acoustic analogy:
    \verbatim
    writeFft    true;

    fftInterpolationScheme linear;  // optional: linear (default), cubicSpline

    observers
    {
        R-A
        {
            position    (0 0 1);
            pRef        2e-5;   // optional, reference pressure for SPL
            fftFreq     1024;   // optional, samples between spectrum updates
        }
    }
    \endverbatim

    The spectrum is computed from the history without its first samples,
    so that the FFT size has small prime factors only (see
    Foam::amplitudeSpectrum). Samples with a variable
    time step are interpolated to the uniform grid with the mean time step
    by the fftInterpolationScheme, constant time step samples are used as
    they are.
    Spectrum columns are frequency [Hz], amplitude [Pa] and SPL [dB], the SPL
    is computed from the rms pressure.

SourceFiles
    SoundObserver.C

SeeAlso
    Foam::amplitudeSpectrum
    Foam::Interpolator

\*---------------------------------------------------------------------------*/

#ifndef SoundObserver_H
#define SoundObserver_H

#include "DynamicList.H"
#include "scalarField.H"
#include "vector.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class SoundObserver Declaration
\*---------------------------------------------------------------------------*/

class SoundObserver
{
    // Private Data

        //- Observer name
        word name_;

        //- Observer position
        vector position_;

        //- Reference pressure for the sound pressure level
        scalar pRef_;

        //- Number of samples between spectrum updates
        label fftFreq_;

        //- Sampled times
        DynamicList<scalar> times_;

        //- Sampled acoustic pressure
        DynamicList<scalar> pressure_;


public:

    // Constructors

        //- Default construct
        SoundObserver();

        //- Construct from components
        SoundObserver
        (
            const word& name,
            const vector& position,
            const scalar pRef,
            const label fftFreq
        );


    // Member Functions

        // Access

            //- Observer name
            const word& name() const noexcept
            {
                return name_;
            }

            //- Observer position
            const vector& position() const noexcept
            {
                return position_;
            }

            //- Reference pressure for the sound pressure level
            scalar pRef() const noexcept
            {
                return pRef_;
            }

            //- Number of samples between spectrum updates
            label fftFreq() const noexcept
            {
                return fftFreq_;
            }

            //- Sampled times
            const UList<scalar>& times() const noexcept
            {
                return times_;
            }

            //- Sampled acoustic pressure
            const UList<scalar>& pressure() const noexcept
            {
                return pressure_;
            }

            //- Last sampled acoustic pressure, zero before the first sample
            scalar apressure() const
            {
                return pressure_.empty() ? 0 : pressure_.back();
            }


        // Edit

            //- Set the observer name
            void name(const word& name)
            {
                name_ = name;
            }

            //- Set the observer position
            void position(const vector& position)
            {
                position_ = position;
            }

            //- Append the acoustic pressure at the given time
            void append(const scalar t, const scalar p);


        // Spectrum

            //- True when the spectrum update is due:
            //- the number of samples is a multiple of fftFreq
            bool fftDue() const;

            //- Amplitude spectrum of the pressure history.
            //  Samples with a variable time step are interpolated with
            //  the given Interpolator scheme. Empty for less than 2 samples
            void fft
            (
                const word& interpolationScheme,
                scalarField& freq,
                scalarField& amp,
                scalarField& spl
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
