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
    Foam::functionObjects::AcousticAnalogy

Description
    Abstract base of the acoustic analogy function objects.

    Reads the common settings and the observers, calls correct() of the
    analogy to compute the acoustic pressure of the observers at the probe
    time steps and writes acousticData/<name>-time.dat and the spectra
    acousticData/fft-<name>-<observer>.dat.

    Common settings:
    \verbatim
    probeFrequency  1;          // optional, time steps between probes
    timeStart       0.001;
    timeEnd         1;
    c0              343;        // sound speed
    rhoInf          1.2;
    writeFft        false;      // optional
    fftInterpolationScheme linear;  // optional, see Foam::SoundObserver

    observers
    {
        R-A
        {
            position    (0 0 1);
        }
    }
    \endverbatim

SourceFiles
    AcousticAnalogy.C

SeeAlso
    Foam::SoundObserver
    Foam::functionObjects::CurleAnalogy
    Foam::functionObjects::FfowcsWilliamsHawkings

\*---------------------------------------------------------------------------*/

#ifndef AcousticAnalogy_H
#define AcousticAnalogy_H

#include "forces.H"
#include "OFstream.H"
#include "SoundObserver.H"
#include "tdPrimitive.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class AcousticAnalogy Declaration
\*---------------------------------------------------------------------------*/

class AcousticAnalogy
:
    public forces
{
protected:

    // Protected Data

        //- Time history output
        autoPtr<OFstream> analogyOutPtr_;

        //- Number of time steps between probes
        label probeFreq_;

        //- Start time of probing
        scalar timeStart_;

        //- End time of probing
        scalar timeEnd_;

        //- Write the observer spectra
        bool writeFft_;

        //- Interpolation of observer samples with a variable time step
        //- to the uniform FFT grid, see Foam::Interpolator
        word fftInterpolationScheme_;

        //- Sonic speed, m/s
        scalar c0_;

        //- List of observers
        List<SoundObserver> observers_;

        //- Time steps since the last probe
        label probeI_;


    // Protected Member Functions

        //- Open the time history output
        void makeFile();

        //- Write the observer spectra that are due
        void writeFft();

        //- Compute and append the acoustic pressure of the observers
        virtual void correct() = 0;

        //- Update the surfaces at the time steps between probes
        virtual void updateSurfaceCf()
        {}


public:

    //- Time derivative helper, see Foam::tdPrimitive
    template<class T>
    using tdPrimitive = Foam::tdPrimitive<T>;


    //- Runtime type information
    TypeName("AcousticAnalogy");


    // Constructors

        //- Construct from given Time and dictionary.
        AcousticAnalogy
        (
            const word& name,
            const Time&,
            const dictionary&
        );

        //- Construct from given objectRegistry and dictionary.
        AcousticAnalogy
        (
            const word& name,
            const objectRegistry&,
            const dictionary&
        );


    //- Destructor
    virtual ~AcousticAnalogy() = default;


    // Member Functions

        // Access

            //- Observers list
            const List<SoundObserver>& getSoundObservers() const noexcept
            {
                return observers_;
            }

            //- Sonic speed
            scalar c0() const noexcept
            {
                return c0_;
            }

            //- Reference density
            scalar rhoRef() const noexcept
            {
                return rhoRef_;
            }


        //- Read the analogy settings
        virtual bool read(const dictionary&);

        //- Execute, does nothing
        virtual bool execute();

        //- Probe the acoustic pressure and write it
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
