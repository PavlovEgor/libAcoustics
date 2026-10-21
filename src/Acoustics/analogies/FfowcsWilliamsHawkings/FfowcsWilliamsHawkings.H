/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::functionObjects::FfowcsWilliamsHawkings

Group
    grpFfowcsWilliamsHawkingsFunctionObjects

Description
    Ffowcs Williams - Hawkings acoustic analogy: acoustic pressure of the
    observers from the flow on the control surfaces.

    The acoustic pressure is the sum of the contributions of all control
    surfaces, so the surfaces are parts of one integration surface.

    \verbatim
    fwh
    {
        type                FfowcsWilliamsHawkings;
        libs                ("libAcoustics.so");

        // AcousticAnalogy settings: timeStart, c0, observers, ...

        patches             ();
        formulationType     Farassat1AFormulation;  // or GTFormulation
        interpolationScheme cellPoint;
        pInf                0;
        U0                  (0 0 0);    // optional, uniform flow for GT
        Ufwh                (0 0 0);    // optional, surface velocity
        nonUniformSurfaceMotion false;  // optional
        fixedResponseDelay  false;      // optional
        responseDelay       0.0001;     // optional, with fixedResponseDelay
        cleanFreq           100;        // optional

        surfaces
        (
            sphere
            {
                type            patch;
                patches         (sphere);
                interpolate     false;
            }
        );
    }
    \endverbatim

SourceFiles
    FfowcsWilliamsHawkings.C
    FfowcsWilliamsHawkingsTemplates.C

SeeAlso
    Foam::functionObjects::AcousticAnalogy
    Foam::functionObjects::fwhFormulation

\*---------------------------------------------------------------------------*/

#ifndef FfowcsWilliamsHawkings_H
#define FfowcsWilliamsHawkings_H

#include "AcousticAnalogy.H"
#include "sampledSurface.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

// Forward Declarations
class fwhFormulation;

/*---------------------------------------------------------------------------*\
                   Class FfowcsWilliamsHawkings Declaration
\*---------------------------------------------------------------------------*/

class FfowcsWilliamsHawkings
:
    public AcousticAnalogy
{
    // Private Data

        //- Name of the FW-H formulation
        word formulationType_;

        //- Fixed delay of the response
        Switch fixedResponseDelay_;

        //- Response delay with the fixed delay
        scalar responseDelay_;

        //- Velocity of the control surfaces
        vector Ufwh_;

        //- Uniform flow velocity, e.g. a wind tunnel, for GTFormulation
        vector U0_;

        //- Surfaces move with the mesh
        bool nonUniformSurfaceMotion_;

        //- Face centres at the previous time step, for the moving surfaces
        List<vectorField> Cf0_;

        //- Face velocities of the control surfaces
        List<vectorField> vS_;

        //- Undisturbed pressure, in Pa
        scalar pInf_;

        //- FW-H formulation
        autoPtr<fwhFormulation> fwhFormulationPtr_;

        //- Time steps between cleaning of the expired sources
        label cleanFreq_;

        //- Interpolation scheme to use for sampling
        word interpolationScheme_;

        //- Control surfaces
        PtrList<sampledSurface> controlSurfaces_;

        //- Dictionary entry of the control surfaces, to detect its change
        string surfacesDescription_;


protected:

    // Protected Member Functions

        //- Compute and append the acoustic pressure of the observers
        virtual void correct();

        //- Update the moving surfaces at the time steps between probes
        virtual void updateSurfaceCf();

        //- Sample or interpolate a volume field on the surface
        template<class Type>
        tmp<Field<Type>> sampleOrInterpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const sampledSurface& surface
        ) const;

        //- Sampled density
        tmp<scalarField> surfaceDensity(const sampledSurface& surface) const;

        //- Sampled velocity
        tmp<vectorField> surfaceVelocity(const sampledSurface& surface) const;

        //- Sampled pressure
        tmp<scalarField> surfacePressure(const sampledSurface& surface) const;

        //- Allocate the surface data and select the formulation
        void initialize();


public:

    //- Runtime type information
    TypeName("FfowcsWilliamsHawkings");


    // Constructors

        //- Construct from given Time and dictionary.
        FfowcsWilliamsHawkings
        (
            const word& name,
            const Time&,
            const dictionary&
        );

        //- Construct from given objectRegistry and dictionary.
        FfowcsWilliamsHawkings
        (
            const word& name,
            const objectRegistry&,
            const dictionary&
        );

        //- No copy construct
        FfowcsWilliamsHawkings(const FfowcsWilliamsHawkings&) = delete;

        //- No copy assignment
        void operator=(const FfowcsWilliamsHawkings&) = delete;


    //- Destructor
    virtual ~FfowcsWilliamsHawkings();


    // Member Functions

        // Access

            //- Control surfaces
            const PtrList<sampledSurface>& controlSurfaces() const noexcept
            {
                return controlSurfaces_;
            }

            //- Control surfaces
            PtrList<sampledSurface>& getSampledSurface() noexcept
            {
                return controlSurfaces_;
            }

            //- Face velocities of the control surfaces
            const List<vectorField>& surfaceVelocities() const noexcept
            {
                return vS_;
            }

            //- Velocity of the control surfaces
            const vector& Ufwh() const noexcept
            {
                return Ufwh_;
            }

            //- Uniform flow velocity
            const vector& U0() const noexcept
            {
                return U0_;
            }

            //- Surfaces move with the mesh
            bool nonUniformSurfaceMotion() const noexcept
            {
                return nonUniformSurfaceMotion_;
            }

            //- Fixed delay of the response
            bool fixedResponseDelay() const noexcept
            {
                return fixedResponseDelay_;
            }

            //- Response delay with the fixed delay
            scalar responseDelay() const noexcept
            {
                return responseDelay_;
            }

            //- Time steps between cleaning of the expired sources
            label cleanFreq() const noexcept
            {
                return cleanFreq_;
            }


        //- Read the FfowcsWilliamsHawkings data
        virtual bool read(const dictionary&);

        //- Mark the control surfaces as needing an update.
        //  Return false if all surfaces were already marked as expired.
        virtual bool expire();

        //- True if any control surface needs an update
        virtual bool needsUpdate() const;

        //- Update the control surfaces, true if any was updated
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FfowcsWilliamsHawkingsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
