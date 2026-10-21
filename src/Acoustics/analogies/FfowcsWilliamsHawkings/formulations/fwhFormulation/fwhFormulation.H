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
    Foam::functionObjects::fwhFormulation

Description
    Abstract base of the Ffowcs Williams - Hawkings formulations.

    Keeps the time histories of the thickness (qds) and loading (fds)
    sources of every observer, control surface and face at the observer
    (retarded) times and sums them up at the current time.

    Formulations are selected by the formulationType entry of
    FfowcsWilliamsHawkings.

SourceFiles
    fwhFormulation.C

\*---------------------------------------------------------------------------*/

#ifndef fwhFormulation_H
#define fwhFormulation_H

#include "DynamicList.H"
#include "Pair.H"
#include "vectorField.H"
#include "autoPtr.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "tdPrimitive.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

// Forward Declarations
class FfowcsWilliamsHawkings;

/*---------------------------------------------------------------------------*\
                       Class fwhFormulation Declaration
\*---------------------------------------------------------------------------*/

class fwhFormulation
{
public:

    // Public Typedefs

        //- Time history of a source at a face: (times, values)
        typedef Pair<DynamicList<scalar>> pointTimeData;

        //- Time histories per observer, control surface and face
        typedef List<List<List<pointTimeData>>> surfaceTimeData;


protected:

    // Protected Data

        //- Analogy that owns the formulation
        const FfowcsWilliamsHawkings& fwh_;

        //- Time steps since the last clean of expired data
        label fwhProbeI_;

        //- Thickness source histories
        surfaceTimeData qds_;

        //- Loading source histories
        surfaceTimeData fds_;

        //- Observer times of the sources of the current time step
        List<List<List<scalar>>> tobs_;

        //- Distances from the faces to the observers
        List<List<List<vector>>> robs_;

        //- Effective distances from the faces to the observers
        List<List<List<scalar>>> magrobs_;

        //- Outward face normals of the control surfaces
        List<tdPrimitive<vector>> ni_;

        //- Orientation of the face normals: 1 or -1
        List<List<scalar>> nl_;

        //- Largest distance to each observer
        List<scalar> rMax_;

        //- Largest delay to each observer
        List<scalar> tauMax_;

        //- Smallest distance to each observer, for the fixed response delay
        List<scalar> rMin_;

        //- Smallest delay to each observer, for the fixed response delay
        List<scalar> tauMin_;


private:

    // Private Member Functions

        //- Allocate the data and compute the geometry
        void initialize();

        //- Distances from the faces to the observers and the delay bounds
        void calcDistances();

        //- Outward face normals
        void calcNormals();


public:

    //- Runtime type information
    TypeName("fwhFormulation");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fwhFormulation,
            analogy,
            (
                const FfowcsWilliamsHawkings& fwh
            ),
            (fwh)
        );


    // Constructors

        //- Construct for the given analogy
        explicit fwhFormulation(const FfowcsWilliamsHawkings& fwh);

        //- No copy construct
        fwhFormulation(const fwhFormulation&) = delete;

        //- No copy assignment
        void operator=(const fwhFormulation&) = delete;


    // Selectors

        //- Return the selected formulation
        static autoPtr<fwhFormulation> New
        (
            const word& type,
            const FfowcsWilliamsHawkings& fwh
        );


    //- Destructor
    virtual ~fwhFormulation() = default;


    // Member Functions

        //- Source value at the time tau, linear interpolation in time.
        //  Zero outside of the stored history
        scalar valueAt
        (
            const surfaceTimeData& data,
            label iObs,
            label iSurf,
            label iFace,
            scalar tau
        ) const;

        //- Index of the last expired sample, -1 if none
        label findExpiredIndex
        (
            const pointTimeData& timeData,
            scalar expiredTime
        ) const;

        //- History without the samples up to expiredIndex
        pointTimeData getNewPointData
        (
            const pointTimeData& timeData,
            label expiredIndex
        ) const;

        //- Thickness source histories
        const surfaceTimeData& getQdsData() const noexcept
        {
            return qds_;
        }

        //- Loading source histories
        const surfaceTimeData& getFdsData() const noexcept
        {
            return fds_;
        }

        //- Append the sources of the surface iSurf at the time ct and
        //- return its contribution to the acoustic pressure of the
        //- observer iObs
        virtual scalar observerAcousticPressure
        (
            const vectorField& Sf,
            const vectorField& uS,
            const scalarField& rhoS,
            const scalarField& pS,
            label iObs,
            label iSurf,
            scalar ct
        ) = 0;

        //- Update the geometry and the observer times
        virtual void update();

        //- Remove the sources that will not reach the observers anymore
        virtual void clearExpiredData();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
