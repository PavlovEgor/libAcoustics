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

Class
    Foam::functionObjects::Farassat1AFormulation

Description
    Farassat formulation 1A of the Ffowcs Williams - Hawkings equation for
    the observers at rest.

    \verbatim
    formulationType Farassat1AFormulation;
    \endverbatim

SourceFiles
    Farassat1AFormulation.C

\*---------------------------------------------------------------------------*/

#ifndef Farassat1AFormulation_H
#define Farassat1AFormulation_H

#include "fwhFormulation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class Farassat1AFormulation Declaration
\*---------------------------------------------------------------------------*/

class Farassat1AFormulation
:
    public fwhFormulation
{
    // Private Data

        //- Normal velocity per observer and surface
        List<List<tdPrimitive<scalar>>> Un_;

        //- Loading vector per observer and surface
        List<List<tdPrimitive<vector>>> L_;

        //- Surface Mach number per observer and surface
        List<List<tdPrimitive<vector>>> M_;

        //- Thickness integral per observer
        scalarList intDotQdS_;

        //- Loading integral per observer
        scalarList intFdS_;


    // Private Member Functions

        //- Append the sources of the surface iSurf at the time ct
        void calculateAcousticPressure
        (
            const vectorField& Sf,
            const vectorField& uS,
            const scalarField& rhoS,
            const scalarField& pS,
            label iObs,
            label iSurf,
            scalar ct
        );


public:

    //- Runtime type information
    TypeName("Farassat1AFormulation");


    // Constructors

        //- Construct for the given analogy
        explicit Farassat1AFormulation(const FfowcsWilliamsHawkings& fwh);


    //- Destructor
    virtual ~Farassat1AFormulation() = default;


    // Member Functions

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
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
