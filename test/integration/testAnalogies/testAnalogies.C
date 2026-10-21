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
    testAnalogies

Description
    Checks of the Curle and FW-H analogies on the test mesh against the
    analytic values: time derivatives at the first and the variable time
    steps, moving surfaces and patches, re-read of the settings and the
    velocity field name.

\*---------------------------------------------------------------------------*/

#include "catch.hpp"

#include "fvMesh.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "IStringStream.H"
#include "mathematicalConstants.H"

#include "tdPrimitive.H"
#include "CurleAnalogy.H"
#include "FfowcsWilliamsHawkings.H"

using namespace Foam;
using namespace functionObjects;

// Time with the given start, step and adjustTimeStep switch
#define CREATE_TIME(start, dt, adjust)                                       \
    dictionary testDict {};                                                  \
    testDict.add<scalar>("deltaT", dt);                                      \
    testDict.add<word>("startFrom", "startTime");                            \
    testDict.add<scalar>("startTime", start);                                \
    testDict.add<scalar>("endTime", 1000);                                   \
    testDict.add<scalar>("writeInterval", 100000);                           \
    testDict.add<word>("adjustTimeStep", adjust ? "yes" : "no");             \
    Foam::Time timeObj(testDict, ".", ".");

// Kinematic pressure, velocity with the given name and the viscosity
#define CREATE_FIELDS(UName)                                                 \
    fvMesh& mesh = meshPtr();                                                \
    mesh.fvSchemes::readOpt(IOobject::MUST_READ);                            \
    mesh.fvSchemes::read();                                                  \
    volScalarField p                                                         \
    (                                                                        \
        IOobject("p", timeObj.timeName(), mesh),                             \
        mesh,                                                                \
        dimensionedScalar(dimPressure/dimDensity, Zero)                      \
    );                                                                       \
    volVectorField U                                                         \
    (                                                                        \
        IOobject(UName, timeObj.timeName(), mesh),                           \
        mesh,                                                                \
        dimensionedVector(dimVelocity, Zero)                                 \
    );                                                                       \
    IOdictionary transportProperties                                         \
    (                                                                        \
        IOobject("transportProperties", timeObj.constant(), mesh)            \
    );                                                                       \
    transportProperties.add("nu", 1e-5);                                     \
    const label capI = mesh.boundaryMesh().findPatchID("testPatch");

const word curleDict
(
    "type CurleAnalogy; patches (testPatch); CofR (0 0 0);"
    "p p; rho rhoInf; rhoInf 1; c0 1; timeStart 0; timeEnd 1000; log false;"
    "observers { A { position (10 0 0); } }"
);

const word fwhDict
(
    "type FfowcsWilliamsHawkings; patches (); CofR (0 0 0);"
    "formulationType Farassat1AFormulation; interpolationScheme cell;"
    "p p; rho rhoInf; rhoInf 1; pInf 0; c0 100;"
    "timeStart 0; timeEnd 1000; log false;"
    "surfaces ( caps { type patch; patches (testPatch); interpolate false; } );"
    "observers { A { position (0 10 0); } }"
);

dictionary makeDict(const string& s)
{
    IStringStream is(s);
    return dictionary(is);
}

// Curle pressure for the compact source at c
scalar curlePressure
(
    const vector& obs,
    const vector& c,
    const vector& F,
    const vector& dotF,
    const scalar c0
)
{
    const vector l = obs - c;
    const scalar r = mag(l);

    return (l & (dotF/c0 + F/r))/(4*constant::mathematical::pi)/r/r;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

TEST_CASE("tdPrimitive derivative at the first step", "[testAnalogies]")
{
    // constant signal starting at t = 1.1, constructed at t = 1
    tdPrimitive<scalar> f(0, 1.0);

    f.value(0) = 5;
    const scalar d1 = f.dot(1.1, 0);
    f.value(0) = 5;
    const scalar d2 = f.dot(1.2, 0);

    INFO("first derivative " << d1 << ", second " << d2);
    REQUIRE(mag(d1) < SMALL);
    REQUIRE(mag(d2) < SMALL);
}


TEST_CASE("Curle pressure at the first probe", "[testAnalogies]")
{
    CREATE_TIME(1.0, 0.01, false)
    #include "createTestMesh.H"
    CREATE_FIELDS("U")

    CurleAnalogy curle("curle", timeObj, makeDict(curleDict));

    // constant force
    p.boundaryFieldRef()[capI][0] = 2;

    ++timeObj;
    curle.write();

    const scalar got = curle.getSoundObservers()[0].apressure();
    const scalar expected =
        curlePressure(vector(10, 0, 0), Zero, curle.forceEff(), Zero, 1);

    INFO("Curle " << got << ", expected " << expected);
    REQUIRE(mag(got - expected) < 1e-6*mag(expected));
}


TEST_CASE("Curle force derivative with variable time step", "[testAnalogies]")
{
    CREATE_TIME(1.0, 0.01, true)
    #include "createTestMesh.H"
    CREATE_FIELDS("U")

    REQUIRE(timeObj.isAdjustTimeStep());

    CurleAnalogy curle("curle", timeObj, makeDict(curleDict));

    // force F ~ t^2, exact derivative 2F/t; steps 0.02 and 0.01 alternate
    scalar maxErr = 0;

    for (label k = 0; k < 40; ++k)
    {
        timeObj.setDeltaT(k % 2 ? 0.01 : 0.02);
        ++timeObj;
        const scalar t = timeObj.value();

        p.boundaryFieldRef()[capI][0] = t*t;
        curle.write();

        const vector F = curle.forceEff();
        const scalar got = curle.getSoundObservers()[0].apressure();
        const scalar expected =
            curlePressure(vector(10, 0, 0), Zero, F, 2*F/t, 1);

        if (k >= 3)
        {
            maxErr = max(maxErr, mag(got - expected)/mag(expected));
        }
    }

    // first order difference: error h/(2t) < 1%
    INFO("max relative error " << maxErr);
    REQUIRE(maxErr < 0.02);
}


TEST_CASE("Velocity of the moving FW-H surface", "[testAnalogies]")
{
    CREATE_TIME(0.0, 0.1, false)
    #include "createTestMesh.H"
    CREATE_FIELDS("U")

    dictionary dict(makeDict(fwhDict));
    dict.add("nonUniformSurfaceMotion", true);

    FfowcsWilliamsHawkings fwh("fwh", timeObj, dict);

    const vector Umesh(0.5, 0, 0);

    for (label k = 0; k < 5; ++k)
    {
        ++timeObj;
        pointField newPoints(mesh.points() + Umesh*timeObj.deltaTValue());
        mesh.movePoints(newPoints);

        fwh.write();

        const vectorField& vS = fwh.surfaceVelocities()[0];

        INFO("step " << k << ": surface velocity x " << vS[0].x());
        REQUIRE(max(mag(vS - Umesh)) < 1e-9);
    }
}


TEST_CASE("Re-read of the FW-H settings", "[testAnalogies]")
{
    CREATE_TIME(0.0, 0.1, false)
    #include "createTestMesh.H"
    CREATE_FIELDS("U")

    const dictionary dict(makeDict(fwhDict));

    FfowcsWilliamsHawkings fwh("fwh", timeObj, dict);

    ++timeObj;
    fwh.write();
    ++timeObj;
    fwh.write();

    // dictionary changed: e.g. timeEnd is extended during the run
    dictionary dict2(dict);
    dict2.set("timeEnd", 2000);
    fwh.read(dict2);

    REQUIRE(fwh.getSoundObservers().size() == 1);

    ++timeObj;
    fwh.write();

    // history is kept
    REQUIRE(fwh.getSoundObservers()[0].pressure().size() == 3);
}


TEST_CASE("FW-H velocity field name", "[testAnalogies]")
{
    CREATE_TIME(0.0, 0.1, false)
    #include "createTestMesh.H"
    CREATE_FIELDS("Uother")

    dictionary dict(makeDict(fwhDict));
    dict.add("U", "Uother");

    FatalError.throwExceptions();

    FfowcsWilliamsHawkings fwh("fwh", timeObj, dict);

    ++timeObj;
    REQUIRE_NOTHROW(fwh.write());
}


TEST_CASE("Curle source centre on the moving mesh", "[testAnalogies]")
{
    CREATE_TIME(1.0, 0.1, false)
    #include "createTestMesh.H"
    CREATE_FIELDS("U")

    CurleAnalogy curle("curle", timeObj, makeDict(curleDict));

    // constant force, the body moves along x
    p.boundaryFieldRef()[capI][0] = 2;
    const vector d(0.5, 0, 0);

    for (label k = 0; k < 5; ++k)
    {
        ++timeObj;
        pointField newPoints(mesh.points() + d);
        mesh.movePoints(newPoints);

        curle.write();

        if (k >= 2)
        {
            const vector c = (k + 1)*d;
            const scalar got = curle.getSoundObservers()[0].apressure();
            const scalar expected =
                curlePressure(vector(10, 0, 0), c, curle.forceEff(), Zero, 1);

            INFO("step " << k << ": Curle " << got << ", expected " << expected);
            REQUIRE(mag(got - expected) < 1e-6*mag(expected));
        }
    }
}


// ************************************************************************* //
