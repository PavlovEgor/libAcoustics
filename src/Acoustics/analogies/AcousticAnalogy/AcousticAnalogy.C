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

#include "AcousticAnalogy.H"
#include "volFields.H"
#include "Time.H"
#include "Interpolator.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(AcousticAnalogy, 0);
}
}

Foam::functionObjects::AcousticAnalogy::AcousticAnalogy
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    forces
    (
        name,
        runTime,
        dict
    ),
    analogyOutPtr_(nullptr),
    probeFreq_(1024),
    timeStart_(-1.0),
    timeEnd_(-1.0),
    writeFft_(true),
    fftInterpolationScheme_("linear"),
    c0_(300.0),
    observers_(0),
    observersChanged_(true),
    probeI_(0)
{
}

Foam::functionObjects::AcousticAnalogy::AcousticAnalogy
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    forces
    (
        name,
        obr,
        dict
    ),
    analogyOutPtr_(nullptr),
    probeFreq_(1024),
    timeStart_(-1.0),
    timeEnd_(-1.0),
    writeFft_(true),
    fftInterpolationScheme_("linear"),
    c0_(300.0),
    observers_(0),
    observersChanged_(true),
    probeI_(0)
{
}

void Foam::functionObjects::AcousticAnalogy::writeFileHeader()
{
    analogyOutPtr_() << "Time" << " ";
    forAll(observers_, iObserver)
    {
        analogyOutPtr_() << observers_[iObserver].name() << "_pFluct ";
    }
    analogyOutPtr_() << endl;
}

void Foam::functionObjects::AcousticAnalogy::makeFile()
{
    // the file is open on the master only
    if (analogyOutPtr_)
    {
        // observers are changed by re-read: header of the new columns
        if (observersChanged_)
        {
            writeFileHeader();
        }

        return;
    }

    fileName ResultsDir;

    if (Pstream::master() && Pstream::parRun())
    {
        ResultsDir = obr_.time().rootPath() + "/" + obr_.time().caseName().path()  + "/acousticData";
        mkDir(ResultsDir);
    }
    else if (!Pstream::parRun())
    {
        ResultsDir = obr_.time().rootPath() + "/" + obr_.time().caseName() + "/acousticData";
        mkDir(ResultsDir);
    }
    else
    {
    }

    // File update
    if (Pstream::master() || !Pstream::parRun())
    {

        analogyOutPtr_.reset
        (
            new OFstream
            (
                ResultsDir + "/" + (name() + "-time.dat")
            )
        );

        writeFileHeader();
    }
}

void Foam::functionObjects::AcousticAnalogy::writeFft()
{
    fileName ResultsDir;

    if (Pstream::master() && Pstream::parRun())
    {
        ResultsDir = obr_.time().rootPath() + "/" + obr_.time().caseName().path()  + "/acousticData";
    }
    else if (!Pstream::parRun())
    {
        ResultsDir = obr_.time().rootPath() + "/" + obr_.time().caseName() + "/acousticData";
    }

    if (Pstream::master() || !Pstream::parRun())
    {
        scalarField freq;
        scalarField amp;
        scalarField spl;

        for (const SoundObserver& obs : observers_)
        {
            if (!obs.fftDue())
            {
                continue;
            }

            obs.fft(fftInterpolationScheme_, freq, amp, spl);

            if (freq.empty())
            {
                continue;
            }

            Log << "Executing fft for obs: " << obs.name() << endl;

            OFstream fftStream
            (
                ResultsDir/("fft-" + name() + "-" + obs.name() + ".dat")
            );

            fftStream << "Freq p\' spl" << nl;

            forAll(freq, k)
            {
                fftStream << freq[k] << ' ' << amp[k] << ' ' << spl[k] << nl;
            }
        }
    }
}

bool Foam::functionObjects::AcousticAnalogy::read(const dictionary& dict)
{
    if (!forces::read(dict))
    {
        return false;
    }
    Info << "Reading analogy settings" << endl;

    probeFreq_ = dict.getOrDefault<label>("probeFrequency", 1);

    dict.readEntry("timeStart", timeStart_);

    dict.readEntry("timeEnd", timeEnd_);

    writeFft_ = dict.getOrDefault("writeFft", false);

    fftInterpolationScheme_ =
        dict.getOrDefault<word>("fftInterpolationScheme", "linear");

    if (!Interpolator::timesConstructorTable(fftInterpolationScheme_))
    {
        FatalIOErrorInLookup
        (
            dict,
            "fftInterpolationScheme",
            fftInterpolationScheme_,
            *Interpolator::timesConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    dict.readEntry("c0", c0_);

    //Ask for rhoRef again, because libforces sometimes do not
    dict.readEntry("rhoInf", rhoRef_);

    // read observers, keep them with their histories if not changed
    {
        const dictionary& obsDict = dict.subDict("observers");
        const wordList names(obsDict.toc());

        List<SoundObserver> observers(names.size());

        forAll(names, i)
        {
            const dictionary& oDict = obsDict.subDict(names[i]);

            observers[i] = SoundObserver
            (
                names[i],
                oDict.get<vector>("position"),
                oDict.getOrDefault<scalar>("pRef", 2e-5),
                oDict.getOrDefault<label>("fftFreq", 1024)
            );
        }

        observersChanged_ = (observers.size() != observers_.size());

        forAll(observers, i)
        {
            if (observersChanged_)
            {
                break;
            }

            const SoundObserver& o = observers[i];
            const SoundObserver& old = observers_[i];

            observersChanged_ =
                o.name() != old.name()
             || o.position() != old.position()
             || o.pRef() != old.pRef()
             || o.fftFreq() != old.fftFreq();
        }

        if (observersChanged_)
        {
            observers_.transfer(observers);
        }
    }

    this->makeFile();

    return true;
}

bool Foam::functionObjects::AcousticAnalogy::execute()
{
    return true;
}

bool Foam::functionObjects::AcousticAnalogy::write()
{
    //use forces library to calculate forces acting on patch
    if (!forces::write())
    {
        return false;
    }

    scalar cTime = obr_.time().value();
    scalar stTime = obr_.time().startTime().value();

    probeI_++;

    if ( mag(probeI_ % probeFreq_) > VSMALL  )
    {
        updateSurfaceCf();
        return true;
    }
    else
    {
        Log << "Starting acoustics probe" << endl;
        probeI_ = 0.0;
    }

    if ( (cTime < timeStart_) || (cTime > timeEnd_))
    {
        return true;
    }

    correct();

    if (Pstream::master() || !Pstream::parRun())
    {
        // time history output
        if (stTime > timeStart_)
        {
            analogyOutPtr_() << (cTime - stTime) << " ";
        }
        else
        {
            analogyOutPtr_() << (cTime - timeStart_) << " ";
        }

        forAll(observers_, iObserver)
        {
            const SoundObserver& obs = observers_[iObserver];
            analogyOutPtr_() << obs.apressure() << " ";
        }

        analogyOutPtr_() << endl;

        //fft output
        if (writeFft_ == true)
        {
            writeFft();
        }

        //output to stdio
        Log << "Acoustic pressure" << endl;
        forAll(observers_, iObserver)
        {
            const SoundObserver& obs = observers_[iObserver];
            Info << iObserver << " Observer: " << obs.name() << " p\' = " << obs.apressure() << endl;
        }
    }

    return true;
}

// ************************************************************************* //
