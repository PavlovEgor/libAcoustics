/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    surfaceNoise

Description
    Utility to perform noise analysis of pressure data in points of surface using the libFFT
    library.

    Control settings are read from the $FOAM_CASE/system/surfaceNoiseDict dictionary,
    or user-specified dictionary using the -dict option. Pressure data prepared by
    surfacePressureSampling utility and read using a simple file reader:

    \heading Usage

    \verbatim

    outputFormat    gmsh; // optionally: gmsh (default), none (Spectrum.dat only)

    inputFileName   "fileName.txt";

    maxFrequency    1000;

    interpolationScheme linear; // optionally: linear (default), cubicSpline

    \endverbatim

    Pressure is interpolated to the uniform grid with the mean sampled time
    step, so frequencies above maxFrequency are not aliased into the
    spectrum. maxFrequency above the Nyquist frequency of the sampled data
    is truncated.

    Output:
    - FFT of the pressure data

SeeAlso
    Interpolator.H
    amplitudeSpectrum.H

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "Pair.H"
#include "OFstream.H"

#include "FileInterface/FileInterface.H"
#include "Interpolator.H"
#include "amplitudeSpectrum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions(true, true);
    argList::noParallel();

    #include "addDictOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    #include "createFields.H"

    Info<< "Reading data file..." << endl;

    // paths are relative to the case directory (-case option)
    fileName pFileName = runTime.path()/"pressureData"/inputFileName;

    FileInterface reader(pFileName);

    // sampled times and pressure per point: pSampled[pointI][timeI]
    scalarList tSampled;
    List<scalarList> pSampled;

    reader.read(tSampled, pSampled);

    //number of sampled (init) times
    label sampledTimePoints = tSampled.size();

    if (sampledTimePoints < 2)
    {
        FatalErrorInFunction
            << "At least 2 time samples are required, got "
            << sampledTimePoints << exit(FatalError);
    }

    //------------------------
    // uniform grid for FFT with the mean sampled time step: resampling
    // keeps the sampled content, frequencies above maxFreq are dropped
    // from the spectrum instead of being aliased into it.
    // FFT size is reduced to small prime factors by dropping the first
    // samples, so for the constant time step the grid stays on the samples

    Info << "Prepare to fft..." << endl;

    amplitudeSpectrum spectrum(tSampled);

    label nFftTimePoints = spectrum.size();
    label startI = spectrum.start();

    scalar startTime = spectrum.t0();

    // duration of the analysed signal
    scalar lastTime = spectrum.period();

    if (startI > 0)
    {
        Info<< "    Dropped first " << startI << " of " << sampledTimePoints
            << " samples (" << 100*(startTime - tSampled[0])
              /(tSampled[sampledTimePoints-1] - tSampled[0])
            << "% of the record) to get FFT size " << nFftTimePoints
            << endl;
    }

    scalar deltaIntT = spectrum.dt();

    // frequencies i/lastTime up to maxFreq, but not above the Nyquist one
    label freqPoints = label(maxFreq*lastTime) + 1;
    const label nyquistPoints = spectrum.nFrequencies();

    if (freqPoints > nyquistPoints)
    {
        WarningInFunction
            << "maxFrequency " << maxFreq
            << " is above the Nyquist frequency " << 0.5/deltaIntT
            << " of the sampled data, spectrum is truncated" << endl;

        freqPoints = nyquistPoints;
    }

    Info<< "    FFT points: " << nFftTimePoints
        << ", time step: " << deltaIntT
        << ", frequency step: " << 1/lastTime
        << ", frequencies: " << freqPoints << endl;

    autoPtr<Interpolator> interpolator
    (
        Interpolator::New
        (
            interpolationScheme,
            tSampled,
            startTime,
            deltaIntT,
            nFftTimePoints
        )
    );

    List<scalar> pInt(nFftTimePoints); // interpolated pressure in one point


    //------------------------
    //fft for all points on surface point-by-point, one plan for all points

    List<scalar> freq(freqPoints,0.0);

    forAll(freq,i)
    {
        freq[i] = spectrum.frequency(i);
    }

    // buffer for all complex pressures: allComplexPressures[freqI][pointI]
    List<List<complex>> allComplexPressures
    (
        freqPoints,
        List<complex>(pSampled.size())
    );

    Info << "Calculate interpolation and FFT..." << endl;

    forAll(pSampled, pointI)
    {
        interpolator->interpolate(pSampled[pointI], pInt);

        spectrum.transform(pInt);

        //write to buffer
        forAll(allComplexPressures, freqI)
        {
            allComplexPressures[freqI][pointI] = spectrum.amplitude(freqI);
        }
    }

    //------------------------

    Info<< "Write gmsh files..." << endl;

    fileName outputDirectory =
        runTime.path()/"complexPressureData"/inputFileName.lessExt();

    // directory must exist before the spectrum file is opened
    mkDir(outputDirectory);

    fileName spectrumFileName = outputDirectory/"Spectrum.dat";

    OFstream spectrumFile(spectrumFileName);

    if (!spectrumFile.good())
    {
        FatalErrorInFunction
            << "Cannot open file " << spectrumFileName << nl
            << exit(FatalError);
    }

    forAll(freq, freqI)
    {
        //write data to msh file
        if (outputFormat == "gmsh")
        {
            fileName freqDirectory = outputDirectory/name(freq[freqI]);

            mkDir(freqDirectory);

            fileName dataFile = freqDirectory/"nodeData.msh";

            FileInterface outputFile(dataFile);

            outputFile.write(freq[freqI],allComplexPressures[freqI]);
        }

        Info << "---------------\n";
        Info << "Frequency: " << freq[freqI] << nl;

        scalar maxAbs = 0.0;

        forAll(allComplexPressures[freqI],k)
        {
            complex currentCA = allComplexPressures[freqI][k];
            scalar currentAbs = Foam::sqrt( pow(currentCA.Re(),2) + pow(currentCA.Im(),2) );

            if (maxAbs < currentAbs)
            {
                maxAbs = currentAbs;
            }
        }

        Info << "Max abs: " << maxAbs << nl;

        spectrumFile << freq[freqI] << ' ' << maxAbs << nl;
    }

    return 0;
}


// ************************************************************************* //
