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

    outputFormat    gmsh; //optionally

    inputFileName   "fileName.txt";

    maxFrequency    1000;

    \endverbatim


    Output:
    - FFT of the pressure data

SeeAlso
    FoamFftwDriver.H

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "Pair.H"
#include "OFstream.H"

#include "FileInterface/FileInterface.H"

#include <fftw3.h>
#include <algorithm>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// linear interpolation from ascending tSampled to the uniform grid t0 + j*dt:
// value[j] = data[lower[j]] + weight[j]*(data[lower[j]+1] - data[lower[j]]);
// the grid is the same for all surface points, so weights are computed once
void uniformInterpolationWeights
(
    const UList<scalar>& tSampled,
    const scalar t0,
    const scalar dt,
    labelList& lower,
    scalarList& weight
)
{
    const label lastI = tSampled.size() - 2;

    label i = 0;

    forAll(lower, j)
    {
        const scalar t = t0 + j*dt;

        while (i < lastI && tSampled[i+1] < t)
        {
            ++i;
        }

        const scalar dtSampled = tSampled[i+1] - tSampled[i];

        lower[j] = i;
        weight[j] =
        (
            dtSampled > 0
          ? min(max((t - tSampled[i])/dtSampled, scalar(0)), scalar(1))
          : scalar(0)
        );
    }
}


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

    fileName pFileName = "pressureData/" + inputFileName;

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

    // sampled signal duration (sampling may start at non-zero time)
    scalar lastTime = tSampled[sampledTimePoints-1] - tSampled[0];

    if (lastTime <= 0)
    {
        FatalErrorInFunction
            << "Sampled times must be ascending, got duration " << lastTime
            << exit(FatalError);
    }

    //------------------------
    // uniform grid for FFT with the mean sampled time step: resampling
    // keeps the sampled content, frequencies above maxFreq are dropped
    // from the spectrum instead of being aliased into it

    Info << "Prepare to fft..." << endl;

    label nFftTimePoints = sampledTimePoints - 1;
    scalar deltaIntT = lastTime/scalar(nFftTimePoints);

    // frequencies i/lastTime up to maxFreq, but not above the Nyquist one
    label freqPoints = label(maxFreq*lastTime) + 1;
    const label nyquistPoints = nFftTimePoints/2 + 1;

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

    labelList intLower(nFftTimePoints);
    scalarList intWeight(nFftTimePoints);

    uniformInterpolationWeights(tSampled, tSampled[0], deltaIntT, intLower, intWeight);

    List<scalar> pInt(nFftTimePoints); // interpolated pressure in one point

    // real-to-complex FFT: one plan for all points, the buffers are reused
    double* fftIn = fftw_alloc_real(nFftTimePoints);
    fftw_complex* fftOut = fftw_alloc_complex(nyquistPoints);

    fftw_plan fftPlan =
        fftw_plan_dft_r2c_1d(nFftTimePoints, fftIn, fftOut, FFTW_ESTIMATE);

    // amplitude normalisation: 2/N, but 1/N at zero and Nyquist frequencies
    scalarList ampNorm(freqPoints, 2.0/nFftTimePoints);

    ampNorm[0] = 1.0/nFftTimePoints;

    if (nFftTimePoints % 2 == 0 && freqPoints == nyquistPoints)
    {
        ampNorm[freqPoints-1] = 1.0/nFftTimePoints;
    }


    //------------------------
    //fft for all points on surface point-by-point

    List<scalar> freq(freqPoints,0.0);

    forAll(freq,i)
    {
        freq[i] = scalar(i)/lastTime;
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
        const scalarList& p = pSampled[pointI];

        forAll(pInt, j)
        {
            const scalar p0 = p[intLower[j]];

            pInt[j] = p0 + intWeight[j]*(p[intLower[j]+1] - p0);
        }

        std::copy(pInt.cbegin(), pInt.cend(), fftIn);

        fftw_execute(fftPlan);

        //write to buffer
        forAll(allComplexPressures, freqI)
        {
            allComplexPressures[freqI][pointI] =
                complex(fftOut[freqI][0], fftOut[freqI][1])*ampNorm[freqI];
        }
    }

    fftw_destroy_plan(fftPlan);
    fftw_free(fftIn);
    fftw_free(fftOut);

    //------------------------

    Info<< "Write gmsh files..." << endl;

    fileName outputDirectory = inputFileName.lessExt();    

    fileName controlFile = "complexPressureData/" + outputDirectory + "/" + "Spectrum.dat";

    autoPtr<OFstream> controlFilePtr;

    if (!bool(controlFilePtr))
    {
        // Open new file at start up
        controlFilePtr.reset
        (
            new OFstream
            (
                controlFile
            )
        );
    }

    forAll(freq, freqI)
    {
        word freqValue = name(freq[freqI]);

        fileName freqDirectory = "complexPressureData/" + outputDirectory + "/" + freqValue;

        mkDir(freqDirectory);

        //write data to msh file
        if (outputFormat == "gmsh")
        {
            fileName dataFile = freqDirectory + "/nodeData.msh";

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

        controlFilePtr() << freq[freqI] << ' ' << maxAbs << nl;
    }

    return 0;
}


// ************************************************************************* //
