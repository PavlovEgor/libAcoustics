    // system/surfaceNoiseDict or the file given by -dict
    const word dictName("surfaceNoiseDict");

    #include "setSystemRunTimeDictionaryIO.H"

    IOdictionary dict(dictIO);

    // upper frequency of frequency band
    scalar maxFreq = dict.lookupOrDefault("maxFrequency", 10000.0);

    // output format of complex amplitudes: gmsh, none (Spectrum.dat only)
    word outputFormat = dict.getOrDefault<word>("outputFormat", "gmsh");

    if (outputFormat != "gmsh" && outputFormat != "none")
    {
        FatalIOErrorInFunction(dict)
            << "Unknown outputFormat " << outputFormat << nl
            << "Valid formats: gmsh none" << nl
            << exit(FatalIOError);
    }

    // output directory (similar with surface name)
    fileName inputFileName = dict.get<word>("inputFileName");

    // interpolation scheme to the uniform FFT grid: linear, cubicSpline
    word interpolationScheme =
        dict.getOrDefault<word>("interpolationScheme", "linear");
