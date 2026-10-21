    // system/surfaceNoiseDict or the file given by -dict
    const word dictName("surfaceNoiseDict");

    #include "setSystemRunTimeDictionaryIO.H"

    IOdictionary dict(dictIO);

    // upper frequency of frequency band
    scalar maxFreq = dict.lookupOrDefault("maxFrequency", 10000.0);

    // output format
    word outputFormat = dict.lookupOrDefault<word>("outputFormat", "gmsh");

    // output directory (similar with surface name)
    fileName inputFileName = dict.get<word>("inputFileName");

    // interpolation scheme to the uniform FFT grid: linear, cubicSpline
    word interpolationScheme =
        dict.getOrDefault<word>("interpolationScheme", "linear");
