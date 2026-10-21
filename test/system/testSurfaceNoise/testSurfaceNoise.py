#!/usr/bin/env python3
"""
System test of the surfaceNoise utility.

Generates synthetic pressure histories, runs surfaceNoise on them with -case
from another directory and compares the written spectra with an independent
numpy implementation of the same algorithm:

    - drop first samples to get FFT size with prime factors 2, 3, 5, 7
    - interpolate to the uniform grid with the mean sampled time step
      (numpy.interp for linear, dense natural cubic spline for cubicSpline)
    - rfft, amplitude normalisation 2/N (1/N at zero and Nyquist frequencies)
    - keep frequencies i/T <= maxFrequency

Also checks analytic tone amplitudes, absence of aliasing and error handling.

Usage: testSurfaceNoise.py [workDir]   (requires numpy and surfaceNoise in PATH,
       or the executable given by the SURFACE_NOISE environment variable)

Exit code is the number of failed checks.
"""

import os
import shutil
import subprocess
import sys

import numpy as np

EXE = os.environ.get("SURFACE_NOISE") or shutil.which("surfaceNoise")

CONTROL_DICT = """\
FoamFile { version 2.0; format ascii; class dictionary; object controlDict; }
application surfaceNoise;
startFrom startTime; startTime 0; stopAt endTime; endTime 1; deltaT 1;
writeControl timeStep; writeInterval 1; writeFormat ascii; writePrecision 15;
"""

failures = []


def check(name, ok, detail=""):
    print("    %-58s %s %s" % (name, "OK" if ok else "FAILED", detail))
    if not ok:
        failures.append(name)


# ---------------------------------------------------------------- reference

def fft_size(n):
    """Largest size <= n with prime factors 2, 3, 5, 7 only"""
    while n > 1:
        m = n
        for p in (2, 3, 5, 7):
            while m % p == 0:
                m //= p
        if m == 1:
            break
        n -= 1
    return n


def natural_spline(t, y, tNew):
    """Natural cubic spline, dense solve: independent of the utility code"""
    n = len(t)
    h = np.diff(t)
    A = np.zeros((n, n))
    rhs = np.zeros((n,) + y.shape[1:])
    A[0, 0] = A[-1, -1] = 1
    for i in range(1, n - 1):
        A[i, i - 1], A[i, i], A[i, i + 1] = h[i - 1], 2*(h[i - 1] + h[i]), h[i]
        rhs[i] = 6*((y[i + 1] - y[i])/h[i] - (y[i] - y[i - 1])/h[i - 1])
    M = np.linalg.solve(A, rhs)
    i = np.clip(np.searchsorted(t, tNew, side="right") - 1, 0, n - 2)
    u = ((tNew - t[i])/h[i])[:, None]
    hh = (h[i]**2/6)[:, None]
    return (
        (1 - u)*y[i] + u*y[i + 1]
      + hh*((1 - u)**3 - (1 - u))*M[i] + hh*(u**3 - u)*M[i + 1]
    )


def reference(t, p, scheme, maxFreq):
    """Spectrum [freqI, pointI], frequencies, duration"""
    M = len(t)
    N = fft_size(M - 1)
    s = M - 1 - N
    T = t[-1] - t[s]
    tGrid = t[s] + np.arange(N)*(T/N)

    if scheme == "linear":
        g = np.stack([np.interp(tGrid, t, p[:, k]) for k in range(p.shape[1])], 1)
    else:
        g = natural_spline(t, p, tGrid)

    F = min(int(maxFreq*T) + 1, N//2 + 1)
    X = np.fft.rfft(g, axis=0)[:F]*2/N
    X[0] /= 2
    if N % 2 == 0 and F == N//2 + 1:
        X[-1] /= 2

    return X, np.arange(F)/T, T, N, s


# ---------------------------------------------------------------- case I/O

def write_case(caseDir, rows, dictEntries, extraLines=(), blankAfter=None):
    shutil.rmtree(caseDir, ignore_errors=True)
    os.makedirs(os.path.join(caseDir, "system"))
    os.makedirs(os.path.join(caseDir, "pressureData"))

    with open(os.path.join(caseDir, "system", "controlDict"), "w") as f:
        f.write(CONTROL_DICT)

    with open(os.path.join(caseDir, "system", "surfaceNoiseDict"), "w") as f:
        f.write("FoamFile { version 2.0; format ascii; class dictionary; "
                "object surfaceNoiseDict; }\n")
        f.write('inputFileName "test.dat";\n')
        for key, value in dictEntries.items():
            f.write("%s %s;\n" % (key, value))

    with open(os.path.join(caseDir, "pressureData", "test.dat"), "w") as f:
        for i, r in enumerate(rows):
            f.write(" ".join("%.17g" % v for v in r) + "\n")
            if i == blankAfter:
                f.write("\n")
        for line in extraLines:
            f.write(line + "\n")


def run(caseDir, workDir, *args):
    """Run surfaceNoise with -case from workDir, return (exit code, log)"""
    res = subprocess.run(
        [EXE, "-case", caseDir] + list(args),
        cwd=workDir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        universal_newlines=True
    )
    with open(os.path.join(caseDir, "log.surfaceNoise"), "w") as f:
        f.write(res.stdout)
    return res.returncode, res.stdout


def read_output(caseDir, nPoints):
    """Spectrum [freqI, pointI] from nodeData.msh files and Spectrum.dat"""
    root = os.path.join(caseDir, "complexPressureData", "test")
    spectrum = np.loadtxt(os.path.join(root, "Spectrum.dat"), ndmin=2)

    X = np.zeros((len(spectrum), nPoints), complex)
    for d in os.listdir(root):
        msh = os.path.join(root, d, "nodeData.msh")
        if not os.path.isfile(msh):
            continue
        k = int(np.argmin(np.abs(spectrum[:, 0] - float(d))))
        lines = open(msh).read().split("\n")
        assert int(lines[8]) == nPoints, "wrong number of nodes in " + msh
        X[k] = [complex(l.split()[1]) for l in lines[10:10 + nPoints]]

    return X, spectrum


# ---------------------------------------------------------------- scenarios

def signals(t, freqs, amps, rng, noise):
    p = np.stack(
        [a*np.sin(2*np.pi*f*t + 0.3*k) for k, (f, a) in enumerate(zip(freqs, amps))],
        1
    )
    return p + noise*rng.random(p.shape)


def compare_with_reference(name, caseDir, t, p, scheme, maxFreq, tol=1e-9):
    X, freq, T, N, s = reference(t, p, scheme, maxFreq)
    out, spectrum = read_output(caseDir, p.shape[1])

    check(name + ": number of frequencies", len(spectrum) == len(freq),
          "(%d, expected %d)" % (len(spectrum), len(freq)))
    if len(spectrum) != len(freq):
        return X, out, T

    check(name + ": frequencies", np.allclose(spectrum[:, 0], freq, rtol=1e-9))

    err = np.abs(out - X).max()
    check(name + ": amplitudes vs numpy", err < tol, "(max diff %.1e)" % err)

    errMax = np.abs(spectrum[:, 1] - np.abs(out).max(1)).max()
    check(name + ": Spectrum.dat max abs", errMax < tol, "(max diff %.1e)" % errMax)

    return X, out, T


def scenario_uniform(work, rng):
    """Constant time step, start time > 0: interpolation must be exact"""
    print("uniform time step, start time 0.05, blank line, damaged last line")

    dt, M = 1e-4, 2500
    t = 0.05 + dt*np.arange(M)
    N = fft_size(M - 1)
    T = N*dt
    bins = [40, 100, 155, 230, 301]
    amps = [1.0, 0.5, 0.25, 0.8, 0.1]
    p = signals(t, [k/T for k in bins], amps, rng, 0.0)

    for scheme in ("linear", "cubicSpline"):
        case = os.path.join(work, "uniform_" + scheme)
        write_case(
            case, np.column_stack([t, p]),
            {"maxFrequency": 4000, "interpolationScheme": scheme},
            extraLines=["%.17g 1 2" % (t[-1] + dt)],
            blankAfter=500
        )

        rc, log = run(case, work)
        check(scheme + ": run", rc == 0, "(exit %d)" % rc)
        if rc:
            continue

        check(scheme + ": damaged last line dropped",
              "Dropped incomplete last line" in log)

        X, out, _ = compare_with_reference(scheme, case, t, p, scheme, 4000)

        # Interpolation is exact on the samples, tones are on bin centres
        got = [abs(out[k, j]) for j, k in enumerate(bins)]
        err = max(abs(g - a) for g, a in zip(got, amps))
        check(scheme + ": tone amplitudes exact", err < 1e-9,
              "(max error %.1e)" % err)


def scenario_adaptive(work, rng):
    """Time step varies +-50%: compare with reference and tone amplitudes"""
    print("adaptive time step +-50%")

    M = 1500
    t = np.concatenate([[0], np.cumsum(1e-4*(1 + 0.5*(2*rng.random(M - 1) - 1)))])
    N = fft_size(M - 1)
    T = t[-1] - t[M - 1 - N]
    dtMean = T/N
    # f*dt = 0.01 .. 0.05 of the mean sampled step
    bins = [int(round(r*T/dtMean)) for r in (0.01, 0.02, 0.03, 0.05)]
    amps = [1.0, 0.7, 0.4, 0.2]
    p = signals(t, [k/T for k in bins], amps, rng, 0.0)

    for scheme, tolAmp in (("linear", 2e-2), ("cubicSpline", 2e-3)):
        case = os.path.join(work, "adaptive_" + scheme)
        write_case(
            case, np.column_stack([t, p]),
            {"maxFrequency": 100000, "interpolationScheme": scheme}
        )

        rc, log = run(case, work)
        check(scheme + ": run", rc == 0, "(exit %d)" % rc)
        if rc:
            continue

        check(scheme + ": maxFrequency above Nyquist truncated",
              "spectrum is truncated" in log)

        X, out, _ = compare_with_reference(scheme, case, t, p, scheme, 1e5)

        err = max(abs(abs(out[k, j]) - a)/a for j, (k, a) in enumerate(zip(bins, amps)))
        check(scheme + ": tone amplitudes (rel. tol %g)" % tolAmp, err < tolAmp,
              "(max rel. error %.1e)" % err)


def scenario_aliasing(work, rng):
    """Tone above maxFrequency must not appear in the spectrum"""
    print("aliasing: tone above maxFrequency")

    dt, M = 1e-4, 2001
    t = dt*np.arange(M)
    T = fft_size(M - 1)*dt
    fIn, fOut = 80/T, 575/T             # 400 Hz and 2875 Hz
    p = signals(t, [fIn, fOut], [0.5, 1.0], rng, 0.0)
    p = np.column_stack([p[:, 0] + p[:, 1]])

    case = os.path.join(work, "aliasing")
    write_case(case, np.column_stack([t, p]), {"maxFrequency": 1000})

    rc, _ = run(case, work)
    check("aliasing: run", rc == 0, "(exit %d)" % rc)
    if rc:
        return

    out, spectrum = read_output(case, 1)
    a = np.abs(out[:, 0])
    k = int(round(fIn*T))
    spurious = np.delete(a, k).max()

    check("aliasing: in-band tone amplitude", abs(a[k] - 0.5) < 1e-9,
          "(%.6f)" % a[k])
    check("aliasing: no aliased tone in band", spurious < 1e-9,
          "(max spurious %.1e)" % spurious)


def scenario_errors(work, rng):
    """Error handling and options"""
    print("options and error handling")

    dt, M = 1e-4, 1001
    t = dt*np.arange(M)
    p = signals(t, [500.0, 700.0], [1.0, 1.0], rng, 0.0)
    rows = np.column_stack([t, p])

    # -dict with another maxFrequency, outputFormat none
    case = os.path.join(work, "options")
    write_case(case, rows, {"maxFrequency": 4000})
    with open(os.path.join(case, "system", "otherDict"), "w") as f:
        f.write("FoamFile { version 2.0; format ascii; class dictionary; "
                "object otherDict; }\n")
        f.write('inputFileName "test.dat";\nmaxFrequency 1000;\n'
                'outputFormat none;\n')

    rc, _ = run(case, work, "-dict", "system/otherDict")
    root = os.path.join(case, "complexPressureData", "test")
    entries = os.listdir(root) if os.path.isdir(root) else []
    nFreq = len(np.loadtxt(os.path.join(root, "Spectrum.dat"))) if rc == 0 else 0
    check("-dict and outputFormat none", rc == 0 and entries == ["Spectrum.dat"]
          and nFreq == int(1000*fft_size(M - 1)*dt) + 1,
          "(exit %d, entries %s, %d frequencies)" % (rc, entries[:3], nFreq))

    def expect_error(name, rows, entries, message, extraLines=()):
        case = os.path.join(work, name)
        write_case(case, rows, entries, extraLines)
        rc, log = run(case, work)
        check(name + ": fatal error", rc == 1 and message in log,
              "(exit %d)" % rc)

    overlap = np.concatenate([rows[:600], rows[400:]])
    expect_error("descendingTimes", overlap, {}, "strictly ascending")

    expect_error("unknownOutputFormat", rows, {"outputFormat": "vtk"},
                 "Unknown outputFormat")

    expect_error("unknownScheme", rows, {"interpolationScheme": "spline"},
                 "Unknown Interpolator type")

    expect_error("brokenLine", rows[:10], {}, "values, expected",
                 ["%.17g 1" % (t[10]), "%.17g 1 2" % (t[11])])

    case = os.path.join(work, "missingFile")
    write_case(case, rows, {})
    os.remove(os.path.join(case, "pressureData", "test.dat"))
    rc, log = run(case, work)
    check("missingFile: fatal error", rc == 1 and "Cannot open file" in log,
          "(exit %d)" % rc)


# ---------------------------------------------------------------- main

def main():
    work = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else "cases")

    if not EXE:
        print("surfaceNoise not found: set SURFACE_NOISE or source OpenFOAM")
        return 1

    print("Using " + EXE)
    os.makedirs(work, exist_ok=True)

    rng = np.random.default_rng(1)

    scenario_uniform(work, rng)
    scenario_adaptive(work, rng)
    scenario_aliasing(work, rng)
    scenario_errors(work, rng)

    print("\n%d checks failed" % len(failures))
    for name in failures:
        print("    " + name)

    return len(failures)


if __name__ == "__main__":
    sys.exit(main())
