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

\*---------------------------------------------------------------------------*/

#include "FileInterface.H"
#include "IFstream.H"
#include "OFstream.H"

#include <cstdlib>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

//- Parse whitespace-separated numbers up to the first non-number,
//- calling action(index, value) for each. Return the number of values.
template<class Action>
Foam::label parseValues(const std::string& line, const Action& action)
{
    const char* p = line.c_str();
    char* end = nullptr;

    Foam::label n = 0;

    while (true)
    {
        const double value = std::strtod(p, &end);

        if (end == p)
        {
            break;
        }

        action(n, value);

        ++n;
        p = end;
    }

    return n;
}


//- True if the line contains only whitespace
bool isBlank(const std::string& line)
{
    return line.find_first_not_of(" \t\r") == std::string::npos;
}


//- Write complex number as Re+Imj
void writeComplex(Foam::Ostream& os, const Foam::complex& z)
{
    os  << z.real() << (z.imag() < 0 ? '-' : '+') << Foam::mag(z.imag())
        << 'j';
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FileInterface::FileInterface(const fileName& name)
:
    name_(name)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FileInterface::read
(
    scalarList& times,
    List<scalarList>& pressure
) const
{
    IFstream is(name_);

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot open file " << name_ << nl
            << exit(FatalError);
    }

    std::string line;

    // First pass: count time steps and values per line,
    // so that the storage is allocated once

    label nTimes = 0;
    label nValues = 0;

    while (is.getLine(line))
    {
        if (isBlank(line))
        {
            continue;
        }

        if (!nTimes)
        {
            nValues = parseValues(line, [](const label, const scalar){});
        }

        ++nTimes;
    }

    if (nValues < 2)
    {
        FatalErrorInFunction
            << "Expected time and at least one pressure value per line in "
            << name_ << nl
            << exit(FatalError);
    }

    times.resize_nocopy(nTimes);
    pressure.resize(nValues - 1);

    for (scalarList& p : pressure)
    {
        p.resize_nocopy(nTimes);
    }

    // Second pass: fill the data transposed to pressure[pointI][timeI]

    is.rewind();

    label lineI = 0;
    label timeI = 0;

    while (timeI < nTimes && is.getLine(line))
    {
        ++lineI;

        if (isBlank(line))
        {
            continue;
        }

        const label n = parseValues
        (
            line,
            [&](const label i, const scalar value)
            {
                if (i == 0)
                {
                    times[timeI] = value;
                }
                else if (i < nValues)
                {
                    pressure[i-1][timeI] = value;
                }
            }
        );

        if (n != nValues)
        {
            // Short last line shows the interrupted write
            if (n < nValues && timeI == nTimes - 1)
            {
                WarningInFunction
                    << "Dropped incomplete last line " << lineI
                    << " of " << name_ << endl;

                break;
            }

            FatalErrorInFunction
                << "Line " << lineI << " of " << name_ << " has " << n
                << " values, expected " << nValues << nl
                << exit(FatalError);
        }

        ++timeI;
    }

    if (timeI < nTimes)
    {
        times.resize(timeI);

        for (scalarList& p : pressure)
        {
            p.resize(timeI);
        }
    }
}


void Foam::FileInterface::write
(
    const scalar frequency,
    const UList<complex>& data
) const
{
    if (name_.ext() != "msh")
    {
        FatalErrorInFunction
            << "Unsupported output file format: " << name_.ext() << nl
            << "Use file extension .msh" << nl
            << exit(FatalError);
    }

    OFstream os(name_);

    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot open file " << name_ << nl
            << exit(FatalError);
    }

    os  << "$NodeData" << nl

        // String tags: number of tags, name of data
        << 1 << nl
        << "node_data" << nl

        // Real tags: number of tags, time (frequency here)
        << 1 << nl
        << frequency << nl

        // Integer tags: number of tags, time step index,
        // number of field components, number of nodes, partition index
        << 4 << nl
        << 0 << nl
        << 1 << nl
        << data.size() << nl
        << 0 << nl;

    forAll(data, pointI)
    {
        os  << pointI + 1 << ' ';
        writeComplex(os, data[pointI]);
        os  << nl;
    }

    os  << "$EndNodeData" << nl;
}


// ************************************************************************* //
