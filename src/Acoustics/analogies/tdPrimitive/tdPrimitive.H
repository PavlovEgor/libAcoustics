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
    Foam::tdPrimitive

Description
    List of values with the time derivative.

    The derivative of each value is computed from its previous values:
    zero at the first step, then the first order difference for the second
    step or with the variable time step, otherwise the second order backward
    difference.

    Usage per time step: set value(i), then call dot(t, i) once.

SourceFiles
    tdPrimitive.H

\*---------------------------------------------------------------------------*/

#ifndef tdPrimitive_H
#define tdPrimitive_H

#include "List.H"
#include "pTraits.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class tdPrimitive Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class tdPrimitive
{
    // Private Data

        //- Time step is variable: use the first order difference
        bool varTimeStep_;

        //- Time of the last derivative
        scalar prevt_;

        //- Last time step
        scalar deltat_;

        //- Current values
        List<T> new_;

        //- Values at the previous derivative
        List<T> old_;

        //- Values before the previous derivative
        List<T> oldOld_;


public:

    // Constructors

        //- Construct with one value and the previous time
        tdPrimitive(T initv, scalar tprev)
        :
            varTimeStep_(false),
            prevt_(tprev),
            deltat_(tprev),
            new_(1, initv),
            old_(0),
            oldOld_(0)
        {}

        //- Default construct with one zero value
        tdPrimitive()
        :
            varTimeStep_(false),
            prevt_(0.0),
            deltat_(GREAT),
            new_(1, pTraits<T>::zero),
            old_(0),
            oldOld_(0)
        {}


    // Member Functions

        //- Number of values
        label size() const
        {
            return new_.size();
        }

        //- Change the number of values
        void resize(label newsz)
        {
            new_.resize(newsz);
        }

        //- Use the first order difference only
        void setTimeStepVariable(bool adjustTimeStep)
        {
            varTimeStep_ = adjustTimeStep;
        }

        //- Current value
        const T& value(label idx) const
        {
            return new_[idx];
        }

        //- Current value
        T& value(label idx)
        {
            return new_[idx];
        }

        //- Time derivative of the current value at the time newt
        T dot(scalar newt, label idx)
        {
            T dvdt(pTraits<T>::zero);

            if (newt > prevt_)
            {
                deltat_ = newt - prevt_;
                prevt_ = newt;
            }

            if (old_.size() <= idx)
            {
                // no previous value: zero derivative
                old_.resize(idx+1);
                old_[idx] = new_[idx];
            }
            else if ((oldOld_.size() <= idx) || varTimeStep_)
            {
                dvdt = (new_[idx] - old_[idx])/deltat_;
                oldOld_.resize(idx+1);
                oldOld_[idx] = old_[idx];
                old_[idx] = new_[idx];
            }
            else
            {
                dvdt =
                    (3.0*new_[idx] - 4.0*old_[idx] + oldOld_[idx])
                   /2.0/deltat_;
                oldOld_[idx] = old_[idx];
                old_[idx] = new_[idx];
            }

            return dvdt;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
