#include "forces.H"
#include "volFieldsFwd.H"
#include "Tuple2.H"
#include "OFstream.H"
#include "SoundObserver.H"
#include "tdPrimitive.H"

#ifndef AcousticAnalogy_H
#define AcousticAnalogy_H

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;


namespace functionObjects
{

class AcousticAnalogy
:
    public forces
{

protected:

        //-
        autoPtr<OFstream> analogyOutPtr_;

        //-
        label probeFreq_;

        //-
        scalar timeStart_;

        //-
        scalar timeEnd_;

        //-
        bool writeFft_;

        //- Interpolation of observer samples with a variable time step
        //- to the uniform FFT grid, see Foam::Interpolator
        word fftInterpolationScheme_;

        //- Sonic speed, m/s
        scalar c0_;

        //- List of observers
        List<SoundObserver> observers_;

        //-
        label probeI_;

protected:

        //-
        void makeFile();

        //-
        void writeFft();

        //-
        virtual void correct();
        
        //-
        virtual void updateSurfaceCf();

public:

    //- Time derivative helper, see Foam::tdPrimitive
    template<class T>
    using tdPrimitive = Foam::tdPrimitive<T>;

        //- Runtime type information
        TypeName("AcousticAnalogy");

        // Constructors

        //- Construct from given Time and dictionary.
        AcousticAnalogy
        (
            const word& name,
            const Time&,
            const dictionary&
        );

        //- Construct from given objectRegistry and dictionary.
        AcousticAnalogy
        (
            const word& name,
            const objectRegistry&,
            const dictionary&
        );

        //- Destructor
        virtual ~AcousticAnalogy();

public:
        //- Observers list
        const List<SoundObserver>& getSoundObservers() const noexcept
        {
            return observers_;
        }

        //- Read the Curle data
        virtual bool read(const dictionary&);

        //-
        virtual bool write();

        //-
        virtual bool execute();
};

} //namespace functionObjects

} //namespace Foam

#endif

//
//END-OF-FILE
//



