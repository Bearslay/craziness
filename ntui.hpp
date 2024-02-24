#ifndef NTypeUI
#define NTypeUI

#include <ncurses.h>
#include <locale.h>
#include <utility>
#include <type_traits>
#include <vector>
#include <unordered_map>

namespace ntui {
    template <typename Type> class MathVector {
        private:
            Type IAxis = 0; 
            Type JAxis = 0;
            Type KAxis = 0;

        public:
            MathVector(Type iAxis = 0, Type jAxis = 0, Type kAxis = 0) {
                // Make sure that Type is only an arithmetic data type
                // TODO: look into seeing if unsigned data types break things
                static_assert(std::is_arithmetic<Type>::value, "Type must be an arithmetic type");
                
                IAxis = iAxis;
                JAxis = jAxis;
                KAxis = kAxis;
            }

            inline const bool operator == (const MathVector<Type> &vector) {return (IAxis == vector.IAxis && JAxis == vector.JAxis && KAxis == vector.KAxis);}

            inline MathVector<Type> operator += (const MathVector<Type> &vector) {
                MathVector<Type> original = MathVector<Type>(IAxis, JAxis, KAxis);

                IAxis += vector.IAxis;
                JAxis += vector.JAxis;
                KAxis += vector.KAxis;

                return original;
            }

            inline MathVector<Type> operator -= (const MathVector<Type> &vector) {
                MathVector<Type> original = MathVector<Type>(IAxis, JAxis, KAxis);

                IAxis -= vector.IAxis;
                JAxis -= vector.JAxis;
                KAxis -= vector.KAxis;

                return original;
            }

            inline MathVector<Type> operator + (const MathVector<Type> &vector) {return MathVector<Type>(IAxis + vector.IAxis, JAxis + vector.JAxis, KAxis + vector.KAxis);}

            inline MathVector<Type> operator - (const MathVector<Type> &vector) {return MathVector<Type>(IAxis - vector.IAxis, JAxis - vector.JAxis, KAxis - vector.KAxis);}

            const Type geti() {return IAxis;}
            const Type getj() {return JAxis;}
            const Type getk() {return KAxis;}

            Type seti(Type iAxis) {
                Type original = IAxis;

                IAxis = iAxis;

                return original;
            }
            Type setj(Type jAxis) {
                Type original = JAxis;
                
                JAxis = jAxis;

                return original;
            }
            Type setk(Type kAxis) {
                Type original = KAxis;
                
                KAxis = kAxis;

                return original;
            }
    };

    class Cell {

    };

    template <typename Arrangement, typename Type> class VectorGrid {
        private:
            std::unordered_map<MathVector<Arrangement>, Type> Contents;
    };
}

#endif
