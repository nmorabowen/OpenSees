// Ladruno (WP-220): the reference load of the path-control integrators.
//
// domainChanged took phat = R(lambda+1), which assumes no unbalance now, so re-issuing the
// integrator left the last step's residual inside phat. The reference load is the difference
// R(lambda+1) - R(lambda). A pattern on a Constant series has no such difference (the load
// does not follow lambda), and decks do drive one: those keep R(lambda+1).
#ifndef LadrunoReferenceLoad_h
#define LadrunoReferenceLoad_h

#include <Vector.h>

inline void ladrunoRemoveUnbalance(Vector &phat, const Vector &unbalanceNow)
{
    Vector difference(phat);
    difference -= unbalanceNow;
    for (int i = 0; i < difference.Size(); i++) {
        if (difference(i) != 0.0) {
            phat = difference;
            return;
        }
    }
}

#endif
