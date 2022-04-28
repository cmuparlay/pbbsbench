#include <algorithm>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "common/geometry.h"
#ifndef WEAVESS_DISTANCE_H
#define WEAVESS_DISTANCE_H

namespace weavess {
    class Distance {
    public:
        template<typename T>
        T compare(const T *a, const T *b, unsigned length) const {
            T result = 0;

            float diff0, diff1, diff2, diff3;
            const T *last = a + length;
            const T *unroll_group = last - 3;

            /* Process 4 items with each loop for efficiency. */
            while (a < unroll_group) {
                diff0 = a[0] - b[0];
                diff1 = a[1] - b[1];
                diff2 = a[2] - b[2];
                diff3 = a[3] - b[3];
                result += diff0 * diff0 + diff1 * diff1 + diff2 * diff2 + diff3 * diff3;
                a += 4;
                b += 4;
            }
            /* Process last 0-3 pixels.  Not needed for standard vector lengths. */
            while (a < last) {
                diff0 = *a++ - *b++;
                result += diff0 * diff0;
            }
            return result;
        }
    };
}

template<class fvec_point>
float distance(fvec_point *p, fvec_point *q, unsigned d){
    weavess::Distance distfunc;
    return distfunc.compare<float>((p->coordinates).begin(), (q->coordinates).begin(), d);

}

#endif //WEAVESS_DISTANCE_H