#ifndef __KHAYYAM_SIMPLEX_H__
#define __KHAYYAM_SIMPLEX_H__

//
// Simplex parameter fitting
//

#ifdef __cplusplus
extern "C" {
#endif

float downhill_simplex (int nvalues, float values[], float deltas[], float maxerror, int nruns, int niterations, float (* distance) (int, const float[], void *), void *data);

#ifdef __cplusplus
}
#endif

#endif
