"""Private bridge between PearsonT3 and Python's Randompack generator."""

from libc.limits cimport INT_MAX
from libc.stdint cimport int32_t, int64_t
from libc.string cimport memcpy
cimport numpy as cnp
import numpy as np

cnp.import_array()

cdef extern from *:
    """
    typedef int (*ac_draw)(void *, int, int, int32_t *);
    void p3_python(int, const double *, const double *, const double *, int, const double *,
                   int, int, int, int, int, ac_draw, void *, double *, double *,
                   double *, int *, double *, double *, double *, double *,
                   double *, double *, int *, int64_t *, int *, int *);
    """
    ctypedef int (*ac_draw)(void *, int, int, int32_t *) noexcept
    void p3_python(int, const double *, const double *, const double *, int, const double *,
                   int, int, int, int, int, ac_draw, void *, double *, double *,
                   double *, int *, double *, double *, double *, double *,
                   double *, double *, int *, int64_t *, int *, int *)


cdef class _DrawStream:
    cdef object rng, remaining, error
    cdef cnp.ndarray buffer
    cdef Py_ssize_t offset

    def __cinit__(self, rng, batches):
        self.rng = rng
        self.remaining = batches
        self.error = None
        self.buffer = np.empty(0, dtype=np.int32)


cdef int _draw(void *context, int n, int upper, int32_t *out) noexcept:
    cdef _DrawStream stream = <_DrawStream><object>context
    cdef Py_ssize_t batches
    try:
        if stream.offset == stream.buffer.size:
            # Both bootstrap levels request the same number and range of starts.
            # Refill in whole samples and never draw beyond the last resample.
            batches = min(stream.remaining, max(1, 65536 // n))
            if batches <= 0:
                raise RuntimeError("unexpected extra PearsonT3 bootstrap draw")
            stream.buffer = stream.rng.int(size=batches * n, a=1, b=upper, dtype=np.int32)
            stream.offset = 0
            stream.remaining -= batches
        memcpy(out, <int32_t *>stream.buffer.data + stream.offset, n * sizeof(int32_t))
        stream.offset += n
        return 0
    except BaseException as exc:
        stream.error = exc
        return -1


def run_t3(const double[::1] t, const double[::1] x, const double[::1] y,
           const double[::1] levels, int b1, int b2, int n_lambda, int block, int linear, rng):
    if not 10 <= t.shape[0] <= INT_MAX // 2 or t.shape[0] != x.shape[0] or t.shape[0] != y.shape[0]:
        raise ValueError("invalid PearsonT3 array lengths")
    if b1 < 2 or b2 < 2 or n_lambda < 1 or not 0 <= block < t.shape[0]:
        raise ValueError("invalid PearsonT3 bootstrap settings")
    if not 1 <= levels.shape[0] <= INT_MAX:
        raise ValueError("invalid number of confidence levels")
    cdef _DrawStream stream = _DrawStream(rng, int(b1) * (int(b2) + 1))
    cdef double corr, se_z
    cdef double tau[2]
    cdef int actual_block, zero_inner, status
    cdef int persistence[2]
    cdef int64_t degenerate
    cdef cnp.ndarray cover = np.empty(n_lambda)
    cdef cnp.ndarray z1 = np.empty(b1)
    cdef cnp.ndarray se2 = np.empty(b1)
    cdef cnp.ndarray grid = np.empty(n_lambda)
    cdef cnp.ndarray ci = np.empty((levels.shape[0], 2))
    cdef cnp.ndarray selected = np.empty(levels.shape[0])
    p3_python(t.shape[0], &t[0], &x[0], &y[0], levels.shape[0], &levels[0], b1, b2, n_lambda,
              block, linear, _draw, <void *>stream, &corr, <double *>ci.data, tau, &actual_block,
              &se_z, <double *>selected.data, <double *>cover.data, <double *>z1.data,
              <double *>se2.data, <double *>grid.data, persistence, &degenerate, &zero_inner, &status)
    if stream.error is not None:
        raise stream.error
    if status:
        return status, tuple(persistence) if status == 2 else None
    return 0, (corr, ci, tuple(tau), actual_block, se_z, selected,
               cover, z1, se2, grid, tuple(persistence), degenerate, zero_inner)
