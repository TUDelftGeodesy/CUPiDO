## construct_st_idealization_covmx.py
#
#   Function to construct the spatio-temporal idealization precsion.
# 
#   (c) Sami Samiei Esfahany, Adriaan van Natijne, Hans van der Marel and Freek van Leijen
#       Delft University of Technology, 2016.
#
#   Version:    1.0 
#   Created:    19 September 2016
#   Modified:   
#



import numpy as _np;
import scipy.io as _sp_io;

def _make_st_cov_ing(IdealNoiseParam, pipjtmtn, DDtable, BMxy, Dates0):
    GPS_flag = 0 in DDtable[:, 0:2];
    if GPS_flag:                    # To account for the imaginary ref point (for GPS)
        BMxy = _np.vstack((_np.ones((1, 2)) *1e20, BMxy));
        DDtable = DDtable.copy();   # Copy DDtable, not to edit the input table.
        DDtable[:, 0:2] += 1;

    pipjtmtn -= 1;                  # Convert to zero indexed.
    p1 = pipjtmtn[0];
    p2 = pipjtmtn[1];
    t1 = pipjtmtn[2];
    t2 = pipjtmtn[3];

    BMx = BMxy[:, 0];
    BMy = BMxy[:, 1];

    DDS = _np.hstack((DDtable[:, 0:2], DDtable[:, 0:2])) -1;
    DDT = _np.hstack((DDtable[:, 2:4], DDtable[:, 2:4])) -1;
    p1ind, p2ind = _np.meshgrid(DDS[:, p1], DDS[:, p2]);
    t1ind, t2ind = _np.meshgrid(DDT[:, t1], DDT[:, t2]);

    Dist2 = _np.abs((BMx[p1ind] -BMx[p2ind]) **2 + (BMy[p1ind] -BMy[p2ind]) **2);
    Tint = _np.abs(Dates0[t1ind].squeeze() -Dates0[t2ind].squeeze());
    T1 = Dates0[t1ind].squeeze();
    T2 = Dates0[t2ind].squeeze();

    Qpijtmn = IdealNoiseParam[0] / 2 * (_np.power(T1,   IdealNoiseParam[2]) +
                                        _np.power(T2,   IdealNoiseParam[2]) -
                                        _np.power(Tint, IdealNoiseParam[2])) * \
                                       _np.exp(-1 / _np.power(IdealNoiseParam[1], 2) * Dist2);

    if GPS_flag:
        temp = p1ind + p2ind;
        tempind = (temp != 0); # Zero indexed.
        Qpijtmn *= tempind;

    return Qpijtmn;

def construct_st_idealization_covmx(IdealNoiseParam, DDtable, BMxy, Dates):
    """Function to evaluate the idealization covariance matrix based on non-stationary (flicker-noise-kind) noise
        in time and stationary gaussian covariance function in space.

    Inputs:
    --------- IdealNoiseParam: 3x1 vector of noise parameters
                               [variance  spatial_corr_range power]
                               for only-temporal compnent (i.e. no spatial correlation spatial_corr_range=eps)
    --------- DDtable: index_table of double differences
                               [bm_index_from bm_index_to project_index_from project_index_to]
                               bm_index of 0 is for imaginary ref point of GPS data outside of AOI!
    --------- BMxy: mx2 matrix of benchmarks coordinates [x y] in meter
    --------- Dates: vector of dates in decimal year

    Output: Qid covariance matrix of idealization noise."""

    # Check if the Dates vector is a floating point (decimal).
    if Dates.dtype.kind != 'f':
        ValueError('Dates should be floating point, decimal year.')

    # Ingerdients ind (results of applying error propagation).
    ing_ind = _np.column_stack((_np.kron([[1, 3], [1, 4], [2, 3], [2, 4]], [[1], [1], [1], [1]]),
                                _np.kron([[1], [1], [1], [1]], [[1, 3], [1, 4], [2, 3], [2, 4]])));
    # Results of applying error propagation.
    signVec = _np.asarray([+1, -1, -1, +1, -1, +1, +1, -1, -1, +1, +1, -1, +1, -1, -1, +1]).T;

    Qid = _np.zeros((DDtable.shape[0], DDtable.shape[0]));
    Qid_series = _np.zeros((ing_ind.shape[0], DDtable.shape[0], DDtable.shape[0]));

    # Create dates wrt random reference date
    Dates0 = Dates - _np.min(Dates);

    k = 0;
    for ind, s in zip(ing_ind, signVec):
        Qpijtmn = _make_st_cov_ing(IdealNoiseParam, ind, DDtable, BMxy, Dates0);
        Qid    += s *Qpijtmn;

        Qid_series[k, ...] = Qid;
        k += 1;

    Qid = _np.triu(Qid, 1).T + _np.triu(Qid);

    return Qid;

