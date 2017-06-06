## construct_t_idealization_covmx.py
#
#   Function to construct the temporal idealization precsion.
# 
#   (c) Sami Samiei Esfahany, Adriaan van Natijne, Hans van der Marel and Freek van Leijen
#       Delft University of Technology, 2016.
#
#   Version:    1.1 
#   Created:    19 September 2016
#   Modified:   v1.1, 29 September 2016, enabling different benchmark classes
#



import numpy as _np;
import scipy.io as _sp_io;

def _make_t_cov_ing(IdealNoiseParam, pipjtmtn, DDtable, BMxy, Dates0):

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

    # Extract the values for the idealisation noise.
    VAR1 = IdealNoiseParam[p1ind, 0];
    VAR1[p1ind != p2ind] = 0;
    POW = IdealNoiseParam[p1ind, 1];

    Dist2 = _np.abs((BMx[p1ind] -BMx[p2ind]) **2 + (BMy[p1ind] -BMy[p2ind]) **2);
    Dist_ind = Dist2 == 0;
    Tint = _np.abs(Dates0[t1ind].squeeze() -Dates0[t2ind].squeeze());
    T1 = Dates0[t1ind].squeeze();
    T2 = Dates0[t2ind].squeeze();

    Qpijtmn = 0.5 *VAR1 *(_np.power(T1,   POW) +
                          _np.power(T2,   POW) -
                          _np.power(Tint, POW)) *Dist_ind.astype(_np.float64);

    return Qpijtmn;

def construct_t_idealization_covmx(IdealNoiseParam, DDtable, BMxy, Types, Dates):
    """Function to evaluate the idealization covariance matrix based on non-stationary (flicker-noise-kind) noise
        in time and stationary gaussian covariance function in space.

    Inputs:
    --------- IdealNoiseParam: Ntypes x 3 array of noise parameters
                               [variance dummy power]
                               for only-temporal compnent (i.e. no spatial correlation spatial_corr_range=eps)
    --------- DDtable: index_table of double differences
                               [bm_index_from bm_index_to project_index_from project_index_to]
                               bm_index of 0 is for imaginary ref point of GPS data outside of AOI!
    --------- BMxy: mx2 matrix of benchmarks coordinates [x y] in meter
    --------- Types: benchmark class
    --------- Dates: vector of dates in decimal year

    Output: Qid covariance matrix of idealization noise."""

    # Check if the Dates vector is a floating point (decimal).
    if Dates.dtype.kind != 'f':
        ValueError('Dates should be floating point, decimal year.');

    # Analyse/process the idealisation parameters
    assert isinstance(IdealNoiseParam, dict), 'IdealNoiseParam should be a dictionary of classes.';
    assert IdealNoiseParam.has_key('DEFAULT'), 'IdealNoiseParam should contain a \'DEFAULT\' class.';
    assert not _np.any([r.has_key('range') for r in IdealNoiseParam.itervalues()]), \
        'Dictionaries of IdealNoiseParams can only be applied for temporal components (range = eps).';

    # Construct a matrix of variances and powers for the Temporal case (range = eps).
    param_default = IdealNoiseParam['DEFAULT'];
    BM_param = _np.asarray([(r['variance'], r['power']) for r in [IdealNoiseParam.get(t, param_default) for t in Types]],
                           dtype=_np.float64);
    del param_default;

    # Ingerdients ind (results of applying error propagation).
    ing_ind = _np.column_stack((_np.kron([[1, 3], [1, 4], [2, 3], [2, 4]], [[1], [1], [1], [1]]),
                                _np.kron([[1], [1], [1], [1]], [[1, 3], [1, 4], [2, 3], [2, 4]])));
    # Results of applying error propagation.
    signVec = _np.asarray([+1, -1, -1, +1, -1, +1, +1, -1, -1, +1, +1, -1, +1, -1, -1, +1]).T;

    Qid = _np.zeros((DDtable.shape[0], DDtable.shape[0]));
    Qid_series = _np.zeros((ing_ind.shape[0], DDtable.shape[0], DDtable.shape[0]));

    # Create a synthetic benchmark, for GPS.
    if 0 in DDtable[:, 0:2]:
        BMxy = _np.vstack((_np.ones((1, 2)) *1e20, BMxy));
        BM_param = _np.vstack((_np.zeros((1, 2)), BM_param));
        DDtable = DDtable.copy();   # Copy DDtable, not to edit the input table.
        DDtable[:, 0:2] += 1;       # Add one, to add zero synthetic benchmark.

    # Create dates wrt random reference date
    Dates0 = Dates - _np.min(Dates);

    k = 0;
    for ind, s in zip(ing_ind, signVec):
        Qpijtmn = _make_t_cov_ing(BM_param, ind, DDtable, BMxy, Dates0);
        Qid    += s *Qpijtmn;

        Qid_series[k, ...] = Qid;
        k += 1;

    Qid = _np.triu(Qid, 1).T + _np.triu(Qid);

    return Qid;
