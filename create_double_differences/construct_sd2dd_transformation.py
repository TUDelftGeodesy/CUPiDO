## construct_sd2dd_transformation.py
#
#   Function to construct double difference levelling observations based on
#   single difference observations.
# 
#   (c) Sami Samiei Esfahany, Adriaan van Natijne, Hans van der Marel and Freek van Leijen
#       Delft University of Technology, 2016.
#
#   Version:    1.1 
#   Created:    19 September 2016
#   Modified:   v1.1,  9 October 2016,  bug fix in case SD observations of a single epoch
#                                       are selected (and therefore no DD can be formed).
#               v1.2, 14 February 2017, the dependent rows of the transformation are excluded
#                                       resulting in a full row-rank Sfinal transformation and 
#                                       independent DDs.
#                                       Removal of maximum loop count of 100.  


import numpy as _np;
import scipy.io as _sp_io;
from scipy.sparse.csgraph import connected_components as _graphconncomp;
from scipy import linalg as _linalg;

def construct_sd2dd_transformation(BENCHMARKS_all, dates, SD_all, BM_ind, SD_list_ind):
    """Function to create double difference (or spatio-temporal) network from a list of single-differencces.
    This function also construct the matrix (Sfinal) which is the sd2dd transformation matrix.

    Inputs:
    ---------BENCHMARKS_all: 4 coloumns cell array of all benchmarks (bm_id,x,y,tech)
    ---------dates: vector of dates of different projects (epochs)
    ---------SD_all: three column matrix of all SD [bm_index_from bm_index_to date_index]
    ---------BM_ind: index of BENCHMARKS within region of interest
    ---------SD_list_ind: index of SD to be considered based on region of interest and period of interest

    Outputs:
    ---------Sfinal: sd2dd transformation matrix
    ---------DD_table: table of double differences (4 columns)
                            [bm_from bm_to date_from date_to]
    ---------SD_table: table of used single differences
                             [bm_from bm_to date]
    """

    # Initial variables.
    BENCHMARKS = BENCHMARKS_all[BM_ind, :];
    SD_list = SD_all[SD_list_ind, :];
    Nbm = BENCHMARKS.shape[0];                          # Benchmark count (length of BENCHMARKS array).
    Nbm_all = BENCHMARKS_all.shape[0];                  # Benchmark count (length of BENCHMARKS array).
    Nepoch = dates.shape[0];                            # Epoch count (length of dates array).

    # Initialisation for the main loop that constructs the final transformation matrix.
    flagSD = _np.zeros((SD_list.shape[0], 1),
                       dtype=_np.bool);                 # Flags, in the first iteration all SDs are unused (flag = 0).
    Sfinal = None;                                      # Final transformation matrix (from SD_list to final DD).
    STfromto_final = None;                              # Final DD from-to table.

    iloop = 0;                                          # Iteration index of the main loop.
                                                        # N.B. Nbm < 100 instead of <= in the loop. Zero-indexed!
    refind = [];                                        # Vector of the selected reference points.

    # Loop, until either: all points are connected, or
    #                     the number of loops exceeds the number of points
    while (~_np.all(flagSD)) and (iloop < Nbm):

        Nconnections = _np.zeros((Nepoch, Nbm_all), dtype=_np.int64);

        for i in xrange(0, Nepoch):
            sd_ind = _np.where(SD_all[:, 2] == i)[0];   # Select indices of all SDs within this epoch.
            sdtemp = SD_all[sd_ind, 0:2];           # Select single differences from this epoch.
                                                        # Convert to zero-indexed.
            Narcs = sdtemp.shape[0];                    # SD count (length of sdtemp array).

            # Create the connectivity matrix.
            A = _np.zeros((Nbm_all, Nbm_all), dtype=_np.bool);

            for j in xrange(0, Narcs):
                A[sdtemp[j, 0], sdtemp[j, 1]] = 1;
                A[sdtemp[j, 1], sdtemp[j, 0]] = 1;

            # Calculate the connectivity between benchmarks.
            S, C = _graphconncomp(A, directed=False);   # Analyse the graph network.

            # Count the connections from/to each benchmark.
            i1, i3, counts = _np.unique(C, return_inverse=True, return_counts=True);
            Nconnections[i, :] = counts[i3];

        temp2 = _np.mean(Nconnections, axis=0);
        temp2 = temp2[BM_ind];                       # This line added in order to select only BM in the AOI.
        temp2[refind] = -1;
        ind = _np.argmax(temp2);                        # Point with maximum connections.
        temp = temp2[ind];

        refind.append(ind);                             # Refind is the ref point ind in the BENCHMARKS, not in the BENCHMARKS_all!

        # Transformation to the reference point (still in SD).
        S_sd2sdref = None;
        Sfrom_to_epoch = None;

        for i in xrange(0, Nepoch):
            sd_ind = _np.where(SD_all[:, 2] == i)[0];   # Select indices of all SDs within this epoch.
            sdtemp = SD_all[sd_ind, 0:2];            # Select single differences from this epoch.
                                                        # Convert to zero-indexed.
            Narcs = sdtemp.shape[0];                    # SD count (length of sdtemp array).

            # Create the connectivity matrix.
            B = _np.zeros((Nbm_all, Nbm_all), dtype=_np.bool);

            for j in xrange(0, Narcs):
                B[sdtemp[j, 0], sdtemp[j, 1]] = 1;
                B[sdtemp[j, 1], sdtemp[j, 0]] = 1;

            # Calculate the connectivity between benchmarks.
            S, C = _graphconncomp(B, directed=False);   # Analyse the graph network.
            bm = _np.where(C == C[BM_ind[refind[iloop]]]); # Find points which are connected to the reference point.
            bm = _np.sort(bm[0]);
            ref_epoch_ind = _np.where(bm == BM_ind[refind[iloop]])[0];

            A = _np.zeros((Narcs, Nbm_all), dtype=_np.float64);
            for j in range(0, Narcs):
                A[j, sdtemp[j, 0]] = -1;
                A[j, sdtemp[j, 1]] =  1;
            A = A[:, bm];
            A = _np.delete(A, ref_epoch_ind, 1);

            S = _np.linalg.pinv(A) if A.size != 0 else A.T; # Pseudo inverse of A. Numpy can not handle empty A.
            S[_np.abs(S) <= 1e-10] = 0;                     # Set small values to zero.

            Sall = _np.zeros((S.shape[0], SD_all.shape[0]), dtype=_np.float64);
            for j in xrange(0, sd_ind.shape[0]):
                Sall[:, sd_ind[j]] = S[:, j];
            temp = bm;
            temp = _np.delete(temp, ref_epoch_ind, 0);

            #Here the BMs that are not in the AOI should be excluded
            #from bm and also from Sall the rows of Sall corresponds to
            #the BMs in 'bm'. The following two lines the exclusion
            # (BM_ind is the index of the BMs in the AOI).
            INDtemp = _np.in1d(temp,BM_ind);
            temp=temp[INDtemp];
            Sall=Sall[INDtemp,:];

            if temp.size != 0:
                temp_rows = _np.column_stack((temp *0 +BM_ind[refind[iloop]],
                                              temp,
                                              temp *0 +i));
                Sfrom_to_epoch = _np.vstack((Sfrom_to_epoch, temp_rows)) if Sfrom_to_epoch is not None \
                                                                         else temp_rows.copy();

            S_sd2sdref = _np.vstack((S_sd2sdref, Sall)) if S_sd2sdref is not None else Sall.copy();

        # Transformation from SD to DD.
        # Sfrom_to_epoch is filled with SDs with a common reference point.
        # A transformation for between-epoch-subtraction will be made based on this matrix.
        if (Sfrom_to_epoch is not None) and (Sfrom_to_epoch.size != 0):

            # Unique rows from Sfrom_to_epoch.
            # This is not implemented as-such in Numpy, so work around the problem.
            # We apply this elegant solution: http://stackoverflow.com/a/16971224 .
            Sfrom_to_epoch_flat = Sfrom_to_epoch[:, 0:2].copy().view(Sfrom_to_epoch.dtype.descr *2);
            SDuniq, indSD1 = _np.unique(Sfrom_to_epoch_flat, return_index=True);
            SDuniq = SDuniq.view(Sfrom_to_epoch.dtype).reshape(-1, 2);

            STfromto = None;
            S_sdref2dd = None;

            for i in xrange(0, indSD1.size):
                ind = _np.where(Sfrom_to_epoch_flat == Sfrom_to_epoch_flat[indSD1[i]])[0];
                tempdates = Sfrom_to_epoch[ind, 2];
                minind = _np.argmin(tempdates);
                datemin = tempdates[minind];

                if ind.size > 1:
                    S = _np.zeros((ind.size, Sfrom_to_epoch.shape[0]), dtype=_np.int32);
                    for j in xrange(0, ind.size):
                        S[j, ind[minind]] = -1;
                        S[j, ind[j]]      =  1;
                    S = _np.delete(S, minind, 0);

                    temp = _np.column_stack((_np.tile(Sfrom_to_epoch[indSD1[i], 0:2], (ind.size -1, 1)),
                                             _np.tile(datemin,                        (ind.size -1, 1)),
                                             tempdates[tempdates != datemin]));
                    STfromto = _np.vstack((STfromto, temp)) if STfromto is not None else temp.copy();
                    S_sdref2dd = _np.vstack((S_sdref2dd, S)) if S_sdref2dd is not None else S.copy();



        # Make the final S tranformation for this loop.
        if isinstance(S_sdref2dd, _np.ndarray) and isinstance(S_sd2sdref, _np.ndarray) \
            and S_sdref2dd.size > 0 and S_sd2sdref.size > 0:
            S_sd2dd = S_sdref2dd.dot(S_sd2sdref);

            Sind = (S_sd2dd != 0);

            k = 1;
            usefullDD = list();

            for v in xrange(0, S_sd2dd.shape[0]):       # Here we want to avoid DDs which are already included (avoid repetition).
                DDft = STfromto[v, 0:2];                # Select the from-to points of each DD.
                ind = _np.where(Sind[v, :] != 0)[0];    # Find the index acossiated to the used SDs for this DD.
                temp = SD_all[ind, 0:2];             # select the used SD from the SD_list. Convert to zero-indexed.

                for j in xrange(0, ind.size):
                  if DDft[1] in temp[j, 0:2] and _np.in1d(ind[j], SD_list_ind):
                    indtemp = _np.where(SD_list_ind == ind[j])[0];
                    if flagSD[indtemp] != 1:
                      usefullDD.append(v);
                
            usefullDD = _np.unique(usefullDD);          # Should be of type set, but set and Numpy don't mix well...
            if usefullDD.size > 0:
                S_sd02dd = S_sd2dd[usefullDD, :];
                STfromto = STfromto[usefullDD, :];

                Sfinal = _np.vstack((Sfinal, S_sd02dd)) if Sfinal is not None else S_sd02dd.copy();
                STfromto_final = _np.vstack((STfromto_final, STfromto)) if STfromto_final is not None else STfromto.copy();

        if Sfinal is not None:
            # Flag unused SDs.
            Sind = (Sfinal != 0);
            for i in xrange(0, Sfinal.shape[0]):
                DDft = STfromto_final[i, 0:2];
                ind = _np.where(Sind[i, :] != 0)[0];
                temp = SD_all[ind, 0:2];

                for j in xrange(0, ind.size):
                    if (len(_np.intersect1d(DDft, temp[j, :]))>0) and _np.in1d(ind[j], SD_list_ind):
                        indtemp = _np.where(SD_list_ind == ind[j])[0];
                        flagSD[indtemp] = 1;

        # Increase the loop number by one.
        iloop += 1;

    if (Sfinal is not None):
        temp = _np.sum(_np.abs(Sfinal), axis=0);
        SD_ind_final = _np.where(temp != 0)[0];
        Sfinal = Sfinal[:,SD_ind_final];

        # Remove dependent rows of the S-transformation
        QQ, RR, PP = _linalg.qr(Sfinal.transpose(),mode='economic',pivoting=True);
        rankS = _np.linalg.matrix_rank(Sfinal);
        ind_indep_obs = sorted(PP[0:rankS]);
        Sfinal = Sfinal[ind_indep_obs,:];
        STfromto_final = STfromto_final[ind_indep_obs,:];

    else:
        STfromto_final = list();
        SD_ind_final = list();
        Sfinal = list();

    return (Sfinal, STfromto_final, SD_ind_final);


# Run when file is called directly.
if __name__ == '__main__':
    """Run an example version of this function."""

    # Import supporting functions.
    from netcdf_reading import read_netCdf as _read_netCdf;

    # Load provided example data.
    sample_data = d = _read_netCdf('./leveling.nc');

    temp_BENCHMARKS = _np.asarray([[r['station_name'].strip(), r['x'], r['y'], r['station_class']] for r in d['Stations']], dtype=_np.object_);
    temp_SD_all = _np.asarray([[r['from_index'], r['to_index'], r['project_index']] for r in d['Observations']], dtype=_np.int32);

    # Execute the function and display results
    sample_out = construct_sd2dd_transformation(temp_BENCHMARKS,
                                                d['Projects']['project_epoch'],
                                                temp_SD_all,
                                                _np.arange(temp_BENCHMARKS.shape[0]),
                                                _np.arange(temp_SD_all.shape[0]));


    # Display results.

    ## Read the reference data
    reference_out = _sp_io.loadmat('result_data_for_DD_making.mat');

    # Compare the Sfinal transformation matrix.
    print('{:=^59s}'.format('Size comparison'));
    print('{:^10s} + {:^10s} + {:^10s} + {:^5}'.format('Variable', 'Matlab', 'Python', 'OK?'));
    for i, k in enumerate(('Sfinal', 'DD_table', 'SD_table')):
        print('{:^10s} | {:^10s} | {:^10s} | {:^5}'.format(k, sample_out[i].shape, reference_out[k].shape,
                                                           sample_out[i].shape == reference_out[k].shape));

    # Visual comparison
    # Sfinal
    print('\n\n{:=^59s}'.format(' Sfinal '));
    print('{:59s}'.format('The same!' if _np.allclose(sample_out[0], reference_out['Sfinal']) else
                          'Matrices are not the same (within tolerance).'));

    # DD_table
    print('\n\n{:=^77}'.format(' DD_table '));
    print('{:^33} | {:^33} | {:^5}'.format('Matlab', 'Python', 'OK?'));
    print((' | '.join(('-+-'.join(('{:-^6}',) * 4),) * 2)).format(*(('from', 'to', 'd_from', 'd_to') * 2)) + ' | -----');
    for r1, r2 in zip(reference_out['DD_table'], sample_out[1]):
        print('{:>6s} > {:6s} : {:>6d} > {:4d}   | {:>6s} > {:6s} : {:>6d} > {:4d}   | {:^5}'.format(
           r1[0][0], r1[1][0], r1[2][0][0], r1[3][0][0],
           r2[0][0], r2[1][0], r2[2],       r2[3],
           (r1[0][0] == r2[0][0] and r1[1][0] == r2[1][0] and
            r1[2][0][0] == r2[2] and r1[3][0][0] == r2[3])));

    # SD_table
    print('\n\n{:=^59}'.format(' SD_table '));
    print('{:^24} | {:^24} | {:^5}'.format('Matlab', 'Python', 'OK?'));
    print((' | '.join(('-+-'.join(('{:-^6}',) *3),) *2)).format(*(('from', 'to', 'date') *2)) + ' | -----');
    for r1, r2 in zip(reference_out['SD_table'], sample_out[2]):
        print('{:>6s} > {:6s} : {:^6d} | {:>6s} > {:6s} : {:^6d} | {:^5}'.format(r1[0][0], r1[1][0], r1[2][0][0],
                                                                                 r2[0][0], r2[1][0], r2[2],
                                                                                 (r1[0][0] == r2[0][0] and
                                                                                  r1[1][0] == r2[1][0] and
                                                                                  r1[2][0][0] == r2[2])));
