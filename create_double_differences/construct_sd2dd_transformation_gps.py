## construct_sd2dd_transformation_gps.py
#
#   Function to construct double difference GPS observations based on
#   single difference observations.
# 
#   (c) Sami Samiei Esfahany, Adriaan van Natijne, Hans van der Marel and Freek van Leijen
#       Delft University of Technology, 2016.
#
#   Version:    1.0 
#   Created:    19 September 2016
#   Modified:   v1.1, 9 October 2016, bug fix in case SD observations of a single epoch
#                                     are selected (and therefore no DD can be formed).
#


import numpy as _np;


def construct_sd2dd_transformation_gps(BENCHMARKS_all, dates, SD_all, BM_ind, SD_list_ind):
    """Function to create double difference (or spatio-temporal) network from a list of single-differencces.
    This function also construct the matrix (Sfinal) which is the sd2dd transformation matrix.

    Inputs:
    ---------BENCHMARKS_all: 4 coloumns cell array of benchmarks (bm_id,x,y,tech)
    ---------dates: vector of dates of different projects (epochs)
    ---------SD_all: three column matrix of SD [bm_index_from bm_index_to date_index]
    ---------BM_ind
    ---------SD_list_ind

    Outputs:
    ---------Sfinal: sd2dd transformation matrix
    ---------DD_table: table of double differences (4 columns)
                            [bm_from bm_to date_from date_to]
    ---------SD_table: table of used single differences
                             [bm_from bm_to date]
    """

    # Select applicable benchmarks.
    BENCHMARKS = BENCHMARKS_all[BM_ind, :];
    SD_list    = SD_all[SD_list_ind, :];

    # Initial variables, unused in Python.
    Sfinal         = None;
    STfromto_final = None;

    if not _np.all(SD_list[:, 0] == SD_list[0, 0]):
        raise ValueError('This function only works if all the SD observations have a same FROM point.');

    Sfrom_to_epoch = SD_list;
    STfromto = None;
    S_sdref2dd = None;

    SDuniq, indSD1 = _np.unique(Sfrom_to_epoch[:, 1], return_index=True);
    # Loop over all unique 'to' stations.
    for i in xrange(0, len(indSD1)):

        ind = _np.where(Sfrom_to_epoch[:, 1] == Sfrom_to_epoch[indSD1[i], 1])[0];
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


    if (S_sdref2dd is not None):
        # Preparing for output (DD_table)
        STfromto_final = STfromto;

        # Preparing for output (SD_table)
        temp = S_sdref2dd.sum(axis=0);
        SD_ind_final = _np.where(temp != 0)[0];

        # Preparing for output (Sfinal)
        Sfinal = S_sdref2dd;
        Sfinal = Sfinal[:, SD_ind_final];

    else:
        STfromto_final = list();
        SD_ind_final = list();
        Sfinal = list();


    return (Sfinal, STfromto_final, SD_ind_final);


# Run when file is called directly.
if __name__ == '__main__':
    """Run an example version of this function."""

    # Import supporting functions.
    from netcdf_reading import read_netCdf as _read_netCdf; # netCdf
    import scipy.io as _sp_io;                              # Reference dataset

    # Load provided example data.
    sample_data = d = _read_netCdf('./gps.nc');

    # Select only observations in a single direction.
    selected_obs = d['Observations']['sensitivity'][:, 0] == 1;

    temp_BENCHMARKS = _np.asarray([[r['station_name'].strip(), r['x'], r['y'], r['station_class']] for r in d['Stations']], dtype=_np.object_);
    temp_SD_all = _np.asarray([[r['from_index'], r['to_index'], r['epoch']] for r in d['Observations'][selected_obs]], dtype=_np.int32);

    # Execute the function and display results
    sample_out = construct_sd2dd_transformation_gps(temp_BENCHMARKS,
                                                    temp_SD_all,
                                                    d['Projects']['project_epoch'],
                                                    _np.ones(temp_BENCHMARKS.shape[0], dtype=_np.bool),
                                                    _np.ones(temp_SD_all.shape[0], dtype=_np.bool));

    # Display results.

    # Some information on the input
    print('{:=^28s}'.format(' Inputs '));
    for k in ('Observations', 'Stations'):
        print('{:^15s} | {:^10s}'.format(k, d[k].shape));
    print('{:^15s} | {:^10d}'.format('Selected Obs', _np.sum(selected_obs)));

    # Compare the Sfinal transformation matrix.
    print('{:=^28s}'.format(' Size comparison '));
    print('{:^15s} + {:^10s}'.format('Variable', 'Python'));
    for i, k in enumerate(('Sfinal', 'STfromto_final', 'SD_ind_final')):
        print('{:^15s} | {:^10s}'.format(k, sample_out[i].shape));

    # Print results
    print('{:=^23s}'.format(' Sfinal '));
    print(sample_out[0]);

    print('{:=^23s}'.format(' STfromto_final '));
    print(sample_out[1]);

    print('{:=^23s}'.format(' SD_ind_final '));
    print(sample_out[2]);
