## csv_idealization_reading.
#
#   This function reads the specified CSV file with model parameters for the
#   idealization precision used by the get_data function.
# 
#   (c) Sami Samiei Esfahany, Adriaan van Natijne, Hans van der Marel and Freek van Leijen
#       Delft University of Technology, 2016.
#
#   Version:    1.0 
#   Created:    19 September 2016
#   Modified:   v1.1, 29 September 2016, allow empty lines in csv file.
#


from csv import reader as _csv;
_float = float;

def read_idealization(f):
    # Create initial structure.
    param = {'ST': {'variance': None,
                    'range':    None,
                    'power':    None},
             'T':  {}};

    # Open file.
    with open(f, 'r') as f_i:
        # Open as CSV.
        f_csv = _csv(f_i);

        # Skip initial header.
        next(f_csv);

        # Process rows
        # type, class, variance, range, power
        # 0   , 1    , 2       , 3    , 4
        for i, r in enumerate(f_csv):
            if len(r) != 5:
                print(('Ignored row in idealization CSV: \'{:s}\'. ' +
                       'Unexpected number of elements on row {:d}.').format(r, i +2));
            elif r[0] == 'T':
                # Check if key already exists.
                if r[1] in param['T']:
                    raise ValueError('Class {} defined twice.'.format(r[1]));

                # Create/add class.
                param['T'][r[1]] = {'variance': _float(r[2]),
                                    'power':    _float(r[4])};
            elif r[0] == 'ST':
                param['ST']['variance'] = _float(r[2]);
                param['ST']['range']    = _float(r[3]);
                param['ST']['power']    = _float(r[4]);
            else:
                raise ValueError('Unrecognised type, should be ST or T.');

    # Check if all values are set.
    assert param['ST']['variance'] is not None, 'Idealization parameters for the ST component should be set.';
    assert len(param['T']) > 0, 'No idealization parameters set for T component.';

    return param;

if __name__ == '__main__':
    d = read_idealization('cupido_example_idealization_precision.csv');

    print d;
