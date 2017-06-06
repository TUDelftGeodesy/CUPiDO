## csv_polygon_reading.
#
#   This function reads the specified CSV file with polygon boundaries used by the cupido function.
#   Furthermore, a function is created that uses the 'contains_points' function from Matplotlib 
#   to determine which benchmarks are within the polygon provided.
# 
#   (c) Sami Samiei Esfahany, Adriaan van Natijne, Hans van der Marel and Freek van Leijen
#       Delft University of Technology, 2016.
#
#   Version:    1.0 
#   Created:    24 September 2016
#   Modified:   
#                                        
#

import numpy as _np;
from csv import reader as _csv;
from matplotlib.path import Path as _path;

def read_polygon(f):
    # Open file.
    with open(f, 'r') as f_i:
        # Open as CSV.
        f_csv = _csv(f_i);

        # Skip initial header.
        next(f_csv);

        # Process rows (x, y)
        polygon = list();
        for i, r in enumerate(f_csv):
            if len(r) != 2:
                print(('Ignored row in polygon CSV: \'{:s}\'. ' +
                       'Unexpected number of elements on row {:d}.').format(r, i + 2));
            else:
                polygon.append(r);

    polygon = _np.asarray(polygon, dtype=_np.float32);
    
    return polygon;


def inside_polygon(pnts, poly):
    poly = _path(poly);
    rtn = poly.contains_points(pnts, radius=1);
    # radius=1 adds 1 meter extra to the polygon to avoid numerical problems.
    
    return rtn;

