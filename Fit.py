"""
Created on Wed Oct 2 14:50:19 2019

Functions for creating models with rotational symmetry and fitting.

Alistair Curd
University of Leeds
2 October 2019

Software Engineering practices applied

Joanna Leng (an EPSRC funded Research Software Engineering Fellow (EP/R025819/1)
University of Leeds
October 2019

---
Copyright 2019 Peckham Lab

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this file except in compliance with the License. You may obtain a copy of the
License at
http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
"""

import os
import sys
import argparse
import datetime

from tkinter import Tk
from tkinter.filedialog import askopenfilename
import numpy as np
import modelling_general as models
from modelling_general import ModelWithFitSettings
import modelstats as stats
from relative_positions import getdistances
import utils
import plotting
import reports
import pandas as pd

class Number:
    """Class object providing only a number, to be passed to rotational
    symmetry models to provide the order of symmetry, without being fitted.
    """
    def __init__(self, number):
        self.number = number


def get_inputs(info):

    #root = Tk()
    Tk().withdraw()
    print('\n\nPlease select input file containing relative positions to assess '
          'for rotational symmetry (.csv or .txt with comma delimiters).')

    infile = askopenfilename()
    #root.destroy()

    print("The file you selected is: ", infile)

    info['in_file_and_path'] = infile
    info['verbose'] = True


def rot_sym_only(separation_values,
                 diameter,
                 broadening,
                 amplitude):
    """Parametric model for distances between localisations on vertices
    of a polygon (order of symmetry = number of vertices). The value of
    the model at a distance is termed the relative position density (RPD)
    at that distance. Broadening is modelled assuming Gaussian imprecision on
    the positions of the vertices.

    Calls sym_order.number from outside, so that the degree of symmetry os
    not fitted when this function is passed to scipy.optimize.curve_fit.

    Args:
        separation_values (numpy array):
            Distances at which density values of the model will be obtained.
        diameter (float):
            Diameter of the circle containing the vertices of the polygon.
        broadening (float):
            Broadening of the peaks located at distances between the vertices.
        amplitude (float):
            Amplitude of the contribution of one inter-vertex distance.

    Returns:
        rpd (numpy array):
            The relative position density given by the model at the
            separation_values.
    """
    vertices = models.generate_polygon_points(sym_order.number,
                                              diameter)
    filter_distance = (2*diameter)
    

    # getdistances includes removel of duplicates 27/11/2019
    relative_positions = getdistances(vertices, filter_distance)
    xy_separations = np.sqrt(relative_positions[:, 0] ** 2
                             + relative_positions[:, 1] ** 2)

    # Initialise array for set of density values.
    rpd = separation_values * 0.

    # Add 2D pair correlations at the distances between vertices.
    for distance in xy_separations:
        rpd = (rpd
               + (amplitude * models.pairwise_correlation_2d(separation_values,
                                                             distance,
                                                             broadening)
                  )
               )

    return rpd


def rot_sym_replocs_substructure_isotropic_bg_with_onset(
        r, dia,
        vertssd, vertsamp,
        replocssd, replocsamp,
        substructsd, substructamp):
    """Parametric model for distances between localisations on vertices
    of a polygon (order of symmetry = number of vertices). The value of
    the model at a distance is termed the relative position density (RPD)
    at that distance. Includes the effect of localisation precision and
    unresolvable substructure at the polygon vertices.

    In this model, we take account of the fact that rotationally
    symmetric structures may be unlikely to be found within one another:
    we allow the isotropic (linearly increasing) background term to
    remain zero until an onset distance (bgonset) is reached.

    Args:
        r (numpy array):
            Distances at which density values of the model will be obtained.
        dia (float):
            Diameter of the circle containing the vertices of the polygon.
        vertsd (float):
            Broadening of the peaks located at
            distances between the vertices.
        vertsamp (float):
            Amplitude of the contribution of one inter-vertex distance.
        replocssd (float):
            Spread representing localisation precision for repeated
            localisations of the same fluorescent molecule.
        replocsamp (float):
            Amplitude of the contribution of repeated localisations
            of the same fluorescent molecule.
        substructsd (float):
            Spread of a contribution resulting from unresolvable
            substructure, or mislocalisations resulting from
            a combination of simultaneous nearby emitters.
        substructamp (float):
            Amplitude of the contribution of unresolvable
            substructure, or mislocalisations resulting from
            a combination of simultaneous nearby emitters.
        bggrad (float):
            Gradient of an isotropic (linearly increasing) background term.
        bgonset (float):
            Onset distance for linearly increasing background term,
            since rotationally symmetric structures may exclude
            one another.
    Returns:
        rpd (numpy array):
            The relative position density given by the model
            at distances r.
    """
    # GETDISTANCES IS CALLED BY relative_positions.py.
    # There is an option to set verbose=True if output to screen if desired.

    # Get RPD arising from rotationally symmetric structure with broadening
    # only from imprecision on single vertex points.
    rpd = rot_sym_only(r, dia, vertssd, vertsamp)

    # Isotropic 2D background after an onset distance.
    # Background is zero before the onset distance.


    # Add pair correlation distribution for repeated localisations.
    rpd = rpd + (replocsamp* models.pairwise_correlation_2d(r,0.,np.sqrt(2) * replocssd))

    # Add pair correlation distribution for unresolvable substructure/
    # mislocalisations of simultaneous nearby emitters.

    rpd = rpd + (substructamp* models.pairwise_correlation_2d(r,0.,np.sqrt(2) *substructsd))


    return rpd


def rot_sym_with_replocs_and_substructure_isotropic_bg_with_onset_vectorargs(
        input_vector):
    """Function to calculate the values given by
    rot_sym_with_replocs_and_substructure_isotropic_bg_with_onset, but using a
    vector input for the parameters, so that the numdifftools package can be
    used to calculate partial derivatives for correct error propagation in the
    model.

    Args:
        input_vector (list or numpy array):
            A concatenation of:
                1. A distance at which density values of the model will be
                obtained (numpy array)

                2. The parameters used by
                rot_sym_with_replocs_and_substructure_isotropic_bg_with_onset.
    Returns:
        rpd (numpy array):
            The relative position density given by the model at the input
            distances (called separation_values_1d).
    """
    (separation_values_1d,
     dia,
     vertssd, vertsamp,
     replocssd, replocsamp,
     substructsd, substructamp) = input_vector

    rpd = rot_sym_replocs_substructure_isotropic_bg_with_onset(
        separation_values_1d,
        dia,
        vertssd, vertsamp,
        replocssd, replocsamp,
        substructsd, substructamp)

    return rpd


def create_default_fitting_params_dicts():
    """Create lists for initial guesses and bounds for parameters during
    fitting. Sets the defaults for scipy.optimise.curve_fit

    Returns:
        lower_bound_dict (dict):
            Dictionary of default lower bound options for parameter values.
        upper_bound_dict (dict):
            Dictionary of default upper bound options for parameter values.
        initial_params_dict (dict):
            Dictionary of default initial parameter value options.
    """
    lower_bound_dict = {'diameter': 0,
                        'vertices_sd': 0,
                        'vertices_amp': 0,
                        'loc_prec_sd': 0,
                        'loc_prec_amp': 0,
                        'substruct_sd': 0,
                        'substruct_amp': 0,
                        #'bg_slope': -10,
                        #'bg_onset': 0.005,
                        }

    upper_bound_dict = {'diameter': binnorm + 10,
                        'vertices_sd': binnorm,
                        'vertices_amp': binnorm,
                        'loc_prec_sd': binnorm,
                        'loc_prec_amp': binnorm,
                        'substruct_sd': binnorm,
                        'substruct_amp': binnorm,
                        #'bg_slope': 100.,
                        #'bg_onset': 2*binnorm,
                        }

    initial_params_dict = {'diameter': 80,
                           'vertices_sd': 4,
                           'vertices_amp': 0.1,
                           'loc_prec_sd': 2.5,
                           'loc_prec_amp': 2,
                           'substruct_sd': 20,
                           'substruct_amp':0.1,
                           #'bg_slope': -10, #hete = -0.005 homo 0.1
                           #'bg_onset': 2*binnorm/3,
                           }

    return lower_bound_dict, upper_bound_dict, initial_params_dict


def set_up_model_replocs_substruct_iso_bg_with_onset_with_fit_settings():
    """Set up the RPD model with fitting settings,
    for a rotationally symmetric model with spread due to repeated
    localisations, spread to unresolvable substructure, and a background
    that increases linearly after an onset distance.
    The fitting settings are to pass to scipy's
    curve_fit, and the vector-input version of the model is for
    differentiation and error propagation with numdifftools.

    Args:
        None

    Returns:
        A ModelWithFitSettings object containing:
            model_rpd (function name):
                Relative position density as a function of separation
                between localisations.
            initial_params (list):
                Starting guesses for the parameter values by
                scipy.optimize.curve_fit
            lower_bounds (list), upper_bounds (list):
                The bounds on allowable parameter values as
                scipy.optimize.curve_fit runs.
    """
    # Generate ModelWithFitSettings object, conatining a model_rpd
    model_replocs_substruct_iso_bg_with_onset_with_fit_settings = (
        ModelWithFitSettings(
            model_rpd=rot_sym_replocs_substructure_isotropic_bg_with_onset
            )
        )

    # Add fitting parameters to ModelWithFitSettings object
    (lower_bound_dict,
     upper_bound_dict,
     initial_params_dict) = create_default_fitting_params_dicts()

    # Can optionally modify these dictionaries here:

    initial_params = [initial_params_dict['diameter'],
                      initial_params_dict['vertices_sd'],
                      initial_params_dict['vertices_amp'],
                      initial_params_dict['loc_prec_sd'],
                      initial_params_dict['loc_prec_amp'],
                      initial_params_dict['substruct_sd'],
                      initial_params_dict['substruct_amp'],
                      #initial_params_dict['bg_slope'],
                      #initial_params_dict['bg_onset']
                      ]

    lower_bounds = [lower_bound_dict['diameter'],
                    lower_bound_dict['vertices_sd'],
                    lower_bound_dict['vertices_amp'],
                    lower_bound_dict['loc_prec_sd'],
                    lower_bound_dict['loc_prec_amp'],
                    lower_bound_dict['substruct_sd'],
                    lower_bound_dict['substruct_amp'],
                    #lower_bound_dict['bg_slope'],
                    #lower_bound_dict['bg_onset']
                    ]

    upper_bounds = [upper_bound_dict['diameter'],
                    upper_bound_dict['vertices_sd'],
                    upper_bound_dict['vertices_amp'],
                    upper_bound_dict['loc_prec_sd'],
                    upper_bound_dict['loc_prec_amp'],
                    upper_bound_dict['substruct_sd'],
                    upper_bound_dict['substruct_amp'],
                    #upper_bound_dict['bg_slope'],
                    #upper_bound_dict['bg_onset']
                    ]

    bounds = (lower_bounds, upper_bounds)

    model_replocs_substruct_iso_bg_with_onset_with_fit_settings.initial_params = (
        initial_params
        )
    model_replocs_substruct_iso_bg_with_onset_with_fit_settings.param_bounds = (
        bounds
        )
    model_replocs_substruct_iso_bg_with_onset_with_fit_settings.vector_input_model = (
        rot_sym_with_replocs_and_substructure_isotropic_bg_with_onset_vectorargs
        )

    return model_replocs_substruct_iso_bg_with_onset_with_fit_settings



# Set up callable symmetry order object, with default value
sym_order = Number(40)
#the number of bins that is normalized
binnorm = 100

# handle the input arguments (flags)
prog = 'Fit_droplets_protein_distribution'
prog_short_name = 'fdpd'
description = ('Fits rotational symmetry models to relative positions '
                   'of proteins at the interface.')

info = {'prog': prog,
            'prog_short_name': prog_short_name,
            'description': description}

start = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')

info['start'] = start

parser = argparse.ArgumentParser(prog, description)

parser.add_argument('-i', '--input_file',
                    dest='input_file',
                    type=argparse.FileType('r'),
                    help='File of localisations which is a .csv (or .txt '
                         'with comma delimiters) or .npy and containing N '
                         'localisations in N rows.',
                    metavar="FILE")
parser.add_argument('-f', '--filter_distance',
                    dest='filter_dist',
                    type=int,
                    # default=150, PROBLEM LINE
                    default=binnorm, help="Filter distance.")
parser.add_argument('-s', '--short_names',action="store_true") 
parser.add_argument('-v', '--verbose', help="Increase output verbosity",action="store_true")
args = parser.parse_args()
print("\nargs: ", args)
if args.verbose:
    print("Verbosity is turned on.\n")
info['filter_dist'] = binnorm
info['verbose'] = args.verbose    
info['short_names'] = args.short_names
if args.input_file is None:
    get_inputs(info)
else:
    info['in_file_and_path'] = args.input_file.name

info['host'], info['ip_address'], info['operating_system'] = utils.find_hostname_and_ip()


utils.secondary_filename_and_path_setup(info)
in_file = info['in_file_and_path']

if not os.path.exists(in_file):
    sys.exit("ERROR; The input file does not exist.")

if in_file[-4:] == '.csv':
    try:
        line = open(in_file).readline()
    except (EOFError, IOError, OSError) as exception:
        print("\n\nCould not open file: ", in_file)
        print("\n\n", type(exception))
        sys.exit("Could not open the input file "+in_file+".\n")
    if (line.__contains__("ch,sum")):
        skip = 1
        try:
            xyz_values = np.loadtxt(in_file, delimiter=',', skiprows=skip)
        except (EOFError, IOError, OSError) as exception:
            print("\n\nCould not read file: ", in_file)
            print("\n\n", type(exception))
            sys.exit("Could not read the input file "+in_file+".\n")
    else:
        xyz_values = 'Ouch'
        print('Sorry, wrong format! This program needs a file output from '
              'relative_positions.py\n')
        sys.exit("The input file "+in_file+" has the wrong format. It needs "
                 "a file output form relative_positions\n")
else:
        xyz_values = 'Ouch'
        print('Sorry, wrong format! This program needs a file output from '
              'relative_positions.py\n')
        sys.exit("The input file "+in_file+" has the wrong format. It needs "
                 "a file output form relative_positions\n")

fitlength = info['filter_dist']
bin_values = xyz_values[:,0]
bin_values = np.append(bin_values,binnorm + 0.00)
xy_histogram = xyz_values[:,1]
   

if info['short_names'] is True:
    try:
        os.makedirs(info['short_results_dir'])
    except OSError:
        print("Unexpected error:", sys.exc_info()[0])
        sys.exit("Could not create directory for the results.")    
else:
    try:
        os.makedirs(info['results_dir'])
    except OSError:
        print("Unexpected error:", sys.exc_info()[0])
        sys.exit("Could not create directory for the results.")

# Define symmetries over which to perform and evaluate fit
symmetries = range(35,36)

    # Prepare to record fit metric (AICc)
aiccs = np.zeros(len(symmetries))


# Define model, including parameter guesses and bounds to pass to
# scipy.optimize.curve_fit
model_with_info = (
        set_up_model_replocs_substruct_iso_bg_with_onset_with_fit_settings()
        )
info['model_name'] = model_with_info.model_rpd.__name__


info['p0'] = model_with_info.initial_params
info['optimisation_bounds'] = model_with_info.param_bounds

curve_values = []
diameter_values = []
table_param_values = []

for i, sym in enumerate(symmetries):

        sym_order.number = sym

        (params_optimised,
         params_covar,                                                                                                                                                                                                                                                                                                                                   
         params_1sd_error) = models.fit_model_to_experiment(xy_histogram,
                                                            model_with_info.model_rpd,
                                                            model_with_info.initial_params,
                                                            model_with_info.param_bounds,
                                                            fitlength=binnorm)
        print("working on...")
        aicc = stats.aic_from_least_sqr_fit(xy_histogram,
                                            model_with_info.model_rpd,
                                            params_optimised,
                                            fitlength=binnorm)[1]
        aiccs[i] = aicc
        x_values = np.linspace(0.5,100,int(binnorm))
        diameter_values.append(params_optimised[0])

        curve_values.append(model_with_info.model_rpd(x_values, *params_optimised))


        if info['verbose']:
            print('\nSymmetry:', sym)
            print("\nParameter    Optimised Value")

        param_string_values = []
        uncertainty_string_values = []

        for count, param_optimised in enumerate(params_optimised):
            uncertainty = params_1sd_error[count]

            value = param_optimised

            value_str, uncertainty_str = utils.plus_and_minus(value, uncertainty)

            param_string_values.append(value_str)
            uncertainty_string_values.append(uncertainty_str)

            if info['verbose']:
                print('  %d            %s +/- %s ' % (count+1, value_str, uncertainty_str))


        params = np.column_stack((params_optimised,
                                  params_1sd_error,
                                  param_string_values,
                                  uncertainty_string_values))


        table_param_values.append(params)
bin_values = bin_values[0:int(binnorm)]
plotting.plot_histogram_with_curves(bin_values,
                                    xy_histogram,
                                    symmetries,
                                    x_values,
                                    curve_values,
                                    info)

# =============================================================================
# for i, sym in enumerate(symmetries):
#     sym_order.number = sym
#     plotting.plot_rot_2d_geometry(sym, diameter_values[i], info)
# =============================================================================


weights = stats.akaike_weights(aiccs)

print('\nSymmetry\tAICc\t\tAkaike weight')

for index, symmetry in enumerate(symmetries):
            print('{0:2d}\t\t{1:.2f} \t{2:.2}'.format(symmetry,
                                                      aiccs[index],
                                                      weights[index]))
nbin = bin_values[0:100]
ncurv = np.hstack(curve_values[0])                                                     
col1 = "Bins_exp"
col2 = "Histogram_exp"
col3 = "Bins_fit"
col4 = "Histogram_fit"
col5 = "Optimised_parameter"
data = pd.DataFrame({col1:nbin,col2:xy_histogram,col3:x_values,col4:ncurv})
data2 = pd.DataFrame({col5:params_optimised})
data.to_csv(info['results_dir'] + r'Fit_result.csv',  index=False)
data2.to_csv(info['results_dir'] + r'Params.csv',  index=False)
reports.write_rot_2d_html_report(info, symmetries, aiccs, weights, table_param_values)

#with open('readme.txt', 'w') as f:
    
#    f.write(f'C value: {popt[2]}\n')
#    f.write('\n')
#    f.close()



    