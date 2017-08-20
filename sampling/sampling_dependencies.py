# System
from sys import exit
from math import log, sqrt, exp
import operator, json, csv
from copy import copy
# Dependencies from 3rd party
try:
    import h5py
except Exception as e:
    print(e);
try:
    import scipy.io
    import numpy as np
    from numpy import histogram, mean, std, loadtxt, savetxt
    import pandas as pd
except Exception as e:
    print(e);
try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
except Exception as e:
    print(e);
# Resources
try:
    from molmass.molmass import Formula
except Exception as e:
    print(e);
# Dependencies from cobra
try:
    from cobra.io.mat import load_matlab_model,save_matlab_model
    # from cobra.mlab import matlab_cobra_struct_to_python_cobra_object
    # from cobra.util.solver import linear_reaction_coefficients
    from cobra.io.sbml import create_cobra_model_from_sbml_file
    from cobra.io.json import load_json_model, save_json_model
    from cobra.flux_analysis import flux_variability_analysis, single_deletion
    from cobra.flux_analysis.parsimonious import optimize_minimal_flux
    #from cobra.flux_analysis.objective import update_objective
except Exception as e:
    print(e);
# Dependencies from optGpSampler
try:
    from optGpSampler.cbModel import CbModel
    from optGpSampler.cbModelSampler import CbModelSampler
    from optGpSampler.cbModelReducer import CbModelReducer
except Exception as e:
    print(e);