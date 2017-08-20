import pytest

# Dependencies from cobra
from cobra.io.json import load_json_model
from cobra.manipulation.modify import convert_to_irreversible

# Other dependencies
import csv,json,sys

# Dependencies
from cobra_sampling.sampling_statistics import cobra_sampling_statistics

from . import data_dir

class test_sampling_statistics():
    
    def __init__(self):
        self.cobra_model = load_json_model(data_dir + '/mini.json')

    def test_descriptive_statistics(self):

        simulation_id_I = 'test_sampling'
        filename_points = simulation_id_I + '_points' + '.json';
        filename_warmup = simulation_id_I + '_warmup' + '.json';
        sampling = cobra_sampling_statistics(
            data_dir_I = data_dir,
            model_I=self.cobra_model);
        sampling.get_points_json(filename_points);
        sampling.get_warmup_json(filename_warmup);
        sampling.calculate_mixFraction();
        assert(len(sampling.points) == 31)
        assert(sampling.mixed_fraction == 1.0) #need to update
        # check if the model contains loops
        sampling.simulate_loops(
            data_fva = data_dir + 'test_loops_fva.json',
            solver_I = 'glpk');
        sampling.find_loops(data_fva = data_dir + 'test_loops_fva.json');
        assert('ENO' in sampling.loops)
        sampling.remove_loopsFromPoints();
        assert(len(sampling.points) == 1)
        assert('EX_glc__D_e' in sampling.points.keys())
        # calculate the flux descriptive statistics
        sampling.descriptive_statistics(points_I='flux');
        assert(sampling.points_statistics['EX_glc__D_e']['n'] == 62)
        assert(sampling.points_statistics['EX_glc__D_e']['ave'] == 1.9176474254840106)
        sampling.convert_points2MetabolitePoints();
        assert('glc__D_e' in sampling.points_metabolite.keys())
        # calculate descriptive stats for metabolites
        sampling.descriptive_statistics(points_I='metabolite');
        assert(sampling.points_statistics['glc__D_e']['n'] == 62)
        assert(sampling.points_statistics['glc__D_e']['ave'] == 0.95882371274200529)
        sampling.convert_points2SubsystemPoints();
        assert('' in sampling.points_subsystem.keys())
        # calculate descriptive stats for subsystems
        sampling.descriptive_statistics(points_I='subsystem');
        assert(sampling.points_statistics['']['n'] == 62)
        assert(sampling.points_statistics['']['ave'] == 1.9176474254840106)
