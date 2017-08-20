import pytest

# Dependencies from cobra
from cobra.io.json import load_json_model
from cobra.manipulation.modify import convert_to_irreversible

# Other dependencies
import csv,json,sys

# Dependencies
from cobra_utilities.optGpSampler_sampling import optGpSampler_sampling
from cobra_utilities.matlab_sampling import matlab_sampling

from . import data_dir

class test_sampling():
    
    def __init__(self):
        self.cobra_model = None

    def test_cobra_model(self):
        cobra_model = load_json_model(data_dir + '/mini.json')
        solution = cobra_model.optimize()
        assert(solution.objective_value == 30.0)
        assert(solution.fluxes['ENO'] == 20.0)
        self.cobra_model = cobra_model
    
    def test_export_sampling_optGpSampler(self):

        sampling = optGpSampler_sampling(data_dir_I = data_dir);
        simulation_id_I = 'test_sampling'
        filename_model = simulation_id_I + '.json';
        filename_script = simulation_id_I + '.py';
        filename_points = simulation_id_I + '_points' + '.json';
        filename_warmup = simulation_id_I + '_warmup' + '.json';

        sampling.export_sampling_optGpSampler(cobra_model=self.cobra_model,
            filename_model=filename_model,
            filename_script=filename_script,
            filename_points=filename_points,
            filename_warmup=filename_warmup,
            solver_id_I = 'glpk',
            n_points_I = 2*len(self.cobra_model.reactions),
            n_steps_I = 5000,
            n_threads_I = 2)
    
    def test_sampling_optGpSampler(self):

        simulation_id_I = 'test_sampling'
        filename_points = simulation_id_I + '_points' + '.json';
        filename_warmup = simulation_id_I + '_warmup' + '.json';
        sampling = optGpSampler_sampling(
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
        sampling.convert_points2MetabolitePoints();
        assert('glc__D_e' in sampling.points_metabolite.keys())
        sampling.convert_points2SubsystemPoints();
        assert('' in sampling.points_subsystem.keys())
    
    def test_export_sampling_matlab(self):

        sampling = matlab_sampling(data_dir_I = data_dir);
        simulation_id_I = 'test_sampling'
        filename_model = simulation_id_I + '.mat';
        filename_script = simulation_id_I + '.m';
        filename_points = simulation_id_I + '_points' + '.mat';

        sampling.export_sampling_matlab(cobra_model=self.cobra_model,
            filename_model=filename_model,
            filename_script=filename_script,
            filename_points=filename_points,
            solver_id_I = 'glpk',
            n_points_I = 2*len(self.cobra_model.reactions),
            n_steps_I = 5000)            
    
    def test_sampling_matlab(self):

        simulation_id_I = 'test_sampling'
        filename_model = simulation_id_I + '.mat';
        filename_points = simulation_id_I + '_points' + '.mat';
        sampling = matlab_sampling(
            data_dir_I = data_dir,
            model_I=self.cobra_model);
        sampling.get_points_matlab(
            matlab_data=filename_points,
            sampler_model_name=filename_model);
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
        sampling.convert_points2MetabolitePoints();
        assert('glc__D_e' in sampling.points_metabolite.keys())
        sampling.convert_points2SubsystemPoints();
        assert('' in sampling.points_subsystem.keys())

    def test_get_points_numpy(self):
        """ """
        pass

    def test_export_points_numpy(self):
        """ """
        pass

    def test_plot_points_histogram(self):
        """ """
        pass

    def test_plot_points_boxAndWhiskers(self):
        """ """
        pass

    def test_remove_noFluxReactionsFromPoints(self):
        """ """
        pass

    def test_remove_points(self):
        """ """
        pass

    def test_remove_points_notInSolutionSpace(self):
        """ """
        pass

    def test_remove_points_notInBounds(self):
        """ """
        pass

    def test_normalize_points2Total(self):
        """ """
        pass

    def test_normalize_points2CarbonInput(self):
        """ """
        pass

    def test_normalize_points2Input(self):
        """ """
        pass
    