try:
    from .sampling import cobra_sampling
    from .sampling_dependencies import *
except ImportError as e:
    from sampling.sampling_dependencies import *
    from sampling.sampling import cobra_sampling

class matlab_sampling(cobra_sampling):
    
    def get_points_matlab(self,matlab_data=None,sampler_model_name='sampler_out'):
        '''load sampling points from MATLAB
        
        Args
            matlab_data (str): name of the matlab data file
            sampler_model_name (str): name of the model file
        '''

        # extract information about the file
        import os, time
        from datetime import datetime
        from stat import ST_SIZE, ST_MTIME
        if matlab_data:
            filename=self.data_dir + '/' + matlab_data;
        else:
            filename=self.data_dir;
        try:
            st = os.stat(filename)
        except IOError:
            print("failed to get information about", filename)
            return;
        else:
            file_size = st[ST_SIZE]
            simulation_dateAndTime_struct = time.localtime(st[ST_MTIME])
            simulation_dateAndTime = datetime.fromtimestamp(time.mktime(simulation_dateAndTime_struct))

        # load model from MATLAB file
        try:
            model = load_matlab_model(filename,sampler_model_name);
        except NotImplementedError as e:
            print(e);
            model_tmp = h5py.File(filename,'r')['sampler_out'];
            #model = matlab_cobra_struct_to_python_cobra_object(matlab_struct)
            model = self.model;

        # load sample points from MATLAB file into numpy array
        try:
            points = scipy.io.loadmat(filename)[sampler_model_name]['points'][0][0];
            mixed_fraction=scipy.io.loadmat(filename)['mixedFrac'][0][0];
        except NotImplementedError as e:
            print(e);
            points = h5py.File(filename,'r')[sampler_model_name]['points'];
            points = np.array(points);
            mixed_fraction=h5py.File(filename,'r')['mixedFrac'][0][0];
        #mat = scipy.io.loadmat('data/EvoWt.mat')
        #points = mat['model_WT_sampler_out']['points'][0][0]

        points_dict = {};
        for i,r in enumerate(model.reactions):
            # convert names:
            r_id_conv = r.id.replace('-','_DASH_');
            r_id_conv = r_id_conv.replace('(','_LPAREN_');
            r_id_conv = r_id_conv.replace(')','_RPAREN_');
            # extract points
            points_dict[r_id_conv]=points[i,:];

        self.points = points_dict;
        self.model = model;
        self.mixed_fraction = mixed_fraction;
        self.simulation_dateAndTime = simulation_dateAndTime;
        
    def export_sampling_matlab(self,cobra_model,fraction_optimal = None,
        filename_model='sample_model.mat',
        filename_script='sample_script.m',
        filename_points='points.mat',
        solver_id_I='glpk',n_points_I = None, n_steps_I = 20000, max_time_I = None):
        '''export model and script for sampling using matlab cobra_toolbox
        
        Args
            cobra_model (cobra.Model)
            fraction_optimal (float): optional fraction of optimal to set the objective value
            filename_model (str): name of the model file
            filename_script (str): name of the script to run the sampling algorithm
            filename_points (str): name of the resulting points file
            solver_id (str): name of the solver to use for sampling
            n_points_I (int): number of sampling points
            n_steps_I (int): number of sampling steps
            max_time_I (float): maximum sampling time
            '''
        ## copy the model:
        #cobra_model_copy = cobra_model.copy();
        # convert numerical input to string
        n_points,n_steps,max_time = '[]','[]','[]';
        if n_points_I:
            n_points = n_points_I;
        if n_steps_I:
            n_steps = n_steps_I;
        if max_time_I:
            max_time = max_time_I;
        # confine the objective to a fraction of maximum optimal
        if fraction_optimal:
            # optimize
            cobra_model.solver = solver_id_I
            cobra_model.optimize();
            objective = [x.id for x in cobra_model.reactions if x.objective_coefficient == 1]
            cobra_model.reactions.get_by_id(objective[0]).upper_bound = fraction_optimal * cobra_model.objective.value;
        # write model to mat
        cobra_model.description = 'model';
        save_matlab_model(cobra_model,self.data_dir + '/' + filename_model);
        ## write model to xml
        #write_sbml_model(cobra_model,'data/sampling/sampler.xml');
        # write the sampling script to file
        #[sampleStructOut, mixedFraction] = gpSampler(sampleStruct, nPoints, bias, maxTime, maxSteps, threads, nPointsCheck)
        mat_script = "% initialize with " + solver_id_I + "\n";
        mat_script +="load('" + self.data_dir + '/' + filename_model + "')\n";
        mat_script +="initCobraToolbox();\n";
        mat_script +="% sample\n";
        mat_script +="[sampler_out, mixedFrac] = gpSampler(" + cobra_model.description
        mat_script +=(", %s, [], %s, %s, [], true);\n" %(n_points,n_steps,max_time));
        mat_script +="[sampler_out, mixedFrac] = gpSampler(sampler_out";
        mat_script +=(", %s, [], %s, %s, [], true);\n" %(n_points,n_steps,max_time));
        mat_script +="save('"+ self.data_dir + '/' + filename_points + "','sampler_out', 'mixedFrac');";
        with open(self.data_dir + '/' + filename_script,'w') as f:
            f.write(mat_script);