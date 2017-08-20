
try:
    from .sampling import cobra_sampling
    from .sampling_dependencies import *
except ImportError as e:
    from sampling.sampling_dependencies import *
    from sampling.sampling import cobra_sampling

class optGpSampler_sampling(cobra_sampling):
    def mix_fraction(self, sample1, sample2, **kwargs):
        """ Compares two sets of sampled points and determines how mixed
        they are.

        Args
         sample1, sample2   Ordered set of points, numpy arrays.  The points must be in
                           the same order otherwise it does not make sense.
        kwargs
         fixed (optional)   The directions which are fixed and are not expected (indices)

        Returns
         mix                the mix fraction.  Goes from 0 to 1 with 1 being
                            completely unmixed and .5 being essentially 
                            perfectly mixed.  


        """
        from numpy import min, isnan, median, outer, ones
    
        if 'fixed' not in kwargs:
            fixed = []
        else:
            fixed = kwargs['fixed']
        

        # ignore NAN rows
        ignore_rows = isnan(min(sample1,1)) | isnan(min(sample2,1))
        if len(fixed) > 0:
            ignore_rows[fixed] = True
        keep_rows = ~ ignore_rows

        sample1_reduced = sample1[keep_rows,:]
        sample2_reduced = sample2[keep_rows,:]

        m1 = median(sample1_reduced, 1)
        LPproblem = median(sample2_reduced, 1)
        n_rxn_reduced, n_points = sample1_reduced.shape

        l1 = (sample1_reduced > (outer(m1, ones([1, n_points]))))
        eq1 = (sample1_reduced == outer(m1, ones([1, n_points])))
        l2 = (sample2_reduced > outer(LPproblem, ones([1, n_points])))
        eq2 = (sample2_reduced == outer(LPproblem, ones([1, n_points])))

        eqtotal = eq1 | eq2

        fr_mix = float(sum(sum((l1 == l2) & (~ eqtotal))))/float(l1.size-sum(sum(eqtotal)))

        return fr_mix

    def calculate_mixFraction(self):
        '''Convert points and warmup points to numpy arrays
        and calculate the mixed fraction
        '''
        samples1_lst = []
        samples2_lst = []
        for k in self.points.keys():
            samples1_lst.append(self.warmup[k]);
            samples2_lst.append(self.points[k]);
        samples1 = np.array(samples1_lst)
        samples2 = np.array(samples2_lst)
        self.mixed_fraction = self.mix_fraction(samples1,samples2);

    def export_sampling_optGpSampler(self,
            cobra_model,
            fraction_optimal = None, 
            filename_model='sample_model.json',
            filename_script='sample_script.py', 
            filename_points='points.json',
            filename_warmup='warmup.json',
            solver_id_I='glpk',
            n_points_I = None, 
            n_steps_I = 20000,
            n_threads_I = 128,
            verbose_I = 1,
            ):
        '''export model and script for sampling using optGpSampler

        Args
            cobra_model (cobra.Model)
            fraction_optimal (float): optional fraction of optimal to set the objective value
            filename_model (str): name of the model file
            filename_script (str): name of the script to run the sampling algorithm
            filename_points (str): name of the resulting points file
            filename_warmup (str): name of the resulting warmup points file
            solver_id (str): name of the solver to use for sampling
            n_points_I (int): number of sampling points
            n_steps_I (int): number of sampling steps
            n_threads_I (int): number of processor threads to dedicate to samples
            verbose_I (int): OptGpSampler option to print out results to the console
        '''

        if n_points_I:
            n_points = n_points_I;
        if n_steps_I:
            n_steps = n_steps_I;
        if n_threads_I:
            n_threads = n_threads_I;
        if solver_id_I=='cglpk':
           solver_id = 'glpk';
        else:
           solver_id = solver_id_I;
        # export the model as json
        save_json_model(cobra_model,self.data_dir + '/' + filename_model);
        # make the sampling script
        script = 'from cobra_utilities.sampling_dependencies import *\n';
        script += 'from cobra_utilities.optGpSampler_sampling import optGpSampler_sampling\n';
        script += 'cobra_model = load_json_model("%s");\n'%filename_model;        
        script += 'sampling = optGpSampler_sampling(data_dir_I = "/");\n';
        script += 'sampling.generate_samples(cobra_model=cobra_model,';
        script += 'filename_points = "%s",'%filename_points;
        script += 'filename_warmup = "%s",'%filename_warmup;
        script += 'solver_id_I = "%s",'%solver_id;
        script += 'n_points_I = %s,'%n_points;
        script += 'n_steps_I = %s,'%n_steps;
        script += 'n_threads_I = %s);\n'%n_threads;
        # export the sampling script
        with open(self.data_dir + '/' + filename_script,'w') as file:
            file.write(script);
            file.close();

    def generate_samples(self,
            cobra_model,
            fraction_optimal = None, 
            filename_points='points.json',
            filename_warmup='warmup.json',
            solver_id_I='glpk',
            n_points_I = None, 
            n_steps_I = 20000,
            n_threads_I = 8,
            verbose_I = 1,
            reduce_model_I = False
            ):
        '''sample the model using optGpSampler
        
        Args
            cobra_model (cobra.Model)
            fraction_optimal (float): optional fraction of optimal to set the objective value
            filename_model (str): name of the model file
            filename_points (str): name of the resulting points file
            filename_warmup (str): name of the resulting warmup points file
            solver_id (str): name of the solver to use for sampling
            n_points_I (int): number of sampling points
            n_steps_I (int): number of sampling steps
            n_threads_I (int): number of processor threads to dedicate to samples
            verbose_I (int): OptGpSampler option to print out results to the console
            reduce_model_I (boolean): OptGpSampler option to reduce model before sampling
        '''

        # confine the objective to a fraction of maximum optimal
        if fraction_optimal:
            # optimize
            cobra_model.optimize(solver_id_I);
            objective = [x.id for x in cobra_model.reactions if x.objective_coefficient == 1]
            cobra_model.reactions.get_by_id(objective[0]).upper_bound = fraction_optimal * cobra_model.objective.value;

        model = cobra_model.to_array_based_model(deepcopy_model=True)
        model = CbModel.convertPyModel(model)

        # reduce the model
        if reduce_model_I:
            rm = CbModelReducer(model)
            rm.fixStoichiometricMatrix(0)
            rm.setTolerance(1e-6)
            rm.setSolverName(solver_id_I)
            rm.setVerbose(verbose_I)
            model = rm.reduce()

        # get warmup points:
        sampler = CbModelSampler(model)
        sampler.setNrSamples(n_points_I)
        sampler.setNrSteps(n_steps_I)
        sampler.setSolverName(solver_id_I)
        sampler.setNrThreads(n_threads_I)
        sampler.setVerbose(1)
        sModel = sampler.sample()
        warmup = sModel.samplePts
        print('done warmup');

        #sample
        sampler = CbModelSampler(model)
        sampler.setNrSamples(n_points_I)
        sampler.setNrSteps(n_steps_I)
        sampler.setSolverName(solver_id_I)
        sampler.setNrThreads(n_threads_I)
        sampler.setWarmupPts(warmup)
        sampler.setVerbose(verbose_I)

        sModel = sampler.sample()
        samples = sModel.samplePts
        print('done sampling');
        LBS = np.tile(sModel.lower_bounds, (n_points_I,1))    
        UBS = np.tile(sModel.upper_bounds, (n_points_I,1))   
        
        max_dev_null = abs(sModel.S * samples).max();
        max_dev_lb = max((LBS.T - samples).max(), 0);
        max_dev_ub = max((samples - UBS.T).max(), 0)
        print("Maximum deviation from the nullspace = " + str(max_dev_null));
        print("Maximum violation of lb = " + str(max_dev_ub));
        print("Maximum violation of ub = " + str(max_dev_lb));

        points_dict = {};
        points_dict = {k:list(samples[i,:]) for i,k in enumerate(model.reactions)};
        warmup_dict = {};
        warmup_dict = {k:list(warmup[i,:]) for i,k in enumerate(model.reactions)};
        self.points = points_dict;
        self.warmup = warmup_dict;
        self.export_points_json(filename_points);
        self.export_warmup_json(filename_warmup);

        #samples = pd.DataFrame(data=samples.T, columns=[i for i in model.reactions])
        #dump(samples, open(filename_points,'w'));