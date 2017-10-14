try:
    from .sampling_dependencies import *
except ImportError as e:
    from sampling.sampling_dependencies import *

class cobra_sampling():

    def __init__(self,data_dir_I=None,model_I=None,loops_I=[]):
        if data_dir_I:self.data_dir =  data_dir_I;
        else: self.data_dir = None;
        if model_I: self.model = model_I;
        else: self.model = None;
        if loops_I: self.loops = loops_I;
        else: self.loops = [];
        self.points = {};
        self.warmup = {};
        self.points_metabolite = {};
        self.points_subsystem = {};
        self.mixed_fraction = None;
        self.points_statistics = {};

    def get_points_numpy(self,numpy_data):
        '''load sampling points from numpy file
        
        Args
            numpy_data (str): name of the file
        '''

        # extract information about the file
        import os, time
        from datetime import datetime
        from stat import ST_SIZE, ST_MTIME
        try:
            st = os.stat(self.data_dir + '/' + numpy_data)
        except IOError:
            print("failed to get information about", self.data_dir + '/' + numpy_data)
            return;
        else:
            file_size = st[ST_SIZE]
            simulation_dateAndTime_struct = time.localtime(st[ST_MTIME])
            simulation_dateAndTime = datetime.fromtimestamp(time.mktime(simulation_dateAndTime_struct))

        # load points from numpy file
        points = loadtxt(numpy_data);

        points_dict = {};
        for i,r in enumerate(self.cobra_model.reactions):
            # extract points
            points_dict[r_id_conv]=points[i,:];

        self.points = points_dict;
        self.simulation_dateAndTime = simulation_dateAndTime;

    def get_points_json(self,json_data=None):
        '''load sampling points from json file
        
        Args
            json_data (str): name of the file
        '''

        # extract information about the file
        import os, time
        from datetime import datetime
        from stat import ST_SIZE, ST_MTIME
        if json_data:
            filename=self.data_dir + '/' + json_data;
        else:
            filename=self.data_dir;

        try:
            st = os.stat(filename)
            file_size = st[ST_SIZE]
            simulation_dateAndTime_struct = time.localtime(st[ST_MTIME])
            simulation_dateAndTime = datetime.fromtimestamp(time.mktime(simulation_dateAndTime_struct))
        except IOError:
            print("failed to get information about", filename)
            return;

        points_dict = json.load(open(filename));
        points_dict = {k.split(':')[0]:v for k,v in points_dict.items()}
        self.points = points_dict;
        self.simulation_dateAndTime = simulation_dateAndTime;

    def get_warmup_json(self,json_data=None):
        '''load sampling warmup points from json file
        
        Args
            json_data (str): name of the file
        '''

        # extract information about the file
        import os, time
        from datetime import datetime
        from stat import ST_SIZE, ST_MTIME
        if json_data:
            filename=self.data_dir + '/' + json_data;
        else:
            filename=self.data_dir;

        try:
            st = os.stat(filename)
            file_size = st[ST_SIZE]
            simulation_dateAndTime_struct = time.localtime(st[ST_MTIME])
            simulation_dateAndTime = datetime.fromtimestamp(time.mktime(simulation_dateAndTime_struct))
        except IOError:
            print("failed to get information about", filename)
            return;

        warmup_dict = json.load(open(filename));
        warmup_dict = {k.split(':')[0]:v for k,v in warmup_dict.items()}
        self.warmup = warmup_dict;

    def export_points_numpy(self,filename):
        '''export sampling points
        
        Args
            filename (str): name of the file
        '''
        savetxt(filename,self.points);

    def export_points_json(self,filename):
        '''export sampling points'''
        with open(filename,'w') as file:
            json.dump(self.points,file, indent=4);

    def export_warmup_json(self,filename):
        '''export sampling warmup points
        
        Args
            filename (str): name of the file
        '''
        with open(filename,'w') as file:
            json.dump(self.warmup,file, indent=4);

    def plot_points_histogram(self,reaction_lst=[]):
        '''plot sampling points as a histogram
        
        Args:
            reaction_lst (list(str)): list of reaction ids
            
        '''
        if not reaction_lst:
            reaction_lst = ['ENO','FBA','FBP','G6PP','GAPD','GLBRAN2',
                        'GLCP','GLCP2','GLDBRAN2','HEX1','PDH','PFK',
                        'PGI','PGK','PGM','PPS','PYK','TPI','ENO_reverse',
                        'FBA_reverse','GAPD_reverse','PGI_reverse',
                        'PGK_reverse','PGM_reverse','TPI_reverse']
        for r in reaction_lst:
            # loop through each reaction in the list
            plt.figure()
            n, bins, patches = plt.hist(self.points[r],50,label = [r])
            plt.legend()
            plt.show()

    def plot_points_boxAndWhiskers(self,reaction_lst=[]):
        '''plot sampling points as box and whiskers plots
        
        Args:
            reaction_lst (list(str)): list of reaction ids
        '''
        if not reaction_lst:
            reaction_lst = ['ENO','FBA','FBP','G6PP','GAPD','GLBRAN2',
                        'GLCP','GLCP2','GLDBRAN2','HEX1','PDH','PFK',
                        'PGI','PGK','PGM','PPS','PYK','TPI','ENO_reverse',
                        'FBA_reverse','GAPD_reverse','PGI_reverse',
                        'PGK_reverse','PGM_reverse','TPI_reverse']
        for r in reaction_lst:
            # loop through each reaction in the list
            plt.figure()
            fig, ax = plt.subplots()
            bp = ax.boxplot(self.points[r], sym='k+',
                            notch=False, bootstrap=False,
                            usermedians=None,
                            conf_intervals=None)
            plt.show()

    def check_loops(self,cobra_model_I=None,solver_I = 'glpk'):
        '''Check if the model contains loops
        
        Args:
            cobra_model_I (cobra.Model)
            solver_I (str): solver name
        '''

        # Change all uptake reactions to 0
        if cobra_model_I: cobra_model = cobra_model_I.copy();
        else: cobra_model = self.model.copy();
        system_boundaries = [x.id for x in cobra_model.reactions if x.boundary == 'system_boundary'];
        for rxn in cobra_model.reactions:
            if rxn.id in system_boundaries:
                cobra_model.reactions.get_by_id(rxn.id).lower_bound = 0.0;
                cobra_model.reactions.get_by_id(rxn.id).upper_bound = 0.0;
        # Set ATPM to 0
        cobra_model.reactions.get_by_id('ATPM').lower_bound = 0.0
        # set the objective function to a default value
        for k,v in linear_reaction_coefficients(cobra_model).items():
            cobra_model.reactions.get_by_id(k.id).lower_bound=0.0
            cobra_model.reactions.get_by_id(k.id).upper_bound=1e6

        loops_bool = True;
        cobra_model.solver = solver_I
        cobra_model.optimize();
        if not cobra_model.objective.value:
            loops_bool = False;

        return loops_bool;

    def simulate_loops(self,cobra_model_I=None,data_fva='loops_fva.json',solver_I='glpk'):
        '''Simulate FVA after closing exchange reactions and setting ATPM to 0
        reactions with flux will be involved in loops
        
        Args:
            cobra_model_I (cobra.Model)
            data_fva (str): name of the output file
            solver_I (str): solver name
        '''
        
        if cobra_model_I: cobra_model = cobra_model_I.copy();
        else: cobra_model = self.model.copy();
        # Change all uptake reactions to 0
        system_boundaries = [x.id for x in cobra_model.reactions if x.boundary == 'system_boundary'];
        for rxn in cobra_model.reactions:
            if rxn.id in system_boundaries:
                cobra_model.reactions.get_by_id(rxn.id).lower_bound = 0.0;
                cobra_model.reactions.get_by_id(rxn.id).upper_bound = 0.0;
        # Set ATPM to 0
        cobra_model.reactions.get_by_id('ATPM').lower_bound = 0.0;
        # set the objective function to a default value
        for k,v in linear_reaction_coefficients(cobra_model).items():
            cobra_model.reactions.get_by_id(k.id).lower_bound=0.0
            cobra_model.reactions.get_by_id(k.id).upper_bound=1e6

        try:
            # calculate the reaction bounds using FVA
            cobra_model.solver = solver_I
            reaction_bounds = flux_variability_analysis(cobra_model, fraction_of_optimum=1.0,
                                              the_reactions=None);
        except Exception as e:
            print(e);
            reaction_bounds = '';

        # Update the data file
        with open(data_fva, 'w') as outfile:
            json_data = dict(zip(list(reaction_bounds.index),reaction_bounds.to_dict('records')))
            json.dump(json_data, outfile, indent=4);

    def find_loops(self,data_fva='loops_fva.json'):
        '''extract out loops from simulate_loops
        
        Args
            data_fva (str): name of the input file
        '''

        data_loops = json.load(open(data_fva))
        rxn_loops = [];
        if data_loops:
            for k,v in data_loops.items():
                if abs(v['minimum'])>1.0 or abs(v['maximum'])>1.0:
                    rxn_loops.append(k);
        #return rxn_loops
        self.loops = rxn_loops;

    def remove_loopsFromPoints(self):
        '''remove reactions with loops from sampling points'''

        points_loopless = {};
        for k,v in self.points.items():
            if k in self.loops: continue
            else: 
                points_loopless[k] = v;

        #return points_loopless_mean;
        self.points = points_loopless;

    def remove_noFluxReactionsFromPoints(self):
        '''remove reactions that carry 0 flux'''

        points_flux = {};
        for k,v in self.points.items():
            # determine the max/min of the data
            max_point = max(v);
            min_point = min(v);
            if max_point == 0.0 and min_point == 0.0: continue;
            else: 
                points_flux[k] = v;

        self.points = points_flux;
        return

    def remove_points(self,rxn_ids_I=[]):
        '''
        remove points
        Args
            rxn_ids_I (list(str)): list of reaction ids to remove
        '''

        for rxn_id in rxn_ids_I:
            self.points.pop(rxn_id);

    # def remove_points_notInSolutionSpace_v1(self):
    #     '''remove points that are not in the solution space
        
    #     '''
    #     pruned_reactions = [];
    #     rxn_ids_noPointsInSolutionSpace = [];
    #     for rxn in self.model.reactions:
    #         points_copy = copy(self.points[rxn.id])
    #         self.points[rxn.id] = self.remove_points_notInBounds(self.points[rxn.id],rxn.lower_bound,rxn.upper_bound);
    #         if len(self.points[rxn.id])<1:
    #             print("no points found in the solution space for rxn_id " + rxn.id + "!");
    #             rxn_ids_noPointsInSolutionSpace.append(rxn.id);
    #         if len(points_copy)!= len(self.points[rxn.id]):
    #             pruned_reactions.append(rxn.id)
    #     return pruned_reactions

    def remove_points_notInSolutionSpace(self,min_points_I=1000):
        '''remove points that are not in the solution space

        Args
            min_points_I (int): minimum number of points
        
        '''
        pruned_reactions = [];
        for rxn in self.model.reactions:
            points_copy = copy(self.points[rxn.id]);
            self.points[rxn.id] = self._remove_points_notInSolutionSpace(points_copy,rxn.lower_bound,rxn.upper_bound,min_points_I);
        return pruned_reactions

    def _remove_points_notInSolutionSpace(self,points_I,lower_bound_I,upper_bound_I,min_points_I):
        '''remove points that are not in the solution space.
        If the minimum number of points is not found, the bounds will be increased
        by +/- (upper_bound_I-lower_bound_I)/4 until the minimum number of points
        is found

        Args
            points_I (list(float)): list of points
            lower_bound_I (float)
            upper_bound_I (float)
            min_points_I (int): minimum number of points
        
        '''
        points_O=[];
        points_O=self.remove_points_notInBounds(points_I,lower_bound_I,upper_bound_I);
        if len(points_O)<min_points_I:
            adjust_bounds = (upper_bound_I-lower_bound_I)/4;
            if adjust_bounds == 0.0:
                adjust_bounds = 1.0;
            lower_bound_new = lower_bound_I-adjust_bounds;
            upper_bound_new = upper_bound_I+adjust_bounds;
            points_O=self._remove_points_notInSolutionSpace(points_I,lower_bound_new,upper_bound_new,min_points_I);
        return points_O;

    def remove_points_notInBounds(self,points_I,lower_bound_I,upper_bound_I):
        '''remove points not in the lower/upper bounds

        Args
            points_I (list(float)): list of points
            lower_bound_I (float)
            upper_bound_I (float)
        '''
        points_O = [p for p in points_I if p >= lower_bound_I and p<= upper_bound_I];
        return points_O;        

    def convert_points2MetabolitePoints(self):
        '''convert the reaction flux to total flux through each metabolite for each sampling point'''
        metabolite_points = {};
        first_loop = True;
        for k,v in self.points.items():
            if first_loop:
                for met in self.model.metabolites:
                    metabolite_points[met.id]=np.zeros_like(v);
                first_loop = False;
            for i,flux in enumerate(v):
                for p in self.model.reactions.get_by_id(k).products:
                    metabolite_points[p.id][i]+=0.5*abs(flux*self.model.reactions.get_by_id(k).get_coefficient(p.id))
                for p in self.model.reactions.get_by_id(k).reactants:
                    metabolite_points[p.id][i]+=0.5*abs(flux*self.model.reactions.get_by_id(k).get_coefficient(p.id))
        self.points_metabolite=metabolite_points;

    def convert_points2SubsystemPoints(self):
        '''convert the reaction flux to total flux through each subsystem for each sampling point'''
        subsystem_points = {};
        subsystems_all = [];
        for r in self.model.reactions:
            subsystems_all.append(r.subsystem);
        subsystems = list(set(subsystems_all)); 
        first_loop = True;
        for k,v in self.points.items():
            if first_loop:       
                for sub in subsystems:
                    subsystem_points[sub]=np.zeros_like(v);
                first_loop = False;
            for i,flux in enumerate(v):
                subsystem_points[self.model.reactions.get_by_id(k).subsystem][i]+=abs(flux);
        self.points_subsystem=subsystem_points;

    def normalize_points2Total(self):
        '''normalize each reaction for a given point to the total
        flux through all reactions for that point'''
        points = self.points;
        points_normalized={};
        n_points = len(points[list(points.keys())[0]]);
        for i in range(n_points):
            total=0.0;
            for k,v in points.items():
                total+=np.abs(v[i]);
            for k,v in points.items():
                if i==0:
                    points_normalized[k]=np.zeros_like(v);
                points_normalized[k][i] = v[i]/total
        self.points = points_normalized;

    def normalize_points2CarbonInput(self):
        '''normalize each reaction for a given point to the total
        carbon input flux for that point'''
        points = self.points;
        system_boundaries = [x.id for x in self.model.reactions if x.boundary == 'system_boundary'];
        points_normalized={};
        n_points = len(points[list(points.keys())[0]]);
        for i in range(n_points):
            total=0.0;
            for k,v in points.items():
                if k in system_boundaries:
                    if self.model.reactions.get_by_id(k).reactants and v[i] < 0:
                        # e.g. glc-D -->
                        mets = self.model.reactions.get_by_id(k).reactants
                        for met in mets:
                            formula_str = met.formula.formula
                            n12C = 0
                            if 'C' not in Formula(formula_str)._elements and 0 in Formula(formula_str)._elements['C']:
                                n12C += Formula(formula_str)._elements['C'][0]; #get the # of Carbons
                            total+=np.abs(v[i])*n12C;
                    elif self.model.reactions.get_by_id(k).products and v[i] > 0:
                        # e.g. --> glc-D
                        mets = self.model.reactions.get_by_id(k).reactants
                        for met in mets:
                            formula_str = met.formula.formula
                            n12C = 0
                            if 'C' not in Formula(formula_str)._elements and 0 in Formula(formula_str)._elements['C']:
                                n12C += Formula(formula_str)._elements['C'][0]; #get the # of Carbons
                            total+=np.abs(v[i])*n12C;
            for k,v in points.items():
                if i==0:
                    points_normalized[k]=np.zeros_like(v);
                points_normalized[k][i] = v[i]/total
        self.points = points_normalized;

    def normalize_points2Input(self, rxn_ids=[]):
        '''normalize each reaction for a given point to the total
        input flux for that point'''
        points = self.points;
        system_boundaries = [x.id for x in self.model.reactions if x.boundary == 'system_boundary'];
        points_normalized={};
        n_points = len(points[list(points.keys())[0]]);
        for i in range(n_points):
            total=0.0;
            for k,v in points.items():
                if k in rxn_ids:
                    total+=np.abs(v[i]);
                elif k in system_boundaries:
                    if self.model.reactions.get_by_id(k).reactants and v[i] < 0:
                        # e.g. glc-D -->
                        total+=np.abs(v[i]);
                    elif self.model.reactions.get_by_id(k).products and v[i] > 0:
                        # e.g. --> glc-D
                        total+=np.abs(v[i]);
            for k,v in points.items():
                if i==0:
                    points_normalized[k]=np.zeros_like(v);
                points_normalized[k][i] = v[i]/total
        self.points = points_normalized;

    def add_data(self,data_dir_I=None,model_I=None,loops_I=[],points_I={}):
        '''add new data'''
        if data_dir_I:self.data_dir =  data_dir_I;
        if model_I: self.model = model_I;
        if loops_I: self.loops = loops_I;
        if points_I: self.points = points_I;

    def remove_data(self):
        '''remove all data'''
        self.data_dir = None;
        self.model = None;
        self.points = {};
        self.mixed_fraction = None;
        self.loops = {};