try:
    from .sampling import cobra_sampling
    from .sampling_dependencies import *
    from .sampling_statisticsDescriptive import sampling_statisticsDescriptive
    from .sampling_statisticsSampledPoints import sampling_statisticsSampledPoints
except ImportError as e:
    from sampling.sampling_dependencies import *
    from sampling.sampling import cobra_sampling
    from sampling.sampling_statisticsDescriptive import sampling_statisticsDescriptive
    from sampling.sampling_statisticsSampledPoints import sampling_statisticsSampledPoints


class cobra_sampling_statistics(cobra_sampling):

    # statistics
    def descriptive_statistics(self,points_I='flux'):
        '''calculate the following:
        1. mean, variance, 95% CI
        2. median, mode, 1st quartile, 3rd quartile, range

        Args
            points_I (st): point type, i.e., 'flux','metabolite','subsystem', default='flux'
        '''

        statisticsDescriptive = sampling_statisticsDescriptive()
        statisticsSampledPoints = sampling_statisticsSampledPoints()
        points_statistics = {};

        if points_I == 'flux': points = self.points;
        elif points_I == 'metabolite': points = self.points_metabolite;
        elif points_I == 'subsystem': points = self.points_subsystem;
        else: 
            print('points not recognized: default "flux" will be used.')
            points = self.points;

        for k,v in points.items():
            # calculate the mean and variance
            n = len(points[k]);
            m,var,lb,ub = statisticsDescriptive.calculate_ave_var(points[k],confidence_I = 0.95);
            # directly calculate the 95% CI
            lb,ub = statisticsSampledPoints.calculate_ciFromPoints(points[k],alpha=0.05)
            #lb,ub = statisticsSampledPoints.bootstrap(points[k]['points'], num_samples=100000, statistic=np.mean, alpha=0.05)
            # calculate the min, max, median, and interquartile ranges
            min,max,median,iq_1,iq_3=statisticsDescriptive.calculate_interquartiles(points[k],iq_range_I = [25,75])
            tmp = {};
            tmp = {
                'n':n,
                'ave':m,
                'var':var,
                'lb':lb,
                'ub':ub,
                'median':median,
                'min':min,
                'max':max,
                'iq_1':iq_1,
                'iq_3':iq_3
                }
            points_statistics[k]=tmp;
        self.points_statistics = points_statistics;

class cobra_sampling_n():

    def __init__(self,data_dir_I=None,model_I=None,loops_I=[],sample_ids_I=[],samplers_I=None,control_I=False):
        #   control_I = True: sample_ids_I[0]=control,sample_ids_I[1:]=perturbation
        #               False: pairwise test is performed on all
        #               controls how the pairwisetests are performed
        if data_dir_I:self.data_dir =  data_dir_I;
        else: self.data_dir = None;
        if model_I: self.model = model_I;
        else: self.model = None;
        if loops_I: self.loops = loops_I;
        else: self.loops = [];
        if sample_ids_I: self.sample_ids = sample_ids_I;
        else: self.sample_ids = []
        if samplers_I: self.samplers = samplers_I;
        else: self.samplers = [];
        if control_I: self.control = control_I;
        else: self.control=control_I;
        self.points = [];
        self.points_metabolites = [];
        self.points_subsystems = [];
        self.data = [];

    def calculate_pairWiseTest(self,
        redundancy_I=False):
        '''Calculate the difference of the mean and median, fold change, and p-value
         between pairwise sampling distribution'''

        statisticsSampledPoints = sampling_statisticsSampledPoints()
        data_O = [];

        test_stat = None
        test_description = 'permutation'

        for i,data_1 in enumerate(self.points):
            sample_id_1 = self.sample_ids[i];
            if self.control and i>0:
                break;
            for j,data_2 in enumerate(self.points):
                sample_id_2 = self.sample_ids[j];
                if i==j: continue;
                if not redundancy_I and i>j: continue;
                for r in self.model.reactions:
                    rxn_id = r.id
                    mean_difference = None
                    median_difference = None
                    pvalue = None
                    fold_change = None;
                    # check point length
                    n_data_1 = len(data_1);
                    n_data_2 = len(data_2);
                    if n_data_1 != n_data_2:
                        print('number of sampled points are not consistent')
                    n_points = n_data_1;
                    if rxn_id in list(data_1.keys()):
                        cond1 = np.array(data_1[rxn_id]);
                    else:
                        cond1 = np.zeros(n_points);
                    if rxn_id in list(data_2.keys()):
                        cond2 = np.array(data_2[rxn_id]);
                    else:
                        cond2 = np.zeros(n_points);
                    mean_difference = cond2.mean()-cond1.mean();
                    median_difference = np.median(cond2) - np.median(cond1) 
                    #pvalue = statisticsSampledPoints.permutation_resampling(cond1,cond2);
                    pvalue = statisticsSampledPoints.calculate_pvalue_permutation(cond1,cond2);
                    fold_change = np.log(np.exp(cond2.mean())/np.exp(cond1.mean()));
                    tmp = {};
                    tmp = {'sample_id_1':sample_id_1,
                           'sample_id_2':sample_id_2,
                           'rxn_id':rxn_id,
                           'pvalue':pvalue,
                           'mean_difference':mean_difference,
                           'median_difference':median_difference,
                           'fold_change':fold_change,
                           'test_stat':test_stat,
                           'test_description':test_description};
                    data_O.append(tmp);
        self.data = data_O;

    def calculate_pairWiseTest_metabolites(self,
                    redundancy_I=False,):
        '''Calculate the difference of the mean and median, fold change, and p-value
         between pairwise sampling distribution for metabolites'''

        data_O = [];

        test_stat = None
        test_description = 'permutation'

        for i,data_1 in enumerate(self.points_metabolites):
            sample_id_1 = self.sample_ids[i];
            if self.control and i>0:
                break;
            for j,data_2 in enumerate(self.points_metabolites):
                sample_id_2 = self.sample_ids[j];
                if i==j: continue;
                if not redundancy_I and i>j: continue;
                for m in self.model.metabolites:
                    met_id = m.id
                    mean_difference = None
                    median_difference = None
                    pvalue = None
                    fold_change = None;
                    # check point length
                    n_data_1 = len(data_1);
                    n_data_2 = len(data_2);
                    if n_data_1 != n_data_2:
                        print('number of sampled points are not consistent')
                    n_points = n_data_1;
                    if met_id in list(data_1.keys()):
                        cond1 = np.array(data_1[met_id]);
                    else:
                        cond1 = np.zeros(n_points);
                    if met_id in list(data_2.keys()):
                        cond2 = np.array(data_2[met_id]);
                    else:
                        cond2 = np.zeros(n_points);
                    mean_difference = cond1.mean() - cond2.mean()
                    median_difference = np.median(cond2) - np.median(cond1) 
                    #pvalue = statisticsSampledPoints.permutation_resampling(cond1,cond2);
                    pvalue = statisticsSampledPoints.calculate_pvalue_permutation(cond1,cond2);
                    fold_change = np.log(np.exp(cond2.mean())/np.exp(cond1.mean()));
                    tmp = {};
                    tmp = {'sample_id_1':sample_id_1,
                           'sample_id_2':sample_id_2,
                           'met_id':met_id,
                           'pvalue':pvalue,
                           'mean_difference':mean_difference,
                           'median_difference':median_difference,
                           'fold_change':fold_change,
                           'test_stat':test_stat,
                           'test_description':test_description};
                    data_O.append(tmp);
        self.data = data_O;

    def calculate_pairWiseTest_subsystems(self,
                    redundancy_I=False,):
        '''Calculate the difference of the mean and median, fold change, and p-value
         between pairwise sampling distribution for metabolites'''

        data_O = [];

        test_stat = None
        test_description = 'permutation'

        for i,data_1 in enumerate(self.points_subsystems):
            sample_id_1 = self.sample_ids[i];
            if self.control and i>0:
                break;
            for j,data_2 in enumerate(self.points_subsystems):
                sample_id_2 = self.sample_ids[j];
                if i==j: continue;
                if not redundancy_I and i>j: continue;
                for sub_id in list(data_1.keys()):
                    mean_difference = None
                    median_difference = None
                    pvalue = None
                    fold_change = None;
                    # check point length
                    n_data_1 = len(data_1);
                    n_data_2 = len(data_2);
                    if n_data_1 != n_data_2:
                        print('number of sampled points are not consistent')
                    n_points = n_data_1;
                    if sub_id in list(data_1.keys()):
                        cond1 = np.array(data_1[sub_id]);
                    else:
                        cond1 = np.zeros(n_points);
                    if sub_id in list(data_2.keys()):
                        cond2 = np.array(data_2[sub_id]);
                    else:
                        cond2 = np.zeros(n_points);
                    mean_difference = cond1.mean() - cond2.mean()
                    median_difference = np.median(cond2) - np.median(cond1) 
                    #pvalue = statisticsSampledPoints.permutation_resampling(cond1,cond2);
                    pvalue = statisticsSampledPoints.calculate_pvalue_permutation(cond1,cond2);
                    fold_change = np.log(np.exp(cond2.mean())/np.exp(cond1.mean()));
                    tmp = {};
                    tmp = {'sample_id_1':sample_id_1,
                           'sample_id_2':sample_id_2,
                           'subsystem_id':sub_id,
                           'pvalue':pvalue,
                           'mean_difference':mean_difference,
                           'median_difference':median_difference,
                           'fold_change':fold_change,
                           'test_stat':test_stat,
                           'test_description':test_description};
                    data_O.append(tmp);
        self.data = data_O;

    def get_points(self,data_points_I=[],remove_loops_I=True,remove_no_flux_I=True,normalize_I=True,remove_points_I=[]):
        '''Get multiple points from sampling'''
        sampling = cobra_sampling();
        for i,sample_id in enumerate(self.sample_ids):
            sampling.add_data(data_dir_I=self.data_dir,model_I=self.model,loops_I=self.loops);
            if self.samplers[i]=='gpSampler':
                sampling.get_points_matlab(data_points_I[i]);
            elif self.samplers[i]=='optGpSampler':
                sampling.get_points_json(data_points_I[i])
            if remove_loops_I: sampling.remove_loopsFromPoints();
            if remove_no_flux_I: sampling.remove_noFluxReactionsFromPoints();
            if normalize_I: sampling.normalize_points2Total();
            if remove_points_I:sampling.remove_points(remove_points_I[i]);
            #sampling.convert_points2MetabolitePoints();
            #sampling.convert_points2SubsystemPoints()
            self.points.append(sampling.points);
            #self.points_metabolites.append(sampling.points_metabolite)
            #self.points_subsystems.append(sampling.points_subsystem)
            sampling.remove_data();     

    def convert_points2MetabolitePoints(self):
        '''Get multiple points from sampling'''
        sampling = cobra_sampling();
        for i,sample_id in enumerate(self.sample_ids):
            sampling.add_data(model_I=self.model,points_I=self.points[i]);
            sampling.convert_points2MetabolitePoints();
            self.points_metabolites.append(sampling.points_metabolite)
            sampling.remove_data();            

    def convert_points2SubsystemPoints(self):
        '''Get multiple points from sampling'''
        sampling = cobra_sampling();
        for i,sample_id in enumerate(self.sample_ids):
            sampling.add_data(model_I=self.model,points_I=self.points[i]);
            sampling.convert_points2SubsystemPoints()
            self.points_subsystems.append(sampling.points_subsystem)
            sampling.remove_data();
