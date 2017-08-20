try:
    from .sampling_dependencies import *
except ImportError as e:
    from sampling.sampling_dependencies import *

class sampling_statisticsDescriptive():

    # calculate the geometric mean and variance:
    def calculate_ave_var_geometric(self,data_I):
        """ calculate the geometric average and var of data
        with 95% confidence intervals
        """

        try:
            data_ave_O = 0.0
            # calculate the average of the sample
            for c in data_I:
                data_ave_O += np.log(c);
            data_ave_O = np.exp(data_ave_O/len(data_I));

            data_var_O = 0.0
            #calculate the variance of the sample
            for c in data_I:
                data_var_O += np.power(np.log(c/data_ave_O),2);
            data_var_O = data_var_O/(len(data_I)-1); #note: we will need to take the exp()
                                                     # to get the final variance
                                                     # but leaving it this way makes the
                                                     # downstream calculations simpler

            #calculate the 95% confidence intervals
            data_se = np.sqrt(data_var_O/len(data_I));
            data_lb_O = np.exp(np.log(data_ave_O) - 1.96*data_se);
            data_ub_O = np.exp(np.log(data_ave_O) + 1.96*data_se);
            
            #correct the variance for use in reporting
            data_var_O = np.exp(data_var_O);

            return data_ave_O, data_var_O, data_lb_O, data_ub_O;
        except Exception as e:
            print(e);
            exit(-1);
    # calculate the mean and variance:
    def calculate_ave_var(self,data_I,confidence_I = 0.95):
        """calculate the average and var of data
        with 95% confidence intervals"""

        try:
            data = np.array(data_I);

            data_ave_O = 0.0
            # calculate the average of the sample
            data_ave_O = np.mean(data);

            data_var_O = 0.0
            #calculate the variance of the sample
            data_var_O = np.var(data);

            #calculate the standard error of the sample
            se = scipy.stats.sem(data)

            #calculate the 95% confidence intervals
            n = len(data);
            h = se * scipy.stats.t._ppf((1+confidence_I)/2., n-1)
            data_lb_O = data_ave_O - h;
            data_ub_O = data_ave_O + h;

            return data_ave_O, data_var_O, data_lb_O, data_ub_O;
        except Exception as e:
            print(e);
            exit(-1);
    # calculate the mean and CV: 
    def calculate_ave_var_cv(self,data_I,confidence_I = 0.95):
        """calculate the average, var, %cv of data
        with 95% confidence intervals"""

        try:
            data = np.array(data_I);

            data_ave_O = 0.0
            # calculate the average of the sample
            data_ave_O = np.mean(data);

            data_var_O = 0.0
            #calculate the variance of the sample
            data_var_O = np.var(data);

            #calculate the standard error of the sample
            se = scipy.stats.sem(data)

            #calculate the CV% of the sample
            data_cv_O = 0.0;
            if data_ave_O !=0.0:
                data_cv_O = np.std(data)/data_ave_O*100;

            #calculate the 95% confidence intervals
            n = len(data);
            h = se * scipy.stats.t._ppf((1+confidence_I)/2., n-1)
            data_lb_O = data_ave_O - h;
            data_ub_O = data_ave_O + h;

            return data_ave_O, data_var_O, data_cv_O, data_lb_O, data_ub_O;
        except Exception as e:
            print(e);
            exit(-1);
    # convert CV to SD:
    def convert_cv2StDev(self,data_ave_I,data_cv_I):
        """convert the %CV to standard deviation
        INPUT
        data_ave_I = float, average/mean
        data_cv_I = float, % coefficient of variation or % relative standard deviation
        OUTPUT
        data_stdev_O = float, standard deviation"""

        try:
            data_stdev_O = 0.0;
            #calculate the SD of the sample
            if not data_ave_I is None and not data_cv_I is None:
                data_stdev_O = data_ave_I*data_cv_I/100.0;

            return data_stdev_O;
        except Exception as e:
            print(e);
            exit(-1);
    def convert_ci2StDev(self,lb_I,ub_I,n_I,alpha_I,distribution_I='lognorm',dist_params_I={},raise_I = False):
        '''
        Convert confidence intervals to standard deviation

        Args:
            lb_I = float, lower bound
            ub_I = float, upper bound
            n_I = int, number of replicates
            alpha_I = float, significance level
            distribution_I = string, scipy distribution function
            dist_params_I = distribution function parameters

        Returns:
            stdev_O = float

        Notes:
            see http://docs.scipy.org/doc/scipy/reference/stats.html
            for a list of distributions and parameters

            alpha must be included in the dist_params_I dictionary
            based on the name given for the particular function
        '''
        
        from scipy import stats

        stdev_O = None;
        try:
            if hasattr(stats, distribution_I):
                distribution_func = getattr(stats, distribution_I);
                stderr = (ub_I-lb_I)/(2*distribution_func.isf(**dist_params_I));
                stdev_O = np.sqrt(n_I)*stderr;
            else: print('scipy.stats does not support ' + distribution_I);
        except Exception as e:
            if raise_I: raise;
            else: print(e);
        return stdev_O;

    # calculate the interquartiles
    def calculate_interquartiles(self,data_I,iq_range_I = [25,75]):
        '''compute the interquartiles and return the min, max, median, iq1 and iq3'''
        try:
            min_O = np.min(data_I);
            max_O = np.max(data_I);
            iq_1_O, iq_2_O = np.percentile(data_I, iq_range_I)
            median_O = np.median(data_I);

            return min_O, max_O, median_O, iq_1_O, iq_2_O;
        except Exception as e:
            print(e)
    # calculate the fold change
    def calculate_foldChange(self,data_1_I,data_2_I,type_I = 'relative', scale_values_I = None,scale_fold_change_I = None):
        """Calculate the fold change between two data value 1 and 2
        Input:
        data_1_I = data value 1
        data_2_I = data value 2
        type_I = 'relative', 'absolute','geometric' (relative is the default);
        scale_values_I = string, factor to scale values by prior to calculating the fold change
        scale_fold_change_I = string factor to scale the fold change by
        """
        # scale the values
        if data_1_I and data_2_I:
            fold_change_O = 0.0;
        elif data_1_I==0.0:
            print("data_1_I cannot be zero!");
            return None;
        else:
            print("bad data provided!");
            return None;
        supported_scale_values = ["log2","log10","ln","abs","exp","exp2","^10","^2"]
        if scale_values_I and scale_values_I in supported_scale_values: #TODO: test
            data_1 = self.scale_values(data_1_I,scale_values_I);
            data_2 = self.scale_values(data_2_I,scale_values_I);
        elif scale_values_I:
            print('scale_values_I not recognized.  No scaling will be applied');
            data_1 = data_1_I;
            data_2 = data_2_I;
        else:
            data_1 = data_1_I;
            data_2 = data_2_I;
        # relative, absolute, or geometric
        if type_I == 'relative':
            fold_change_O = data_2/data_1;
        elif type_I == 'absolute':
            fold_change_O = np.abs(data_2/data_1);
        elif type_I == 'geometric':
            fold_change_O = np.log(np.exp(data_2)/np.exp(data_1));
        else:
            fold_change_O = data_2/data_1;
            print("type_I not recognized.  Relative type will be applied as default.");
        # scale the fold change
        supported_scale_fold_change = ["log2","log10","ln","abs","sqrt"]
        if scale_fold_change_I and scale_fold_change_I == "log2":
            fold_change_O = self.scale_values(fold_change_O,scale_fold_change_I);
        elif scale_fold_change_I:
            print("scale_fold_change_I not recognized.  No scaling will be applied.");
            fold_change_O = fold_change_O; 
        else:
            fold_change_O = fold_change_O;
        #check for nan
        if np.isnan(fold_change_O):
            fold_change_O = 0.0;
        return fold_change_O;
    def scale_values(self,data_1_I,scale_I=None):
        '''Scale values
        INPUT:
        data_1_I = float, value to scale
        scale_I = string, type of scale
        OUTPUT:
        data_O = float, scaled version of data_1
        '''
        data_O = None;
        if scale_I == "log2":
            data_O = np.log2(data_1_I);
        elif scale_I == "log10":
            data_O = np.log10(data_1_I);
        elif scale_I == "ln":
            data_O = np.log(data_1_I);
        elif scale_I == "abs":
            data_O = np.abs(data_1_I);
        elif scale_I == "exp":
            data_O = np.exp(data_1_I);
        elif scale_I == "exp2":
            data_O = np.exp2(data_1_I);
        elif scale_I == "^10":
            data_O = np.power(data_1_I,10);
        elif scale_I == "^2":
            data_O = np.power(data_1_I,2);
        elif scale_I == "sqrt":
            data_O = np.sqrt(data_1_I);
        elif scale_I:
            print("scale_I not recognized.  No scaling will be applied.");
            data_O = data_1_I;
        else:
            print("scale_I not recognized.  No scaling will be applied.");
            data_O = data_1_I;
        return data_O;
    # calculate the difference
    def calculate_difference(self,data_1_I,data_2_I,type_I = 'relative',scale_values_I = None,scale_difference_I = None):
        """Calculate the differencefold change between two data value 1 and 2
        Input:
        data_1_I = data value 1
        data_2_I = data value 2
        type_I = 'relative', 'absolute','geometric' (relative is the default);
        scale_values_I = string, factor to scale values by prior to calculating the difference
        scale_difference_I = string factor to scale the fold change by
        """
        if data_1_I and data_2_I:
            difference_O = np.array(0.0);
            data_1_I = np.array(data_1_I);
            data_2_I = np.array(data_2_I);
        else:
            print("bad data provided!");
            return None;
        if scale_values_I and scale_values_I == "log2":
            data_1 = np.log2(data_1_I);
            data_2 = np.log2(data_2_I);
        elif scale_values_I and scale_values_I == "log10":
            data_1 = np.log10(data_1_I);
            data_2 = np.log10(data_2_I);
        elif scale_values_I and scale_values_I == "ln":
            data_1 = np.log(data_1_I);
            data_2 = np.log(data_2_I);
        elif scale_values_I and scale_values_I == "abs":
            data_1 = np.abs(data_1_I);
            data_2 = np.abs(data_2_I);
        elif scale_values_I and scale_values_I == "exp":
            data_1 = np.exp(data_1_I);
            data_2 = np.exp(data_2_I);
        elif scale_values_I and scale_values_I == "exp2":
            data_1 = np.exp2(data_1_I);
            data_2 = np.exp2(data_2_I);
        elif scale_values_I and scale_values_I == "^10":
            data_1 = np.power(data_1_I,10);
            data_2 = np.power(data_2_I,10);
        elif scale_values_I and scale_values_I == "^2":
            data_1 = np.power(data_1_I,2);
            data_2 = np.power(data_2_I,2);
        elif scale_values_I:
            print("scale_values_I not recognized.  No scaling will be applied.");
            data_1 = data_1_I;
            data_2 = data_2_I;
        else:
            data_1 = data_1_I;
            data_2 = data_2_I;
        # relative, absolute, or geometric
        if type_I == 'relative':
            difference_O = data_2-data_1;
        elif type_I == 'absolute':
            difference_O = np.abs(data_2-data_1);
        elif type_I == 'geometric': #i.e. distance
            difference_O = np.log(np.exp(data_2)+np.exp(data_1));
        elif type_I == 'euclidean': #i.e. distance
            difference_O = np.sqrt(np.power(data_2,2)+np.power(data_1,2));
        else:
            difference_O = data_2-data_1;
            print("type_I not recognized.  Relative type will be applied as default.");
        
        if scale_difference_I and scale_difference_I == "log2":
            difference_O = np.log2(difference_O);
        elif scale_difference_I and scale_difference_I == "log10":
            difference_O = np.log10(difference_O);
        elif scale_difference_I and scale_difference_I == "ln":
            difference_O = np.log(difference_O);
        elif scale_difference_I and scale_difference_I == "abs":
            difference_O = np.abs(difference_O);
        elif scale_difference_I and scale_difference_I == "sqrt":
            difference_O = np.sqrt(difference_O);
        elif scale_difference_I:
            print("scale_difference_I not recognized.  No scaling will be applied.");
            difference_O = difference_O; 
        else:
            difference_O = difference_O;

        if type(difference_O)==type(np.array(0.0)) and \
            len(difference_O)==1:
            return difference_O[0];
        else:
            return difference_O;

    def calculate_descriptiveStats(self,
                data_I,
                confidence_I = 0.95,
                iq_range_I = [25,75]):
        '''calculate the mean, var, cv, 95% CI,
        min, max, median, and IQ ranges for the data
        INPUT:
        data_I = array of data points
        OUTPUT:
        descriptiveStats = {} with fields
            'n',
            'mean'
            'var'
            'cv'
            'ci'
            'lb'
            'ub'
            'min'
            'max'
            'median'
            'iq_1'
            'iq_3',
            'iq'
        '''
        descriptiveStats_O = None;
        try:
            mean,var,cv,lb,ub = self.calculate_ave_var_cv(data_I,confidence_I =confidence_I);
            min, max, median, iq_1, iq_3 = self.calculate_interquartiles(data_I,iq_range_I=iq_range_I);
            descriptiveStats_O = {
                'n':len(data_I),
                'mean':mean,
                'var':var,
                'cv':cv,
                'ci':confidence_I,
                'lb':lb,
                'ub':ub,
                'min':min,
                'max':max,
                'median':median,
                'iq_1':iq_1,
                'iq_3':iq_3,
                'iq':iq_range_I};
        except Exception as e:
            print(e);
            descriptiveStats_O = {
                'n':len(data_I),
                'mean':None,
                'var':None,
                'cv':None,
                'ci':confidence_I,
                'lb':None,
                'ub':None,
                'min':None,
                'max':None,
                'median':None,
                'iq_1':None,
                'iq_3':None,
                'iq':iq_range_I};
        return descriptiveStats_O;
