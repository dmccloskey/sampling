try:
    from .sampling_dependencies import *
except ImportError as e:
    from sampling.sampling_dependencies import *

class sampling_statisticsSampledPoints():
    # calculate the confidence intervals
    def calculate_ciFromPoints(self,data_I, alpha=0.05):
        """Calculate the confidence intervals from sampled points"""
        data_sorted = np.sort(data_I)
        n = len(data_sorted)
        lb = data_sorted[int((alpha/2.0)*n)]
        ub = data_sorted[int((1-alpha/2.0)*n)]
        return lb,ub
    def bootstrap(self,data, num_samples=100000, statistic=np.mean, alpha=0.05):
        """Returns bootstrap estimate of 100.0*(1-alpha) CI for statistic."""
        n = len(data)
        idx = npr.randint(0, n, (num_samples, n))
        samples = data[idx]
        stat = np.sort(statistic(samples, 1))
        return (stat[int((alpha/2.0)*num_samples)],
                stat[int((1-alpha/2.0)*num_samples)])
    # calculate the p-value difference
    def permutation_resampling(self,case, control, num_samples=50, statistic=np.mean):
        '''calculate the pvalue of two data sets using a resampling approach'''

        observed_diff = abs(statistic(case) - statistic(control))
        num_case = len(case)

        combined = np.concatenate([case, control])
        diffs = []
        for i in range(num_samples):
            xs = npr.permutation(combined)
            diff = np.mean(xs[:num_case]) - np.mean(xs[num_case:])
            diffs.append(diff)

        pval = (np.sum(diffs > observed_diff) +
                np.sum(diffs < -observed_diff))/float(num_samples)
        return pval, observed_diff, diffs
    def calculate_pvalue_permutation(self,data_1_I,data_2_I,n_permutations_I=10,n_resamples_I=10):    
        '''calculate the pvalue of two data by determining
        the lack of overlap between sample points using a permutation test.
        If the sample points of the data sets is not equal,
        a subset of samples of matching length is resampled from the larger data set'''
        
        data_1 = None; #sample set with fewer points
        data_2 = None; #sample set with more points
        n_resamples = 0;

        # check the length of data_1 and data_2
        if len(data_1_I)>len(data_2_I):
            data_1=np.array(data_2_I);
            data_2=np.array(data_1_I);
            n_resamples=n_resamples_I;
        elif len(data_1_I)<len(data_2_I):
            data_1=np.array(data_1_I);
            data_2=np.array(data_2_I);
            n_resamples=n_resamples_I;
        else:
            data_1=np.array(data_1_I);
            data_2=np.array(data_2_I);

        n_samples_min = len(data_1);

        vals = []
        for i in range(0,n_permutations_I):
            if n_resamples==0:
                cond1 = np.random.permutation(data_1)
                cond2 = np.random.permutation(data_2)
                z = cond1 - cond2
                x = len(z[z>0]) + 1
                y = len(z[z<0]) + 1
                k = min(x,y)
                vals.append(k)
            else:
                cond1 = np.random.permutation(data_1)
                cond2 = np.random.permutation(data_2)
                for resample in range(n_resamples):
                    cond2_int = np.random.randint(0,n_samples_min);
                    z = cond1 - cond2[cond2_int]
                    x = len(z[z>0]) + 1
                    y = len(z[z<0]) + 1
                    k = min(x,y)
                    vals.append(k)
        p = np.mean(vals)/len(data_1)*2
        return p;
    # sample points from a known distribution
    def sample_pointsFromDistribution(self,n_points_I,distribution_I='lognormal',dist_params_I={},seed_I=1,raise_I = False):
        '''
        Sample points from a known distribution

        Args:
            n_points_I = int, number of points to sample
            distribution_I = function/object, distribution function
            dist_params_I = {}, of parameters for the distribution function
            seed_I = int, seed to initialize the random number generator

        Returns:
            points_O = ndarray of points of length n_points_I

        Example:
            n_points_I = 50
            distribution_I = lognormal
            distribution_params_I = {'mean':0.0,'sigma'=1.0}

        Notes:
            see http://docs.scipy.org/doc/numpy/reference/routines.random.html
            for a list of supported distributions and associated parameters
        '''
        import numpy.random as npr

        dist_params_I['size'] = n_points_I;
        points_O = None;
        try:
            if seed_I:
                npr.seed(seed_I);
            if hasattr(npr, distribution_I):
                distribution_func = getattr(npr, distribution_I);
                points_O = distribution_func(**dist_params_I);
            else: print('numpy.random does not support ' + distribution_I);
        except Exception as e:
            if raise_I: raise;
            else: print(e);
        return points_O;
    def sample_pointsFromPoints(self,n_points_I,points_I,resample_params_I={},seed_I=1,raise_I = False):
        '''
        Sample points from a distribution of points

        Args:
            n_points_I = int, number of points to sample
            points_I = the set of known points
            resample_params_I = {}, of parameters for the distribution function

        Returns:
            points_O = ndarray of points of length n_points_I

        Example:
            n_points_I = 50
            points_I = [...]
            resample_params_I = {'replace':True,'p'=None}

        Notes:
            see http://docs.scipy.org/doc/numpy/reference/generated/numpy.random.choice.html#numpy.random.choice
            for a list of associated parameters
        '''
        import numpy.random as npr

        resample_params_I['size'] = n_points_I;
        points_O = None;
        try:
            if seed_I:
                npr.seed(seed_I);
            points_O = npr.choice(points_I,**resample_params_I);
        except Exception as e:
            if raise_I: raise;
            else: print(e);
        return points_O;