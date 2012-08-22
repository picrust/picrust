#!/usr/bin/env python
# File created on April 10  2012
from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2012, The PICRUST project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"
from collections import defaultdict
from numpy import  array,ravel,around
from copy import copy
from cogent.maths.stats.test import correlation 
from cogent.maths.stats.distribution import tprob,t_high
from biom.table import table_factory,DenseOTUTable


def merge_by_column_union(old_table,new_table,new_field_modifier):
    """Merge tables using observations"""
    #print "Old new table:",new_table
    obs_ids,obs_meta,sample_ids,sample_meta,data =\
      new_table.ObservationIds, new_table.ObservationMetadata,new_table.SampleIds,new_table.SampleMetadata,new_table._data
    new_sample_ids = []
    for sample_id in sample_ids:
        new_sample_ids.append(sample_id +'_'+ new_field_modifier)
    #Note that sample metadata is a list of dicts, so renaming should be OK
    new_table = table_factory(data,new_sample_ids,obs_ids,sample_meta,obs_meta,constructor=DenseOTUTable)
    #print "\nnew new table:",new_table
    #print "\nold table:",old_table
    merged_table =  old_table.merge(new_table,Sample='union',Observation='intersection')
    #print "\nmerged table:",merged_table
    #raise ValueError("")
    return merged_table

def update_pooled_data(obs_table,exp_table,tags,pooled_observations,\
   pooled_expectations,new_unique_id,verbose=False):
    """Update pooled data using an observed and expected table, tags, and a dict of pooled_observations"""

    #Update the global pooled obs,exp tables with the current results
    #Format results for printing
    #NOTE: currently merged values actually add value into a single column rather than creating new columns for each  
    #print "Non-pooled fields:", "_".join(non_pooled_fields)
    #print "pooled_observations:",pooled_observations
    #print "pooled_expectations:",pooled_expectations
    for tag in tags:
        #pooled entries will either be empty or be a .BIOM table
        if pooled_observations.get(tag,None) is None:
            pooled_observations[tag] = obs_table
            pooled_expectations[tag] = exp_table
        else:    
            #Entry should already be a table.  So we want to update it by merging in 
            #the current result
            #if verbose:
            # 
            #    print "Merging observations with existing biom table for ", tag
            old_obs_table = pooled_observations[tag]
            
            pooled_observations[tag]=\
                merge_by_column_union(old_obs_table,obs_table,new_unique_id)
            
            #pooled_observations[tag] =\
            #    old_obs_table.merge(obs_table,Sample='union',Observation='union')
            
            #if verbose:
            #    print "Merging observations with existing biom table for ", tag
            
            old_exp_table = pooled_expectations[tag]
            #pooled_expectations[tag] =\
            #    old_exp_table.merge(exp_table,Sample='union',Observation='union')
            pooled_expectations[tag]=\
                merge_by_column_union(old_exp_table,exp_table,new_unique_id)
    
    return pooled_observations,pooled_expectations

def unzip(l):
    """Unzip a list into a tuple"""
    return tuple(zip(*l))

def run_accuracy_calculations_on_biom_table(obs_table,exp_table,metadata,verbose=False):
    """Run correlation and AUC measures for paired obs,exp .biom Table objects"""

    scatter_data_points, correlations=\
        evaluate_test_dataset(obs_table,exp_table)

    #For AUC, format = [(all_obs_points,all_exp_points)]
    new_trial = unzip(scatter_data_points)

    new_scatter_lines = format_scatter_data(scatter_data_points,metadata)
    #scatter_lines.extend(new_scatter_lines)

    new_correlation_lines = format_correlation_data(correlations,metadata)
    #correlation_lines.extend(new_correlation_lines)

    if verbose:
        for l in new_correlation_lines:
            print l
    return new_scatter_lines,new_correlation_lines,new_trial

def run_accuracy_calculations_on_pooled_data(pooled_observations,pooled_expectations,roc_success_criteria=['binary','exact'],verbose=False):
    """Run pearson, spearman calculations on pooled observation & expectation datai
    
    pooled_observations -- a dict of observations, with keys set to a descriptive tag and values equal to the observed .biom Table object

    pooled_expectations -- a dict of expectations, with keys set to a descriptive tag and values equal to the observed .biom Table object
    success_criterion -- criterion for success in ROC trial.  Either 'binary' or 'exact'
    verbose -- if set to True, print verbose output
    """
    all_scatter_lines = []
    all_correlation_lines = []
    trials = defaultdict(list)
    for tag in sorted(pooled_observations.keys()):
        if verbose:
            print "calculating scatter,correlations,trials for tag:",tag

        metadata= tag.split('\t')
        obs_table = pooled_observations[tag]
        exp_table = pooled_expectations[tag]
        scatter_lines,correlation_lines,new_trial =\
           run_accuracy_calculations_on_biom_table(obs_table,\
           exp_table,metadata,verbose=verbose)
        #trials["_".join(map(str,[holdout_method,prediction_method]))].append(new_trial)
        trials["_".join(map(str,[metadata[0],metadata[1]]))].append(new_trial)
        #field[0] = holdout_method, field[1] = prediction_method
        all_scatter_lines.extend(scatter_lines)
        all_correlation_lines.extend(correlation_lines)

    if verbose:
        print "Running ROC analysis..."


    # Now that we have all trials calculated, we can produce AUC results
    # (ROC stands for receiver-operating characteristics)

    #TODO: calculate the 'binary' and 'exact' ROC curves for each dataset
    all_roc_result_lines={}
    all_roc_auc_lines={}
    for success_criterion in roc_success_criteria:
        if verbose:
            print "Calculating ROC graph points, AUC for criterion:",success_criterion
        roc_result_lines, roc_auc_lines = run_and_format_roc_analysis(trials,success_criterion=success_criterion,verbose=verbose)
        if verbose:
            for l in roc_auc_lines:
                print l
        all_roc_result_lines[success_criterion]=roc_result_lines
        all_roc_auc_lines[success_criterion]=roc_auc_lines

    return all_scatter_lines, all_correlation_lines, all_roc_result_lines, all_roc_auc_lines

def evaluate_test_dataset(observed_table,expected_table):
    """ evaluate the correlation between an observed and expected
    biom table.

    Returns data points for a scatter plot of observed v. expected values,
    and a dict of correlations keyed by method (each containing the r value,
    then the probability) 
    """
    # identify the overlapping otus that can be used to predict metagenomes 
    overlapping_ids = list(set(observed_table.ObservationIds) &
                            set(expected_table.ObservationIds))
    #print "dir(observed_table):\n",dir(observed_table)
    #print overlapping_ids

    if len(overlapping_ids) < 1:
        print "obs ids:",observed_table.ObservationIds[0:10]
        print "exp ids:",expected_table.ObservationIds[0:10]
        
        raise ValueError,\
         "No ids are in common  between the observed and expected tables, so no evaluations can be performed."

    # create lists to contain filtered data - we're going to need the data in 
    # numpy arrays, so it makes sense to compute this way rather than filtering
    # the tables
    obs_data = []
    exp_data = []

    # build lists of filtered data
    for obs_id in overlapping_ids:
        obs_data.append(observed_table.observationData(obs_id))
        exp_data.append(expected_table.observationData(obs_id))

    #print obs_data
    #print exp_data
    flat_obs_data = ravel(array(obs_data))
    flat_exp_data = ravel(array(exp_data))
    #print flat_obs_data
    #print flat_exp_data

    #GET THE SCATTER PLOT POINTS
    scatter_data_points =\
      zip(flat_obs_data,flat_exp_data)
    
    correlations = {}
    if len(scatter_data_points) <= 2:
        #can't validly calc correlation
        correlations["pearson"] = (None,None)
        correlations["spearman"] = (None,None)
        return scatter_data_points,correlations    
    # CALCULATE CORRELATIONS

    
    pearson_r,pearson_t_prob =\
      correlation(flat_obs_data,flat_exp_data)
    
    correlations["pearson"] = (pearson_r,pearson_t_prob)
    pearson_r2 = pearson_r**2
    correlations["pearson_r2"] = [pearson_r2]
    spearman_r,spearman_t_prob =\
      spearman_correlation(flat_obs_data,flat_exp_data)
    
    correlations["spearman"] = (spearman_r,spearman_t_prob)
    spearman_r2 = spearman_r**2
    correlations["spearman_r2"] = [spearman_r2]
   
    return scatter_data_points,correlations     



##########################
#  Formatting functions  #
##########################
def run_and_format_roc_analysis(trials, success_criterion='binary',verbose=False):
    """Run a receiver-operating characteristics(ROC) analysis and format results
    trials:  a dict of trials, keyed by relevant metadata strings.  Each trial
      must be a list of observed,expected values  (so the length of each trial is always two)
    
    success_criterion -- method for determining if a quantitative prediciton is a 'success' or 'failure' for purposes of constructing ROC curves
    """  
    
    roc_result_lines = []
    roc_auc_lines = []
    for roc_analysis_type in trials.keys():
        roc_trial_data = trials[roc_analysis_type]
        roc_points, roc_auc, = roc_analysis(roc_trial_data,\
          success_criterion=success_criterion,verbose=verbose)
        
        new_roc_result_lines = format_scatter_data(roc_points,metadata=[roc_analysis_type])
        roc_result_lines.extend(new_roc_result_lines)

        new_roc_auc_line = '\t'.join(map(str,[roc_analysis_type] + [roc_auc]))+"\n"
        roc_auc_lines.append(new_roc_auc_line)

    return roc_result_lines, roc_auc_lines



def format_roc_data(roc_points,metadata=[],delimiter="\t"):
    """Convert evaluation data to delimited lines
    metadata -- fixed strings that will be added to their own columns
    in the order supplied (i.e. which organism/method/distance the results
    for which the results were calculated)

    """
    lines =[]
    fixed_fields = metadata
    for x,y in scatter_points:
        data_fields = map(str,[x,y])
        lines.append('\t'.join(fixed_fields + data_fields)+"\n")
    return lines



def format_scatter_data(scatter_points,metadata=[],delimiter="\t"):
    """Convert evaluation data to delimited lines
    metadata -- fixed strings that will be added to their own columns
    in the order supplied (i.e. which organism/method/distance the results
    for which the results were calculated)

    """
    lines =[]
    fixed_fields = metadata
    for x,y in scatter_points:
        data_fields = map(str,[x,y])
        lines.append('\t'.join(fixed_fields + data_fields)+"\n")
    return lines


def format_correlation_data(correlations,metadata=[],delimiter="\t"):
    """Convert a dict of correlation data to formatted lines
    correlations -- a dict of correlation data, keyed by method 
    metadata-- a list of strings representing metadata for the correlation
    (e.g. which organism, what prediction method etc)
    delimiter -- delimiter to separate output.
    
    """
    correlation_lines = []
    for corr_type in sorted(correlations.keys()):

        new_correlation_fields =\
          metadata + [corr_type] + map(str,[d for d in correlations[corr_type]])

        new_correlation_line = "\t".join(new_correlation_fields)+"\n"
        correlation_lines.append(new_correlation_line)

    return correlation_lines




#########################
# Evaluation functions  #
#########################

def convert_vals_to_spearman_ranks(vals):
    """Return the rank of each val in list for Spearman rank correlation
    vals -- a list of floats

    Ranks will be returned in the same order as the entered values.
    Calculating rank:  if there are no ties, rank is just position in the ordered
    list.  If there are ties, then the rank is equal to the mean of all tied ranks
    """
    
    #In order to find ranks, we need to 
    #sort values from largest to smallest
    ordered_vals = sorted(vals,reverse=True)
    

    #Next we need to generate a mapping between 
    #each value and its rank using the sorted
    #list, so that we can apply it to the original 
    #list
    
    val_to_rank = {}
    for i,val1 in enumerate(ordered_vals):
        
        #Get indices (numbering from 1 not 0) of all occurences of value
        ranks= [j+1 for j,val2 in enumerate(ordered_vals) if val2 == val1]
        
        #The overall rank is always the mean of the list of ranks
        rank = float(sum([r for r in ranks]))/float(len(ranks))
        
        val_to_rank[val1] = rank

    result = []
    for val in vals:
        result.append(val_to_rank[val])

    return result

def calc_spearman_t(r,n,tails='two-tailed',eps=1e-10):
    """Calculate the t value, and probability for Spearman correlation
    
    Equation from Wikipedia article on Spearman rank correlation:
    http://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient
    
    TODO:  verify with additional / published sources.
    TODO:  want to find the correct handling for the 
    case of perfect correlation
    """
    #Right now I subtract an arbitrarily tiny
    #epsilon if r is exactly 1.0 to prevent divide-by-zero
    #errors.  

    if r == 1.0:
        r = r-eps
    
    t = r*(((n-2)/(1.0-r**2))**0.5)
    
    if tails == 'two-tailed':
        prob = tprob(t,n-2)
    elif tails == 'high':
        prob = t_high(t,n-2)
    elif tails == 'low':
        prob = t_low(t,n-2)
    else:
        raise RuntimeError("Valid prob. methods are 'two-tailed','high',and 'low'")

    return prob

        
def spearman_correlation(x_array,y_array,tails = 'two-tailed'):
    """calculate the Spearman rank correlation for x and y
    
    x_array -- a 1D NumPy array
    y_array -- a 1D NumPy array
    """
    
    #Convert absolute values to ranks
    x_ranks = convert_vals_to_spearman_ranks(x_array)
    y_ranks = convert_vals_to_spearman_ranks(y_array)
    
    #Now we get r by performing Pearson correlation
    #on the rank data.
    r, pearson_prob = correlation(x_ranks, y_ranks)
    
    #However, the conversion to ranks affects the prob
    #so we need the corrected version of the t statistic
    #not the generic version used by Pearson correlation
    spearman_t_prob =\
      calc_spearman_t(r,n=len(x_array),tails=tails)

    
    #return r,spearman_t_prob
    return r,spearman_t_prob


def calculate_accuracy_stats_from_observations(obs,exp,success_criterion='binary',\
  verbose=False):
    """Return statistics derived from the confusion matrix
    obs -- a list of floats representing observed values
    exp -- a list of floats for expected values (same order as obs)
    verbose -- print verbose output 
    """
    
    tp,fp,fn,tn =\
      confusion_matrix_from_data(obs,exp,success_criterion=success_criterion,verbose=verbose)
    
    result =\
      calculate_accuracy_stats_from_confusion_matrix(tp,fp,fn,tn,verbose=verbose)
    return result
    
def calculate_accuracy_stats_from_confusion_matrix(tp,fp,fn,tn,allow_zero_results=True,verbose=False):
    """Calculate accuracy,sensitivity,specificity,etc from the number of true positives,false positives,false negatives, and true negatives
    
    tp -- number of true positives
    fp -- number of false positives
    fn -- number of false negatives
    tn -- number of true negatives
    verbose -- print verbose output

    allow_zero_results --
    """
    if allow_zero_results:
        #Add a tiny pseudocount to all categories
        small_num = 1e-10
        tp += small_num
        tn += small_num
        fn += small_num
        fp += small_num
    
    n = sum([tp,fp,tn,fn])
    results ={} 
    results['positive_predictive_value'] = tp/(tp+fp) # == precision
    results['sensitivity'] =  tp/(tp+fn) # == TPR, hit rate, recall
    results['specificity'] =  tn/(tn+fp) # == 1- FPR, TNR
    results['false_positive_rate'] = fp/(fp+tn)
    results['accuracy'] = (tp + tn)/n  

    return results
    
def confusion_matrix_from_data(obs,exp,success_criterion='binary',verbose=False):
    """Return the number of true_positive,false_positive,false_negatives,false_positives from paired lists of observed and expected values


    obs -- a list of floats representing observed values
    exp -- a list of floats for expected values (same order as obs)
    verbose -- print verbose output 
    """
    tp_indices,fp_indices,fn_indices,tn_indices =\
      confusion_matrix_results_by_index(obs,exp,\
        success_criterion=success_criterion,verbose=verbose)
    #print obs,exp
    #print tp_indices,fp_indices,fn_indices,tn_indices 
    tp = float(len(tp_indices))
    fp = float(len(fp_indices))
    fn = float(len(fn_indices))
    tn = float(len(tn_indices)) 
    #print tp,fp,fn,tn 
    return tp,fp,fn,tn

def confusion_matrix_results_by_index(obs,exp,success_criterion='binary',verbose=False):
    """Return indices in obs that are true positives, false positives, false negatives or true negatives

    
    obs -- a list of floats representing observed values
    exp -- a list of floats for expected values (same order as obs)
    verbose -- print verbose output 
    success_criterion -- choose either 'binary' or 'exact'.  If binary,
      all non-zero values are considered equal. If 'exact', the exp value is compared
      to the observation.  If 'int_exact', floats are rounded to integers, then 'exact' comparison is performed.
    All values >= 1 are considered presence and scored
    identically.  0s are scored as absence.
    """
    
    if success_criterion == 'int_exact':
        obs = around(obs)
        exp = around(exp)
        #Like exact, but first round numbers such that e.g. 0.7 == 1.0
        success_criterion = 'exact'

    if success_criterion == 'exact':
        
        tp_indices =\
            [i for i,f in enumerate(obs) if exp[i] == obs[i] and exp[i] != 0]
        fp_indices =\
            [i for i,f in enumerate(obs) if obs[i] > exp[i]]
        fn_indices =\
            [i for i,f in enumerate(obs) if obs[i] < exp[i]]
        tn_indices =\
            [i for i,f in enumerate(obs) if exp[i] == obs[i] and exp[i] == 0]
    
    elif success_criterion == 'binary':
        #Binary method (all non-zero values equal success 
        tp_indices =\
            [i for i,f in enumerate(obs) if exp[i] >= 1 and obs[i]>=1]
    
        fp_indices =\
            [i for i,f in enumerate(obs) if exp[i] == 0 and obs[i]>=1]

        fn_indices =\
            [i for i,f in enumerate(obs) if exp[i] >= 1 and obs[i]==0]

        
        tn_indices =\
            [i for i,f in enumerate(obs) if exp[i] == 0 and obs[i]==0]
    
    else: 
        raise ValueError('Specify either "binary" or "exact" as success criterion')
    return tp_indices,fp_indices,fn_indices,tn_indices

def roc_analysis(trials,success_criterion='binary',verbose= False):
    """Perform ROC analysis for a set of trials, each a list of obs,exp values"""
    points = roc_points(trials,success_criterion=success_criterion,verbose=verbose)
    #These points are now (FPR,TPR) tuples 
    #print "ROC AUC points (FPR,TPR):", sorted(points)
    points = average_y_values(points) 
    area_under_the_curve = roc_auc(points)
    return points,area_under_the_curve


def roc_points(trials,success_criterion='binary',verbose=False):
    """get points for a roc curve from lists of paired observed and expected values
    trials -- a list of (obs,exp) tuples, where obs and exp are lists of observed vs. expected data values
    
    Returns a list of points on the ROC curve in x,y format, with the False Positive Rate (FPR, 1-specificity) on the x axis, and the True Positive Rate (sensitivity) on the y axis."""
    points = []
    #print trials
    for obs,exp in trials:
        curr_result = calculate_accuracy_stats_from_observations(obs,exp,success_criterion=success_criterion,verbose=verbose)
        x = curr_result["false_positive_rate"]
        y = curr_result["sensitivity"]
        points.append((x,y))
    #print points
    return points

def roc_auc(points, add_endpoints = True):
    """Get the ROC Area Under the Curve given a list of FPR,TPR values for trials
    points -- a list of (FPR, TPR) tuples.  That is, a list of tuples each providing a false_positive_rate, sensitivity.
    
    add_endpoints -- if True, add values for 0,0 and 1.0,1.0 if *the corresponding x values* are not present.  This allows natural calculation of the AUC, even if there is only one data point present.
    
    This is the average prob. of ranking a random positive above a random negative.
    """
    predict_none = (0.0,0.0)
    predict_all = (1.0,1.0)
    if add_endpoints:
        if min([p[0] for p in points]) > 0.0: 
            #print "Adding %s" % str(predict_none)
            points.append(predict_none)
        if max([p[0] for p in points]) < 1.0: 
            #print "Adding %s" % str(predict_all)
            points.append(predict_all)
   
    #Average identical x values with
    #different y values to produce a
    #single curve (when multiple or 
    #stochastic measurements are available)
    points = average_y_values(points) 

   
    #gini_coefficient will order points for us
    #G = gini_coefficient(points)
    #print "G:",G
    #area_under_the_curve = (G+1.0)/2.0 
    
    
    area_under_the_curve =\
      trapezoidal_approx_to_auc(points)
    #print "AUC:",area_under_the_curve
    
    return area_under_the_curve

def average_y_values(points):
    """Given a list of points, return a list averaging all y values with identical x values
    
    points -- a list of (x,y) tuples
    
    Example:  [1,3],[1,5] --> (1,4)
    """

    tprs = defaultdict(list)
    for fpr, tpr in points:
        tprs[fpr].append(tpr)

    for y in tprs.keys():
        tprs[y] = sum(tprs[y])/len(tprs[y])
     
    result = zip(tprs.keys(),tprs.values())
    return result


def trapezoidal_area(x1,y1,x2,y2):
    """Calculate the area of a trapezoid given x,y coordinates
    
    
    """
    
    #formula = (a + b)/2 *h
    height = max(x2,x1)-min(x2,x1)
    a = y1
    b = y2
    average_side = (a+b)/2.0
    area =  average_side * height
    return area

def trapezoidal_approx_to_auc(points):
    """Return the trapezoidal approximation to the AUC, given a list of (TPR,FPR) points
    
    The area under the curve is calculated using a trapezoidal approximation.
    Effectively, this results in linear interpolation between the different
    True Positive Rate and False Positive Rate values.
    """
    
    points = sorted(points)
    indices_to_consider = [(i,i-1) for i in range(1,len(points))]
    X,Y = zip(*points)
    area = sum([float(trapezoidal_area(X[k],Y[k],X[k_minus_one],\
           Y[k_minus_one])) for k,k_minus_one in indices_to_consider])
    return area
      


def gini_coefficient(points):
    """Calculate the Gini coefficient using trapezoidal approximations
    
    trials -- a list of (TPR,FPR) tuples, where TPR  and FPR are lists of true positive rate vs. false_positive_rate values

    """
    points = sorted(points)
    indices_to_consider = [(i,i-1) for i in range(1,len(points))]
    #print indices_to_consider
    #print points
    X,Y = zip(*points)
    print X,Y
    gini = 1.0 - \
        sum([float((X[k]-X[k_minus_one])*(Y[k]+Y[k_minus_one])) for k,k_minus_one in indices_to_consider])
                
    return gini
