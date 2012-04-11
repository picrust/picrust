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

from numpy import  array,ravel
from copy import copy
from cogent.maths.stats.test import correlation 
from cogent.maths.stats.distribution import tprob,t_high
from biom.table import table_factory




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

    scatter_data_points =\
      zip(flat_obs_data,flat_exp_data)
    
    correlations = {}

    
    pearson_r,pearson_t_prob =\
      correlation(flat_obs_data,flat_exp_data)
    
    correlations["pearson"] = (pearson_r,pearson_t_prob)
    
    spearman_r,spearman_t_prob =\
      spearman_correlation(flat_obs_data,flat_exp_data)
    
    correlations["spearman"] = (spearman_r,spearman_t_prob)
    return scatter_data_points,correlations     


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
    print r,n,tails
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

