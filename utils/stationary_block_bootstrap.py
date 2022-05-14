#Copyright (c) 2011, Josh Hemann  (hemann @ colorado . edu)
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in the
#      documentation and/or other materials provided with the distribution.
#    * Neither the name of the code's author, Josh Hemann, nor the
#      names of its contributors, may be used to endorse or promote products
#      derived from this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
#DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
#ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Modified to support for Python 3.8 (by Taylor Bolt)

import numpy as np
from collections import defaultdict




def resample(_data, p, seed=None):
    """
    Performs a stationary block bootstrap resampling of elements from a time 
    series, or rows of an input matrix representing a multivariate time series
	
    Inputs:
        data - An MxN numerical array of data to be resampled. It is assumed
               that rows correspond to measurement observations and columns to 
	       items being measured. 1D arrays are allowed.
        p    - A scalar float in (0,1] defining the parameter for the geometric
               distribution used to draw block lengths. The expected block 
	       length is then 1/p				

    Keywords:
        seed - A scalar integer defining the seed to use for the underlying
	   random number generator.
	
    Return:
        A three element list containing
	- A resampled version, or "replicate" of data
	- A length M array of integers serving as indices into the rows
	  of data that were resampled, where M is the number of rows in 
	  data. Thus, if indices[i] = k then row i of the replicate data 
	  contains row k of the original data.
	- A dictionary containing a mapping between indices into data and
	  indices into the replicate data. Specifically, the keys are the
	  indices for unique numbers in data and the associated dict values 
	  are the indices into the replicate data where the numbers are.

    Example:			
        In [1]: import numpy as np
        In [2]: x = np.random.randint(0,20, 10)
        In [3]: import stationary_block_bootstrap as sbb
        In [4]: x_star, x_indices, x_indice_dict = sbb.resample(x, 0.333333)
        In [5]: x
        Out[5]: array([19,  2,  9,  9,  9, 10,  2,  2,  0, 11])
        In [6]: x_star
        Out[6]: array([19, 11,  2,  0, 11, 19,  2,  2, 19,  2])
        In [7]: x_indices
        Out[7]: array([0, 9, 7, 8, 9, 0, 6, 7, 0, 1])

        So, x_indices[1] = 9 means that the 1th element of x_star corresponds 
        to the 9th element of x
		
        In [8]: x_indice_dict
        Out[8]: {0: [0, 5, 8], 7: [2, 6, 7, 9], 8: [3], 9: [1, 4]}
	
        So, x[0] = 19 occurs in position 0, 5, and 8 in x_star. Likewise, 
        x[9] = 11 occurs in positions 1 and 4 in x_star
	
"""
    num_obs = np.shape(_data)[0]
    num_dims = np.ndim(_data)
    assert num_dims == 1 or num_dims == 2, "Input data must be a 1 or 2D array"
    #There is a probably smarter way to wrap the series without doubling
    #the data in memory; the approach below is easy but wasteful
    if num_dims == 1:
        wrapped_data = np.concatenate((_data, _data)) 
    elif num_dims == 2:
        wrapped_data = np.row_stack((_data, _data)) 
	
    assert p > 0 and p <=1, "p must be in (0,1]"
	
    if seed is not None:
        np.random.seed(seed=seed)

    #Make the random variables used in the resampling ahead of time. Could be
    #problematic for memory if num_obs is huge, but doing it this way cuts down 
    #on the function calls in the for loop below...
    choices = np.random.randint(0, num_obs, num_obs)
    unifs = np.random.uniform(0, 1, num_obs)
	
    #Let x and x* be length-n vectors with x*[0] = x[0]. 
    #Then x*[1] = x[1] with probability 1-p. With probability p, x*[1] will
    #equal a random i for x[i]. The expected run length is 1/p = "block length"
    indices = -np.ones(num_obs, dtype=int)
    indices[0] = 0
		
    for i in range(1, num_obs):
        if (unifs[i] > p): 
            indices[i] = indices[i-1] + 1 
        else:
            indices[i] = choices[i]

    if num_dims == 1:		
        resampled_data = wrapped_data[indices]   
        index_to_data_map = dict((x, i) for i, x in enumerate(wrapped_data))
        bootstrap_indices = map(index_to_data_map.get, resampled_data)
    elif num_dims == 2:
        #Mapping is the same for each column with respect to which rows are
        #resampled, so just consider one variable when mapping indices to data...
        resampled_data = wrapped_data[indices, :]   
        index_to_data_map = dict((x, i) for i, x in enumerate(wrapped_data[:,0]))
        bootstrap_indices = map(index_to_data_map.get, resampled_data[:,0])
		
    bootstrap_indices = [index % num_obs for index in bootstrap_indices]
	
    #The code below is making a dictionary mapping of observations resampled to
    #where in the array of indices that observation shows up. Some observations
    #are resampled multiple times, others not at all, in any given replicate
    #data set. The straight-forward code is
    # try:
    #   items = dict[key]
    # except KeyError:
    #   dict[key] = items = [ ]
    #   items.append(value)
    
    index_occurences = defaultdict(list)
    for pos, index in enumerate(bootstrap_indices):
        index_occurences[index].append(pos)
	
    index_dict = dict(index_occurences)

    #Need to make the indices we save be bounded by 0:num_obs. For example, 
    #data[0,:] = data[num_obs,:]  and data[1,:] = data[num_obs+1,*] etc     
    #because we wrapped the data. But, with respect to the data arrays used 
    #elsewhere, an index of num_obs+1 is out of bounds, so num_obs should be 
    #converted to 0, num_obs+1 to 1, etc...   

    return [resampled_data, indices % num_obs, index_dict] 
    #end resample() 
	
	
if __name__ == "__main__":
    x = np.random.randint(0,20, 10)
    x_star, x_indices, x_indice_dict = resample(x, 0.333333) #expected block length = 3
    print('Original series: ', x)
    print('Resampled series: ', x_star)
    print('Indices into original series: ', x_indices)
    print('Dictionary where key=index into original data, value=index into\n' \
          + '  resampled data where that value occurs: ''', x_indice_dict)
    print('\n\n')

    y = np.arange(18).reshape((6,3))
    y_star, y_indices, y_indice_dict = resample(y, 0.5) #expected block length = 2
    print('Original MxN series of M observations over N variables:\n', y)
    print('Resampled series:\n', y_star)
    print('Indices into rows of original series: ', y_indices)
    print('Dictionary where key=row in original data, value=rows in\n' \
          + '  resampled data where that observation occurs: ', y_indice_dict)

	
