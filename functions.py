from CB_02704 import *
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy as sc

# NormalizeGenotypeArray function takes in a genotype_array array and normalizes it (both mean and std dev)
def StandardizeGenotypeArray(genotype_array):
    p = genotype_array.mean(axis=1, keepdims = True)/2
    # print(p)
    genotype_array_standardized = (genotype_array-2*p)/np.sqrt(2*p*(1-p))
    genotype_array_standardized.set_fill_value(0)
    genotype_array_standardized = genotype_array_standardized.filled()
    return genotype_array_standardized

# FilterOutLowVarianceFeatures function takes in a genotype array (np masked array), and filters out the
# the features which have high variance and only keeps top_x_features 
def FilterLowVarianceFeatures(genotype_array, top_x_features):
    p = genotype_array.mean(axis=1).data/2
    # get the variance values assuming hardy-weinberg principle
    variance_values = 2*p*(1-p)
    # sort the variance values
    variance_indices_sorted = np.argsort(variance_values)[::-1] # the sorted values' indices in descending order

    # now just keep top_x_features number of indices
    genotype_array_filtered = genotype_array[variance_indices_sorted[:top_x_features], :]

    return genotype_array_filtered 

# PCA function takes in the genotype data (numpy masked array with snps as rownames and 
# individuals as colnames) and computes the PCA for that data
def PCA(genotype_array, top_x_features): 
    # first, subset the genotype_array to top_x_features number of features, so that we can do pca fast
    # the top features are chosen according to their variance
    genotype_array_filtered = FilterLowVarianceFeatures(genotype_array, top_x_features)

    # first mean centre (the features, which are snps in this case should have mean zero)
    # will use the StandardizeGenotype function from the assignments, hence also correcting variance
    genotype_array_standardized = StandardizeGenotypeArray(genotype_array_filtered)

    # then compute covariance matrix (it should have dimensions of features x features)
    cov = (1/genotype_array.shape[0])*(genotype_array_standardized @ genotype_array_standardized.T)
    
    # then compute eigenvalues and eigenvectors of the matrix
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    # eigenvalues, eigenvectors = sc.linalg.eigh(cov)


    # then perform projection of the datapoints in the direction of the two top eigenvectors
    # get the coordinate of each individual in the direction of the top eigenvector
    x_values = eigenvectors[:,-1] @ genotype_array_standardized
    # get the coordinate of each individual in the direction of the second top eigenvector
    y_values = eigenvectors[:,-2] @ genotype_array_standardized

    return x_values, y_values, eigenvalues, cov

# PlotVarianceOfFeatures function takes in a genotype array (np masked array), and plots the variance of
# the different snps. This is just to get an idea of how many snps may be less important while doing pca
def PlotVarianceOfFeatures(genotype_array, title):
    p = genotype_array.mean(axis=1).data/2
    # get the variance values assuming hardy-weinberg principle
    variance_values = 2*p*(1-p)
    # sort the variance values
    # variance_values_sorted = np.sort(variance_values)[::-1] # the sorted values in descending order
    # sns.lineplot(variance_values_sorted)

    sns.histplot(variance_values)
    plt.title(title)
    plt.xlabel('Variance of SNP')
    plt.show()
    return

# PlotPCAResults function takes in the arrays for xvalues and yvalues of a set of points, and plots them
def PlotPCA(x_values, y_values, population_list, title):
    df = pd.DataFrame()
    df['PC1'] = x_values
    df['PC2'] = y_values
    df['population'] = population_list
    sns.scatterplot(df, x = 'PC1', y = 'PC2', hue='population')
    plt.title(title)
    plt.show()
    return

def CalculateCorrelationBetweenSNPs(snp1_index, snp2_index, pop_gen, masked=True):
    
    # first find the indices where both snps are not none
    if masked == True:
        ppl_indices_with_data_from_both_snps = ~((pop_gen[snp1_index].mask) | (pop_gen[snp2_index].mask))
    else:
        ppl_indices_with_data_from_both_snps = [True]*pop_gen.shape[1]

    # copies of (variant) allele at any index is basically the row in pop_gen at that index (but needs to 
    # be filtered to have only those individuals which are present in both the alleles)
    ga = pop_gen[snp1_index][ppl_indices_with_data_from_both_snps] # ga is of type masked numpy array, and has dimension 1.
    gb = pop_gen[snp2_index][ppl_indices_with_data_from_both_snps]

    # estimate the covariance between ga and gb
    cov = CovarianceBetweenArrays(ga, gb)

    # estimate the variance of ga
    var_ga = VarianceOfArray(ga)
    var_gb = VarianceOfArray(gb)

    r2 = (cov*cov)/(var_ga*var_gb)

    return r2

# CovarianceBetweenArrays function takes as input two numpy arrays and return the covariance between them
# taking N as the normalization number instead of N-1
def CovarianceBetweenArrays (arr1, arr2):

    assert arr1.shape[0] == arr2.shape[0], 'Arrays must be of same length'

    arr1_mean = arr1.mean()
    arr2_mean = arr2.mean()

    # mean corrected arrays
    arr1_mean_corrected = arr1 - arr1_mean
    arr2_mean_corrected = arr2 - arr2_mean

    numerator = np.dot(arr1_mean_corrected, arr2_mean_corrected)
    denominator = arr1.shape[0]

    return numerator/denominator

# VarianceOfArray function takes as input a 1d np array and returns its variance
def VarianceOfArray(arr):
    
    arr_mean = arr.mean()

    # mean corrected array
    arr_mean_corrected = arr - arr_mean

    numerator = np.sum(arr_mean_corrected*arr_mean_corrected)
    denominator = arr.shape[0]

    return numerator/denominator

# AverageR2BetweenAdjacentSNPs function takes in a pop_gen array, and computes the correlation
# between every pair of adjacent snps and returns the average value. 
def AverageR2BetweenAdjacentSNPs(pop_gen, masked = True): # masked tells if the input array is masked type or not

    # num snps variable will help in the number of iterations in the loop below
    num_snps = pop_gen.shape[0]

    # make r2 list which will store the r2 values (num_snps-1 in count)
    r2_list = []

    # make list variables which will store the ga and gb value for each 
    for snp1_index in range(num_snps-1):
        snp2_index = snp1_index + 1

        r2 = CalculateCorrelationBetweenSNPs(snp1_index, snp2_index, pop_gen, masked)

        r2_list.append(r2) 

        # print the index numbers to keep track of progress
        # if snp1_index%10000 == 0: print(f'snp_index {snp1_index} out of {num_snps} reached.')

    r2_avg = np.array(r2_list).mean()
    return r2_avg

# DataframeFromDict function to make the plot of the data using dicts mapping values to list (plots variance too)
def DataframeFromDict(input_dict, label): # label is the identifier for which population we are processing
    # first, we need to convert the dict into a df in long form
    correlation_values = []
    key_values = []
    for key in input_dict.keys():
        # create a list of same value having keys and append the values for the key to the ongoing list 
        # of correlation values
        correlation_values = correlation_values + input_dict[key]
        key_values = key_values + [key]*len(input_dict[key])
    
    # now we have the two cols ready for the df
    df = pd.DataFrame()
    df['distance'] = key_values
    df['correlation_values'] = correlation_values
    df['label'] = label

    return df

# AvgCorrelationBetweenPairOfSetOfSNPs function takes in two lists of snp indices and computes the correlation
# between every pair between the two lists and returns the average of that
def AvgCorrelationBetweenPairOfSetOfSNPs(indices_list1, indices_list2, pop_gen):
    correlation_sum = 0
    # loop through each snp in indices_list1 and indices_list2 and compute correlation between them
    for snp1_index in indices_list1:
        for snp2_index in indices_list2:
            correlation_sum += CalculateCorrelationBetweenSNPs(snp1_index, snp2_index, pop_gen)
    # now take average
    avg_correlation = correlation_sum/(len(indices_list1)*len(indices_list2))
    return avg_correlation