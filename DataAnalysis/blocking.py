import numpy as np

def blocking(sample):
    ''' perform blocking transformation on an array of sampled energies.
    In case there is an odd number of energies, the last one is discarded.'''

    x1 = sample[0:-1:2]
    x2 = sample[1::2]
    return 0.5*(x1+x2)

def block_variance(sample, talk=False):
    '''estimate the variance of the mean, using the blocking method.'''
    var_array = np.empty(0)
    nk=10
    sample_mean = np.mean(sample)
    while(nk>2):
        blocked_sample = blocking(sample)
        nk = blocked_sample.shape[0]
        if talk:
            print(f"Blocking a sample of size {nk}")

        samplevar = np.var(blocked_sample, axis=0)
        old_est = var_est
        var_est = samplevar/(nk-1)
        var_array = np.append(var_array, var_est, axis=0)
        sample = blocked_sample
        #exit if change is less than the estimated uncertainty (aka standard deviation)
        uncertainty = var_est*np.sqrt(2/(nk-1))
        if np.abs(old_est-var_est)<uncertainty:
            break

    return var_est, sample_mean



