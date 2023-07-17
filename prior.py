import scipy.stats

def normal_Prior(center,width,random):
    return scipy.stats.norm(center, width).ppf(random)

def flat_prior(low, up,random):
    return low + (up - low) * random



