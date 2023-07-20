import scipy.stats

def normal_Prior(center,width,random):
    """Generate a normal prior distribution for a given parameter,
    
    Args:
        center (float): center value of this gaussian distribution.
        width (float): width of this gaussian distribution, this is the 1-sigma width.
        random (float): random number generated to do inference, this is follow the 
        definition of baysian workflow of UltraNest, here default to be cube[i]
        
    Returns:
        ppf (float): ppf of this distribution function
        
    """
    return scipy.stats.norm(center, width).ppf(random)

def flat_prior(low, up,random):
    """Generate a flat prior distribution for a given parameter,
    
    Args:
        low (float): lower bound of this flat distribution.
        up (float): upper bound of this flat distribution.
        random (float): random number generated to do inference, this is follow the 
        definition of baysian workflow of UltraNest, here default to be cube[i]
        
    Returns:
        ppf (float): ppf of this distribution function
        
    """
    return low + (up - low) * random



