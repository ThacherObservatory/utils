from numpy import sqrt,float64
from scipy.special import erfinv,erf

def conf_sig(sigma=None,conf=None,verbose=False):
    """
    Routine to calculate either the confidence interval of a 
    Gaussian distribution given the number of standard deviations
    OR
    the equivalent standard deviation given a confidence interval
    
    examples
    -------
    In [1]: conf_sig(sigma=3,verbose=True)
    Confidence interval = 0.997300, Sigma = 3.000000
    Out[1]: 0.99730020393673979

    In [1]: sigma = conf_sig(conf=0.99)
    In [2]: sigma
    Out[2]: 2.5758293035489004
    """
        
    if sigma:
        # Calculate confidence interval for a given sigma
        val = erf(float64(sigma)/sqrt(2))
        conf = val
    elif conf:
        # Calculate a sigma value for a given confidence interval
        val = erfinv(float64(conf))*sqrt(2)
        sigma = val
    else:
        # Either sigma or conf has to be input 
       return 'You must input either a sigma or a confidence interval'

    out = 'Confidence interval = %f, Sigma = %f' % (conf,sigma)

    if verbose:
        print out 

    return val
