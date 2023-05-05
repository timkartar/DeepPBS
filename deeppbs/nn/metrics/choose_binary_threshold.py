# third party modules
import numpy as np

def chooseBinaryThreshold(y_gt, probs, metric_fn, 
        score='metric_value',
        criteria='max',
        beta=2,
        n_samples=25,
        minimize_threshold=False,
        **kwargs
    ):
    """ Choose a threshold value which meets the following criteria:
        y_gt: ...
        probs: ...
        score (string): determine what we are going to evaluate
            metric_value - the metric its self
            F-beta - the F-beta score of the metric and threshold with beta weighting the metric.
                     This is useful if we want to choose higher or lower thresholds while still 
                     preferring a good metric score.
        criteria (string):
            min - minimizes the score
            max - maximizes the score
        beta (float): the beta value to use when score=F-beta
    """
    # sample thresholds
    thresholds = np.linspace(0, 1, n_samples+2)[1:-1] # skip 0 and 1 values
    
    # choose what we are actually evaluating
    m = lambda t: metric_fn(y_gt, probs >= t, **kwargs)
    values = np.array(list(map(m, thresholds)))
    if score == 'f-beta':
        if minimize_threshold:
            t = 1 - thresholds
        values = (1+beta**2)*(t*values)/((beta**2)*t + values)
    
    # choose how to evaluate
    if criteria == 'max':
        idx = np.argmax(values)
    elif criteria == 'min':
        idx = np.argmin(values)
    
    return thresholds[idx], values[idx]
