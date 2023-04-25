import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln
from scipy.integrate import quad

def hist_fit2(data_, bins, fit_function, theta_0):

    N = np.sum(data_)

    def binned_fit_function(i, theta):
        return quad(fit_function(*theta), bins[i], bins[i+1], epsabs=1e-3, epsrel=1e-3)[0]

    def get_log_L(theta):
        log_L = 0
        for i in range(len(data_)):
            if binned_fit_function(i, theta) == 0.0:
                log_p_i = -300
            else:
                log_p_i = np.log(binned_fit_function(i, theta))
            log_L += data_[i] * log_p_i - gammaln(data_[i]+1) 
        log_L += N*np.log(N)+N
        return -log_L
    
    reps = 1
    theta_min=theta_0
    for _ in range(reps):
        res = minimize(get_log_L, x0=theta_min, method='SLSQP', bounds=[(0,10),(0,20),(0,20)])
        print(res.x)
        theta_min = res.x
    return theta_min


def hist_fit(data, bins_range, n_bins, fit_function, theta_0):

    data_, bins = np.histogram(data, bins=np.linspace(bins_range[0], bins_range[1], n_bins))

    N = np.sum(data_)

    def binned_fit_function(i, theta):
        return quad(fit_function(theta), bins[i], bins[i+1])[0]

    def get_log_L(theta):
        log_L = 0
        for i in range(len(data_)):
            if binned_fit_function(i, theta) == 0.0:
                log_p_i = -300
            else:
                log_p_i = np.log(binned_fit_function(i, theta))
            log_L += data_[i] * log_p_i - gammaln(data_[i]+1) 
        log_L += N*np.log(N)+N
        return -log_L

    theta_min = minimize(get_log_L, x0=theta_0)
    return theta_min


