from cmath import inf
import numpy as np
from scipy import stats
from scipy.special import gamma, binom

# six commonly used moment statistics for inference
# this is observed moments


def statisData(data):
    data_mean = np.mean(data, axis=0)
    data_var = np.var(data, axis=0)
    data_cv2 = data_var / (data_mean ** 2)
    data_fano = data_var / data_mean
    data_sk = stats.skew(data) + 1  # no idea why 1 is added ??
    data_kt = stats.kurtosis(data, fisher=False)
    data_bc = ((data_sk - 1) ** 2 + 1) / data_kt

    return np.array([data_mean, data_cv2, data_fano, data_sk, data_kt, data_bc])


"""
    Inputs:
    param: A structure contains model parameters (possibly dict or object)
    param.kon and param.ron: OFF dwell time distribution f_off(t) = ron^(kon) * t^(kon - 1) * e^( - ron * t) / gamma(kon)
    param.koff and param.roff: ON dwell time distribution f_on(t) = roff^(koff) * t^(koff - 1) * e^( - roff * t) / gamma(koff)
    param.mu: Transcriptional rate
    param.delta: Degradation rate.

    k ?? do not know what this is yet

    Output: 
    [mean,noise,sk,kt,bc,ff]
"""
def statisGTM(param, k):

    # parameters setting
    kon = param['kon']
    ron = param['ron']
    koff = param['koff']
    roff = param['roff']
    mu = param['mu']
    delta = param['delta']

    Laplace_s = np.array(range(0, k + 1))  # similar to Laplace_s = 0:k; ??
    Lfon = (ron / (Laplace_s + ron)) ** kon  # Laplace fon(x)
    Lfoff = (roff / (Laplace_s + roff)) ** koff  # Laplace foff(x)

    mean_tau_off = kon / ron
    mean_tau_on = koff / roff

    c = 1 / (mean_tau_off + mean_tau_on)

    # Lfon[1:] ~ Lfon(2:end)
    # LFon = [mean_tau_off (1 - Lfon(2:end)). / Laplace_s(2:end)];
    LFon = np.insert((1 - Lfon[1:]) / Laplace_s[1:], 0,
                     mean_tau_off)  # Laplace Fon(x)
    # LFoff = [mean_tau_on (1 - Lfoff(2:end)). / Laplace_s(2:end)]; % Laplace Foff(x)
    LFoff = np.insert((1 - Lfoff[1:]) / Laplace_s[1:], 0,
                      mean_tau_on)  # Laplace Foff(x)

    Ck = np.zeros(k)
    Ck[0] = (1 - Lfoff[1]) / (mean_tau_off + mean_tau_on)
    for iter in range(2, k + 1):
        i = np.array(range(1, iter + 1))
        Ck_coef = np.array(
            (mu ** i[0:-1])
            * Lfon[iter - i[0:-1]]
            * Ck[iter - 1 - i[0:-1]]
            / (((1 - Lfoff[iter - i[0:-1]] * Lfon[iter - i[0:-1]]) * gamma(i[0:-1] + 1)) ))
        Ck_coef_last_element = c * mu ** (iter - 1) / gamma(iter + 1)
        Ck_coef = np.append(Ck_coef, Ck_coef_last_element)

        Ck_sum = np.zeros(iter)
        for iter_i in range(1, iter + 1):
            for iter_j in range(0, iter_i + 1):
                Ck_sum[iter_i - 1] = Ck_sum[iter_i - 1] + \
                    binom(iter_i, iter_j) * \
                    (- 1) ** (iter_i - iter_j) * Lfoff[iter - iter_j]
        Ck[iter - 1] = np.sum(Ck_coef * Ck_sum)

    bk = np.zeros(k)
    for iter in range(1, k + 1):
        i = np.array(range(1, iter + 1))
        bk_coef_first_element = c * mu ** (iter) / gamma(iter + 1)
        bk_coef = np.array(
            (mu / iter)
            * (mu ** (iter - i[0:-1]))
            * Lfon[i[0:-1]]
            * Ck[i[0:-1] - 1]
            / (gamma(iter - i[0:-1]) * (1 - Lfon[i[0:-1]] * Lfoff[i[0:-1]]) )
        )
        bk_coef = np.insert(bk_coef, 0, bk_coef_first_element)
        bk_sum = np.zeros(iter)
        for iter_i in range(0, iter):
            for iter_j in range(0, iter - iter_i):
                bk_sum[iter_i] = bk_sum[iter_i] + binom(iter - iter_i - 1, iter_j) * (- 1) ** (
                    iter - iter_i - iter_j - 1) * LFoff[iter - iter_j - 1]
        temp = np.sum(bk_coef * bk_sum)
        
        # if temp==0. or temp == np.inf or temp == np.nan:
        #     temp = np.finfo(float).eps

        bk[iter - 1] = temp

    b = bk


    m = b[0]
    v = 2 * b[1] + b[0] - b[0] ** 2
    cv2 = v / (m ** 2)
    fano = v / m
    # print(bk)
    sk = (6 * b[2] + 6 * b[1] + b[0] - 3 * b[0] * (2 * b[1] + b[0]) + 2 * b[0] ** 3) / ((2 * b[1] + b[0] - b[0] ** 2) ** 1.5 ) + 1
    kt = (24 * b[3] + 36 * b[2] + 14 * b[1] + b[0] - 4 * b[0] * (6 * b[2] + 6 * b[1] + b[0])
          + 6 * b[0] ** 2 * (2 * b[1] + b[0]) - 3 * b[0] ** 4) / ((2 * b[1] + b[0] - b[0] ** 2) ** 2 )
    bc = ((sk - 1) ** 2 + 1) / (kt )
    statis = np.array([m, cv2, fano, sk, kt, bc])

    return statis





"""This code simulates scRNA-seq data under the GTM.
Inputs:
   param: A structure contains model parameters 
   param.kon and param.ron: OFF dwell time distribution f_off(t) = ron^(kon) * t^(kon-1) * e^(-ron * t) / gamma(kon)
   param.koff and param.roff: ON dwell time distribution f_on(t) = roff^(koff) * t^(koff-1) * e^(-roff * t) / gamma(koff)
   param.mu: Transcriptional rate
   param.delta: Degradation rate
   param.x0; The initial point of the simulation
   param.tottime: The total time of the simulation

Outputs:
   x - Time series of [OFF ON mRNA]
   t - Corresponds to time of x
"""
# def simulGTM(param):
#     koff = param['koff']
#     roff = param['roff']
#     kon = param['kon']
#     ron = param['ron']
#     mu = param['mu']
#     delta = param['delta']
#     x = param['x0']
#     tottime = param['tottime']

#     r_mu = np.matrix('-1. 1. 0.; 1. -1. 0; 0 0 1; 0 0 -1', dtype=np.float64)
    
#     t = 0

# while 
