import numpy as np
from scipy import stats
import time

"""Sequential Monte Carlo Sampler for approximate Bayesian computaion

Inputs:
   N - the number of particles
   prior - prior distribution sampler
   f - function that generates simulated data give a parameters set
   rho - discrepancy metric, treated as a function of simulated data only
   epsilon - a sequence of discrepancy acceptance thresholds
   T - number of rounds
   proposal - proposal kernel process, generates new trials
   proposal_pdf - the PDF of the proposal kernel
   gene - sequence number of the gene.

Outputs:
   result - the sequence of particles
"""
def ABCSMCSampler(N,prior,f,rho,epsilon,T,proposal,proposal_pdf,gene):

    # initialize
    result, flag = ABCRejectionSampler(N,prior,f,rho,epsilon,T,gene)
    # result, flag = np.load('data/prior/gene_Atf4.npy', allow_pickle='True')

    # start_posterior = time.time()
    # print(f"Start Processing posterior: Gene = {gene}")

    
    W = (1 / N) * np.ones((N,T))
    param = {}

    # sequential sampling
    W = (1 / N) * np.ones((N,T))
    
    if flag == True: 
        for t in range(2, T + 1):
            # temp_start = time.time()
            # print(f"Start Processing Posterior: Gene = {gene}, T = {t}")


            temp_dist_array = np.array([d[0]['dist'] for d in result[:,t-2]])        
            epsilon = np.percentile(temp_dist_array, 50)

            # generate new generation of particles
            for i in range(1, N+1):

                # rejections steps
                r = np.inf
                total_time = 0

                while total_time < 10 and r > epsilon:
                    start = time.time()

                    j = np.random.choice(N, size = 1, replace = True, p = W[:, t-2])[0]
                    param_temp = result[j,t-2]

                    # generate new particle based on proposal kernel
                    param_proposal = proposal(
                        np.log(
                            [
                                param_temp[0]['kon'],
                                param_temp[0]['ron'],
                                param_temp[0]['koff'],
                                param_temp[0]['roff'],
                                param_temp[0]['mu']
                            ]
                        )
                    )
                    
                    param['kon'] = param_proposal[0]
                    param['ron'] = param_proposal[1]
                    param['koff'] = param_proposal[2]
                    param['roff'] = param_proposal[3]
                    param['mu'] = param_proposal[4]
                    param['delta'] = 1

                    result[i-1,t-1][0] = {}
                    static_temp = f(param)

                    if (len(static_temp[static_temp<0]) > 0) or (np.sum(static_temp) == 0.):
                        static_temp = np.abs(static_temp)
                        
                    r = rho(static_temp)
                    
                    result[i-1,t-1][0]['kon'] = param['kon']
                    result[i-1,t-1][0]['ron'] = param['ron']
                    result[i-1,t-1][0]['koff'] = param['koff']
                    result[i-1,t-1][0]['roff'] = param['roff']
                    result[i-1,t-1][0]['mu'] = param['mu']
                    result[i-1,t-1][0]['delta'] = 1
                    result[i-1,t-1][0]['dist'] = r

                    end = time.time()
                    elapsedTime = end - start
                    total_time = total_time + elapsedTime
                
                if total_time > 10:
                    flag = False
                    break

                # recompute particle weight using optimal backward kernel
                back_K = 0

                # test = 0
                for j in range (1 , N + 1):
                    # print('--------------------------')
                    # print(result[i-1,t-1][0])
                    # print(result[j-1,t-2][0])
                    # print('--------------------------')

                    # if(test > 5): break
                    # test+=1

                    back_K = back_K + W[j-1,t-2] * proposal_pdf(
                        result[i-1,t-1][0]['kon'],
                        result[j-1,t-2][0]['kon'], 
                        result[i-1,t-1][0]['ron'],
                        result[j-1,t-2][0]['ron'], 
                        result[i-1,t-1][0]['koff'],
                        result[j-1,t-2][0]['koff'], 
                        result[i-1,t-1][0]['roff'],
                        result[j-1,t-2][0]['roff'], 
                        result[i-1,t-1][0]['mu'],
                        result[j-1,t-2][0]['mu']
                    )
                W[i-1,t-1] = ((1/5)*(1/5)*(1/30)*1/(4*result[i-1,t-1][0]['ron'] * result[i-1,t-1][0]['roff']))/back_K
            
            if flag == False:
                break

            # resample
            if t < T: #need to validate
                result_rs  = result[:,t-1]
                W[:,t-1] = W[:,t-1]/sum(W[:,t-1])
                J =  np.random.choice(N, size = N, replace = True, p = W[:, t-1])
                result[:,t-1] = result_rs[J]
            
            # re-set weights
            W[:,t-1] = 1./N

            # temp_end = time.time()
            # temp_time = temp_end - temp_start
            # print(f"End Processing Posterior: Gene = {gene}, T = {t}, time = {temp_time}")
            # print()
    
    # end_posterior = time.time()
    # posterior_time = end_posterior - start_posterior
    # print(f"Finished Processing posterior: Time = {posterior_time}")
    # print()
    return (result, flag)   


"""Rejection Sampler for approximate Bayesian computaion

Inputs:
   N - the number of ABC posterior samples
   prior - function that generates iid samples from the parameter joint
       prior distribution
   f - function that computes statictis given a parameters set
   rho - discrepancy metric, treated as a function of simulated data only
   epsilon - the discrepancy acceptance threshold
   T - number of rounds
   gene - sequence number of the gene.

Outputs:
   result - a matrix of ABC prior samples
"""
def ABCRejectionSampler(N,prior,f,rho,epsilon,T,gene):
    # start_prior = time.time()
    # print(f"Start Processing Prior: Gene = {gene}")
    result = np.empty([N,T,1], dtype=object)
    result0 = []
    total_time = 0.
    flag = True


    size_val = 0
    # test = 1

    while total_time < 5*N and size_val < 5*N:
        start = time.time()
        param = {}

        # generate trial from the prior
        # e.g = theta_trial = [1.30695591  0.10584372  1.34860353  3.42510636 96.17662629  1.]
        theta_trial = prior()
        param['kon'] = theta_trial[0]
        param['ron'] = theta_trial[1]
        param['koff'] = theta_trial[2]
        param['roff'] = theta_trial[3]
        param['mu'] = theta_trial[4]
        param['delta'] = theta_trial[5]

        # compute theorical statictis of parameters
        static_theo = f(param)

        if len(static_theo[static_theo<0]) > 0:
            continue

        

        dist = rho(static_theo)

        param['dist'] = dist


        # accept or reject
        if dist <= epsilon:

            # print('--------------------------')

            # print(static_theo)
            # print(dist)


            # print('--------------------------')

            # if(test > 5): break
            # test+=1

            # print(size_val, dist)
            result_obj = {}
            result_obj['kon'] = param['kon']
            result_obj['ron'] = param['ron']
            result_obj['koff'] = param['koff']
            result_obj['roff'] = param['roff']
            result_obj['mu'] = param['mu']
            result_obj['delta'] = param['delta']
            result_obj['dist'] = param['dist']
            
            result0.append(result_obj)
            size_val += 1
            
        
        end = time.time()
        elapsedTime = end - start
        total_time = total_time + elapsedTime

    if total_time > 5*N:
        flag = False
        print(f'Gene {gene}:wrong!!\n')
        # time.sleep(5)
    else:
        index0 = np.argsort([d['dist'] for d in result0])
        result0 = np.take(result0, index0[:N])
        for index2 in range (1, N + 1):
            result[index2-1,0] = result0[index2-1]
    
    # np.save(f'data/prior/gene_{gene}.npy', np.array([result, flag], dtype=object))

    # end_prior = time.time()
    # prior_time = end_prior - start_prior
    # print(f"Finished Processing Prior: Time = {prior_time}")
    # print()

    return (result, flag)
