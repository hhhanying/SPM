import numpy as np
import pymc3 as pm
import theano

def train_supervised(X, Y, b, alpha, mu_Mu, sigma2_Mu, alpha_Lambda, beta_Lambda, T, ntrace, nchain, nskip):
    '''
    Given training set and hyperparameters, estimate the parameters of the model.

    X, Y: training set

    b, alpha: hyperparameters for the prior of a and rho
    
    mu_Mu, sigma2_Mu: the prior of Mu is N(mu_Mu, sigma2_Mu/Lambda)
    alpha_Lambda, beta_Lambda: the prior of Lambda is Gamma(alpha_Lambda, beta_Lambda)

    T: transformation matrix. A list of np arrays.

    ntrace, nchain, nskip: parameters of sampling

    output: a directionary.
    a,rho, Mu Lambda: estimate of parameters
    G: estimate of memberships of training set
    '''

    N = len(Y) # sample size
    dg = len(alpha) # dim(untransformed membership)
    K = len(T[0])
    d = len(X[0]) # dim(x)
    T_arr = np.array(T) # nlabel * K * dg

    model = pm.Model()
    with model:
        Lambda_ = pm.Gamma("Lambda", alpha = alpha_Lambda, beta = beta_Lambda, shape = (K, d)) # notice np.random.gamma require shape and scale
        Mu_ = pm.Normal("Mu", mu = mu_Mu, sigma = pm.math.sqrt(sigma2_Mu / Lambda_), shape = (K, d)) # notice sigma rather than sigma^2
        S_ = 1 / Lambda_
        Tau_ = Mu_ / S_       

        a_ = pm.Exponential("a", b)  # mean = 1/b
        rho_ = pm.Dirichlet("rho", a = alpha)
        G_ = pm.Dirichlet("G", a = a_ * rho_, shape = (N, dg))        
        T_ = theano.shared(T_arr)
    
        for i in range(N):
            u_ = pm.math.dot(G_[i: (i + 1)], T_[Y[i]].T) # 1 * K
            Taux_ = pm.math.dot(u_, Tau_)
            Lambdax_ = pm.math.dot(u_, Lambda_)
            Sx_ = 1 / Lambdax_
            Mux_ = Sx_ * Taux_

            X_ = pm.Normal("x" + str(i), mu = Mux_, sigma = pm.math.sqrt(Sx_), observed = X[i])

        trace = pm.sample(ntrace, chains = nchain)

    nsample = ntrace * nchain
    nsave = nsample // nskip
    index_save = [i * nskip - 1 for i in range(1, 1 + nsave)]
    estimate = {}
    
    for para in ["a", 'rho', 'Mu', "Lambda", "G"]:
        estimate[para] = trace[para][index_save].mean(axis = 0)
    return estimate
