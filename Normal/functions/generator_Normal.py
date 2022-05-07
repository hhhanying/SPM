import numpy as np

def document_generator(a, rho, T, Lambda, Tau, N, nrare = 0, rare_rate = 0, prob_rate = 0):
    '''
    a, rho: corpus-level parameters
    T: transformation matrix. ntopic * K * dg
    Lambda, Tau: topics. K*d matrix. Lambda are positive.
    N: the number of documents.
    
    Lambda = 1/sigma^2
    Tau = mu/sigma^2
    or
    sigma = 1/Lambda
    mu = Tau/Lambda

    x|u ~ normal(lambda_x, tau_x), where lambda_x = sum(u_i * lambda_i) and tau_x = sum(u_i * tau_i)
    
    output: 
    X: N*d, X[i] = document[i]
    Y: Y[i] = label[i]
    G: membership
    U: transformed membership
    '''

    nlabel = len(T) # number of classes
    d = len(Tau[0]) # dim(x)
    
    Y = np.random.choice(list(range(nlabel)), N) # labels

    if prob_rate: # if we want an unbalanced dataset
        if rare_rate:
            nrare = int(nlabel * rare_rate)
        p = np.array([prob_rate] * nrare + [1] * (nlabel - nrare))
        p = p/p.sum()
        Y = np.random.choice(list(range(nlabel)), N, p = p) # labels
        
    G = np.random.dirichlet(a * rho, N)
    U = np.array([np.dot(T[Y[i]], G[i]) for i in range(N)])



    return X, Y, G, U