import numpy as np

def document_generator(a, rho, T, Lambda, Tau, N, w = None, seed = None):
    '''
    a, rho: corpus-level parameters
    T: transformation matrix. ntopic * K * dg
    Lambda, Tau: topics. K * d matrix. Lambda are positive.
    N: the number of documents.
    w: probability of labels
    seed: random seed number
    
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
    if not (seed is None): # if random seed is indicated, set seed
        np.random.seed(seed)

    nlabel = len(T) # number of classes
    d = len(Tau[0]) # dim(x)
    
    if w is None: # is w is not given, set w to be uniform
        w = np.ones(nlabel) /  nlabel 

    Y = np.random.choice(a = list(range(nlabel)), size = N, replace = True, p = w) # draw labels
    G = np.random.dirichlet(alpha = a * rho, size = N) # draw memberships
    U = np.array([np.dot(T[Y[i]], G[i]) for i in range(N)]) # tranform memberships, U is an N * K matrix

    LambdaX = np.dot(U, Lambda)
    TauX = np.dot(U, Tau)
    SigmaX = np.sqrt(1 / LambdaX) 
    MuX = TauX / LambdaX

    X = np.random.normal(loc = MuX, scale = SigmaX, size = (N, d))

    return X, Y, G, U