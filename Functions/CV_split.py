import numpy as np

def CV_split(N, fold_number, no_fold, seed):
    '''
    This funuction is used to split dataset into training set and test set in CV.
    N: the sample size
    fold_number: how many folds
    no_fold: the order of current fold, starting from 1
    seed: random seed number to get consistent results in different fold
    Output: train_index and test_index
    '''
    np.random.seed(seed) # set seed

    fold_assign = np.random.choice(a = list(range(1, 1 + fold_number)), size = N, replace = True)
    
    train_index = [i for i in range(N) if fold_assign[i] == no_fold]
    test_index = [i for i in range(N) if fold_assign[i] != no_fold]
    
    return train_index, test_index