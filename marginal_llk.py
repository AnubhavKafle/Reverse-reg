
###################################################################
#
#Functions for optimizing sigma_beta for regularization parameter#
#
##################################################################
def det_lemma(Y, sigmax, sigmabeta):
    ngenes = Y.shape[0]
    sigmabeta2 = sigmabeta * sigmabeta
    sigmax2 = sigmax * sigmax
    Yt = Y.T # shape N x G
    A = np.dot(Yt, Yt.T) * sigmabeta2 / sigmax2
    A[np.diag_indices(A.shape[0])] += 1
    logdetA = np.linalg.slogdet(A)
    logdet = 2 * ngenes * (np.log(sigmax) - np.log(sigmabeta)) + logdetA[1]
    return logdet


def logml(X, Y, S, U, sigmabeta, sigmax):
    nsnps = X.shape[0]
    nsamples = X.shape[1]
    ngenes = Y.shape[0]
    sigmabeta2 = sigmabeta * sigmabeta
    sigmax2    = sigmax    * sigmax
    
    logdetA = det_lemma(Y , sigmax, sigmabeta) #A = np.dot(Yt.T, Yt)
    
    Smod = np.diag(np.square(S) / (np.square(S) + sigmax2 / sigmabeta2))
    W = np.dot(U, np.dot(Smod, U.T))
        
    partres = - 0.5 * ngenes * np.log(2 * np.pi * sigmabeta2) - 0.5 * logdetA
    
    snp_term = 0.5 * np.diag(np.dot(X, np.dot(W, X.T))) / sigmax2

    lml = nsnps * (partres) + np.sum(snp_term)
    
    return lml 

def gradient_lml (s,u,sigma_beta,sigma_geno,ngene,gt):  # in log(sigma beta) space
    sigmabeta2 = sigma_beta * sigma_beta
    sigmax2 = sigma_geno * sigma_geno 
    factor = sigmax2 / sigmabeta2
    nsamples = u.shape[0]
    nsnps = gt.shape[0]
    
    Smod = np.square(s) * sigmabeta2 / sigmax2
    second_term = (ngene - nsamples) + np.sum(1 / (Smod + 1))
    
    diagw = np.square(s) / np.square(np.square(s) + factor)
    Wsvd = np.dot(u, np.dot(np.diag(diagw), u.T))
    third_term  = np.diag(np.dot(gt, np.dot(Wsvd,  gt.T))) / sigmabeta2
    
    gradient =  nsnps * (- ngene + second_term ) + np.sum(third_term) 
    
    return gradient 
