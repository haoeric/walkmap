walkmap <- function(X, pca = TRUE, perplexity=30){
    
    X = as.matrix(X)
    n <- nrow(X)
    
    ## Apply PCA to tidy the data
    if (pca) {
        cat('Apply PCA to tidy the data and initialize the embedding...')
        pca_result <- prcomp(X, retx=TRUE)
        X <- pca_result$x
    }
    
    ## P calculation, diagnol random walk probability matrix
    x2p = calx2p(X, perplexity, 1e-5)
    P <- x2p$P
    P = (P + t(P)) / 2
    for(i in 1:nrow(P)){
        for (j in 1:ncol(P)){
            if(i>=j) P[i,j] <- 0
        }
    }
    
    ## reweight probability with random walk
    newP <- matrix(0, n, n)    # initialize new P matrix
    thisWalk <- rep(0,n)       # initialize first walk vector
    pb <- txtProgressBar(min = 0, max = n, style = 3) # create progress bar
    for(i in 1:n){
        setTxtProgressBar(pb, i)
        thisWalk[i] <- 1
        while(sum(thisWalk) > 0){
            nextWalk <- randomWalk(thisWalk, P)
            thisWalk <- nextWalk
            newP[,i] <- newP[,i] + nextWalk
        }
    }
    close(pb)
    
    ## transform new pairwise probability to distance
    newP <- t(newP)
    ditsp <- -log10(newP)
    

}


randomWalk <- function(v, m){
    o <- v %*% m
    return(t(o)[,1])
}


## conditional probability between point i and j.
## for Pj|i, neighbors are picked in proportion to their probability density(perplexity) 
## under a gaussian centered at point i 

calx2p <- function(X, perplexity = 30, tol = 1e-4){
    cat("Calculate the pairwise P value...")
    n <- nrow(X)
    D <- as.matrix(dist(X))
    P <- matrix(0, n, n )		
    beta <- rep(1, n)           # beta = 0.5/sigma^2
    ## test with variate perplexity
    if(length(perplexity) < n)
        perplexity <- rep(perplexity, n)
    logU <- log2(perplexity)    # logU equal to entropy H
    
    for (i in 1:n){
        Di <- D[i, -i]
        beta[i] <- ifelse(i>1, beta[i-1], beta[i])
        betamin <- -Inf
        betamax <- Inf
        tries <- 0
        while(tries < 1000){
            hbeta <- calHbeta(Di, beta[i])
            H <- hbeta$H
            thisP <- hbeta$P
            Hdiff <- H - logU[i]
            
            if (abs(Hdiff) <= tol) break
            if (Hdiff > 0){
                betamin <- beta[i]
                if (is.infinite(betamax)){
                    beta[i] <- beta[i] * 2
                }else {
                    beta[i] <- (beta[i] + betamax)/2
                }
            } else{
                betamax <- beta[i]
                if (is.infinite(betamin)){
                    beta[i] <- beta[i]/ 2
                }else {
                    beta[i] <- (beta[i]+betamin)/2
                }
            }
            tries <- tries + 1
        }
        P[i,-i] <- thisP	
    }	
    cat("DONE!\n")
    r <- list()
    r$P <- P
    r$sigma <- sqrt(0.5/beta)
    r 
}



## calculate entropy H and probability P with beta
calHbeta <- function(D, beta){
    P = exp(-D^2 * beta)
    eps = 2^(-52) 
    sumP = sum(P)
    if (sumP == 0){
        P = D * 0
        H = 0
    } else {
        P = P/sumP
        P[P < eps] <- eps
        H = -sum(P*log2(P))
    }
    P[P <= eps] <- 0
    r = list()
    r$P = P
    r$H = H
    r
}