# All codes reference the cluster methods studied in my scientific initiation project

# sigma 2 estimation

sigma2est1 = function(X, frac = .5) {
  
  X <- as.matrix(X)
  h = 1:ncol(X)
  xmin <- X[,which(h%%2==1)]
  xmax <- X[,which(h%%2==0)]
  
  n = nrow(xmin)
  m = floor(n*frac)
  idx1 = sample(1:n, m, replace = T)
  idx2 = sample(1:n, m, replace = T)
  tmpmin = (xmin[idx1,, drop = F] - xmin[idx2,, drop = F])^2
  tmpmax = (xmax[idx1,, drop = F] - xmax[idx2,, drop = F])^2
  dist = apply(tmpmin, 1, sum)+apply(tmpmax, 1, sum)
  mean(quantile(dist[dist != 0], probs = c(.9, .1)))
}

KM = function(X, ngrupos, maxiter = 100) {
  
  euclidean.distance = function(A, B){
    n1 = nrow(A)
    n2 = nrow(B)
    indic = as.logical( rep(0:1,(ncol(A)/2)) )
    AL = A[,!indic]
    BL = B[,!indic]
    AS = A[,indic]
    BS = B[,indic]
    dist = matrix(NA, n1, n2)
    for (i in 1:n1)
      for (j in 1:n2)
        dist[i,j] = t(AL[i,]-BL[j,])%*%(AL[i,]-BL[j,])+t(AS[i,]-BS[j,])%*%(AS[i,]-BS[j,])
    dist
  }
  
  X = as.matrix(X)
  n = nrow(X)
  K = ngrupos
  
  # InicializaÃ§Ã£o
  
  smp = sample.int(n, K)
  gtmp = X[smp,]
  
  d2 = euclidean.distance(X, gtmp)
  P = apply(d2, 1, which.min)
  J = 0
  for (k in 1:K)
    J = J+sum(d2[P == k, k])
  # Etapa iterativa
  t = 0
  repeat {
    t = t+1
    # CÃ¡lculo dos protÃ³tipos
    g = gtmp
    for (k in 1:K) {
      if (sum(P == k) > 0) {
        if (sum(P == k) == 1)
          g[k,] = X[P == k,]
        else
          g[k,] = apply(X[P == k,], 2, mean)
      }
    }
    gtmp = g
    # CÃ¡lculo das distÃ¢ncias
    d2 = euclidean.distance(X, g)
    # CÃ¡lculo da partiÃ§Ã£o
    Pnew = apply(d2, 1, which.min)
    # CritÃ©rio de parada
    if (all(Pnew == P) || t > maxiter) break
    else {
      P = Pnew
      Jnew = 0
      for (k in 1:K)
        Jnew = Jnew+sum(d2[P == k, k])
    }
    J = c(J, Jnew)
  }
  (result = list(partition = P, prototypes = g, iterations = t, criterion = J))
}



KKMEC_1 <- function(X,ngrupos,sig2,maxinter = 50){
  
  
  # Função RBF de uma componente
  
  RBF1 <- function(ai,bi,al,bl,sig2){
    
    ai <- as.vector(ai)
    bi <- as.vector(bi)
    al <- as.vector(al)
    bl <- as.vector(bl)
    
    norma <- sum(  (ai-al)^2 + (bi-bl)^2 )
    const <- 1/(2*sig2)
    return( exp(-const * norma) )
  }
  
  
  matriz.kernel <- function(Ai,Bi,Ga,Gb,sig2){
    
    n <- nrow(Ai); K <- nrow(Ga); p <- ncol(Ai)
    
    KMs <- array(NA,c(n,K,p))
    for(j in 1:p)
      for(k in 1:K)
        for(i in 1:n)
          KMs[i,k,j] <- RBF1(Ai[i,j],Bi[i,j],Ga[k,j],Gb[k,j],sig2)   
    
    return(KMs)
  }
  
  
  # Calculando a distancia para cada grupo individual( para todas as variaveis em cada grupo)
  
  dist2<- function(MJ,Part,n,p,K){
    
    # MJ é um array com p matrizes n x K;
    # Part é um vetor de partição;
    # n,p e K descrevem o número de observações, variáveis e Grupos, respectivamente.
    
    dist.square <- array(NA,c(n,p,K))
    
    for(k in 1:K){
      
      P <- which( Part == k)
      card <- length(P)
      
      if(card == 0)
        next
      
      for(j in 1:p){
        
        Const <- MJ[P,P,j]
        for(i in 1:n)
          dist.square[i,j,k] <- 1 -  (2 * sum( MJ[i,P,j]) )/ card +  (sum(Const))/ (card^2)
        
      }
      
    }
    
    dist.final <- matrix(NA,n,p)
    for(k in 1:K)
      dist.final[,k] <- apply(dist.square[,,k],1, sum)
    
    return(dist.final)
  }
  
  
  
  ############################################################################################
  
  X <- as.matrix(X)
  n <- nrow(X); p <- (ncol(X)/2); K <- ngrupos
  
  h <- c(1:ncol(X))
  A <- X[,h[h%%2==1]]
  B <- X[,h[h%%2==0]]
  
  
  # Selecionando uma partição aleatória
  P <- sample(K,n,replace = TRUE)
  
  # Matriz de Kernel
  KM <- matriz.kernel(A,B,A,B,sig2) 
  
  # Calculando a Distancia
  
  d2 <- dist2(KM,P,n,p,K)
  P <- apply(d2, 1, which.min)
  
  J <- 0
  for(k in 1:K)
    J <- J + sum( d2[P == k,k] )
  
  
  inte <- 0
  repeat{
    
    inte <- inte + 1
    
    
    # Calculando a Distancia
    d2 <- dist2(KM,P,n,p,K)
    Pn <- apply(d2, 1, which.min)
    
    
    if(all(Pn == P) || inte >= maxinter) break
    else{
      P <- Pn
      
      Jn <- 0
      for(k in 1:K)
        Jn <- Jn + sum(d2[P == k,k])
    }
    
    J <- c(J,Jn)
    
    
  }
  
  Jsort = sort(J,decreasing = TRUE)
  
  if( all(Jsort == J))
    return( list(Partição = P, Iterações = inte, Critério = J) )
  else
    stop("Error: the objective function did not decrease")
}

KKMEC_2 <- function(X,ngrupos,sig2,maxinter = 50){
  
  RBF2 <- function(ai,bi,al,bl,sig2){
    
    ai <- as.vector(ai)
    bi <- as.vector(bi)
    al <- as.vector(al)
    bl <- as.vector(bl)
    
    norma.inf <- sum(  (ai-al)^2 )
    norma.sup <- sum(  (bi-bl)^2 )
    const <- 1/(2*sig2)
    
    return( exp(-const * norma.inf) +  exp(-const * norma.sup) )
  }
  
  
  matriz.kernel <- function(Ai,Bi,Ga,Gb,sig2){
    
    n <- nrow(Ai); K <- nrow(Ga); p <- ncol(Ai)
    
    KMs <- array(NA,c(n,K,p))
    for(j in 1:p)
      for(k in 1:K)
        for(i in 1:n)
          KMs[i,k,j] <- RBF2(Ai[i,j],Bi[i,j],Ga[k,j],Gb[k,j],sig2)   
    
    return(KMs)
  }
  
  
  # Calculando a distancia para cada grupo individual( para todas as variaveis em cada grupo)
  
  dist2<- function(MJ,Part,n,p,K){
    
    # MJ é um array com p matrizes n x K;
    # Part é um vetor de partição;
    # n,p e K descrevem o número de observações, variáveis e Grupos, respectivamente.
    
    dist.square <- array(NA,c(n,p,K))
    
    for(k in 1:K){
      
      P <- which( Part == k)
      card <- length(P)
      
      if(card == 0)
        next
      
      for(j in 1:p){
        
        Const <- MJ[P,P,j]
        for(i in 1:n)
          dist.square[i,j,k] <- 2 -  (2 * sum( MJ[i,P,j]) )/ card +  (sum(Const))/ (card^2)
        
      }
      
    }
    
    dist.final <- matrix(NA,n,p)
    for(k in 1:K)
      dist.final[,k] <- apply(dist.square[,,k],1, sum)
    
    return(dist.final)
  }
  
  #######################################################################
  
  X <- as.matrix(X)
  n <- nrow(X); p <- (ncol(X)/2); K <- ngrupos
  
  h <- c(1:ncol(X))
  A <- X[,h[h%%2==1]]
  B <- X[,h[h%%2==0]]
  
  # Selecionando uma partição aleatória
  P <- sample(K,n,replace = TRUE)
  
  # Matriz de Kernel
  KM <- matriz.kernel(A,B,A,B,sig2) 
  
  # Calculando a Distancia
  
  d2 <- dist2(KM,P,n,p,K)
  P <- apply(d2, 1, which.min)
  
  J <- 0
  for(k in 1:K)
    J <- J + sum( d2[P == k,k] )
  
  
  inte <- 0
  repeat{
    
    inte <- inte + 1
    
    
    # Calculando a Distancia
    d2 <- dist2(KM,P,n,p,K)
    Pn <- apply(d2, 1, which.min)
    
    
    if(all(Pn == P) || inte >= maxinter) break
    else{
      P <- Pn
      
      Jn <- 0
      for(k in 1:K)
        Jn <- Jn + sum(d2[P == k,k])
    }
    
    J <- c(J,Jn)
    
    
  }
  
  Jsort = sort(J,decreasing = TRUE)
  
  if( all(Jsort == J))
    return( list(Partição = P, Iterações = inte, Critério = J) )
  else
    stop("Error: the objective function did not decrease")
  
}



KKMEC_GP1 <- function(X,ngrupos,sig2,maxinter = 50){
  
  RBF1 <- function(ai,bi,al,bl,sig2){
    
    ai <- as.vector(ai)
    bi <- as.vector(bi)
    al <- as.vector(al)
    bl <- as.vector(bl)
    
    norma <- sum(  (ai-al)^2 + (bi-bl)^2 )
    const <- 1/(2*sig2)
    return( exp(-const * norma) )
  }
  
  matriz.kernel <- function(Ai,Bi,Ga,Gb,sig2){
    
    n <- nrow(Ai); K <- nrow(Ga); p <- ncol(Ai)
    
    KMs <- array(NA,c(n,K,p))
    for(j in 1:p)
      for(k in 1:K)
        for(i in 1:n)
          KMs[i,k,j] <- RBF1(Ai[i,j],Bi[i,j],Ga[k,j],Gb[k,j],sig2)   
    
    return(KMs)
  }
  
  
  # Cálculo da distância
  
  dist2<- function(MJ,Part,n,p,K){
    
    # MJ é um array com p matrizes n x n;
    
    # Part é um vetor de partição;
    
    # n,p e K descrevem o número de observações, variáveis e Grupos, respectivamente;
    
    dist.square <- array(NA,c(n,p,K))
    
    for(k in 1:K){
      
      P <- which( Part == k)
      card <- length(P)
      
      for(j in 1:p){
        
        Const <- sum( MJ[P,P,j] )
        for(i in 1:n)
          dist.square[i,j,k] <-  1 -  (2 * sum( MJ[i,P,j]) )/ card +  (Const)/ (card^2) 
        
      }
      
    }
    
    return(dist.square)
  }
  
  
  update.weights <- function(MK,Part,K,p){
    
    # MK é uma array com K matrizes n x p referente a distância
    # Part é uma partição 
    
    Soma <- matrix(NA,K,p)
    for(k in 1:K){
      
      P <- which(Part == k)
      mtk <- MK[,,k][P,]
      
      if(length(P) == 1)
        Soma[k,] <- mtk
      else
        Soma[k,] <- apply( mtk, 2,sum)
    }
    
    Sum.final <- apply(Soma,2,sum)
    Const <- prod(Sum.final) ^ (1/p)
    Pesos <- Const/Sum.final
    return(Pesos)
  }
  
  # Estimando sig2
  
  
  #####################################################################################  
  
  X <- as.matrix(X)
  n <- nrow(X); p <- (ncol(X)/2); K <- ngrupos
  
  h <- c(1:ncol(X))
  A <- X[,h[h%%2==1]]
  B <- X[,h[h%%2==0]]
  
  
  
  # Selecionando uma partição aleatória
  P <- sample(K,n,replace = TRUE)
  
  # Matriz de Kernel (não mudará ao longo do algoritmo)
  KM <- matriz.kernel(A,B,A,B,sig2) 
  
  #Inicializando a Matriz de Pesos
  wgh <- rep(1,p)
  
  # Calculando a Distancia Para cada váriavel
  d2k <- dist2(KM,P,n,p,K)
  
  # Calculando a distância final (Somando para a todas as variáveis)
  
  d2 <- matrix(NA,n,K)
  for(k in 1:K)
    d2[,k] <- apply(d2k[,,k],1 ,sum)
  
  # Obtendo a partição
  P <- apply(d2, 1, which.min)
  
  J <- 0
  for(k in 1:K)
    J <- J + sum( d2[P == k,k] )
  
  inte <- 0
  repeat{
    
    inte <- inte + 1
    # Calculando a distância para cada grupo e para cada variável
    d2k <- dist2(KM,P,n,p,K)
    
    # Encontrando os pesos das variáveis
    wgh <- update.weights(d2k,P,K,p)
    
    # Calculando a distância final (Somando para a todas as variáveis e já ponderada pelos pesos)
    
    d2 <- matrix(NA,n,K)
    for(k in 1:K)
      for(i in 1:n)
        d2[i,k] <- wgh %*% d2k[i,,k]
    
    # Selecionando a partição
    Pn <- apply(d2, 1, which.min)
    
    
    if(all(Pn == P) || inte >= maxinter) break
    else{
      P <- Pn
      
      Jn <- 0
      for(k in 1:K)
        Jn <- Jn + sum(d2[P == k,k])
    }
    
    J <- c(J,Jn)
    
    
  }
  
  name.var <- colnames(X)
  name.var <- name.var[h[h%%2==1]]
  names(wgh) <- name.var
  
  
  Jsort = sort(J,decreasing = TRUE)
  
  if( all(Jsort == J))
    return( list(Partição = P, Iterações = inte, Critério = J) )
  else
    stop("Error: the objective function did not decrease")
  
}


KKMEC_GP2 <- function(X,ngrupos,sig2,maxinter = 50){
  
  RBF2 <- function(ai,bi,al,bl,sig2){
    
    ai <- as.vector(ai)
    bi <- as.vector(bi)
    al <- as.vector(al)
    bl <- as.vector(bl)
    
    norma.inf <- sum(  (ai-al)^2 )
    norma.sup <- sum(  (bi-bl)^2 )
    const <- 1/(2*sig2)
    
    return( exp(-const * norma.inf) +  exp(-const * norma.sup) )
  }
  
  matriz.kernel <- function(Ai,Bi,Ga,Gb,sig2){
    
    n <- nrow(Ai); K <- nrow(Ga); p <- ncol(Ai)
    
    KMs <- array(NA,c(n,K,p))
    for(j in 1:p)
      for(k in 1:K)
        for(i in 1:n)
          KMs[i,k,j] <- RBF2(Ai[i,j],Bi[i,j],Ga[k,j],Gb[k,j],sig2)   
    
    return(KMs)
  }
  
  
  # Cálculo da distância
  
  dist2<- function(MJ,Part,n,p,K){
    
    # MJ é um array com p matrizes n x n;
    
    # Part é um vetor de partição;
    
    # n,p e K descrevem o número de observações, variáveis e Grupos, respectivamente;
    
    dist.square <- array(NA,c(n,p,K))
    
    for(k in 1:K){
      
      P <- which( Part == k)
      card <- length(P)
      
      for(j in 1:p){
        
        Const <- sum( MJ[P,P,j] )
        for(i in 1:n)
          dist.square[i,j,k] <-  2 -  (2 * sum( MJ[i,P,j]) )/ card +  (Const)/ (card^2) 
        
      }
      
    }
    
    return(dist.square)
  }
  
  
  update.weights <- function(MK,Part,K,p){
    
    # MK é uma array com K matrizes n x p referente a distância
    # Part é uma partição 
    
    Soma <- matrix(NA,K,p)
    for(k in 1:K){
      
      P <- which(Part == k)
      mtk <- MK[,,k][P,]
      
      if(length(P) == 1)
        Soma[k,] <- mtk
      else
        Soma[k,] <- apply( mtk, 2,sum)
    }
    
    Sum.final <- apply(Soma,2,sum)
    Const <- prod(Sum.final) ^ (1/p)
    Pesos <- Const/Sum.final
    return(Pesos)
  }
  
  
  ############################################################################
  
  X <- as.matrix(X)
  n <- nrow(X); p <- (ncol(X)/2); K <- ngrupos
  
  h <- c(1:ncol(X))
  A <- X[,h[h%%2==1]]
  B <- X[,h[h%%2==0]]
  
  # Selecionando uma partição aleatória
  P <- sample(K,n,replace = TRUE)
  
  # Matriz de Kernel (não mudará ao longo do algoritmo)
  KM <- matriz.kernel(A,B,A,B,sig2) 
  
  #Inicializando a Matriz de Pesos
  wgh <- rep(1,p)
  
  # Calculando a Distancia Para cada váriavel
  d2k <- dist2(KM,P,n,p,K)
  
  # Calculando a distância final (Somando para a todas as variáveis)
  
  d2 <- matrix(NA,n,K)
  for(k in 1:K)
    d2[,k] <- apply(d2k[,,k],1 ,sum)
  
  # Obtendo a partição
  P <- apply(d2, 1, which.min)
  
  J <- 0
  for(k in 1:K)
    J <- J + sum( d2[P == k,k] )
  
  inte <- 0
  repeat{
    
    inte <- inte + 1
    # Calculando a distância para cada grupo e para cada variável
    d2k <- dist2(KM,P,n,p,K)
    
    # Encontrando os pesos das variáveis
    wgh <- update.weights(d2k,P,K,p)
    
    # Calculando a distância final (Somando para a todas as variáveis e já ponderada pelos pesos)
    
    d2 <- matrix(NA,n,K)
    for(k in 1:K)
      for(i in 1:n)
        d2[i,k] <- wgh %*% d2k[i,,k]
    
    # Selecionando a partição
    Pn <- apply(d2, 1, which.min)
    
    
    if(all(Pn == P) || inte >= maxinter) break
    else{
      P <- Pn
      
      Jn <- 0
      for(k in 1:K)
        Jn <- Jn + sum(d2[P == k,k])
    }
    
    J <- c(J,Jn)
    
    
  }
  
  name.var <- colnames(X)
  name.var <- name.var[h[h%%2==1]]
  names(wgh) <- name.var
  
  
  Jsort = sort(J,decreasing = TRUE)
  
  if( all(Jsort == J))
    return( list(Partição = P, Iterações = inte, Critério = J) )
  else
    stop("Error: the objective function did not decrease")

}





KKMEC_GS1 <- function(X,ngrupos,sig2,beta,maxinter = 50){
  
  RBF1 <- function(ai,bi,al,bl,sig2){
    
    ai <- as.vector(ai)
    bi <- as.vector(bi)
    al <- as.vector(al)
    bl <- as.vector(bl)
    
    norma <- sum(  (ai-al)^2 + (bi-bl)^2 )
    const <- 1/(2*sig2)
    return( exp(-const * norma) )
  }
  
  matriz.kernel <- function(Ai,Bi,Ga,Gb,sig2){
    
    n <- nrow(Ai); K <- nrow(Ga); p <- ncol(Ai)
    
    KMs <- array(NA,c(n,K,p))
    for(j in 1:p)
      for(k in 1:K)
        for(i in 1:n)
          KMs[i,k,j] <- RBF1(Ai[i,j],Bi[i,j],Ga[k,j],Gb[k,j],sig2)   
    
    return(KMs)
  }
  
  
  # Cálculo da distância
  
  dist2<- function(MJ,Part,n,p,K){
    
    # MJ é um array com p matrizes n x n;
    
    # Part é um vetor de partição;
    
    # n,p e K descrevem o número de observações, variáveis e Grupos, respectivamente;
    
    dist.square <- array(NA,c(n,p,K))
    
    for(k in 1:K){
      
      P <- which( Part == k)
      card <- length(P)
      
      for(j in 1:p){
        
        Const <- sum( MJ[P,P,j] )
        for(i in 1:n)
          dist.square[i,j,k] <-  1 -  (2 * sum( MJ[i,P,j]) )/ card +  (Const)/ (card^2) 
        
      }
      
    }
    
    return(dist.square)
  }
  
  
  update.weights <- function(MK,beta,Part,K,p){
    
    # MK é uma arry com K matrizes n x p referente a distâncias
    # Part é uma partição
    
    Soma <- matrix(NA,K,p); Pesos <- NULL
    for(k in 1:K){
      
      P <- which(Part == k)
      mtk <- MK[,,k][P,]
      
      if(length(P) == 1)
        Soma[k,] <- mtk
      else
        Soma[k,] <- apply(mtk, 2, sum)
    }   
    
    sum.final <- apply(Soma,2, sum)
    
    for(j in 1:p){
      
      ratio.sum <- (sum.final[j]/sum.final)^( 1/(beta-1) )
      Pesos[j] <- (sum(ratio.sum))^(-1)
      
    }
    
    return(Pesos)
  }
  
  
  ################################################################################
  
  X <- as.matrix(X)
  n <- nrow(X); p <- (ncol(X)/2); K <- ngrupos
  
  h <- c(1:ncol(X))
  A <- X[,h[h%%2==1]]
  B <- X[,h[h%%2==0]]
  
  
  # Selecionando uma partição aleatória
  P <- sample(K,n,replace = TRUE)
  
  # Matriz de Kernel (não mudará ao longo do algoritmo)
  KM <- matriz.kernel(A,B,A,B,sig2) 
  
  #Inicializando a Matriz de Pesos
  wgh <- rep(1/p,p)
  
  # Calculando a Distancia Para cada váriavel
  d2k <- dist2(KM,P,n,p,K)
  
  # Calculando a distância final (Somando para a todas as variáveis)
  
  d2 <- matrix(NA,n,K)
  for(k in 1:K)
    for(i in 1:n)
      d2[i,k] <- ( wgh ^ beta ) %*% d2k[i,,k]
  
  # Obtendo a partição
  P <- apply(d2, 1, which.min)
  
  J <- 0
  for(k in 1:K)
    J <- J + sum( d2[P == k,k] )
  
  inte <- 0
  repeat{
    
    inte <- inte + 1
    # Calculando a distância para cada grupo e para cada variável
    d2k <- dist2(KM,P,n,p,K)
    
    # Encontrando os pesos das variáveis
    wgh <- update.weights(d2k,beta,P,K,p)
    
    # Calculando a distância final (Somando para a todas as variáveis e já ponderada pelos pesos)
    
    d2 <- matrix(NA,n,K)
    for(k in 1:K)
      for(i in 1:n)
        d2[i,k] <- ( wgh ^ beta ) %*% d2k[i,,k]
    
    # Selecionando a partição
    Pn <- apply(d2, 1, which.min)
    
    
    if(all(Pn == P) || inte >= maxinter) break
    else{
      P <- Pn
      
      Jn <- 0
      for(k in 1:K)
        Jn <- Jn + sum(d2[P == k,k])
    }
    
    J <- c(J,Jn)
    
    
  }
  
  name.var <- colnames(X)
  name.var <- name.var[h[h%%2==1]]
  names(wgh) <- name.var
  
  Jsort = sort(J,decreasing = TRUE)
  
  if( all(Jsort == J))
    return( list(Partição = P, Iterações = inte, Critério = J) )
  else
    stop("Error: the objective function did not decrease")

}





KKMEC_GS2 <- function(X,ngrupos,beta,sig2,maxinter = 50){
  
  RBF2 <- function(ai,bi,al,bl,sig2){
    
    ai <- as.vector(ai)
    bi <- as.vector(bi)
    al <- as.vector(al)
    bl <- as.vector(bl)
    
    norma.inf <- sum(  (ai-al)^2 )
    norma.sup <- sum(  (bi-bl)^2 )
    const <- 1/(2*sig2)
    
    return( exp(-const * norma.inf) +  exp(-const * norma.sup) )
  }
  
  matriz.kernel <- function(Ai,Bi,Ga,Gb,sig2){
    
    n <- nrow(Ai); K <- nrow(Ga); p <- ncol(Ai)
    
    KMs <- array(NA,c(n,K,p))
    for(j in 1:p)
      for(k in 1:K)
        for(i in 1:n)
          KMs[i,k,j] <- RBF2(Ai[i,j],Bi[i,j],Ga[k,j],Gb[k,j],sig2)   
    
    return(KMs)
  }
  
  
  # Cálculo da distância
  
  dist2<- function(MJ,Part,n,p,K){
    
    # MJ é um array com p matrizes n x n;
    
    # Part é um vetor de partição;
    
    # n,p e K descrevem o número de observações, variáveis e Grupos, respectivamente;
    
    dist.square <- array(NA,c(n,p,K))
    
    for(k in 1:K){
      
      P <- which( Part == k)
      card <- length(P)
      
      for(j in 1:p){
        
        Const <- sum( MJ[P,P,j] )
        for(i in 1:n)
          dist.square[i,j,k] <-  2 -  (2 * sum( MJ[i,P,j]) )/ card +  (Const)/ (card^2) 
        
      }
      
    }
    
    return(dist.square)
  }
  
  
  update.weights <- function(MK,beta,Part,K,p){
    
    # MK é uma arry com K matrizes n x p referente a distâncias
    # Part é uma partição
    
    Soma <- matrix(NA,K,p); Pesos <- NULL
    for(k in 1:K){
      
      P <- which(Part == k)
      mtk <- MK[,,k][P,]
      
      if(length(P) == 1)
        Soma[k,] <- mtk
      else
        Soma[k,] <- apply(mtk, 2, sum)
    }   
    
    sum.final <- apply(Soma,2, sum)
    
    for(j in 1:p){
      
      ratio.sum <- (sum.final[j]/sum.final)^( 1/(beta-1) )
      Pesos[j] <- (sum(ratio.sum))^(-1)
      
    }
    
    return(Pesos)
  }
  
  
  #############################################################################
  
  X <- as.matrix(X)
  n <- nrow(X); p <- (ncol(X)/2); K <- ngrupos
  
  h <- c(1:ncol(X))
  A <- X[,h[h%%2==1]]
  B <- X[,h[h%%2==0]]
  
  # Selecionando uma partição aleatória
  P <- sample(K,n,replace = TRUE)
  
  # Matriz de Kernel (não mudará ao longo do algoritmo)
  KM <- matriz.kernel(A,B,A,B,sig2) 
  
  #Inicializando a Matriz de Pesos
  wgh <- rep(1/p,p)
  
  # Calculando a Distancia Para cada váriavel
  d2k <- dist2(KM,P,n,p,K)
  
  # Calculando a distância final (Somando para a todas as variáveis)
  
  d2 <- matrix(NA,n,K)
  for(k in 1:K)
    for(i in 1:n)
      d2[i,k] <- ( wgh ^ beta ) %*% d2k[i,,k]
  
  # Obtendo a partição
  P <- apply(d2, 1, which.min)
  
  J <- 0
  for(k in 1:K)
    J <- J + sum( d2[P == k,k] )
  
  inte <- 0
  repeat{
    
    inte <- inte + 1
    # Calculando a distância para cada grupo e para cada variável
    d2k <- dist2(KM,P,n,p,K)
    
    # Encontrando os pesos das variáveis
    wgh <- update.weights(d2k,beta,P,K,p)
    
    # Calculando a distância final (Somando para a todas as variáveis e já ponderada pelos pesos)
    
    d2 <- matrix(NA,n,K)
    for(k in 1:K)
      for(i in 1:n)
        d2[i,k] <- ( wgh ^ beta ) %*% d2k[i,,k]
    
    # Selecionando a partição
    Pn <- apply(d2, 1, which.min)
    
    
    if(all(Pn == P) || inte >= maxinter) break
    else{
      P <- Pn
      
      Jn <- 0
      for(k in 1:K)
        Jn <- Jn + sum(d2[P == k,k])
    }
    
    J <- c(J,Jn)
    
    
  }
  
  name.var <- colnames(X)
  name.var <- name.var[h[h%%2==1]]
  names(wgh) <- name.var
  
  
  Jsort = sort(J,decreasing = TRUE)
  
  if( all(Jsort == J))
    return( list(Partição = P, Iterações = inte, Critério = J) )
  else
    stop("Error: the objective function did not decrease")
  
}





KKMEC_LP1 <- function(X,ngrupos,sig2,maxinter = 50){
  
  RBF1 <- function(ai,bi,al,bl,sig2){
    
    ai <- as.vector(ai)
    bi <- as.vector(bi)
    al <- as.vector(al)
    bl <- as.vector(bl)
    
    norma <- sum(  (ai-al)^2 + (bi-bl)^2 )
    const <- 1/(2*sig2)
    return( exp(-const * norma) )
  }
  
  matriz.kernel <- function(Ai,Bi,Ga,Gb,sig2){
    
    n <- nrow(Ai); K <- nrow(Ga); p <- ncol(Ai)
    
    KMs <- array(NA,c(n,K,p))
    for(j in 1:p)
      for(k in 1:K)
        for(i in 1:n)
          KMs[i,k,j] <- RBF1(Ai[i,j],Bi[i,j],Ga[k,j],Gb[k,j],sig2)   
    
    return(KMs)
  }
  
  
  # Cálculo da distância
  
  dist2<- function(MJ,Part,n,p,K){
    
    # MJ é um array com p matrizes n x n;
    
    # Part é um vetor de partição;
    
    # n,p e K descrevem o número de observações, variáveis e Grupos, respectivamente;
    
    dist.square <- array(NA,c(n,p,K))
    
    for(k in 1:K){
      
      P <- which( Part == k)
      card <- length(P)
      
      for(j in 1:p){
        
        Const <- sum( MJ[P,P,j] )
        for(i in 1:n)
          dist.square[i,j,k] <-  1 -  (2 * sum( MJ[i,P,j]) )/ card +  (Const)/ (card^2) 
        
      }
      
    }
    
    return(dist.square)
  }
  
  
  update.weights <- function(MK,Pesos,Part,K,p){
    
    # MK é um array de distâncias c(n,p,K)
    # Part é uma partição
    
    for(k in 1:K){
      
      if(sum(Part == k) == 0 || sum(Part == k) == 1 ){
        Pesos[k,] <- rep(1,p)
      }
      
      else{
        
        P <- which(Part == k)
        mtk <- MK[,,k][P,]
        Soma <- apply( mtk, 2,sum)
        
        Pesos[k,] <-( (prod(Soma) ) ^ (1/p) )/ Soma
      }
      
    }
    
    return(Pesos)
  }
  
  
  #################################################################################
  
  X <- as.matrix(X)
  n <- nrow(X); p <- (ncol(X)/2); K <- ngrupos
  
  h <- c(1:ncol(X))
  A <- X[,h[h%%2==1]]
  B <- X[,h[h%%2==0]]
  
  
  # Selecionando uma partição aleatória
  P <- sample(K,n,replace = TRUE)
  
  # Matriz de Kernel (não mudará ao longo do algoritmo)
  KM <- matriz.kernel(A,B,A,B,sig2) 
  
  #Inicializando a Matriz de Pesos
  wgh <- matrix(1,K,p)
  
  # Calculando a Distancia Para cada váriavel
  d2k <- dist2(KM,P,n,p,K)
  
  # Calculando a distância final (Somando para a todas as variáveis)
  
  d2 <- matrix(NA,n,K)
  for(k in 1:K)
    d2[,k] <- apply(d2k[,,k],1 ,sum)
  
  # Obtendo a partição
  P <- apply(d2, 1, which.min)
  
  J <- 0
  for(k in 1:K)
    J <- J + sum( d2[P == k,k] )
  
  inte <- 0
  repeat{
    
    inte <- inte + 1
    # Calculando a distância para cada grupo e para cada variável
    d2k <- dist2(KM,P,n,p,K)
    
    # Encontrando os pesos das variáveis
    wgh <- update.weights(d2k,wgh,P,K,p)
    
    # Calculando a distância final (Somando para a todas as variáveis e já ponderada pelos pesos)
    
    d2 <- matrix(NA,n,K)
    for(k in 1:K)
      for(i in 1:n)
        d2[i,k] <- wgh[k,] %*% d2k[i,,k]
    
    # Selecionando a partição
    Pn <- apply(d2, 1, which.min)
    
    
    if(all(Pn == P) || inte >= maxinter) break
    else{
      P <- Pn
      
      Jn <- 0
      for(k in 1:K)
        Jn <- Jn + sum(d2[P == k,k])
    }
    
    J <- c(J,Jn)
    
    
  }
  
  name.var <- colnames(X)
  name.var <- name.var[h[h%%2==1]]
  colnames(wgh) <- name.var
  
  
  Jsort = sort(J,decreasing = TRUE)
  
  if( all(Jsort == J))
    return( list(Partição = P, Iterações = inte, Critério = J) )
  else
    stop("Error: the objective function did not decrease")
  
}




KKMEC_LP2 <- function(X,ngrupos,sig2,maxinter = 50){
  
  RBF2 <- function(ai,bi,al,bl,sig2){
    
    ai <- as.vector(ai)
    bi <- as.vector(bi)
    al <- as.vector(al)
    bl <- as.vector(bl)
    
    norma.inf <- sum(  (ai-al)^2 )
    norma.sup <- sum(  (bi-bl)^2 )
    const <- 1/(2*sig2)
    
    return( exp(-const * norma.inf) +  exp(-const * norma.sup) )
  }
  
  
  matriz.kernel <- function(Ai,Bi,Ga,Gb,sig2){
    
    n <- nrow(Ai); K <- nrow(Ga); p <- ncol(Ai)
    
    KMs <- array(NA,c(n,K,p))
    for(j in 1:p)
      for(k in 1:K)
        for(i in 1:n)
          KMs[i,k,j] <- RBF2(Ai[i,j],Bi[i,j],Ga[k,j],Gb[k,j],sig2)   
    
    return(KMs)
  }
  
  
  # Cálculo da distância
  
  dist2<- function(MJ,Part,n,p,K){
    
    # MJ é um array com p matrizes n x n;
    
    # Part é um vetor de partição;
    
    # n,p e K descrevem o número de observações, variáveis e Grupos, respectivamente;
    
    dist.square <- array(NA,c(n,p,K))
    
    for(k in 1:K){
      
      P <- which( Part == k)
      card <- length(P)
      
      for(j in 1:p){
        
        Const <- sum( MJ[P,P,j] )
        for(i in 1:n)
          dist.square[i,j,k] <-  2 -  (2 * sum( MJ[i,P,j]) )/ card +  (Const)/ (card^2) 
        
      }
      
    }
    
    return(dist.square)
  }
  
  
  update.weights <- function(MK,Pesos,Part,K,p){
    
    # MK é um array de distâncias c(n,p,K)
    # Part é uma partição
    
    for(k in 1:K){
      
      if(sum(Part == k) == 0 || sum(Part == k) == 1 ){
        Pesos[k,] <- rep(1,p)
      }
      
      else{
        
        P <- which(Part == k)
        mtk <- MK[,,k][P,]
        Soma <- apply( mtk, 2,sum)
        
        Pesos[k,] <-( (prod(Soma) ) ^ (1/p) )/ Soma
      }
      
    }
    
    return(Pesos)
  }
  
  
  #####################################################################
  
  X <- as.matrix(X)
  n <- nrow(X); p <- (ncol(X)/2); K <- ngrupos
  
  h <- c(1:ncol(X))
  A <- X[,h[h%%2==1]]
  B <- X[,h[h%%2==0]]
  
  # Selecionando uma partição aleatória
  P <- sample(K,n,replace = TRUE)
  
  # Matriz de Kernel (não mudará ao longo do algoritmo)
  KM <- matriz.kernel(A,B,A,B,sig2) 
  
  #Inicializando a Matriz de Pesos
  wgh <- matrix(1,K,p)
  
  # Calculando a Distancia Para cada váriavel
  d2k <- dist2(KM,P,n,p,K)
  
  # Calculando a distância final (Somando para a todas as variáveis)
  
  d2 <- matrix(NA,n,K)
  for(k in 1:K)
    d2[,k] <- apply(d2k[,,k],1 ,sum)
  
  # Obtendo a partição
  P <- apply(d2, 1, which.min)
  
  J <- 0
  for(k in 1:K)
    J <- J + sum( d2[P == k,k] )
  
  inte <- 0
  repeat{
    
    inte <- inte + 1
    # Calculando a distância para cada grupo e para cada variável
    d2k <- dist2(KM,P,n,p,K)
    
    # Encontrando os pesos das variáveis
    wgh <- update.weights(d2k,wgh,P,K,p)
    
    # Calculando a distância final (Somando para a todas as variáveis e já ponderada pelos pesos)
    
    d2 <- matrix(NA,n,K)
    for(k in 1:K)
      for(i in 1:n)
        d2[i,k] <- wgh[k,] %*% d2k[i,,k]
    
    # Selecionando a partição
    Pn <- apply(d2, 1, which.min)
    
    
    if(all(Pn == P) || inte >= maxinter) break
    else{
      P <- Pn
      
      Jn <- 0
      for(k in 1:K)
        Jn <- Jn + sum(d2[P == k,k])
    }
    
    J <- c(J,Jn)
    
    
  }
  
  name.var <- colnames(X)
  name.var <- name.var[h[h%%2==1]]
  colnames(wgh) <- name.var
  
  
  Jsort = sort(J,decreasing = TRUE)
  
  if( all(Jsort == J))
    return( list(Partição = P, Iterações = inte, Critério = J) )
  else
    stop("Error: the objective function did not decrease")
  
}


KKMEC_LS1 <- function(X,ngrupos,beta,sig2,maxinter = 50){
  
  RBF1 <- function(ai,bi,al,bl,sig2){
    
    ai <- as.vector(ai)
    bi <- as.vector(bi)
    al <- as.vector(al)
    bl <- as.vector(bl)
    
    norma <- sum(  (ai-al)^2 + (bi-bl)^2 )
    const <- 1/(2*sig2)
    return( exp(-const * norma) )
  }
  
  matriz.kernel <- function(Ai,Bi,Ga,Gb,sig2){
    
    n <- nrow(Ai); K <- nrow(Ga); p <- ncol(Ai)
    
    KMs <- array(NA,c(n,K,p))
    for(j in 1:p)
      for(k in 1:K)
        for(i in 1:n)
          KMs[i,k,j] <- RBF1(Ai[i,j],Bi[i,j],Ga[k,j],Gb[k,j],sig2)   
    
    return(KMs)
  }
  
  
  # Cálculo da distância
  
  dist2<- function(MJ,Part,n,p,K){
    
    # MJ é um array com p matrizes n x n;
    
    # Part é um vetor de partição;
    
    # n,p e K descrevem o número de observações, variáveis e Grupos, respectivamente;
    
    dist.square <- array(NA,c(n,p,K))
    
    for(k in 1:K){
      
      P <- which( Part == k)
      card <- length(P)
      
      for(j in 1:p){
        
        Const <- sum( MJ[P,P,j] )
        for(i in 1:n)
          dist.square[i,j,k] <-  1 -  (2 * sum( MJ[i,P,j]) )/ card +  (Const)/ (card^2) 
        
      }
      
    }
    
    return(dist.square)
  }
  
  
  update.weights <- function(MK,Pesos,beta,Part,K,p){
    
    # MK é um array de K matrizes n x p, referentes a distância
    
    for(k in 1:K){
      
      if(sum(Part == k) == 0 || sum(Part == k) == 1)
        Pesos[k,] <- rep((1/p),p)
      
      else{
        
        P <- which(Part == k)
        mtk <- MK[,,k][P,]
        Soma <- apply( mtk, 2,sum)
        
        for(j in 1:p){
          
          ratio.sum <- (Soma[j]/Soma)^( 1/(beta-1) )
          Pesos[k,j] <- (sum(ratio.sum))^ (-1)
          
        }
      }
      
    }
    return(Pesos)
  }
  
  
  
  X <- as.matrix(X)
  n <- nrow(X); p <- (ncol(X)/2); K <- ngrupos
  
  h <- c(1:ncol(X))
  A <- X[,h[h%%2==1]]
  B <- X[,h[h%%2==0]]
  
  # Selecionando uma partição aleatória
  P <- sample(K,n,replace = TRUE)
  
  # Matriz de Kernel (não mudará ao longo do algoritmo)
  KM <- matriz.kernel(A,B,A,B,sig2) 
  
  #Inicializando a Matriz de Pesos
  wgh <- matrix(1/p,K,p)
  
  # Calculando a Distancia Para cada váriavel
  d2k <- dist2(KM,P,n,p,K)
  
  # Calculando a distância final (Somando para a todas as variáveis)
  
  d2 <- matrix(NA,n,K)
  for(k in 1:K)
    for(i in 1:n)
      d2[i,k] <- (wgh[k,] ^ beta) %*% d2k[i,,k]
  
  # Obtendo a partição
  P <- apply(d2, 1, which.min)
  
  J <- 0
  for(k in 1:K)
    J <- J + sum( d2[P == k,k] )
  
  inte <- 0
  repeat{
    
    inte <- inte + 1
    # Calculando a distância para cada grupo e para cada variável
    d2k <- dist2(KM,P,n,p,K)
    
    # Encontrando os pesos das variáveis
    wgh <- update.weights(d2k,wgh,beta,P,K,p)
    
    # Calculando a distância final (Somando para a todas as variáveis e já ponderada pelos pesos)
    
    d2 <- matrix(NA,n,K)
    for(k in 1:K)
      for(i in 1:n)
        d2[i,k] <- (wgh[k,] ^ beta) %*% d2k[i,,k]
    
    # Selecionando a partição
    Pn <- apply(d2, 1, which.min)
    
    
    if(all(Pn == P) || inte >= maxinter) break
    else{
      P <- Pn
      
      Jn <- 0
      for(k in 1:K)
        Jn <- Jn + sum(d2[P == k,k])
    }
    
    J <- c(J,Jn)
    
    
  }
  
  name.var <- colnames(X)
  name.var <- name.var[h[h%%2==1]]
  colnames(wgh) <- name.var
  
  
  Jsort = sort(J,decreasing = TRUE)
  
  if( all(Jsort == J))
    return( list(Partição = P, Iterações = inte, Critério = J) )
  else
    stop("Error: the objective function did not decrease")
  
}





KKMEC_LS2 <- function(X,ngrupos,beta,sig2,maxinter = 50){
  
  RBF2 <- function(ai,bi,al,bl,sig2){
    
    ai <- as.vector(ai)
    bi <- as.vector(bi)
    al <- as.vector(al)
    bl <- as.vector(bl)
    
    norma.inf <- sum(  (ai-al)^2 )
    norma.sup <- sum(  (bi-bl)^2 )
    const <- 1/(2*sig2)
    
    return( exp(-const * norma.inf) +  exp(-const * norma.sup) )
  }
  
  
  matriz.kernel <- function(Ai,Bi,Ga,Gb,sig2){
    
    n <- nrow(Ai); K <- nrow(Ga); p <- ncol(Ai)
    
    KMs <- array(NA,c(n,K,p))
    for(j in 1:p)
      for(k in 1:K)
        for(i in 1:n)
          KMs[i,k,j] <- RBF2(Ai[i,j],Bi[i,j],Ga[k,j],Gb[k,j],sig2)   
    
    return(KMs)
  }
  
  
  # Cálculo da distância
  
  dist2<- function(MJ,Part,n,p,K){
    
    # MJ é um array com p matrizes n x n;
    
    # Part é um vetor de partição;
    
    # n,p e K descrevem o número de observações, variáveis e Grupos, respectivamente;
    
    dist.square <- array(NA,c(n,p,K))
    
    for(k in 1:K){
      
      P <- which( Part == k)
      card <- length(P)
      
      for(j in 1:p){
        
        Const <- sum( MJ[P,P,j] )
        for(i in 1:n)
          dist.square[i,j,k] <-  2 -  (2 * sum( MJ[i,P,j]) )/ card +  (Const)/ (card^2) 
        
      }
      
    }
    
    return(dist.square)
  }
  
  
  update.weights <- function(MK,Pesos,beta,Part,K,p){
    
    # MK é um array de K matrizes n x p, referentes a distância
    
    for(k in 1:K){
      
      if(sum(Part == k) == 0 || sum(Part == k) == 1)
        Pesos[k,] <- rep((1/p),p)
      
      else{
        
        P <- which(Part == k)
        mtk <- MK[,,k][P,]
        Soma <- apply( mtk, 2,sum)
        
        for(j in 1:p){
          
          ratio.sum <- (Soma[j]/Soma)^( 1/(beta-1) )
          Pesos[k,j] <- (sum(ratio.sum))^ (-1)
          
        }
      }
      
    }
    return(Pesos)
  }
  
  ###################################################################
  
  X <- as.matrix(X)
  n <- nrow(X); p <- (ncol(X)/2); K <- ngrupos
  
  h <- c(1:ncol(X))
  A <- X[,h[h%%2==1]]
  B <- X[,h[h%%2==0]]
  
  # Selecionando uma partição aleatória
  P <- sample(K,n,replace = TRUE)
  
  # Matriz de Kernel (não mudará ao longo do algoritmo)
  KM <- matriz.kernel(A,B,A,B,sig2) 
  
  #Inicializando a Matriz de Pesos
  wgh <- matrix(1/p,K,p)
  
  # Calculando a Distancia Para cada váriavel
  d2k <- dist2(KM,P,n,p,K)
  
  # Calculando a distância final (Somando para a todas as variáveis)
  
  d2 <- matrix(NA,n,K)
  for(k in 1:K)
    for(i in 1:n)
      d2[i,k] <- (wgh[k,] ^ beta) %*% d2k[i,,k]
  
  # Obtendo a partição
  P <- apply(d2, 1, which.min)
  
  J <- 0
  for(k in 1:K)
    J <- J + sum( d2[P == k,k] )
  
  inte <- 0
  repeat{
    
    inte <- inte + 1
    # Calculando a distância para cada grupo e para cada variável
    d2k <- dist2(KM,P,n,p,K)
    
    # Encontrando os pesos das variáveis
    wgh <- update.weights(d2k,wgh,beta,P,K,p)
    
    # Calculando a distância final (Somando para a todas as variáveis e já ponderada pelos pesos)
    
    d2 <- matrix(NA,n,K)
    for(k in 1:K)
      for(i in 1:n)
        d2[i,k] <- (wgh[k,] ^ beta) %*% d2k[i,,k]
    
    # Selecionando a partição
    Pn <- apply(d2, 1, which.min)
    
    
    if(all(Pn == P) || inte >= maxinter) break
    else{
      P <- Pn
      
      Jn <- 0
      for(k in 1:K)
        Jn <- Jn + sum(d2[P == k,k])
    }
    
    J <- c(J,Jn)
    
    
  }
  
  name.var <- colnames(X)
  name.var <- name.var[h[h%%2==1]]
  colnames(wgh) <- name.var
  
  
  Jsort = sort(J,decreasing = TRUE)
  
  if( all(Jsort == J))
    return( list(Partição = P, Iterações = inte, Critério = J) )
  else
    stop("Error: the objective function did not decrease")
  
}




KKMEE_1 <- function(X,ngrupos,sig2,maxinter = 50){
  
  RBF1 <- function(ai,bi,al,bl,sig2){
    
    ai <- as.vector(ai)
    bi <- as.vector(bi)
    al <- as.vector(al)
    bl <- as.vector(bl)
    
    norma <- sum(  (ai-al)^2 + (bi-bl)^2 )
    const <- 1/(2*sig2)
    return( exp(-const * norma) )
  }
  
  matriz.kernel <- function(Ai,Bi,Ga,Gb,sig2){
    
    n <- nrow(Ai); K <- nrow(Ga); p <- ncol(Ai)
    KMs <- array(NA,c(n,p,K))
    
    for(k in 1:K){
      
      KM <- matrix(0,n,p)
      for(i in 1:n)
        for (j in 1:p) 
          KM[i,j] <- RBF1(Ai[i,j],Bi[i,j],Ga[k,j],Gb[k,j],sig2)
        
        KMs[,,k] <- KM
        
    }
    return(KMs)
  }
  
  dist.square <- function(KM,n,p,K){
    
    dist.indv <- array(NA,c(n,p,K))
    
    for(k in 1:K)
      for(j in 1:p)
        dist.indv[,j,k] <- 2*(1 - KM[,j,k])
      
      
      dist.final <- matrix(NA,n,K)
      
      for (k2 in 1:K) 
        dist.final[,k2] <- apply(dist.indv[,,k2], 1,sum)
      
      return(dist.final)
  }
  
  
  
  #########################################################
  
  X <- as.matrix(X)
  n <- nrow(X); p <- (ncol(X)/2); K <- ngrupos
  
  xx <- unique(X); nn <- nrow(xx)
  is <- sample.int(nn,K)
  
  h <- c(1:ncol(X))
  A <- X[,h[h%%2==1]]; aa <- xx[,h[h%%2==1]]
  B <- X[,h[h%%2==0]]; bb <- xx[,h[h%%2==0]]
  
  ga <- aa[is, ]
  gb <- bb[is, ]
  
  KM <- matriz.kernel(A,B,ga,gb,sig2) 
  
  d2 <- dist.square(KM,n,p,K)
  P <- apply(d2,1, which.min)
  
  J <- 0
  for(k in 1:K)
    J <- J + sum( d2[P == k,k] )
  
  
  inte <- 0
  repeat{
    
    inte <- inte + 1
    ginf <- ga
    gsup <- gb
    
    for(k in 1:K){
      
      if( sum(P == k) == 0)
        next
      
      if( sum(P == k) == 1){
        ginf[k,] <- A[P == k,]
        gsup[k,] <- B[P == k,]
      }  
      
      
      if(sum(P == k) > 1 ){
        
        for(j in 1:p) {
          
          ginf[k,j] <- weighted.mean(A[P == k,j], KM[P == k,j,k])
          gsup[k,j] <- weighted.mean(B[P == k,j], KM[P == k,j,k])
        }
      }
    }
    
    ga <- ginf
    gb <- gsup
    
    # Calculando a matriz de kernel
    KM <- matriz.kernel(A,B,ga,gb,sig2)
    
    
    #calculando as distancias 
    d2 <- dist.square(KM,n,p,K)
    
    # Encontrando a melhor partição
    Pn <- apply(d2, 1, which.min)
    
    
    if(all(Pn == P) || inte >= maxinter) break
    else{
      P <- Pn
      
      Jn <- 0
      for(k in 1:K)
        Jn <- Jn + sum(d2[P == k,k])
    }
    
    J <- c(J,Jn)
    
    
  }
  
  
  Jsort = sort(J,decreasing = TRUE)
  
  if( all(Jsort == J))
    return( list(Partição = P, Iterações = inte, Critério = J) )
  else
    stop("Error: the objective function did not decrease")
  
}




KKMEE_GP1 <- function(X,ngrupos,sig2,maxinter = 50){
  
  RBF1 <- function(ai,bi,al,bl,sig2){
    
    ai <- as.vector(ai)
    bi <- as.vector(bi)
    al <- as.vector(al)
    bl <- as.vector(bl)
    
    norma <- sum(  (ai-al)^2 + (bi-bl)^2 )
    const <- 1/(2*sig2)
    return( exp(-const * norma) )
  }
  
  matriz.kernel <- function(Ai,Bi,Ga,Gb,sig2){
    
    n <- nrow(Ai); K <- nrow(Ga); p <- ncol(Ai)
    KMs <- array(NA,c(n,p,K))
    
    for(k in 1:K){
      
      KM <- matrix(0,n,p)
      for(i in 1:n)
        for (j in 1:p) 
          KM[i,j] <- RBF1(Ai[i,j],Bi[i,j],Ga[k,j],Gb[k,j],sig2)
        
        KMs[,,k] <- KM
        
    }
    return(KMs)
  }
  
  
  
  update.weights <- function(MK,Part,K,p){
    
    Soma <- matrix(NA,K,p)
    for(k in 1:K){
      
      P <- which(Part == k)
      mtk <- MK[,,k][P,]
      
      if(length(P) == 1)
        Soma[k,] <- (1 - mtk)
      else
        Soma[k,] <- apply( (1 - mtk), 2,sum)
    }
    
    Sum.final <- apply(Soma,2,sum)
    Const <- prod(Sum.final) ^ (1/p)
    Pesos <- Const/Sum.final
    return(Pesos)
  }
  
  d2.pond <- function(KM,lamb,n,p,K){
    
    dist.lamb <- array(NA,c(n,p,K))
    
    for(j in 1:p)
      dist.lamb[,j,] <- 2 * lamb[j] *(1 - KM[,j,])
    
    
    dist.final <- matrix(NA,n,K)
    
    for (k2 in 1:K) 
      dist.final[,k2] <- apply(dist.lamb[,,k2], 1,sum)
    
    return(dist.final)
  }
  
  
  ##################################################################
  
  X <- as.matrix(X)
  n <- nrow(X); p <- (ncol(X)/2); K <- ngrupos
  
  xx <- unique(X); nn <- nrow(xx)
  is <- sample.int(nn,K)
  
  h <- c(1:ncol(X))
  A <- X[,h[h%%2==1]]; aa <- xx[,h[h%%2==1]]
  B <- X[,h[h%%2==0]]; bb <- xx[,h[h%%2==0]]
  
  ga <- aa[is, ]
  gb <- bb[is, ]
  
  KM <- matriz.kernel(A,B,ga,gb,sig2) 
  wgh <- matrix(1,K,p)
  d2 <- d2.pond(KM,wgh,n,p,K)
  P <- apply(d2,1, which.min)
  
  J <- 0
  for(k in 1:K)
    J <- J + sum( d2[P == k,k] )
  
  
  inte <- 0
  repeat{
    
    inte <- inte + 1
    ginf <- ga
    gsup <- gb
    
    for(k in 1:K){
      
      if( sum(P == k) == 0)
        next
      
      if( sum(P == k) == 1){
        ginf[k,] <- A[P == k,]
        gsup[k,] <- B[P == k,]
      }  
      
      
      if(sum(P == k) > 1 ){
        
        for(j in 1:p) {
          
          ginf[k,j] <- weighted.mean(A[P == k,j], KM[P == k,j,k])
          gsup[k,j] <- weighted.mean(B[P == k,j], KM[P == k,j,k])
        }
      }
    }
    
    ga <- ginf
    gb <- gsup
    
    
    # Calculando a matrix de kernel
    KM <- matriz.kernel(A,B,ga,gb,sig2)
    
    # cálculo dos pesos
    weights <- update.weights(KM,P,K,p)
    
    #calculando as distancias ponderada pelos pesos
    d2 <- d2.pond(KM,weights,n,p,K)
    
    # Encontrando a melhor partição
    Pn <- apply(d2, 1, which.min)
    
    
    if(all(Pn == P) || inte >= maxinter) break
    else{
      P <- Pn
      
      Jn <- 0
      for(k in 1:K)
        Jn <- Jn + sum(d2[P == k,k])
    }
    
    J <- c(J,Jn)
  
  }
  
  name.var <- colnames(X)
  name.var <- name.var[h[h%%2==1]]
  names(weights) <- name.var
  
  
  Jsort = sort(J,decreasing = TRUE)
  
  if( all(Jsort == J))
    return( list(Partição = P,Pesos = weights, Iterações = inte, Critério = J) )
  else
    stop("Error: the objective function did not decrease")
  
}




KKMEE_GS1 <- function(X,ngrupos,beta,sig2,maxinter = 50){
  
  ################################## FUNÇÕES AUXILIARES #######################################
  
  RBF1 <- function(ai,bi,al,bl,sig2){
    
    ai <- as.vector(ai)
    bi <- as.vector(bi)
    al <- as.vector(al)
    bl <- as.vector(bl)
    
    norma <- sum(  (ai-al)^2 + (bi-bl)^2 )
    const <- 1/(2*sig2)
    return( exp(-const * norma) )
  }
  
  matriz.kernel <- function(Ai,Bi,Ga,Gb,sig2){
    
    n <- nrow(Ai); K <- nrow(Ga); p <- ncol(Ai)
    KMs <- array(NA,c(n,p,K))
    
    for(k in 1:K){
      
      KM <- matrix(0,n,p)
      for(i in 1:n)
        for (j in 1:p) 
          KM[i,j] <- RBF1(Ai[i,j],Bi[i,j],Ga[k,j],Gb[k,j],sig2)
        
        KMs[,,k] <- KM
        
    }
    return(KMs)
  }
  
  
  
  update.weights <- function(MK,beta,Part,K,p){
    
    # FUNÇÃO AUXILIAR (SUBSTITUIR ZEROS)
    
    replace.zero <- function(vector,increase){
      
      P <- which(vector == 0)
      
      if( length(P) == 0)
        return(vector)
      else
        vector[P] <- increase; return(vector)
      
    }
    
    ############################# FUNÇÃO PRINCIPAL ################################
    
    Soma <- matrix(NA,K,p); Pesos <- NULL
    for(k in 1:K){
      
      P <- which(Part == k)
      mtk <- MK[,,k][P,]
      
      if(length(P) == 1)
        Soma[k,] <- (1 - mtk)
      else
        Soma[k,] <- apply((1 - mtk), 2, sum)
    }   
    
    sum.final <- apply(Soma,2, sum)
    sum.final <- replace.zero(sum.final, (10^-3))
    
    for(j in 1:p){
      
      ratio.sum <- (sum.final[j]/sum.final)^( 1/(beta-1) )
      Pesos[j] <- (sum(ratio.sum))^(-1)
      
    }
    
    return(Pesos)
  }
  
  
  
  d2.pond <- function(KM,lamb,beta,n,p,K){
    
    dist.lamb <- array(NA,c(n,p,K))
    
    for(j in 1:p)
      dist.lamb[,j,] <- 2 * (lamb[j])^ (beta) *(1 - KM[,j,])
    
    
    dist.final <- matrix(NA,n,K)
    
    for (k2 in 1:K) 
      dist.final[,k2] <- apply(dist.lamb[,,k2], 1,sum)
    
    return(dist.final)
  }
  
  
  ################################### FUNÇÃO PRINCIPAL ########################################
  
  X <- as.matrix(X)
  n <- nrow(X); p <- (ncol(X)/2); K <- ngrupos
  
  xx <- unique(X); nn <- nrow(xx)
  is <- sample.int(nn,K)
  
  h <- c(1:ncol(X))
  A <- X[,h[h%%2==1]]; aa <- xx[,h[h%%2==1]]
  B <- X[,h[h%%2==0]]; bb <- xx[,h[h%%2==0]]
  
  ga <- aa[is, ]
  gb <- bb[is, ]
  
  KM <- matriz.kernel(A,B,ga,gb,sig2) 
  wgh <- matrix( (1/p),K,p)
  d2 <- d2.pond(KM,wgh,beta,n,p,K)
  P <- apply(d2,1, which.min)
  
  J <- 0
  for(k in 1:K)
    J <- J + sum( d2[P == k,k] )
  
  
  inte <- 0
  repeat{
    
    inte <- inte + 1
    ginf <- ga
    gsup <- gb
    
    for(k in 1:K){
      
      if( sum(P == k) == 0)
        next
      
      if( sum(P == k) == 1){
        ginf[k,] <- A[P == k,]
        gsup[k,] <- B[P == k,]
      }  
      
      
      if(sum(P == k) > 1 ){
        
        for(j in 1:p) {
          
          ginf[k,j] <- weighted.mean(A[P == k,j], KM[P == k,j,k])
          gsup[k,j] <- weighted.mean(B[P == k,j], KM[P == k,j,k])
        }
      }
    }
    
    ga <- ginf
    gb <- gsup
    
    KM <- matriz.kernel(A,B,ga,gb,sig2)
    
    # cálculo dos pesos
    
    weights <- update.weights(KM,beta,P,K,p)
    
    
    #calculando as distancias ponderada pelos pesos
    d2 <- d2.pond(KM,weights,beta,n,p,K)
    
    # Encontrando a melhor partição
    Pn <- apply(d2, 1, which.min)
    
    
    if(all(Pn == P) || inte >= maxinter) break
    else{
      P <- Pn
      
      Jn <- 0
      for(k in 1:K)
        Jn <- Jn + sum(d2[P == k,k])
    }
    
    J <- c(J,Jn)
    
    
  }
  
  name.var <- colnames(X)
  name.var <- name.var[h[h%%2==1]]
  names(weights) <- name.var
  
  
  Jsort = sort(J,decreasing = TRUE)
  
  if( all(Jsort == J))
    return( list(Partição = P,Pesos=weights, Iterações = inte, Critério = J) )
  else
    stop("Error: the objective function did not decrease")
  
}



KKMEE_LP1 <- function(X,ngrupos,sig2,maxinter = 50){
  
  ################################## FUNÇÕES AUXILIARES #######################################
  
  RBF1 <- function(ai,bi,al,bl,sig2){
    
    ai <- as.vector(ai)
    bi <- as.vector(bi)
    al <- as.vector(al)
    bl <- as.vector(bl)
    
    norma <- sum(  (ai-al)^2 + (bi-bl)^2 )
    const <- 1/(2*sig2)
    return( exp(-const * norma) )
  }
  
  matriz.kernel <- function(Ai,Bi,Ga,Gb,sig2){
    
    n <- nrow(Ai); K <- nrow(Ga); p <- ncol(Ai)
    KMs <- array(NA,c(n,p,K))
    
    for(k in 1:K)
      for(i in 1:n)
        for (j in 1:p) 
          KMs[i,j,k] <- RBF1(Ai[i,j],Bi[i,j],Ga[k,j],Gb[k,j],sig2)
    
    return(KMs)
  }
  
  
  
  update.weights <- function(MK,Pesos,Part,K,p){
    
    # FUNÇÃO AUXILIAR (SUBSTITUIR ZEROS)
    
    replace.zero <- function(vector,increase){
      
      P <- which(vector == 0)
      
      if( length(P) == 0)
        return(vector)
      else
        vector[P] <- increase; return(vector)
      
    }
    
    ############################# FUNÇÃO PRINCIPAL ################################ 
    
    for(k in 1:K){
      
      if(sum(Part == k) == 0 || sum(Part == k) == 1 ){
        Pesos[k,] <- rep(1,p)
      }
      
      else{
        
        P <- which(Part == k)
        mtk <- MK[,,k][P,]
        Soma <- apply( (1 - mtk), 2,sum)
        Soma <- replace.zero(Soma,(10^-3))
        Pesos[k,] <-( (prod(Soma) ) ^ (1/p) )/ Soma
      }
      
    }
    
    return(Pesos)
  }
  
  
  
  d2.pond <- function(KM,lamb,n,p,K){
    
    dist.lamb <- array(NA,c(n,p,K))
    
    for(k in 1:K)
      for(j in 1:p)
        dist.lamb[,j,k] <- 2 * lamb[k,j] *(1 - KM[,j,k])
      
      
      dist.final <- matrix(NA,n,K)
      
      for (k2 in 1:K) 
        dist.final[,k2] <- apply(dist.lamb[,,k2], 1,sum)
      
      return(dist.final)
  }
  
  
  ################################### FUNÇÃO PRINCIPAL ########################################
  
  X <- as.matrix(X)
  n <- nrow(X); p <- (ncol(X)/2); K <- ngrupos
  
  xx <- unique(X); nn <- nrow(xx)
  is <- sample.int(nn,K)
  
  h <- c(1:ncol(X))
  A <- X[,h[h%%2==1]]; aa <- xx[,h[h%%2==1]]
  B <- X[,h[h%%2==0]]; bb <- xx[,h[h%%2==0]]
  
  ga <- aa[is, ]
  gb <- bb[is, ]
  
  KM <- matriz.kernel(A,B,ga,gb,sig2) 
  wgh <- matrix(1,K,p)
  d2 <- d2.pond(KM,wgh,n,p,K)
  P <- apply(d2,1, which.min)
  
  J <- 0
  for(k in 1:K)
    J <- J + sum( d2[P == k,k] )
  
  
  inte <- 0
  repeat{
    
    inte <- inte + 1
    ginf <- ga
    gsup <- gb
    
    for(k in 1:K){
      
      if( sum(P == k) == 0)
        next
      
      if( sum(P == k) == 1){
        ginf[k,] <- A[P == k,]
        gsup[k,] <- B[P == k,]
      }  
      
      
      if(sum(P == k) > 1 ){
        
        for(j in 1:p) {
          
          ginf[k,j] <- weighted.mean(A[P == k,j], KM[P == k,j,k])
          gsup[k,j] <- weighted.mean(B[P == k,j], KM[P == k,j,k])
        }
      }
    }
    
    ga <- ginf
    gb <- gsup
    
    KM <- matriz.kernel(A,B,ga,gb,sig2)
    
    # cálculo dos pesos
    
    weights <- update.weights(KM,wgh,P,K,p)
    
    #calculando as distancias ponderada pelos pesos
    d2 <- d2.pond(KM,weights,n,p,K)
    
    # Encontrando a melhor partição
    Pn <- apply(d2, 1, which.min)
    
    if(all(Pn == P) || inte >= maxinter) break
    else{
      P <- Pn
      
      Jn <- 0
      for(k in 1:K)
        Jn <- Jn + sum(d2[P == k,k])
    }
    
    J <- c(J,Jn)
    
    
  }
  
  name.var <- colnames(X)
  name.var <- name.var[h[h%%2==1]]
  colnames(weights) <- name.var
  
  
  Jsort = sort(J,decreasing = TRUE)
  
  if( all(Jsort == J))
    return( list(Partição = P, Pesos = weights, Iterações = inte, Critério = J) )
  else
    stop("Error: the objective function did not decrease")
  
}



KKMEE_LS1 <- function(X,ngrupos,beta,sig2,maxinter = 50){
  
  ################################## FUNÇÕES AUXILIARES #######################################
  
  
  RBF1 <- function(ai,bi,al,bl,sig2){
    
    ai <- as.vector(ai)
    bi <- as.vector(bi)
    al <- as.vector(al)
    bl <- as.vector(bl)
    
    norma <- sum(  (ai-al)^2 + (bi-bl)^2 )
    const <- 1/(2*sig2)
    return( exp(-const * norma) )
  }
  
  matriz.kernel <- function(Ai,Bi,Ga,Gb,sig2){
    
    n <- nrow(Ai); K <- nrow(Ga); p <- ncol(Ai)
    KMs <- array(NA,c(n,p,K))
    
    for(k in 1:K){
      
      KM <- matrix(0,n,p)
      for(i in 1:n)
        for (j in 1:p) 
          KM[i,j] <- RBF1(Ai[i,j],Bi[i,j],Ga[k,j],Gb[k,j],sig2)
        
        KMs[,,k] <- KM
        
    }
    return(KMs)
  }
  
  
  
  update.weights <- function(MK,Pesos,beta,Part,K,p){
    
    # FUNÇÃO AUXILIAR (SUBSTITUIR ZEROS)
    
    replace.zero <- function(vector,increase){
      
      P <- which(vector == 0)
      
      if( length(P) == 0)
        return(vector)
      else
        vector[P] <- increase; return(vector)
      
    }
    
    ############################# FUNÇÃO PRINCIPAL ################################ 
    
    for(k in 1:K){
      
      if(sum(Part == k) == 0 || sum(Part == k) == 1)
        Pesos[k,] <- rep((1/p),p)
      
      else{
        
        P <- which(Part == k)
        mtk <- MK[,,k][P,]
        Soma <- apply( (1 - mtk), 2,sum)
        Soma <- replace.zero(Soma,(10^-3))
        
        for(j in 1:p){
          
          ratio.sum <- (Soma[j]/Soma)^( 1/(beta-1) )
          Pesos[k,j] <- (sum(ratio.sum))^ (-1)
          
        }
      }
      
    }
    return(Pesos)
  }
  
  
  
  d2.pond <- function(KM,lamb,beta,n,p,K){
    
    dist.lamb <- array(NA,c(n,p,K))
    
    for(k in 1:K)
      for(j in 1:p)
        dist.lamb[,j,k] <- 2 * (lamb[k,j])^ (beta) *(1 - KM[,j,k])
      
      
      dist.final <- matrix(NA,n,K)
      
      for (k2 in 1:K) 
        dist.final[,k2] <- apply(dist.lamb[,,k2], 1,sum)
      
      return(dist.final)
  }
  
  
  ################################### FUNÇÃO PRINCIPAL ########################################
  
  X <- as.matrix(X)
  n <- nrow(X); p <- (ncol(X)/2); K <- ngrupos
  
  xx <- unique(X); nn <- nrow(xx)
  is <- sample.int(nn,K)
  
  h <- c(1:ncol(X))
  A <- X[,h[h%%2==1]]; aa <- xx[,h[h%%2==1]]
  B <- X[,h[h%%2==0]]; bb <- xx[,h[h%%2==0]]
  
  ga <- aa[is, ]
  gb <- bb[is, ]
  
  KM <- matriz.kernel(A,B,ga,gb,sig2) 
  wgh <- matrix( (1/p),K,p)
  d2 <- d2.pond(KM,wgh,beta,n,p,K)
  P <- apply(d2,1, which.min)
  
  J <- 0
  for(k in 1:K)
    J <- J + sum( d2[P == k,k] )
  
  
  inte <- 0
  repeat{
    
    inte <- inte + 1
    ginf <- ga
    gsup <- gb
    
    for(k in 1:K){
      
      
      if( sum(P == k) == 0)
        next
      
      
      if( sum(P == k) == 1){
        ginf[k,] <- A[P == k,]
        gsup[k,] <- B[P == k,]
      }  
      
      
      if(sum(P == k) > 1 ){
        
        for(j in 1:p) {
          
          ginf[k,j] <- weighted.mean(A[P == k,j], KM[P == k,j,k])
          gsup[k,j] <- weighted.mean(B[P == k,j], KM[P == k,j,k])
        }
      }
    }
    
    ga <- ginf
    gb <- gsup
    
    KM <- matriz.kernel(A,B,ga,gb,sig2)
    
    # cálculo dos pesos
    
    weights <- update.weights(KM,wgh,beta,P,K,p)
    
    
    #calculando as distancias ponderada pelos pesos
    d2 <- d2.pond(KM,weights,beta,n,p,K)
    
    # Encontrando a melhor partição
    Pn <- apply(d2, 1, which.min)
    
    
    if(all(Pn == P) || inte >= maxinter) break
    else{
      P <- Pn
      
      Jn <- 0
      for(k in 1:K)
        Jn <- Jn + sum(d2[P == k,k])
    }
    
    J <- c(J,Jn)
    
    
  }
  
  name.var <- colnames(X)
  name.var <- name.var[h[h%%2==1]]
  colnames(weights) <- name.var
  
  
  Jsort = sort(J,decreasing = TRUE)
  
  if( all(Jsort == J))
    return( list(Partição = P, Pesos = weights, Iterações = inte, Critério = J) )
  else
    stop("Error: the objective function did not decrease")
  
}
