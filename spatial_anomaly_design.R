library(acebayes)
library(FNN)

# prop=0.1
# scenario=3
# seed=1285

global_i <- 0
global_metrics<- list()
global_accuracy <- c()
global_anom_proportion <- c()
global_groundtruth_anom <- c()
global_y_anoms <- c()
global_ground_truth <- c()
global_y <- c()

d_phi <- 2
e_phi <- 3
c_tau2 <- 0.00001
a_sig2 <- 3
b_sig2 <- 1
k <- 2
r <- 10
n0 <- r ^ k
x0 <- seq(from = 0, to = 1, length.out = r)
d0 <- as.matrix(expand.grid(x0, x0))
delta <- 0.25
n <- 6


rho <- function(X1, X2, phi) {
  k <- ncol(X1)
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  A <- matrix(0, nrow = n1, ncol = n2)
  for(i in 1:k) {
    A <- A - phi * (matrix(rep(X1[, i], n2), nrow =n1) -
                      matrix(rep(X2[, i], each = n1), nrow = n1)) ^ 2
  }
  exp(A)
}

# generate anomalies
induce_anomalies <- function(n, B, prop = 0.001){
  set.seed(3)
  anom_idx <- sample(0:1, n*B, replace = T, prob = c(1-prop,prop))
  dim(anom_idx) <- c(n, B)
  anom_idx
}

accuracy <- function(matrix1, matrix2) {
  # Ensure both matrices have the same dimensions
  if (length(matrix1) != length(matrix2)) {
    stop("Matrices must have the same dimensions.")
  }
  
  # Calculate the accuracy
  num_correct <- sum(matrix1 == matrix2)
  total_elements <- length(matrix1)
  accuracy <- num_correct / total_elements
  
  return(accuracy)
}

confusion_matrix <- function(actual, predicted){
  cm <- matrix(0, nrow = 2, ncol = 2)
  for (i in 1:length(actual)) {
    true_class <- actual[i]
    pred_class <- predicted[i]
    cm[true_class + 1, pred_class + 1] <- cm[true_class + 1, pred_class + 1] + 1
  }
  
  return(cm)
}

calculate_metrics <- function(cm){
  true_positive <- cm[2, 2]
  true_negative <- cm[1, 1]
  false_positive <- cm[1, 2]
  false_negative <- cm[2, 1]
  total <- sum(cm)
  
  specificity <- true_negative / (true_negative + false_positive)
  accuracy <- (true_positive + true_negative) / total
  mcc <- ((cm[2, 2] * cm[1, 1]) -
                     (cm[1, 2] * cm[2, 1])) /
    sqrt((cm[2, 2] + cm[1, 2]) *
           as.numeric((cm[2, 2] + cm[2, 1]) *
                        (cm[1, 1] + cm[1, 2]) *
                        (cm[1, 1] + cm[2, 1])))
  metrics <- c(specificity, accuracy, mcc)
  return(metrics)
}

detect_anomalies <- function(y, d, y_train, threshold = 3, k = 3) {
  n <- length(y)
  anomalies <- numeric(n)
  
  # Calculate weighted means using inverse Euclidean distances
  weighted_means <- numeric(n)
  sds <- apply(y_train, 1, sd)
  
  for (i in 1:n) {
    # Find indices and distances of k-nearest neighbors
    nn_info <- get.knnx(d, d[i,, drop = FALSE], k = k)
    nn_indices <- nn_info$nn.index[-(1)]
    nn_distances <- nn_info$nn.dist[-(1)]
    nn_distances <- ifelse(nn_distances < 0.01, 0.01, nn_distances) # apply floor
    
    # Calculate weights as inverse distances
    weights <- 1 / nn_distances
    total_weight <- sum(weights)
    
    # Calculate the weighted mean for the ith point
    weighted_means[i] <- mean(y[nn_indices]*(weights / total_weight))
  }
  
  # Detect anomalies
  for (i in 1:n) {
    if (abs(y[i] - weighted_means[i]) > threshold * sds[i]) {
      anomalies[i] <- 1
    } else {
      anomalies[i] <- 0
    }
  }
  
  return(anomalies)
}

anomaly_score <- function(y, d, y_train, ground_truth){
  
  # generate anomalies on y
  y_anoms <- ifelse(ground_truth==0, y, y+rnorm(length(y), 5, 10))
  
  # run anomaly detection algorithm 
  predicted_anoms <- detect_anomalies(y_anoms, d, y_train)
  
  # get 'clean' version of data (drop anomalies)
  y_clean <- ifelse(predicted_anoms==0, y_anoms, NA)
  
  list("predicted_anoms" = predicted_anoms, "y_clean" = y_clean)
  
}

rmse_loss <- function(actual_values, predicted_values) {
  if (length(actual_values) != length(predicted_values)) {
    stop("Input vectors must have the same length.")
  }
  squared_errors <- (actual_values - predicted_values)^2
  mean_squared_error <- mean(squared_errors, na.rm = TRUE)
  
  return(sqrt(mean_squared_error))
}

utilpred <- function(d, B, n_train=100) {
  set.seed(seed+global_i)
  n <- dim(d)[1]
  
  phi <- rep(1,B)
  tau2 <- rep(1e-5,B)
  sig2 <- 1 / rgamma(n = B, shape = 0.5 * a_sig2, rate = 0.5 * b_sig2)
  
  # training data
  y_train <- matrix(NA, nrow = n, ncol = n_train)
  for(b in 1:n_train) {
    C <- rho(d, d, phi[b])
    Sigmatilde <- C + tau2[b] * diag(n)
    cSigmatilde <- chol(Sigmatilde)
    y_train[, b] <- sqrt(sig2[b]) * t(rnorm(n) %*% cSigmatilde)
  }
  
  # n x B anomaly ground truth
  ground_truth <- induce_anomalies(n, B, prop=prop)
  global_groundtruth_anom <<- sum(ground_truth)/length(ground_truth)
  
  # for generated 'real' data
  ytilde <- matrix(NA, nrow = n + n0, ncol = B)
  postpredmean <- matrix(NA, nrow = n0, ncol = B)
  y0 <- matrix(NA, nrow = n0, ncol = B)
  predicted_anom <- matrix(NA, nrow = n, ncol = B)
  accuracy <- rep(0,B)
  
  for(b in 1:B) {
    C <- rho(d, d, phi[b])
    S <- rho(d, d0, phi[b])
    C0 <- rho(d0, d0, phi[b])
    Sigmatilde <- rbind(cbind(C, S), cbind(t(S), C0)) + tau2[b] * diag(n + n0)
    cSigmatilde <- chol(Sigmatilde)
    ytilde[, b] <- sqrt(sig2[b]) * t(rnorm(n + n0) %*% cSigmatilde)
    
    # suppose that y is anomalous, get detection score, and 'clean' observations
    y <- ytilde[ (1:n), b]
    out <- anomaly_score(y, d, y_train, ground_truth[,b])
    predicted_anom[,b] <- out$predicted_anom
    y_clean <- out$y_clean
    y_clean[is.na(y_clean)] <- 0 # fill missing with 0
    postpredmean[, b] <- t(S) %*% solve(C + tau2[b] * diag(n)) %*% y_clean
    
    y0[,b] <- ytilde[ - (1:n), b]
    accuracy[b] <- 1/rmse_loss(postpredmean[, b], y0[, b])
  }
  
  # get accuracy of anomaly classification
  cm <- confusion_matrix(ground_truth, predicted_anom)
  metrics <- calculate_metrics(cm)
  
  global_i <<- global_i+1
  global_metrics[[global_i]] <<- metrics
  global_accuracy <<- c(global_accuracy, mean(accuracy, na.rm =TRUE))
  
  utility <- accuracy * metrics[1]
  utility[!is.na(utility)]
}

set.seed(seed)
start.d <- randomLHS(n = n, k = k)
ex44 <- ace(utilpred, start.d, c(1500,1000), lower = 0, upper = 1, N1=30, N2 = 0)

# save
file.name=paste('spatial_', scenario, '_', seed, 'anom_new.Rdata', sep='')
save.image(file= file.name)
