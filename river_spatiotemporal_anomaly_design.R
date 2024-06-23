library(oddstream)
library(dplyr)
library(acebayes)
library(SSN)

# prop=0.0015
# lambda=1.75
# scenario=4

seed.network=179
set.seed(seed.network)

global_i <- 0
global_metrics<- list()
global_accuracy <- c()
global_anom_proportion <- c()
global_groundtruth_anom <- c()
global_y_anoms <- c()
global_ground_truth <- c()
global_y <- c()

# network 
n.segments = 300
n.sites = 14
ssn.obj.path = paste('./nsegments',n.segments, 'seed',seed.network,'.ssn',sep='')

generate.network <- function(ssn.obj.path, n.segments){
  if(dir.exists(ssn.obj.path)) {
    # read in empty network object
    n1 <- importSSN(ssn.obj.path, predpts = 'preds')
  } else {
    # create network
    n1 <- createSSN(c(n.segments), 
                    obsDesign = systematicDesign(spacing=5),
                    predDesign = systematicDesign(spacing=0.5),
                    path=ssn.obj.path, importToR = TRUE, treeFunction = iterativeTreeLayout)
    # create distances file
    createDistMat(n1, o.write=TRUE, predpts = 'preds', amongpreds = TRUE)
  }
  return(n1)
}

get.dist.matrix <- function(ssn.obj.path, netID = 1, PredID = "preds"){
  dist.junc <- readRDS(paste(ssn.obj.path, '/distance/', PredID, '/dist.net', netID, '.RData', sep=''))
  dist.junc
}

# fix network, fix locations on it
ssn <- generate.network(ssn.obj.path, n.segments)
rawDFobs <- getSSNdata.frame(ssn, Name = "Obs")
rawDFpred <- getSSNdata.frame(ssn, Name = "preds")

## Extract the geographic coordinates from the SpatialStreamNetwork object and add to data.frames
obs_data_coord <- data.frame(ssn@obspoints@SSNPoints[[1]]@point.coords)
obs_data_coord$pid<- as.numeric(rownames(obs_data_coord))
rawDFobs<- rawDFobs %>% left_join(obs_data_coord, by = c("pid"),
                                  keep = FALSE)
rawDFobs$point <- "predict" ## Create label for observed points

pred_data_coord <- data.frame(ssn@predpoints@SSNPoints[[1]]@point.coords)
pred_data_coord$pid<- as.numeric(rownames(pred_data_coord))
rawDFpred<- rawDFpred %>% left_join(pred_data_coord, by = "pid",
                                    keep = FALSE)
rawDFpred$point <- "design" ## Create label for prediction points

# fix covariate values and add to dataframes
set.seed(seed)
rawDFpred[,"X1"] <- rnorm(length(rawDFpred[,1]))
rawDFobs[,"X1"] <- rnorm(length(rawDFobs[,1]))

## Ensure the rownames still match the pid values used in the SpatialStreamNetwork object
rownames(rawDFobs)<- as.character(rawDFobs$pid)
rownames(rawDFpred)<- as.character(rawDFpred$pid)

## Put the new covariates back in the SpatialStreamNetwork object
ssn <- putSSNdata.frame(rawDFobs, ssn, Name = 'Obs')
ssn <- putSSNdata.frame(rawDFpred, ssn , Name = 'preds')

get.full.dist.matrix <- function(ssn.obj.path, netID = 1, PredID = "preds"){
  
  # get distance matrix, assume created
  obs2obs <- readRDS(paste(ssn.obj.path, '/distance/obs/dist.net',netID,'.RData', sep=''))
  pred2pred <- readRDS(paste(ssn.obj.path, '/distance/', PredID, '/dist.net', netID, '.RData', sep=''))
  pred2obs <- readRDS(paste(ssn.obj.path, '/distance/', PredID, '/dist.net',netID,'.a.RData', sep=''))
  obs2pred <- readRDS(paste(ssn.obj.path, '/distance/', PredID, '/dist.net',netID,'.b.RData', sep=''))
  dist.junc <- cbind(rbind(obs2obs, obs2pred), rbind(pred2obs, pred2pred))
  dist.junc
}

exp.tailup <- function(dist.hydro, weight, parsil, range)
{
  parsil*exp(-3*dist.hydro/range)*weight
}

get.covariance.matrix <- function(d.pids, dist.hydro, weight, parsil, range)
{
  d.dist.hydro = dist.hydro[d.pids, d.pids]
  d.weight = weight[d.pids, d.pids]
  exp.tailup(d.dist.hydro, d.weight, parsil, range)
}

rmse_loss <- function(actual_values, predicted_values) {
  if (length(actual_values) != length(predicted_values)) {
    stop("Input vectors must have the same length.")
  }
  squared_errors <- (actual_values - predicted_values)^2
  mean_squared_error <- mean(squared_errors, na.rm = TRUE)
  
  return(sqrt(mean_squared_error))
}

# combine dataframes - entire discrete design space
all_space = rbind(rawDFobs, rawDFpred)
design_space = rawDFpred

# get distance matrix for all space
addfunccol <- "addfunccol"
dist.junc <- get.full.dist.matrix(ssn.obj.path)
b.mat <- pmin(dist.junc, base::t(dist.junc))
flow.con.mat <- 1 - (b.mat > 0) * 1
afv <- all_space[c('locID', addfunccol)]
n.all.sites=nrow(all_space)
w.matrix <- sqrt(pmin(outer(afv[, addfunccol],rep(1, times = n.all.sites)),
                      base::t(outer(afv[, addfunccol],rep(1, times = n.all.sites) ))) /
                   pmax(outer(afv[, addfunccol],rep(1, times = n.all.sites)),
                        base::t(outer(afv[, addfunccol], rep(1, times = n.all.sites))))) * flow.con.mat
dist.hydro <- (dist.junc + t(dist.junc))

rownames(dist.hydro) <- colnames(dist.hydro) <- as.character(all_space$pid)
rownames(w.matrix) <- colnames(w.matrix) <- as.character(all_space$pid)

design.space.flow.con.mat = flow.con.mat[as.character(rawDFpred$pid), as.character(rawDFpred$pid)]
rownames(design.space.flow.con.mat) <- colnames(design.space.flow.con.mat) <- as.character(rawDFpred$pid)

# for each rid, select max upDist
up_preds_df = design_space[design_space$shreve == 1,]
up_preds_df$rid <- as.factor(up_preds_df$rid)
up_preds_df = do.call(rbind, lapply(split(up_preds_df,up_preds_df$rid), function(x) {return(x[which.max(x$upDist),])}))

# each unique path starting point pid
up_pids = as.character(up_preds_df$pid)

t = 50 # time points
n_d = 14 # number of sensors in design

# example design
set.seed(seed)

# pre-determined locations (pids) to predict in krieging
pids_pred <- as.character(rawDFobs$pid)
n_pred <- nrow(rawDFobs)

# randomly select paths for searching
PATH_UP_PIDS = sample(up_pids, n_d, replace=TRUE)

# get whole path for PATH_UP_PIDS
PATH_PIDS <- list()
for(i in 1:length(PATH_UP_PIDS)){
  # get all flow connected indexes, for series of points
  PATH_PIDS[i] = list(names(which(design.space.flow.con.mat[PATH_UP_PIDS[i],]==1)))
}

get_pid_from_d_ <- function(d){
  pids_d <- c()
  for(i in 1:n_d){
    candidate_up_dist = design_space[PATH_PIDS[[i]],]$upDist
    clost_upDist_idx = which.min(abs(candidate_up_dist - as.numeric(d[i,])))
    pid <- as.character(design_space[PATH_PIDS[[i]],][clost_upDist_idx,]$pid)
    pids_d <- c(pids_d, pid)
  }
  as.character(pids_d)
}

get_pid_from_d <- function(d){
  pids_d <- c()
  for(i in 1:n_d){
    candidate_up_dist <- design_space[PATH_PIDS[[i]],]$upDist
    distances <- abs(candidate_up_dist - as.numeric(d[i,]))
    sorted_indices <- order(distances)  # Order the indices by their distance to the target
    
    for (idx in sorted_indices) {
      pid <- as.character(design_space[PATH_PIDS[[i]],][idx,]$pid)
      if (!pid %in% pids_d) {
        pids_d <- c(pids_d, pid)
        break  # Break the inner loop once we find a unique pid
      }
    }
  }
  as.character(pids_d)
}

limits <- function(d, i, j) {
  # filter out any existing pid's in design from candidate selection
  pids_d <- get_pid_from_d(d)
  permissable_pids <- PATH_PIDS[[i]][!(PATH_PIDS[[i]] %in% pids_d)]
  design_space[permissable_pids,]$upDist
}

start.d <- data.frame(matrix(nrow = n_d, ncol = 1))
colnames(start.d) <- NULL
rownames(start.d) <- 1:n_d

for(i in 1:n_d){
  new_pt <- sample(limits(start.d,i), 1)
  start.d[i,] <- new_pt
}
get_pid_from_d(start.d)

# generate persistent anomalies
sensor_anom <- function(y, prop = 0.012, lambda = 1.7, seed = 1){
  set.seed(seed+1456)
  anom_idx <- sample(0:1, length(y), replace = T, prob = c(1-prop,prop))
  anom_idx <- sample(anom_idx)
  
  for(i in 1:length(anom_idx)){
    if(anom_idx[i] == 1) {
      len <- rpois(1,lambda)
      anom_idx[i:min(i+len, length(anom_idx))] <- 1
      
      i = i+len
    }
  }
  dim(anom_idx) <- dim(y)
  anom_idx
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

anomaly_score <- function(Y_d, B){
  
  # transform into long timeseries for oddstream input
  y <- matrix(NA, t*B, n_d)
  
  for (k in 1:B){
    y_k = Y_d[,k]
    dim(y_k) = c(n_d, t) # N x T
    # stack k iterations  vertically to synthesize longer timeframe
    y[((k - 1) * t + 1):(k * t),] <- t(y_k) # T*B x N
  }
  
  # define hyperparameters
  train_test_idx <- 50
  window_length <- 50
  n <- n_d
  
  # split into train/test data
  y_train <- y[1:train_test_idx, ]
  y_test <- y[-(1:train_test_idx), ]
  
  # get anomaly indexes
  ground_truth <- sensor_anom(y_test, prop=prop, lambda=lambda)
  ground_truth[1:window_length, ] <- 0 # a quirk of this method is the first window is disregarded so just n anoms here in the first place
  
  global_groundtruth_anom <<- sum(ground_truth)/length(ground_truth)
  
  # generate anomalies on y
  y_anoms <- ifelse(ground_truth==0, y_test, y_test+rnorm(length(y_test), 5, 10))
  
  if(global_i==0){ # just store one realisation for plotting
    global_y_anoms <<- y_anoms
    global_ground_truth <<- ground_truth
    global_y <<- y
  }
  
  # run anomaly detection algorithm (oddstream)
  output <- find_odd_streams(y_train, y_anoms, trials = 10, window_length=window_length)
  predicted_matrix <- output$out_marix
  
  # parse anomaly_idxs to same size as oddstream output
  n_summary <- nrow(ground_truth) %/% window_length - 1
  
  ground_truth_matrix <- matrix(NA, n_summary, n)
  
  # Iterate through the rows and aggregate by w rows
  for (i in 1:n_summary) {
    # Extract a subset of rows to aggregate
    subset_rows <- ((i * window_length)+1) : ((i +1)* window_length)
    
    # Compute the max for each column in the subset
    max_values <- apply(ground_truth[subset_rows, ], 2, max)
    
    # Store the max values in the summary matrix
    ground_truth_matrix[i, ] <- max_values
  }
  
  global_anom_proportion <<- c(global_anom_proportion, sum(ground_truth_matrix)/length(ground_truth_matrix))
  
  # the anomaly detection algorithm disregards the first window in test, and any trailing excess points
  output_idxs <- (window_length+1):((n_summary +1)* window_length) 
  
  # repeat predicted window anomalies for time scale, B
  predicted_anoms <- predicted_matrix[rep(1:nrow(predicted_matrix), each = window_length), ]
  y_anoms_trimmed <- y_anoms[output_idxs,]
  y_anoms_trimmed[predicted_anoms == 1] <- NA
  
  y_clean <- matrix(NA, nrow=B*t, ncol=n)
  y_clean[1:(train_test_idx+output_idxs[1]-1), ] = y[1:(train_test_idx+output_idxs[1]-1), ]
  y_clean[output_idxs+train_test_idx,] = y_anoms_trimmed
  
  Y_clean_reversed <- matrix(NA, n_d * t, B)
  for (k in 1:B) {
    y_k_transposed <- y_clean[((k - 1) * t + 1):(k * t), ]
    y_k <- t(y_k_transposed)
    Y_clean_reversed[, k] <- as.vector(y_k)
  }
  
  list("ground_truth_matrix" = ground_truth_matrix, "predicted_matrix" = predicted_matrix, "y_clean" = Y_clean_reversed)
}

utilpred <- function(d, B) {
  set.seed(seed+global_i)
  
  nugget <- runif(B, min = 0.00001, max = 0.0001)
  psill <- rep(1.5,B)
  range <- rep(2,B)
  phi <- rep(1,B)
  betas <- rep(c(1,2), each=B)
  dim(betas) <- c(B,2)
  
  # find closest pid on dense grid, given design (upDist)
  pids_d <- get_pid_from_d(d)
  
  # establish ordering and index, given new design
  pids_all <- c(pids_pred, pids_d)
  n_all = n_d + n_pred
  
  # filter design space for relevant pids, sort by pid_all
  df_s <- all_space[pids_all, , drop = FALSE]  # Use 'drop = FALSE' to keep it as a dataframe
  rownames(df_s) <- pids_all  
  
  # extend in time
  df_st <- do.call("rbind", replicate(t, df_s, simplify = FALSE)) # replicating the df
  df_st$date <- rep(1:t, each = (nrow(df_st)/t)) # set date variable
  rownames(df_st) = 1:nrow(df_st) # reindex for unique reference
  
  # index values in df_st for preds and design
  idx_pred_st = which(df_st$pid %in% pids_pred)
  idx_d_st = which(df_st$pid %in% pids_d)
  
  # covariates as matrix form for kriging calculation
  X_pred = df_st$X1[idx_pred_st]
  X_pred = cbind(rep(1, length(X_pred)), X_pred)
  X_d = df_st$X1[idx_d_st]
  X_d = cbind(rep(1, length(X_d)), X_d)
  
  # spatial covariance
  C_tu <- array(NA, c(n_all, n_all, B))
  dimnames(C_tu) <- list(pids_all, pids_all, NULL)
  inv_d =  array(NA, c(n_d,n_d,B))
  
  for(i in 1:B){
    C_tu[,,i] = get.covariance.matrix(pids_all, dist.hydro, w.matrix, psill[i], range[i]) + (nugget[i] * rnorm(n_all) * diag(n_all))
    inv_d[,,i]  <- solve(C_tu[pids_d, pids_d,i])
  }
  
  # transition matrix
  PHI <- array(0, dim = c(n_all, n_all, B))
  
  for(i in 1:B){
    PHI[,,i] = diag(phi[i], n_all, n_all)
  }
  
  # ytilde sampled at design (observation), and prediction sites
  Y_d <- matrix(NA, (nrow(X_d)/t)*t, B)
  Y_pred <- matrix(NA, (nrow(X_pred)/t)*t, B)
  
  # krieging prediction
  Y_postpred <- matrix(NA, (nrow(X_pred)/t)*t, B)
  
  for(k in 1:B){
    start_row_d <- 1
    end_row_d <- n_d
    start_row_pred <- 1
    end_row_pred <- n_pred
    
    # calculate y_1 and X_1
    mu_d_t = X_d[start_row_d:end_row_d, ] %*% betas[k,]
    mu_pred_t = X_pred[start_row_pred:end_row_pred, ] %*% betas[k,]
    
    # estimates
    Y_tilde <- rnorm(c(mu_d_t,mu_pred_t), C_tu[c(pids_d,pids_pred),c(pids_d,pids_pred),k])
    Y_d[start_row_d:end_row_d,k] <- Y_tilde[1:n_d]
    Y_pred[start_row_pred:end_row_pred,k] <- Y_tilde[-(1:n_d)]
    
    for(i in 2:t){
      # store info for t-1
      X_d_t_1 = X_d[start_row_d:end_row_d, ]
      X_pred_t_1 = X_pred[start_row_pred:end_row_pred, ]
      Y_d_t_1 = Y_d[start_row_d:end_row_d,k]
      Y_pred_t_1 = Y_pred[start_row_pred:end_row_pred,k]
      
      # update index to timestep t
      start_row_d <- (i - 1) * n_d + 1
      end_row_d <- i * n_d
      start_row_pred <- (i - 1) * n_pred + 1
      end_row_pred <- i * n_pred
      
      # calculate mean
      mu_d_t = X_d[start_row_d:end_row_d, ] %*% betas[k,] + PHI[1:n_d, 1:n_d,k ] %*%(Y_d_t_1 - X_d_t_1%*% betas[k,])
      mu_pred_t = X_pred[start_row_pred:end_row_pred, ] %*% betas[k,] + PHI[1:n_pred, 1:n_pred,k ]%*%(Y_pred_t_1 - X_pred_t_1%*% betas[k,])
      
      # estimates
      Y_tilde <- rnorm(c(mu_d_t,mu_pred_t), C_tu[c(pids_d,pids_pred),c(pids_d,pids_pred),k])
      Y_d[start_row_d:end_row_d,k] <- Y_tilde[1:n_d]
      Y_pred[start_row_pred:end_row_pred,k] <- Y_tilde[-(1:n_d)]
  
    }
  }
  
  # suppose that y is anomalous, get detection score, and 'clean' observations
  out <- anomaly_score(Y_d, B)
  Y_d_clean <- out$y_clean

  # utility to maximise, inverse of mse loss
  inv_mse <- rep(0,B)
  
  for(k in 1:B){
    start_row_d <- 1
    end_row_d <- n_d
    start_row_pred <- 1
    end_row_pred <- n_pred
    
    # calculate y_1 and X_1
    mu_d_t = X_d[start_row_d:end_row_d, ] %*% betas[k,]
    mu_pred_t = X_pred[start_row_pred:end_row_pred, ] %*% betas[k,]
    
    # posterior predictive mean
    Y_minus_mu <- (Y_d_clean[start_row_d:end_row_d,k] - mu_d_t)
    Y_minus_mu[is.na(Y_minus_mu)] <- 0
    Y_postpred[start_row_pred:end_row_pred,k] = mu_pred_t + C_tu[pids_pred,pids_d,k] %*% 
      inv_d[,,k] %*% Y_minus_mu
    
    for(i in 2:t){
      # store info for t-1
      X_d_t_1 = X_d[start_row_d:end_row_d, ]
      X_pred_t_1 = X_pred[start_row_pred:end_row_pred, ]
      Y_d_t_1 = Y_d[start_row_d:end_row_d,k]
      Y_pred_t_1 = Y_pred[start_row_pred:end_row_pred,k]
      
      # update index to timestep, t
      start_row_d <- (i - 1) * n_d + 1
      end_row_d <- i * n_d
      start_row_pred <- (i - 1) * n_pred + 1
      end_row_pred <- i * n_pred
      
      # calculate mean
      mu_d_t = X_d[start_row_d:end_row_d, ] %*% betas[k,] + PHI[1:n_d, 1:n_d,k ] %*%(Y_d_t_1 - X_d_t_1%*% betas[k,])
      mu_pred_t = X_pred[start_row_pred:end_row_pred, ] %*% betas[k,] + PHI[1:n_pred, 1:n_pred,k ]%*%(Y_pred_t_1 - X_pred_t_1%*% betas[k,])
      
      # posterior predictive mean
      Y_minus_mu <- (Y_d_clean[start_row_d:end_row_d,k] - mu_d_t)
      Y_minus_mu[is.na(Y_minus_mu)] <- 0
      Y_postpred[start_row_pred:end_row_pred,k] = mu_pred_t + C_tu[pids_pred,pids_d,k] %*% 
        inv_d[,,k] %*% Y_minus_mu
    }

    # utility based on accuracy of prediction
    inv_mse[k] <- 1/rmse_loss(Y_postpred[,k], Y_pred[, k])
  }
  
  # get accuracy of anomaly classification
  cm <- confusion_matrix(out$"ground_truth_matrix", out$"predicted_matrix")
  metrics <- calculate_metrics(cm)
  
  global_i <<- global_i+1
  global_metrics[[global_i]] <<- metrics
  global_accuracy <<- c(global_accuracy, mean(inv_mse, na.rm =TRUE))
  
  utility <- inv_mse * metrics[3]
  utility[!is.na(utility)]
}

colnames(start.d) <- 'd'
start.d$d <- as.numeric(start.d$d)

upper <- matrix(nrow=n_d, ncol=1)
for(i in 1:n_d){
  upper[i,] <- max(design_space[PATH_PIDS[[i]],]$upDist)
}

# assume start.d is a data.frame
ex44 <- ace(utilpred, start.d, c(1500,1000), limits=limits, N1=40, N2 = 0, lower=0, upper=upper)

# save
file.name=paste('river', scenario, '_', seed,'_mcc.Rdata', sep='')
save.image(file= file.name)
