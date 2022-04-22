# function to run DABOM on multiple cores equal to number of chains using
# dclone and jags.parfit.
run_dabom_parallel <- function(model, data, jags_params, inits,
                      n.chains, n.adapt, n.burn, n.iter, thin, filename){
  require(rjags)
  require(dclone)  
  
  n.cores <- n.chains
  cl <- makePSOCKcluster(n.cores)
  tmp <- clusterEvalQ(cl, library(dclone))
  
  timer <- proc.time()
  
  dabom_mod <- try(jags.parfit(cl = cl, 
                               data = jags_data,
                               params = jags_params,
                               model = final_mod_file,
                               inits = init_fnc,
                               n.chains = n.chains, 
                               n.adapt = n.adapt, 
                               n.update = n.burn,
                               n.iter = n.iter, 
                               thin = thin),
                   silent = TRUE)
  
  time.taken <- proc.time() - timer
  run_time <- round(time.taken[3] / 60)
  
  if(inherits(dabom_mod, "try-error")){
    
    stopCluster(cl)
    
    stop('DABOM failed to run after ', run_time, ' minutes with the following error message: \n', dabom_mod[1])
  }
  
  stopCluster(cl)
  
  dabom_output <- list('dabom_mod' = dabom_mod,
                       'jags_data' = data,
                       'run_time' = run_time)
  
  save(dabom_output,
       file = filename)
  
  cat('DABOM ran successfully and took', run_time, 'minutes to complete on', n.cores,
      'parallel cores.', '\n')  
  
  return(dabom_output)
}