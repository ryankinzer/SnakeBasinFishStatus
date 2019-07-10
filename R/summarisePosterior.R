summarisePosterior <- function(data, summary_var, group_var, pt_est_nm = NULL, cred_int_prob = 0.95){
  
  summary_var = enquo(summary_var)
  group_var = enquo(group_var)
  
  df <- data %>%
    filter(!is.na(!!summary_var))
  
  credInt <- df %>%
    select(iter, !!group_var, !!summary_var) %>%
    spread(!!group_var, !!summary_var) %>%
    select(-iter) %>%
    coda::as.mcmc() %>%
    coda::HPDinterval(prob = cred_int_prob) %>%
    as.data.frame() %>%
    mutate(TRT = rownames(.)) %>%
    rename(lowerCI = lower,
           upperCI = upper) %>%
    tbl_df() %>%
    select(TRT, everything())
  
  
  post_summ = df %>%
    group_by(!!group_var) %>%
    summarise(mean = mean(!!summary_var),
              median = median(!!summary_var),
              mode = estMode(!!summary_var),
              sd = sd(!!summary_var),
              cv = sd / mean) %>%
    mutate_at(vars(mean, median, mode, sd),
              funs(ifelse(. < 0, 0, .))) %>%
    left_join(credInt,
              by = 'TRT') %>%
    mutate_at(vars(mean, median, mode),
              funs(round)) %>%
    mutate_at(vars(sd, lowerCI, upperCI),
              funs(round),
              digits = 1) %>%
    mutate_at(vars(cv),
              funs(round),
              digits = 3)
  
  # if(!is.null(pt_est_nm) & pt_est_nm %in% c('mean', 'median', 'mode')) {
  #   names(post_summ)[match(pt_est_nm, names(post_summ))] = estimate
  #   post_summ = post_summ %>%
  #     select(area, estimate, sd:upperCI)
  # }
  
  return(post_summ)
}
