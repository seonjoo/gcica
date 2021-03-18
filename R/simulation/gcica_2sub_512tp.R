# .libPaths('/ifs/scratch/msph/LeeLab/software/R/hpc') # set the R library in C2B2
library(dplyr)
library(coloredICA)
library(doParallel)
library(ggplot2)
registerDoParallel(cores = 12)
source('../gcica_R_translate_single_group.R')

n_sub = 2
S = 10
M = 3
N = 512

start_time = Sys.time()

male_simulation = lapply(1:S, function(i){
  A = rerow(matrix(runif(M^2)-0.5,M,M))
  W = solve(A)

  male = list(
    rbind(
      arima.sim(list(order=c(1,0,0), ar=-0.8),N),
      arima.sim(list(order=c(2,0,0), ar=c(0.9, -0.2)),N),
      arima.sim(list(order=c(2,0,0), ar=c(1.6,-0.64)),N)),
    rbind(
      arima.sim(list(order=c(1,0,0), ar=-0.8),N),
      arima.sim(list(order=c(2,0,0), ar=c(0.9, -0.2)),N),
      arima.sim(list(order=c(2,0,0), ar=c(1.6,-0.64)),N))
  )

  male_input = lapply(1:n_sub, function(i){A %*% male[[i]]})

  start_sim_time = Sys.time()
  gcica = gcica_bss_dwst_single(male_input, M = M, iter_print = 50)
  total_sim_time = Sys.time() - start_sim_time

  male_cor_list = lapply(1:n_sub, function(i){
    lapply(1:M, function(m){
      lapply(1:M, function(j){
        cor.test(male[[i]][m,], gcica$S[[i]][j,], method = 'spearman')$estimate %>% abs()
      }) %>% unlist() %>% max()
    })
  })

  result = new.env()
  result$male_cor_list = male_cor_list
  result$wlik = gcica$wlik
  result$amari = gcica$amari
  result$time_cost = total_sim_time
  result = as.list(result)
  return(result)
})

print(Sys.time() - start_time)

male_cor_df = do.call("rbind", lapply(1:S, function(s){
  do.call("rbind", lapply(1:n_sub, function(j){
    do.call("rbind", lapply(1:M, function(m){
      c(n_subject = n_sub, n_source = M, n_time_point = N,
        simulation = s, subject = j, source = m,
        cor = male_simulation[[s]]$male_cor_list[[j]][[m]],
        wlik = male_simulation[[s]]$wlik, amari = male_simulation[[s]]$amari,
        time_cost = male_simulation[[s]]$time_cost)
    }))
  }))
})) %>%
  as.data.frame() %>%
  mutate(
    simulation = as.factor(simulation),
    subject = paste('subject', subject),
    source = as.factor(source)
  )

write.csv(male_cor_df, paste(n_sub, 'sub_', M, 'comp', N, 'tp.csv', sep = ''), row.names=FALSE)

# png(paste(n_sub, 'sub_', N, 'tp.png', sep = ''))
# male_cor_df %>%
#   ggplot(aes(x = source, y = cor)) +
#   geom_boxplot() + facet_wrap(~subject)
# dev.off()
# save.image('/ifs/scratch/msph/LeeLab/qz2392/')

