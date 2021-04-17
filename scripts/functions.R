
# create table
make_MCMC_summary = function(MCMCmodel){
  # saved MCMC model
  summary_output = summary(MCMCmodel)
  
  mod_table = bind_rows(as.data.frame(summary_output$solutions),
                        as.data.frame(summary_output$Gcovariances),
                        as.data.frame(summary_output$Rcovariances))
  mod_table = mod_table %>% 
    round(2)
  DIC = c(round(summary_output$DIC,2), rep("", (dim(mod_table)[1]-1)))
  lambda = MCMCmodel$VCV[,'phylo']/
    (MCMCmodel$VCV[,'phylo']+MCMCmodel$VCV[,'units']+MCMCmodel$VCV[,'MSW93_Binomial'])
  Lambda = c(round(mean(lambda),2), round(HPDinterval(lambda)[1], 2),
             round(HPDinterval(lambda)[2],2),
             rep("", (dim(mod_table)[1]-3)))
  names(DIC)=names(mod_table)
  names(Lambda)=names(mod_table)
  mod_table = rbind(mod_table, DIC, Lambda)
  rownames(mod_table)[(nrow(mod_table)-1):nrow(mod_table)] = c("DIC", "Lambda")
  mod_table = mod_table %>% 
    add_column(Predictor = row.names(mod_table), .before=1)
  mod_table = flextable(mod_table)
  return(mod_table)
}
