rm(list = ls())
source("Scripts/Utilities.r")
graphics.off()


simulation_id = "simulation_1"
dir_path = paste0("Simulations/",simulation_id)

nchains = 3
chains = NULL

princ_name = '/chain_data_fitRDDlambda_2_'

for ( i in 1:nchains)
  chains = rbind(chains, samples_from_chain(paste0(dir_path,princ_name,i,'.csv')))


### Prepare for plot
#to_plot = chains%>% select(starts_with("beta") | starts_with("iter") | starts_with("chain"))
to_plot = chains%>% select(starts_with("sigma") |starts_with("tau") | starts_with("rho" )| starts_with("l_rho")| starts_with("iter") | starts_with("chain"))

to_plot  = to_plot %>% gather(key, value, -chain, -iter) %>% filter(iter > 100)
x11() 
to_plot %>% 
  ggplot(aes(x=iter, y=value, color=factor(chain))) + 
  facet_wrap(~key, scales = 'free_y') +
  geom_line()