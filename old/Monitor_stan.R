rm(list = ls())
graphics.off()
samples_from_chain = function(filename) {
  stanfit = rstan::read_stan_csv(filename)
  samples = data.frame(stanfit@sim$samples)
  samples %>% 
    filter(lp__ != 0) %>%  # Remove non-sampled iterations
    mutate(
      iter = 1:nrow(.),
      chain = gsub('.csv', '', last(strsplit(filename, '_')[[1]]))  # Get chain number [filename_1.csv]
    )
}

library(tidyverse)


#to_plot$beta_sigma.1.1 = exp(to_plot$beta_sigma.1.1)
#to_plot$beta_theta.1.1 = pi/2 * rstanarm::invlogit(to_plot$beta_theta.1.1)
#to_plot$beta_lambda_1.1.1 = exp(to_plot$beta_lambda_1.1.1)
#to_plot$beta_lambda_2.1.1 = exp(to_plot$beta_lambda_2.1.1)



chains = rbind(
  samples_from_chain('chain_data_nngp_complex_no_nuis_1.csv'),
  samples_from_chain('chain_data_nngp_complex_no_nuis_2.csv'),
  samples_from_chain('chain_data_nngp_complex_no_nuis_3.csv'),
  samples_from_chain('chain_data_nngp_complex_no_nuis_4.csv'),
  samples_from_chain('chain_data_nngp_complex_no_nuis_5.csv'),
  samples_from_chain('chain_data_nngp_complex_no_nuis_6.csv')
)

### Prepare for plot
to_plot = chains[,c(10:18,20,21)] 




graphics.off()
chains = rbind(
  samples_from_chain('chain_data_easy_1.csv'),
  samples_from_chain('chain_data_easy_2.csv'),
  samples_from_chain('chain_data_easy_3.csv'),
  samples_from_chain('chain_data_easy_4.csv'),
  samples_from_chain('chain_data_easy_5.csv'),
  samples_from_chain('chain_data_easy_6.csv')
)



to_plot = chains[,c(31:60,62,63)]
to_plot  = to_plot %>% gather(key, value, -chain, -iter) %>% filter(iter > 0)

# Visualize it
graphics.off()
x11()
to_plot %>% 
  ggplot(aes(x=iter, y=value, color=factor(chain))) + 
  facet_wrap(~key, scales = 'free_y') +
  geom_line()