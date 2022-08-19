### Analysis on the posterior distributions of spatial parameters
rm(list = ls())
graphics.off()
source("Scripts/Utilities.r")

simulation_id = "simulation_1"
dir_path = paste0("Simulations/",simulation_id)

fit = readRDS(paste0(dir_path,"/fitted.rds"))
fitting_data = readRDS(paste0(dir_path,"/fitting_input.rds"))
real_pars = readRDS(paste0(dir_path,"/real_pars.rds"))
ord = readRDS((paste0(dir_path,"/ord.rds")))

real_pars = real_pars[ord,]
nchains = 5
chains = NULL
real_pars$phi = sqrt(real_pars$lambda_1/real_pars$lambda_2)
for ( i in 1:nchains)
chains = rbind(chains, samples_from_chain(paste0(dir_path,'/chain_data_nngp_',i,'.csv')))


warmup = 200

chains = chains %>% select(starts_with("beta") | starts_with("iter")) %>% filter(iter > warmup) %>% select(starts_with("beta"))

distribution_spatial_pars = from_chains_to_pars_distr(chains,fitting_data)

est_fun = median
est_spatial_pars = from_chains_to_pars_est(chains,fitting_data,est_fun)

p = multiple_heatmaps(est_spatial_pars)
ggsave(filename = paste0(dir_path,"/plots/est_pars_linmodel.pdf"),plot = p,width = 30,height = 15,dpi = "retina")


error_funs = list()


error_funs[["mu"]] = abs_err

for ( name in names(est_spatial_pars)[-(1:3)]) error_funs[[name]] = relative_err

log_rescale = c(FALSE,FALSE,TRUE,TRUE,FALSE,FALSE)

errs = compute_errors(data = est_spatial_pars,real = real_pars, error_funs = error_funs,log_rescale = log_rescale)
p = multiple_heatmaps(errs)
ggsave(filename = paste0(dir_path,"/plots/errs_linmodel.pdf"),plot = p,width = 30,height = 15,dpi = "retina")



p = draw_pars_histograms(errs[3:8])
ggsave(filename = paste0(dir_path,"/plots/errs_linmodel_hists.pdf"),plot = p,width = 30,height = 15,dpi = "retina")
