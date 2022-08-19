rm(list = ls())
graphics.off()
source("Scripts/Utilities.r")
simulation_id = "simulation_7"
dir_path = paste0("Simulations/",simulation_id)
load(paste0(dir_path,"/simulated_process_list.Rdata"))
#for ( N in c(1000,2000,5000,8000,10000)){
for ( N in c(8000,10000)){
for ( K in c(1,4,8,16,32)){
#  for ( K in c(16,32)){
    if( K==32 & N ==1000 | (N == 2000 & K == 32 ))break
    ret = list()
    RDD_estimate = list()
    for ( j in 1:10 ){
      sim_const = simulated_process_list[[j]]
      if(N < nrow(sim_const)){
        set.seed(14061997)
        simulated_process = sim_const[sample(1:nrow(sim_const),N),]
      }
      else simulated_process = sim_const
    B = 120
    if (K == 1) B = 1
    struct = c("K-Bessel", "Nugget Effect")
    param = c("P=2")
    lower = c("P=2")
    upper = c("P=2")
    dirs = seq(0,150, by = 30)
    estimates = RDD_fit(data = simulated_process$process,coords = simulated_process[,1:2],center_grid = simulated_process[,1:2],
                        K = K, B = B, clusters = 12,dir_path = dir_path , struct = struct,dirs = dirs,
                        param = param,lower = lower, upper = upper, threshold = 15)
    ret[[j]] = RDD_predict(estimates,grid = sim_const[,1:2])
    ret[[j]] = simplify2array(ret[[j]])
    RDD_estimate[[j]] = extract_spatial_predictions_median(ret[[j]])
    print(paste0(N,"_",K,"_",j))
    }
    saveRDS(ret,paste0(dir_path,"/RDD_fitting_mr_",K,N,".rds"))
    saveRDS(RDD_estimate,paste0(dir_path,"/RDD_estimate_mr_",K,N,".rds"))
  
}

}