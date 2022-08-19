rm(list = ls())
graphics.off()
source("Scripts/Utilities.r")
colorado = read.csv("Scripts/Colorado_data.csv",)

spatial_data = colorado %>% select(Longitude, Latitude, logPrecip)

p = ggplot(spatial_data) + geom_point(aes(x = Longitude, y = Latitude, color = logPrecip),size=3) + scale_color_viridis(option = "viridis", direction = -1)+ theme_pubclean(base_size = 30)+ 
  theme(legend.text=element_text(size=15), legend.key.size = unit(1.5, 'cm'))+ coord_fixed()

ggsave(filename = "Colorado/Coloradodata.pdf",plot = p, width = 30,height = 15,dpi = "retina")
## Vediamo la superficie di Kriging
library(gstat)
library(sp)

x_grid = seq(min(spatial_data$Longitude), max(spatial_data$Longitude), length.out = 100)
y_grid = seq(min(spatial_data$Latitude), max(spatial_data$Latitude), length.out = 100)

predgrid = expand.grid(x_grid,y_grid)
grid_plain = predgrid
names(predgrid) = c("Lon","Lat")
coordinates(spatial_data) = c("Longitude","Latitude")
coordinates(predgrid) = c("Lon","Lat")

model = vgm(model = "Mat",kappa = 5/2,nugget = 0.05)
bin = variogram(logPrecip ~ 1,data = spatial_data)
x11()
model = fit.variogram(bin,model,fit.sills = c(T,T))
plot(bin,model)

ord_kriging =  krige(logPrecip ~ 1, locations = spatial_data, newdata = predgrid, model = model)
p = plot_ly(x = grid_plain$Var1, y = grid_plain$Var2, z = ord_kriging$var1.pred,type = "scatter3d", size = 0.1, color = ord_kriging$var1.pred) 
print(p)

pred_data = data.frame(lon = grid_plain$Var1,lat = grid_plain$Var2,precip = ord_kriging$var1.pred)
x11()
ggplot(pred_data) + geom_tile(aes(x = lon, y = lat, fill = precip)) + scale_fill_viridis(direction = -1)+coord_fixed() + theme_pubclean(base_size = 20)

# Facciamo RDD
B = 360
K = 6
cl = makePSOCKcluster(12)
registerDoParallel(cl)
dir_path = "Colorado"
writeLines(c(""), paste0(dir_path,"/log.txt"))
N_obs = nrow(spatial_data)
N_grid = nrow(grid_plain)

colorado = read.csv("Scripts/Colorado_data.csv",)
spatial_data = colorado %>% select(Longitude, Latitude, logPrecip)
ret = foreach (i=1:B,.packages = c("RGeostats","radiant.data")) %dopar% {
  
  centers = grid_plain[sample(1:N_grid, size = K),]
  
  assignment_obs = assign(spatial_data[,1:2],centers,K,N_obs)
  
  assignment_pred = assign(grid_plain,centers, K, N_grid)
  
  sink(paste0(dir_path,"/log.txt"),append = T)
  
  while(min(table(assignment_obs)) < 10){
    print("assign again")
    centers = grid_plain[sample(1:N_grid, size = K),]
    
    assignment_obs = assign(spatial_data[,1:2],centers,K,N_obs)
    
    assignment_pred = assign(grid_plain,centers, K, N_grid)
  }
  
  mat_param = t(matrix(NA,N_grid + N_obs,5))
  for ( j in 1:K)
  { 
    idxs_obs = which(assignment_obs == j)
    idxs_pred = which(assignment_pred == j)
    db = db.create(spatial_data[idxs_obs,],flag.grid = F,ndim = 2,autoname = F)
    dirs = seq(0,150, by = 30)
    vario = vario.calc(db,dir = dirs)
    fit = model.auto(vario,struct = c("K-Bessel"),param = c("P=5/2","A=60"),lower = c("P=5/2","A = 10"),upper = c("P=5/2","A = 90"),verbose = 0) 
    mat_param[2:5,idxs_pred] = unlist(get_param(fit[1]))
    mat_param[1,idxs_pred] = mean(spatial_data[idxs_obs,3])
    mat_param[2:5,N_grid + idxs_obs] = unlist(get_param(fit[1]))
    mat_param[1,N_grid + idxs_obs] = mean(spatial_data[idxs_obs,3])
  }
  print(i)
  sink()
  t(mat_param)
}
stopCluster(cl)

ret = simplify2array(ret)
RDD_estimate = extract_spatial_predictions_median(ret)

saveRDS(ret,paste0(dir_path,"/RDD_fitting",K,".rds"))
saveRDS(RDD_estimate,paste0(dir_path,"/RDD_estimate",K,".rds"))

saveRDS(grid_plain,paste0(dir_path,"/predgrid.rds"))
