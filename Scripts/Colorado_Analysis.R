rm(list = ls())
graphics.off()
source("Scripts/Utilities.r")
colorado = read.csv("Colorado/Colorado_data.csv")

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
B = 120
for(K in c(2,4,6,8)){
cl = makePSOCKcluster(12)
registerDoParallel(cl)
dir_path = "Colorado"
writeLines(c(""), paste0(dir_path,"/log.txt"))
N_obs = nrow(spatial_data)
N_grid = nrow(grid_plain)

colorado = read.csv("Colorado/Colorado_data.csv",)
spatial_data = colorado %>% select(Longitude, Latitude, logPrecip)


struct = c("K-Bessel", "Nugget Effect")
param = c("P=2")
lower = c("P=2")
upper = c("P=2")
dirs = seq(0,150, by = 30)
estimates = RDD_fit(data = colorado$logPrecip,
                    coords = colorado[,c("Longitude","Latitude")],
                    center_grid = grid_plain,
                    K = K, B = B, 
                    clusters = 12,
                    dir_path = dir_path, 
                    struct = struct,
                    dirs = dirs,
                    param = param,lower = lower, upper = upper, threshold = 15)


ret = RDD_predict(estimates, grid_plain) 

ret = simplify2array(ret)

RDD_estimate = extract_spatial_predictions_median(ret)

saveRDS(ret,paste0(dir_path,"/RDD_fitting",K,".rds"))
saveRDS(RDD_estimate,paste0(dir_path,"/RDD_estimate",K,".rds"))
saveRDS(grid_plain,paste0(dir_path,"/predgrid.rds"))
}


n_knot_grid = 4
library(LocallyStationaryModels)
d = cbind(colorado$Longitude,colorado$Latitude)
a = find_anchorpoints.lsm(d,n = n_knot_grid,plot_output = F)
z = as.matrix(colorado$logPrecip)
vario = variogram.lsm(z = z,
                      d = d,
                      a = a$anchorpoints,
                      epsilon = 0.25,
                      n_angles=6,
                      n_intervals=16,
                      dim=1,
                      kernel_id="gaussian")
solu = findsolutions.lsm(vario,lower.delta = 0.01, "exponentialnugget", c(0.01,0.01,0.4,3,0.3), 
                         upper.bound = c(0.5,0.5,pi/2,10,1),
                         remove_not_convergent = T)



pars = smooth.lsm(solu,newpoints = as.matrix(grid_plain))

pred = predict.lsm(solu,
                   as.matrix(grid_plain),
                   plot_output = F,predict_y = F)

pars_df = cbind(grid_plain,pars,pred$smoothed_means[,1])
Fouedjofit_model = pars_df
Fouedjofit_anchors = solu


saveRDS(file = paste0(dir_path,"/Fouedjo_anchors.rds"),object = Fouedjofit_anchors)
saveRDS(file = paste0(dir_path,"/Fouedjo_results.rds"),object = Fouedjofit_model)


names(Fouedjofit_model) = c("x_1","x_2",
                            "lambda[1]",
                            "lambda[2]",
                            "theta[deg]",
                            "sigma",
                            "tau",
                            "mu")
Fouedjofit_model$`theta[deg]` = 180/pi*Fouedjofit_model$`theta[deg]`
Fouedjofit_model$phi = Fouedjofit_model$`lambda[1]`/ Fouedjofit_model$`lambda[2]`
Fouedjofit_model$`lambda[1]` = Fouedjofit_model$`lambda[1]`^2
Fouedjofit_model$`lambda[2]` = Fouedjofit_model$`lambda[2]`^2

Fouedjofit_model = Fouedjofit_model[,c(1,2,8,6,3,4,5,7,9)]


saveRDS(file = paste0(dir_path,"/Fouedjo_results.rds"),object = Fouedjofit_model)

