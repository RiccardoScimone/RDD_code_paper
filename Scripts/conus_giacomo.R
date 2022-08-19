rm(list = ls())
source("Scripts/Utilities.r")
library(data.table)
download = FALSE
dir_path = "GHCN/"
if(download)
{
  download.file("https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/ghcnd-countries.txt", destfile = "ghcnd-countries.txt")
  download.file("https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/ghcnd-inventory.txt", destfile = "ghcnd-inventory.txt")
  download.file("https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt", destfile = "ghcnd-stations.txt")
  download.file("https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt", destfile = "readme.txt")
  download.file("https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/ghcnd_all.tar.gz", destfile = "ghcnd-all.tar.gz")
  stnscsv = paste0(dir_path,"stations.csv")
  typedcols = c( "A11", "F9", "F10", "F7", "X1","A2",
                 "X1","A30", "X1", "A3", "X1", "A3", "X1", "A5" )
  stns = read.fortran(paste0(dir_path,"ghcnd-stations.txt"),
                      typedcols,
                      comment.char="")
  hdrs = c("ID", "LAT", "LON", "ELEV", "ST", "NAME","GSN", "HCN", "WMOID")
  names(stns) = hdrs
  write.csv(stns,stnscsv)

  inventorycsv = paste0(dir_path,"inventory.csv")
  invcols = c( "A11", "X1", "F8", "X1", "F9", "X1","A4",
               "X1","I4", "X1", "I4" )
  inv = read.fortran(paste0(dir_path,"ghcnd-inventory.txt"),
                     invcols,
                     comment.char="")
  invhdrs = c("ID", "LAT", "LON", "ELEM" , "FIRST", "LAST")
  names(inv) = invhdrs
  write.csv(inv,inventorycsv)
}


### filter us stations, we are only interested in precipitation
islands = c("AK","HI","PI","UM")
inv = fread(inventorycsv)
stns = fread(stnscsv)
us_inv_prcp = inv %>% filter(substr(ID,1,2) == "US", ELEM == "PRCP", FIRST <= 2017, LAST >= 2018)
us_st_prcp = stns %>% filter(substr(ID,1,2) == "US")
us_infos = us_inv_prcp %>% left_join(us_st_prcp %>% select(ID,ELEV,ST,NAME),by = c("ID" = "ID")) %>% filter(!(ST %in% islands))

numFiles = length(us_infos$ID)
dirname = paste0(dir_path,"ghcnd_all/")

library(foreach)
library(doParallel)
cl = makePSOCKcluster(12)
registerDoParallel(cl)
writeLines(c(""), paste0(dir_path,"/log.txt"))
foreach (i = 1:numFiles,.packages = c("tidyverse","data.table")) %dopar% {
  sink(paste0(dir_path,"/log.txt"),append = T)
  infile = paste0(dirname, us_infos$ID[i], ".dly")
  outfile = paste0(dirname, us_infos$ID[i], ".csv")
  cols = c( "A11", "I4", "I2", "A4",
             rep( c( "I5", "A1", "A1", "A1"), 31) )
  df = read.fortran(infile, cols, na.strings="-9999")
  tmp = c("Val","xxM","xxQ","xxS") # xx so we can ditch them later
  vhdrs = paste(   rep(tmp,31),   rep(1:31,each=4), sep="")
  hdrs = c("ID", "year", "month", "element", vhdrs)
  names(df) = hdrs
  df = df %>% filter(year <= 2018 & year >= 2017 & element == "PRCP" )
  df_out = dplyr::select(df, -matches("xx*")) # get rid of M, Q, S
  fwrite(df_out, outfile)
  print(i)
  sink()
}
stopCluster(cl)

station_list = list()
good = 0

for ( i in 1:numFiles) ## Recupero i dati di Risser e Calder (Spererei)
{
  filename = paste0(dirname, us_infos$ID[i], ".csv")
  df = fread(filename)
  df = df %>% filter((year == 2017 & month >= 10) | (year == 2018 & month <= 9))
  if(nrow(df) == 12 & sum(is.na(df)) <= 7){
    good = good + 1
    station_list[[good]] = df
    print(good)
  }
}


full_df = bind_rows(station_list)
full_df[is.na(full_df)] = 0
full_df = full_df %>% pivot_longer(cols = 5:35,names_to = "DAYS", values_to = "PRECIPITATION")
full_df$element = NULL
full_df = full_df %>% group_by(ID) %>% summarise(logprep = log(mean(0.1*PRECIPITATION)))
fwrite(full_df,"average_dataset.csv")

## Insert latitude and longitude
full_df = full_df %>% left_join(us_infos%>%select(ID,LAT,LON,ELEV))
full_df = full_df %>% filter(logprep != -Inf)
fwrite(full_df,"GHCN/average_dataset.csv")

rm(list = ls())
source("Scripts/Utilities.r")
library(data.table)
library(ggmap)
library(maps)
library(mapdata)
library(geoR)
download = FALSE
dir_path = "GHCN/"
full_df = fread("GHCN/average_dataset_risser.csv") 
spatial_data = as.matrix(full_df%>% dplyr::select(LON,LAT,logprep) )
full_df$log_precipitation = trim(full_df$logprep,-1.7,2.2)
states = map_data("state")
p = ggplot() + geom_polygon(data = states, aes(x=long, y = lat, group = group),fill = NA, color = "blue") +  coord_fixed(1.3) +  theme_pubclean(base_size = 30)+ theme(legend.text=element_text(size=15), legend.key.size = unit(1.5, 'cm'))+guides(fill=FALSE)
x11()
pioggia = p + geom_point(data = full_df, aes(x = LON, y=LAT, col = log_precipitation)) + scale_color_viridis(direction = -1,breaks = seq(-1.5,2, by = 0.5), limits = c(-1.7,2.2))+  theme_pubclean(base_size = 30)+ theme(legend.text=element_text(size=15), legend.key.size = unit(1.5, 'cm'))
plot(pioggia)
borders = map_data("usa")
grid_ini = expand.grid(seq(-124,-67,by = 0.3),seq(25,49.5,by = 0.3))
grid = polygrid(seq(-124,-67,by = 0.3),seq(25,49.5,by = 0.3),borders,vec.inout = T)
check = pioggia + geom_point(data = grid$xypoly, aes(x = x, y=y), col = "red")
grid2 = polygrid(seq(-124,-67,by = 0.3),seq(25,49.5,by = 0.3),rbind(c(-125,49),c(-101.5,29.5), c(-72,41)),vec.inout = T)
idxs = grid$vec.inout | grid2$vec.inout

grid_plain = grid_ini[idxs,]
names(grid_plain) = c("long","lat")
check = pioggia + geom_point(data = grid_plain, aes(x = long, y=lat), col = "red")
x11()
plot(check)
graphics.off()
grid_plain = as.matrix(grid_plain)
B = 120
K = 8
dir_path = "GHCN"
struct = c("K-Bessel")
param = c("P=2")
lower = c("P=2")
upper = c("P=2")
dirs = seq(0,150, by = 30)

estimates = RDD_fit(data = full_df$logprep,coords = spatial_data[,1:2],center_grid = grid_plain,
                    K = K, B = B, clusters = 12,dir_path = dir_path , struct = struct,dirs = dirs,
                    param = param,lower = lower, upper = upper, threshold = 20)


ret = RDD_predict(RDD_output = estimates, grid = grid_plain)
ret = simplify2array(ret)
RDD_estimate = extract_spatial_predictions_median(ret)

saveRDS(ret,paste0(dir_path,"/RDD_fitting",K,".rds"))
saveRDS(RDD_estimate,paste0(dir_path,"/RDD_estimate",K,".rds"))

saveRDS(grid_plain,paste0(dir_path,"/predgrid.rds"))
