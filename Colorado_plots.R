rm(list = ls())
graphics.off()
source("Scripts/Utilities.r")
colorado = read.csv("Colorado/Colorado_data.csv")

spatial_data = colorado %>% select(Longitude, Latitude, logPrecip)

p = ggplot(spatial_data) + 
  geom_point(aes(x = Longitude, y = Latitude, color = logPrecip),size=4) + 
  scale_color_viridis(option = "viridis", direction = -1)+ 
  theme_pubclean(base_size = 40)+ 
  theme(legend.text=element_text(size=40), legend.key.size = unit(2.5, 'cm'))+
  labs(col = "Logarithm of Yearly precipitation")+ ggtitle("(a)") + 
  coord_fixed()

ggsave(filename = "Colorado/Coloradodata.pdf",
       plot = p, width = 20,height = 15,dpi = "retina")


knitr::plot_crop("Colorado/Coloradodata.pdf")

##orography

grid = read.csv("Colorado/Colorado_grid.csv")

p = ggplot(grid) + 
  geom_tile(aes(x = longitude, y = latitude, fill = elevation)) + 
  scale_fill_viridis(option = "turbo")+ 
  theme_pubclean(base_size = 40)+ 
  theme(legend.text=element_text(size=40), legend.key.size = unit(2.5, 'cm'))+
  labs(col = "Elevation")+ ggtitle("(b)") + 
  coord_fixed()

ggsave(filename = "Colorado/Coloradoelev.pdf",
       plot = p, width = 20,height = 15,dpi = "retina")


knitr::plot_crop("Colorado/Coloradoelev.pdf")


##parameters estimation by RDD

est_rdd = read_rds("Colorado/RDD_estimate6.rds")
x_grid = seq(min(spatial_data$Longitude), max(spatial_data$Longitude), length.out = 100)
y_grid = seq(min(spatial_data$Latitude), max(spatial_data$Latitude), length.out = 100)

predgrid = expand.grid(x_grid,y_grid)

est_rdd$x_1 = predgrid$Var1
est_rdd$x_2 = predgrid$Var2

pars_rdd_toplot = est_rdd %>% select(x_1,x_2,mu,sigma,tau = nugget)

rdd_pars_plot = multiple_heatmaps(pars_rdd_toplot)





