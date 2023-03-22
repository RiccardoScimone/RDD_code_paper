rm(list = ls())
graphics.off()
source("Scripts/Utilities.r")
colorado = read.csv("Colorado/Colorado_data.csv")
est_rdd = read_rds("Colorado/RDD_estimate8.rds")
spatial_data = colorado %>% select(Longitude, Latitude, logPrecip)
x_grid = seq(min(spatial_data$Longitude), max(spatial_data$Longitude), length.out = 100)
y_grid = seq(min(spatial_data$Latitude), max(spatial_data$Latitude), length.out = 100)

predgrid = expand.grid(x_grid,y_grid)


p = ggplot(spatial_data) + 
  geom_point(aes(x = Longitude, y = Latitude, color = logPrecip),size=4) + 
  scale_color_viridis(option = "viridis", direction = -1)+ 
  theme_pubclean(base_size = 50)+ 
  theme(legend.text=element_text(size=40), legend.key.size = unit(2.5, 'cm'),plot.title = element_text(hjust = 0.5))+
  labs(col = "Logarithm of Yearly precipitation")+ ggtitle("(b)") + 
  coord_fixed()

ggsave(filename = "Colorado/Coloradodata.pdf",
       plot = p, width = 20,height = 15,dpi = "retina")


knitr::plot_crop("Colorado/Coloradodata.pdf")

##orography

grid = read.csv("Colorado/Colorado_grid.csv") %>% filter(longitude >= x_grid[1], longitude <= x_grid[100],
                                                         latitude >= y_grid[1],  latitude <= y_grid[100])
p = ggplot(grid) + 
  geom_tile(aes(x = longitude, y = latitude, fill = elevation)) + 
  scale_fill_viridis(option = "turbo")+ 
  theme_pubclean(base_size = 50)+ 
  theme(legend.text=element_text(size=40), legend.key.size = unit(2.5, 'cm'),plot.title = element_text(hjust = 0.5))+
  labs(col = "Elevation")+ggtitle("(a)") + 
  coord_fixed()

ggsave(filename = "Colorado/Coloradoelev.pdf",
       plot = p, width = 20,height = 15,dpi = "retina")


knitr::plot_crop("Colorado/Coloradoelev.pdf")


##parameters estimation by RDD



est_rdd$longitude = predgrid$Var1
est_rdd$latitude = predgrid$Var2

pars_rdd_toplot = est_rdd 


rdd_colorado_mu = ggplot(pars_rdd_toplot) + 
  geom_tile(aes(x = longitude, y = latitude, fill = mu)) + 
  scale_fill_viridis(option = "viridis", direction = -1)+ 
  theme_pubclean(base_size = 50)+ 
  theme(legend.text=element_text(size=40), 
        legend.key.size = unit(2.5, 'cm'),plot.title = element_text(hjust = 0.5))+
  labs(fill = expression(mu~" "))+ggtitle("(c)") + 
  coord_fixed()


rdd_colorado_sigma = ggplot(pars_rdd_toplot) + 
  geom_tile(aes(x = longitude, y = latitude, fill = sigma)) + 
  scale_fill_viridis(option = "viridis", direction = -1)+ 
  theme_pubclean(base_size = 50)+ 
  theme(legend.text=element_text(size=40), 
        legend.key.size = unit(2.5, 'cm'),plot.title = element_text(hjust = 0.5))+
  labs(fill = expression(sigma~" "))+ggtitle("(d)") + 
  coord_fixed()


rdd_colorado_tau = ggplot(pars_rdd_toplot) + 
  geom_tile(aes(x = longitude, y = latitude, fill = nugget)) + 
  scale_fill_viridis(option = "viridis", direction = -1)+ 
  theme_pubclean(base_size = 50)+ 
  theme(legend.text=element_text(size=40), 
        legend.key.size = unit(2.5, 'cm'),plot.title = element_text(hjust = 0.5))+
  labs(fill = expression(tau~" "))+
  coord_fixed()

ggsave(filename = "Colorado/Coloradomurdd.pdf",
       plot = rdd_colorado_mu, width = 20,height = 15,dpi = "retina")


knitr::plot_crop("Colorado/Coloradomurdd.pdf")

ggsave(filename = "Colorado/Coloradosigmardd.pdf",
       plot = rdd_colorado_sigma, width = 20,height = 15,dpi = "retina")


knitr::plot_crop("Colorado/Coloradosigmardd.pdf")


ggsave(filename = "Colorado/Coloradotaurdd.pdf",
       plot = rdd_colorado_tau, width = 20,height = 15,dpi = "retina")


knitr::plot_crop("Colorado/Coloradotaurdd.pdf")

p_ell = add_ellipses(est_rdd %>% rename(x_1 = longitude,x_2 = latitude),p,scale = 0.5) + ggtitle("RDD estimates")+ 
  scale_x_continuous(breaks = seq(-108,-102, by = 2)) +scale_y_continuous(breaks = 37:41)

ggsave(filename = "Colorado/Coloradoellipsesrdd.pdf",
       plot = p_ell, width = 20,height = 15,dpi = "retina")


knitr::plot_crop("Colorado/Coloradoellipsesrdd.pdf")

Pac_est = read_rds("Colorado/Paciorek_estimates.rds")

p_ell = add_ellipses(Pac_est%>%
                       filter(x_1 >= x_grid[1], x_1 <= x_grid[100],
                              x_2 >= y_grid[1],  x_2 <= y_grid[100]),p,scale = 0.35) + ggtitle("estimates by model (3.5)")+ 
  scale_x_continuous(breaks = seq(-108,-102, by = 2)) +scale_y_continuous(breaks = 37:41)

ggsave(filename = "Colorado/Coloradoellipsespac.pdf",
       plot = p_ell, width = 20,height = 15,dpi = "retina")


knitr::plot_crop("Colorado/Coloradoellipsesrdd.pdf")
## Confronto stime di teta

rdd_colorado_teta = ggplot(pars_rdd_toplot) + 
  geom_tile(aes(x = longitude, y = latitude, fill = theta_deg)) + 
  scale_fill_viridis(option = "viridis", direction = -1)+ 
  theme_pubclean(base_size = 50)+ 
  theme(legend.text=element_text(size=40), 
        legend.key.size = unit(2.5, 'cm'),plot.title = element_text(hjust = 0.5))+
  labs(fill = expression(theta~"(degrees)"))+
  coord_fixed()


ggsave(filename = "Colorado/Coloradothetardd.pdf",
       plot = rdd_colorado_teta, width = 20,height = 15,dpi = "retina")


knitr::plot_crop("Colorado/Coloradothetardd.pdf")


Pac_colorado_teta  = ggplot(Pac_est %>% rename(longitude = x_1,latitude = x_2)%>%
                              filter(longitude >= x_grid[1], longitude <= x_grid[100],
                                     latitude >= y_grid[1],  latitude <= y_grid[100]) ) + 
  geom_tile(aes(x = longitude, y = latitude, fill = theta_deg)) + 
  scale_fill_viridis(option = "viridis", direction = -1)+ 
  theme_pubclean(base_size = 50)+ 
  theme(legend.text=element_text(size=40), 
        legend.key.size = unit(2.5, 'cm'),plot.title = element_text(hjust = 0.5))+
  labs(fill = expression(theta~"(degrees)"))+
  coord_fixed()


ggsave(filename = "Colorado/ColoradothetaPac.pdf",
       plot = Pac_colorado_teta, width = 20,height = 15,dpi = "retina")


knitr::plot_crop("Colorado/ColoradothetaPac.pdf")


