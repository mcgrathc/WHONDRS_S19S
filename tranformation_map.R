#### Plotting ICR sample data on a map
# Kings Creek (N04D) in Kansas had a latitude listed as ~69...which is not in Kansas. The other King's Creek lat ~39, so we
# adjusted accordingly

# Load in libraries
library(reshape2)
library(ggplot2)
library(ggpubr)
library(stringr)
library(dplyr)

# ################## #
#### Load in data ####
# ################## #

setwd("C:\\Users\\mcgr323\\OneDrive - PNNL\\Documents\\WHONDRS\\Transformation\\mapping")

# Load in ICR data
data = read.csv("S19S_Transformation_Mean_Difference_2.csv", stringsAsFactors = F)

# Load in meta-data
S19S.lat.long = read.csv("S19S_metadata_april.csv", stringsAsFactors = F)

# Ordering meta data
S19S.lat.long = S19S.lat.long[order(S19S.lat.long$ID),]

# Subsetting only to lat.long
S19S.lat.long = S19S.lat.long[,c("ID", "Latitude", "Longitude")]

# Renaming first column
colnames(S19S.lat.long)[1] = "Site"

trans_ll  <-as.data.frame(cbind(S19S.lat.long, data))

# Contiguous US boundaries
north = 49.3457868
south = 24.7433195
east = -66.9513812
west = -124.7844079

missi = -89.253333

# Filtering data based on US boundaries
trans_ll.us = trans_ll[which(trans_ll$Longitude < east & trans_ll$Longitude > west 
                             & trans_ll$Latitude > south & trans_ll$Latitude < north),]

### Performing East/West analyses
# Idenifying east-west of the Mississippi River
trans_ll.us$East.West = "East"
trans_ll.us$East.West[which(trans_ll.us$Longitude < missi)] = "West"


trans_ll.us <- na.omit(trans_ll.us)
##Amino acid transformation##
# Limits on ggplot scale
limitsaa = c(min(trans_ll.us$amino.acid), 
           max(trans_ll.us$amino.acid))
#plot the data
aa_plot <- ggplot(data = trans_ll.us)+
  borders("state", colour = "black", fill = "white")+ 
  geom_point(aes(x = Longitude, y = Latitude, color = amino.acid), size = 3.5) +
  scale_color_gradient2(limits = limitsaa, low = "dodgerblue2",mid ="goldenrod2",
                        high = "firebrick2", midpoint = (max(limitsaa)+min(limitsaa))/2, name = "Difference (%)")+
  ggtitle("Amino Acid transformations difference by site")+
  #coord_sf(xlim = c(-124, -70), ylim = c(25, 49), expand = FALSE)+
  theme_bw() + theme(axis.line=element_blank(),
                                     axis.text.x=element_blank(),
                                     axis.text.y=element_blank(),
                                     axis.ticks=element_blank(),
                                     axis.title.x=element_blank(),
                                     axis.title.y=element_blank(),
                                     panel.background=element_blank(),
                                     panel.border=element_blank(),
                                     panel.grid.major=element_blank(),
                                     panel.grid.minor=element_blank(),
                                     plot.background=element_blank())

##N transformation##
# Limits on ggplot scale
limitsn = c(min(trans_ll.us$N), 
           max(trans_ll.us$N))
#plot the data
N_plot <- ggplot(data = trans_ll.us)+
  borders("state", colour = "black", fill = "white")+ 
  geom_point(aes(x = Longitude, y = Latitude, color = N), size = 3.5) +
  scale_color_gradient2(limits = limitsn, low = "dodgerblue2",mid ="goldenrod2",
                        high = "firebrick2", midpoint = (max(limitsn)+min(limitsn))/2, name = "Difference (%)")+
  ggtitle("Nitrogen transformation difference by site")+
  #coord_sf(xlim = c(-124, -70), ylim = c(25, 49), expand = FALSE)+
  theme_bw() + theme(axis.line=element_blank(),
                     axis.text.x=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks=element_blank(),
                     axis.title.x=element_blank(),
                     axis.title.y=element_blank(),
                     panel.background=element_blank(),
                     panel.border=element_blank(),
                     panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),
                     plot.background=element_blank())
##P transformation##
# Limits on ggplot scale
limits =c(min(trans_ll.us$P), 
          max(trans_ll.us$P))
#plot the data
P_plot <- ggplot(data = trans_ll.us)+
  borders("state", colour = "black", fill = "white")+ 
  geom_point(aes(x = Longitude, y = Latitude, color = P), size = 3.5) +
  scale_color_gradient2(limits = limitsp, low = "dodgerblue2",mid ="goldenrod2",
                        high = "firebrick2", midpoint = (max(limitsp)+min(limitsp))/2, name = "Difference (%)")+
  ggtitle("Phosphorus transformation difference by site")+
  #coord_sf(xlim = c(-124, -70), ylim = c(25, 49), expand = FALSE)+
  theme_bw() + theme(axis.line=element_blank(),
                     axis.text.x=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks=element_blank(),
                     axis.title.x=element_blank(),
                     axis.title.y=element_blank(),
                     panel.background=element_blank(),
                     panel.border=element_blank(),
                     panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),
                     plot.background=element_blank())

##S transformation##
# Limits on ggplot scale
limits = c(min(trans_ll.us$S), 
           max(trans_ll.us$S))
#plot the data
S_plot <- ggplot(data = trans_ll.us)+
  borders("state", colour = "black", fill = "white")+ 
  geom_point(aes(x = Longitude, y = Latitude, color = S), size = 3.5) +
  scale_color_gradient2(limits = limitss, low = "dodgerblue2",mid ="goldenrod2",
                        high = "firebrick2", midpoint = (max(limitss)+min(limitss))/2, name = "Difference (%)")+
  ggtitle("Sulfur transformation difference by site")+
  #coord_sf(xlim = c(-124, -70), ylim = c(25, 49), expand = FALSE)+
  theme_bw() + theme(axis.line=element_blank(),
                     axis.text.x=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks=element_blank(),
                     axis.title.x=element_blank(),
                     axis.title.y=element_blank(),
                     panel.background=element_blank(),
                     panel.border=element_blank(),
                     panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),
                     plot.background=element_blank())

##CHO transformation##
# Limits on ggplot scale
limitscho = c(min(trans_ll.us$CHO), 
           max(trans_ll.us$CHO))
#plot the data
CHO_plot <- ggplot(data = trans_ll.us)+
  borders("state", colour = "black", fill = "white")+ 
  geom_point(aes(x = Longitude, y = Latitude, color = CHO), size = 3.5) +
  scale_color_gradient2(limits = limitscho, low = "dodgerblue2",mid ="goldenrod2",
                        high = "firebrick2", midpoint = (max(limitscho)+min(limitscho))/2, name = "Difference (%)")+
  ggtitle("CHO transformation difference by site")+
  #coord_sf(xlim = c(-124, -70), ylim = c(25, 49), expand = FALSE)+
  theme_bw() + theme(axis.line=element_blank(),
                     axis.text.x=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks=element_blank(),
                     axis.title.x=element_blank(),
                     axis.title.y=element_blank(),
                     panel.background=element_blank(),
                     panel.border=element_blank(),
                     panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),
                     plot.background=element_blank())

##OH transformation##
# Limits on ggplot scale
limitsoh = c(min(trans_ll.us$OH), 
           max(trans_ll.us$OH))
#plot the data
OH_plot <- ggplot(data = trans_ll.us)+
  borders("state", colour = "black", fill = "white")+ 
  geom_point(aes(x = Longitude, y = Latitude, color = OH), size = 3.5) +
  scale_color_gradient2(limits = limitsoh, low = "dodgerblue2",mid ="goldenrod2",
                        high = "firebrick2", midpoint = (max(limitsoh)+min(limitsoh))/2, name = "Difference (%)")+
  ggtitle("OH transformation difference by site")+
  #coord_sf(xlim = c(-124, -70), ylim = c(25, 49), expand = FALSE)+
  theme_bw() + theme(axis.line=element_blank(),
                     axis.text.x=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks=element_blank(),
                     axis.title.x=element_blank(),
                     axis.title.y=element_blank(),
                     panel.background=element_blank(),
                     panel.border=element_blank(),
                     panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(),
                     plot.background=element_blank())

ggsave(aa_plot,filename="Amino_difference.png")

##Create bivariate plots for lat/ long separately with associated regression models ##
##### Install and load packages

trans_ll.us <- trans_ll.us %>% mutate(group = case_when(
              CHO > 5 & Longitude > -123 ~ 1,
              CHO > 5 & Longitude < -123 ~ 2,
              CHO < 5  ~ 2 ))
##phosphorus 
p_lat_biv <- ggplot(trans_ll.us, aes(y=P, x=Latitude)) +
  geom_point(aes(colour = factor(group)), size=4) +
  scale_colour_manual(values = c("orange","grey"))+
  geom_smooth(method=lm) +
  labs(y = "Relative abundance difference (%)", 
  title = "Phosphorus Latitude") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 40, label.y = 0.35) +
  stat_regline_equation(label.x = 40, label.y = 0.3)+
  theme_bw() + theme(plot.background=element_blank())+
  theme(legend.position = "none") 
#  ggsave(aa_lat_biv,filename="Amino_Lat_Regression.png")

p_long_biv <- ggplot(trans_ll.us, aes(y=P, x=Longitude)) +
  geom_point(size=4) +
  geom_point(aes(colour = factor(group)), size=4) +
  scale_colour_manual(values = c("orange","grey"))+
  geom_smooth(method=lm, color = "red") +
  labs(y = "Relative abundance difference (%)", 
       title = "Phosphorus Longitude") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = -90, label.y = 0.35) +
  stat_regline_equation(label.x = -90, label.y = 0.3)+
  theme_bw() + theme(plot.background=element_blank())+
  theme(legend.position = "none") 

library("gridExtra")
grid.arrange(p_lat_biv,p_long_biv,
             ncol = 2, nrow = 1)

##Amino acids
aa_lat_biv <- ggplot(trans_ll.us, aes(y=amino.acid, x=Latitude)) +
  geom_point(aes(colour = factor(group)), size=4) +
  scale_colour_manual(values = c("orange","grey"))+
  geom_smooth(method=lm) +
  labs(y = "Relative abundance difference (%)", 
       title = "Amino acids Latitude") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 40, label.y = 0.4) +
  stat_regline_equation(label.x = 40, label.y = 0.35)+
  theme_bw() + theme(plot.background=element_blank())+
  theme(legend.position = "none") 

aa_long_biv <- ggplot(trans_ll.us, aes(y=amino.acid, x=Longitude)) +
  geom_point(aes(colour = factor(group)), size=4) +
  scale_colour_manual(values = c("orange","grey"))+
  geom_smooth(method=lm, color = "red") +
  labs(y = "Relative abundance difference (%)", 
       title = "Amino acids Longitude") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = -90, label.y = 0.4) +
  stat_regline_equation(label.x = -90, label.y = 0.35)+
  theme_bw() + theme(plot.background=element_blank())+
  theme(legend.position = "none") 

  grid.arrange(aa_lat_biv,aa_long_biv,
               ncol = 2, nrow = 1)

  ##OH 
  oh_lat_biv <- ggplot(trans_ll.us, aes(y=OH, x=Latitude)) +
    geom_point(aes(colour = factor(group)), size=4) +
    scale_colour_manual(values = c("orange","grey"))+
    geom_smooth(method=lm) +
    labs(y = "Relative abundance difference (%)", 
         title = "OH Latitude") +
    stat_cor( aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
             label.x = 40, label.y = 1.25) +
    stat_regline_equation(label.x = 40, label.y = 1)+
    theme_bw() + theme(plot.background=element_blank())+
    theme(legend.position = "none") 
  
  oh_long_biv <- ggplot(trans_ll.us, aes(y=OH, x=Longitude)) +
    geom_point(aes(colour = factor(group)), size=4) +
    scale_colour_manual(values = c("orange","grey"))+
    geom_smooth(method=lm, color = "red") +
    labs(y = "Relative abundance difference (%)", 
         title = "OH Longitude") +
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
      label.x = -90, label.y = 1.25) +
    stat_regline_equation(label.x = -90, label.y = 1)+
    theme_bw() + theme(plot.background=element_blank())+
    theme(legend.position = "none") 
  
  library("gridExtra")
  grid.arrange(oh_lat_biv,oh_long_biv,
               ncol = 2, nrow = 1)
  
  ##CHO 
  cho_lat_biv <- ggplot(trans_ll.us, aes(y=CHO, x=Latitude)) +
    geom_point(aes(colour = factor(group)), size=4) +
    scale_colour_manual(values = c("orange","grey"))+
    geom_smooth(method=lm) +
    labs(y = "Relative abundance difference (%)", 
         title = "CHO Latitude") +
    stat_cor( aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
              label.x = 40, label.y = 15) +
    stat_regline_equation(label.x = 40, label.y = 14)+
    theme_bw() + theme(plot.background=element_blank())+
    theme(legend.position = "none") 
  
  cho_long_biv <- ggplot(trans_ll.us, aes(y=CHO, x=Longitude)) +
    geom_point(aes(colour = factor(group)), size=4) +
    scale_colour_manual(values = c("orange","grey"))+
    geom_smooth(method=lm, color = "red") +
    labs(y = "Relative abundance difference (%)", 
         title = "CHO Longitude") +
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
             label.x = -90, label.y = 15) +
    stat_regline_equation(label.x = -90, label.y = 14)+
    theme_bw() + theme(plot.background=element_blank())+
    theme(legend.position = "none") 
  
  library("gridExtra")
  grid.arrange(cho_lat_biv,cho_long_biv,
               ncol = 2, nrow = 1)
  
  ##N
  N_lat_biv <- ggplot(trans_ll.us, aes(y=N, x=Latitude)) +
    geom_point(aes(colour = factor(group)), size=4) +
    scale_colour_manual(values = c("orange","grey"))+
    geom_smooth(method=lm) +
    labs(y = "Relative abundance difference (%)", 
         title = "Nitrogen Latitude") +
    stat_cor( aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
              label.x = 40, label.y = 8) +
    stat_regline_equation(label.x = 40, label.y = 7.5)+
    theme_bw() + theme(plot.background=element_blank())+
    theme(legend.position = "none") 
  
  N_long_biv <- ggplot(trans_ll.us, aes(y=N, x=Longitude)) +
    geom_point(aes(colour = factor(group)), size=4) +
    scale_colour_manual(values = c("orange","grey"))+
    geom_smooth(method=lm, color = "red") +
    labs(y = "Relative abundance difference (%)", 
         title = "Nitrogen Longitude") +
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
             label.x = -90, label.y = 8) +
    stat_regline_equation(label.x = -90, label.y = 7.5)+
    theme_bw() + theme(plot.background=element_blank())+
    theme(legend.position = "none") 
  
  library("gridExtra")
  grid.arrange(N_lat_biv,N_long_biv,
               ncol = 2, nrow = 1)
  
  ##S
  S_lat_biv <- ggplot(trans_ll.us, aes(y=S, x=Latitude)) +
    geom_point(aes(colour = factor(group)), size=4) +
    scale_colour_manual(values = c("orange","grey"))+
    geom_smooth(method=lm) +
    labs(y = "Relative abundance difference (%)", 
         title = "Sulfur Latitude") +
    stat_cor( aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
              label.x = 40, label.y = 1.5) +
    stat_regline_equation(label.x = 40, label.y = 1.25)+
    theme_bw() + theme(plot.background=element_blank())+
    theme(legend.position = "none") 
  
  S_long_biv <- ggplot(trans_ll.us, aes(y=S, x=Longitude)) +
    geom_point(aes(colour = factor(group)), size=4) +
    scale_colour_manual(values = c("orange","grey"))+
    geom_smooth(method=lm, color = "red") +
    labs(y = "Relative abundance difference (%)", 
         title = "Sulfur Longitude") +
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
             label.x = -90, label.y = 1.5) +
    stat_regline_equation(label.x = -90, label.y = 1.25)+
    theme_bw() + theme(plot.background=element_blank())+
    theme(legend.position = "none") 
  
  library("gridExtra")
  grid.arrange(S_lat_biv,S_long_biv,
               ncol = 2, nrow = 1)