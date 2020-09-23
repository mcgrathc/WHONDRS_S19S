### Script to examine specific transformations

#### Data Pre-processing ####
# Load in libraries
library(reshape2) # For the melting capability
library(ggplot2) # For for prettiers graphs
library(ggthemes) # For different color palettes
library(ggpubr) # For combining plots (amongst other functions)
library(dplyr)

# Load in data
setwd("C:/Users/mcgr323/OneDrive - PNNL/Documents/WHONDRS/Transformation/")
results <- "C:/Users/mcgr323/OneDrive - PNNL/Documents/WHONDRS/Transformation/"
trans = read.csv("S19S_Wat-Sed_8.12_Trans_Profiles.csv", row.names = 1)

# Cleaning up the transformation profile
trans = trans[,-which(colnames(trans) %in% "Mass")]

# Creating factor sheet
#create dataframe of unique sample names 
unique_sample_names <- as.data.frame(unique(names(trans)))

#spilt sample names by "_"
temp <- strsplit(unique_sample_names[,1], '_')

#create dataframe of separated values  
max_length <- max(sapply(temp,length))
temp2 <- as.data.frame(t(sapply(temp, function(x){
  c(x, rep(NA, max_length - length(x)))
})))

#make new columns for spatial location, sample type, and site ID 
factors <- temp2 %>% 
  mutate(Location = case_when(temp2[,4] == "ICR.1" ~ 1,
                                      temp2[,4] == "ICR.2" ~ 2,
                                      temp2[,4] == "ICR.3" ~ 3,
                                      temp2[,6] == "ICR.U" ~ 4, #upstream
                                      temp2[,6] == "ICR.M" ~ 5, #midstream
                                      temp2[,6] == "ICR.D" ~ 6),#downstream
         Type = ifelse(temp2[,4] == "Sed", "sediment", "surface water"),
         ID = temp2[,3])

factors <- cbind(unique_sample_names[,1], factors[,8:10])


# Creating ggplot themes
hori_x_theme = theme_bw()+
  theme(text = element_text(size = 14),
        axis.title = element_text(colour = "black", size = 15),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        panel.border = element_rect(size = 1, colour = "black"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")


#### Analysis Component ####
# Converting to rel. abund.
trans = as.data.frame(apply(trans, 2, function(x) x/sum(x)))
trans = trans*100

# Looking at all amino acids
amino = c("alanine", "arginine", "asparagine", "aspartic acid",
          "cyteine", "glutamic acid", "glutamine", "glycine",
          "histidine", "isoleucine", "lysine",
          "methoionine", "phenylalanine", "proline",
          "serine", "theronine", "tryptophan",
          "tyrosine", "valine")

amino.trans = as.data.frame(t(trans[grep(paste(amino, collapse = "|"), row.names(trans), ignore.case = T),]))
amino.trans$Type = factors$Type
amino.trans$ID = factors$ID
amino.trans$Location = factors$Location

amino.trans = melt(amino.trans, id.vars = c("Type", "ID"))

amino_mean <- amino.trans %>%
  group_by(ID, Type) %>% 
  summarise_each(funs(mean))

# # Bulk amino acid comparison
# amino.trans = as.data.frame(t(trans[grep(paste(amino, collapse = "|"), row.names(trans), ignore.case = T),]))
# amino.trans = as.data.frame(rowSums(amino.trans)); colnames(amino.trans) = "value"
# amino.trans$Type = factors$Type
# amino.trans$ID = factors$ID
# amino.trans$Location = factors$Location

# Looking at all N-based transformations
# N-based transformations (with AAs)
N.trans = as.data.frame(t(trans[grep("N|ine|phan|urea|biotinyl|co-enzyme|uracil|amin|adenyl|Aspart", 
                               row.names(trans), ignore.case = T),]))
N.trans$Type = factors$Type
N.trans$ID = factors$ID
N.trans$Location = factors$Location

to.remove = grep("N/A|glucose-N-|Na_", colnames(N.trans))
N.trans = N.trans[-to.remove,]
CHO.remove = row.names(N.trans)

N.trans = melt(N.trans, id.vars = c("Type", "ID"))

N_mean <- N.trans %>%
  group_by(ID, Type) %>% 
  summarise_each(funs(mean))

# S-transformations
S.trans = as.data.frame(t(trans[grep("S|Cysteine|Cystine|glutathi|Methionine|co-enzyme|biotinyl|sulfate", 
                                     row.names(trans), ignore.case = T),]))
S.trans$Type = factors$Type
S.trans$ID = factors$ID
S.trans$Location = factors$Location

to.remove = grep("Serine", colnames(S.trans))
S.trans = S.trans[-to.remove,]
CHO.remove = c(CHO.remove, colnames(S.trans))

S.trans = melt(S.trans, id.vars = c("Type", "ID"))

S_mean <- S.trans %>%
  group_by(ID, Type) %>% 
  summarise_each(funs(mean))

# P-transformations
P.trans = as.data.frame(t(trans[grep("P|co-enzyme|phosphate|Phosphate|adenylate", 
                                     row.names(trans), ignore.case = T),]))
P.trans$Type = factors$Type
P.trans$ID = factors$ID
P.trans$Location = factors$Location
to.remove = grep("2ndIP|Phenylalanine|Proline", colnames(P.trans))
P.trans = P.trans[-to.remove,]
CHO.remove = c(CHO.remove, colnames(P.trans))

P.trans = melt(P.trans, id.vars = c("Type", "ID"))

P_mean <- P.trans %>%
  group_by(ID, Type) %>% 
  summarise_each(funs(mean))

# Non-NSP transformations
CHO = trans
to.remove = which(row.names(CHO) %in% CHO.remove)
to.remove = c(to.remove, grep("Cl|Na|N/A", row.names(CHO))) # Removing transformations with chloride, sodium, or that are ambiguous
CHO = CHO[-to.remove,]
CHO.trans = as.data.frame(t(CHO[grep("C|A|a|e", row.names(CHO)),])) # Selecting those transformtions which are carbon containing
CHO.trans$Type = factors$Type
CHO.trans$ID = factors$ID
CHO.trans$Location = factors$Location

CHO.trans = melt(CHO.trans, id.vars = c("Type", "ID"))

CHO_mean <- CHO.trans %>%
  group_by(ID, Type) %>% 
  summarise_each(funs(mean))

# Hydrogen transformations
OH = trans
OH = OH[-to.remove,]
OH.trans = as.data.frame(t(OH[-grep("C|A|a|e|P", row.names(OH)),])) # Selecting those transformtions which are O/H, but not C containing
OH.trans$Type = factors$Type
OH.trans$ID = factors$ID
OH.trans$Location = factors$Location

OH.trans = melt(OH.trans, id.vars = c("Type", "ID"))

OH_mean <- OH.trans %>%
  group_by(ID, Type) %>% 
  summarise_each(funs(mean))


site.graph <- function(df, na.rm = TRUE, ...){
  site_list <- unique(df$ID)
  setwd("C:/Users/mcgr323/OneDrive - PNNL/Documents/WHONDRS/Transformation/")
  # create for loop to produce ggplot2 graphs 
  for (i in seq_along(site_list)) { 
    # create plot for each county in df 
    site_plots <- ggplot(data = df, aes(x = Type, y = value))+
      geom_boxplot(aes(color = Type))+
      labs(x="Sample Type", y = "Relative Abundance")+
      scale_color_stata()+
      ggtitle("Site", paste( site_list[i]))+
      hori_x_theme
    ggsave(site_plots,filename=paste("Site",site_list[i],".png",sep=""))
  } 
}


# # Combining plots
# ggarrange(CHO.plot, N.plot, S.plot, P.plot)
# 
# ### Calculating stats
# stats = wilcox.test(CHO~Type, data = CHO)
# stats = data.frame(Comparison = "CHO", MWU.stat = stats$statistic, p.value = stats$p.value)
# 
# temp = wilcox.test(N~Type, data = N)
# stats = rbind(stats, 
#               data.frame(Comparison = "N", MWU.stat = temp$statistic, p.value = temp$p.value))
# 
# temp = wilcox.test(S~Type, data = S)
# stats = rbind(stats, 
#               data.frame(Comparison = "S", MWU.stat = temp$statistic, p.value = temp$p.value))
# 
# temp = wilcox.test(P~Type, data = P)
# stats = rbind(stats, 
#               data.frame(Comparison = "P", MWU.stat = temp$statistic, p.value = temp$p.value))
# 
# stats$p.value = p.adjust(stats$p.value, method = "fdr")
# stats$p.value = round(stats$p.value, digits = 4)
