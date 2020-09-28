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

amino.trans = as.data.frame(trans[grep(paste(amino, collapse = "|"), row.names(trans), ignore.case = T),])
amino_acid = data.frame(Sample = colnames(amino.trans), value = colSums(amino.trans), Type = factors$Type, ID = factors$ID)
amino_acid$ID= as.numeric(as.character(amino_acid$ID))

for (i in amino_acid$ID) { 
  sub_data <- subset(amino_acid, ID == i)
  AA.plot = ggplot(data = sub_data, aes(x = Type, y = value))+
    geom_boxplot(aes(color = Type))+
    xlab(NULL)+
    scale_y_continuous(limits=c(0,1))+
    ylab("Amino Acid Trans. (%)")+
    hori_x_theme+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(AA.plot,filename=paste("Amino_Site",i,".png",sep=""))
}

# Looking at all N-based transformations
# N-based transformations (with AAs)
N.trans = as.data.frame(trans[grep("N|ine|phan|urea|biotinyl|co-enzyme|uracil|amin|adenyl|Aspart", 
                               row.names(trans), ignore.case = T),])

to.remove = grep("N/A|glucose-N-|Na_", row.names(N.trans))
N.trans = N.trans[-to.remove,]
CHO.remove = row.names(N.trans)

N = data.frame(Sample = colnames(N.trans), value = colSums(N.trans), Type = factors$Type, ID = factors$ID)
N$ID= as.numeric(as.character(N$ID))

for (i in N$ID) { 
  sub_data <- subset(N, ID == i)
  N.plot = ggplot(data = sub_data, aes(x = Type, y = value))+
    geom_boxplot(aes(color = Type))+
    xlab(NULL)+
    scale_y_continuous(limits=c(4,28))+
    ylab("N Trans. (%)")+
    hori_x_theme+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(N.plot,filename=paste("N_Site",i,".png",sep=""))
}

# S-transformations
S.trans = as.data.frame(trans[grep("S|Cysteine|Cystine|glutathi|Methionine|co-enzyme|biotinyl|sulfate", 
                                     row.names(trans), ignore.case = T),])
to.remove = grep("Serine", row.names(S.trans))
S.trans = S.trans[-to.remove,]
CHO.remove = c(CHO.remove, row.names(S.trans))

S = data.frame(Sample = colnames(S.trans), value = colSums(S.trans), Type = factors$Type, ID = factors$ID)
S$ID= as.numeric(as.character(S$ID))

for (i in S$ID) { 
  sub_data <- subset(S, ID == i)
  S.plot = ggplot(data = sub_data, aes(x = Type, y = value))+
    geom_boxplot(aes(color = Type))+
    xlab(NULL)+
    scale_y_continuous(limits=c(0,5))+
    ylab("S Trans. (%)")+
    hori_x_theme+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(S.plot,filename=paste("S_Site",i,".png",sep=""))
}

# P-transformations
P.trans = as.data.frame(trans[grep("P|co-enzyme|phosphate|Phosphate|adenylate", 
                                     row.names(trans), ignore.case = T),])

to.remove = grep("2ndIP|Phenylalanine|Proline", row.names(P.trans))
P.trans = P.trans[-to.remove,]
CHO.remove = c(CHO.remove, row.names(P.trans))

P = data.frame(Sample = colnames(P.trans), value = colSums(P.trans), Type = factors$Type, ID = factors$ID)
P$ID= as.numeric(as.character(P$ID))

for (i in P$ID) { 
  sub_data <- subset(P, ID == i)
  P.plot = ggplot(data = sub_data, aes(x = Type, y = value))+
    geom_boxplot(aes(color = Type))+
    xlab(NULL)+
    scale_y_continuous(limits=c(0,3))+
    ylab("P Trans. (%)")+
    hori_x_theme+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(P.plot,filename=paste("P_Site",i,".png",sep=""))
}

# Non-NSP transformations
CHO = trans
to.remove = which(row.names(CHO) %in% CHO.remove)
to.remove = c(to.remove, grep("Cl|Na|N/A", row.names(CHO))) # Removing transformations with chloride, sodium, or that are ambiguous
CHO = CHO[-to.remove,]
CHO.trans = as.data.frame(CHO[grep("C|A|a|e", row.names(CHO)),]) # Selecting those transformtions which are carbon containing
CHO = data.frame(Sample = colnames(CHO.trans), value = colSums(CHO.trans), Type = factors$Type, ID = factors$ID)
CHO$ID= as.numeric(as.character(CHO$ID))

for (i in CHO$ID) { 
  sub_data <- subset(CHO, ID == i)
  CHO.plot = ggplot(data = sub_data, aes(x = Type, y = value))+
    geom_boxplot(aes(color = Type))+
    xlab(NULL)+
    scale_y_continuous(limits=c(50,85))+
    ylab("CHO Trans. (%)")+
    hori_x_theme+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(CHO.plot,filename=paste("CHO_Site",i,".png",sep=""))
}


# Hydrogen transformations
OH = trans
OH = OH[-to.remove,]
OH.trans = as.data.frame(OH[-grep("C|A|a|e|P", row.names(OH)),]) # Selecting those transformtions which are O/H, but not C containing
OH = data.frame(Sample = colnames(OH.trans), value = colSums(OH.trans), Type = factors$Type, ID = factors$ID)
OH$ID= as.numeric(as.character(OH$ID))

for (i in OH$ID) { 
  sub_data <- subset(OH, ID == i)
OH.plot = ggplot(data = sub_data, aes(x = Type, y = value))+
  geom_boxplot(aes(color = Type))+
  xlab(NULL)+
  scale_y_continuous(limits=c(9,30))+
  ylab("O/H Trans. (%)")+
  hori_x_theme+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(OH.plot,filename=paste("OH_Site",i,".png",sep=""))
}


site.graph <- function(df, na.rm = TRUE, ...){
  site_list <- unique(df$ID)
  setwd("C:/Users/mcgr323/OneDrive - PNNL/Documents/WHONDRS/Transformation/")
  # create for loop to produce ggplot2 graphs 
  for (i in seq_along(site_list)) { 
    sub_data <- subset(df, ID == "i")
    # create plot for each county in df 
    site_plots <- ggplot(data = df, aes(x = Type, y = value))+
      geom_boxplot(aes(color = Type))+
      labs(x="Type", y = "Relative Abundance (%)")+
      scale_color_stata()+
      ggtitle("Site", paste( site_list[i]))+
      hori_x_theme
    ggsave(site_plots,filename=paste("Site",site_list[i],".png",sep=""))
  } 
}


  site_list <- unique(df$ID)
  setwd("C:/Users/mcgr323/OneDrive - PNNL/Documents/WHONDRS/Transformation/")
  # create for loop to produce ggplot2 graphs 
  for (i in seq_along(site_list)) { 
    # create plot for each county in df 
    sub_data <- subset(df, ID == "i")
    site_plots <- ggplot(data = sub_data, aes(x = Type, y = value))+
      geom_boxplot(aes(color = Type))+
      labs(x="Sample Type", y = "Relative Abundance")+
      scale_color_stata()+
      ggtitle("Site", paste( site_list[i]))+
      hori_x_theme
    ggsave(site_plots,filename=paste("Site",site_list[i],".png",sep=""))
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
