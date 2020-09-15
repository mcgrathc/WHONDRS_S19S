### Script to examine specific transformations

#### Data Pre-processing ####
# Load in libraries
library(reshape2) # For the melting capability
library(ggplot2) # For for prettiers graphs
library(ggthemes) # For different color palettes
library(ggpubr) # For combining plots (amongst other functions)

# Load in data
setwd("C:/Users/mcgr323/OneDrive - PNNL/Documents/WHONDRS/Transformation/")
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

amino.trans = melt(amino.trans, id.vars = "Type")

ggplot(data = amino.trans, aes(x = Type, y = value, group = Type))+
  geom_boxplot(aes(color = Type))+
  facet_wrap(~ variable, ncol = 3)+
  scale_color_stata()+
  hori_x_theme

# Bulk amino acid comparison
amino.trans = as.data.frame(t(trans[grep(paste(amino, collapse = "|"), row.names(trans), ignore.case = T),]))
amino.trans = as.data.frame(rowSums(amino.trans)); colnames(amino.trans) = "value"
amino.trans$Type = factors$Type
amino.trans$ID = factors$ID
amino.trans$Location = factors$Location

ggplot(data = amino.trans, aes(x = Type, y = value, group = Type))+
  geom_boxplot(aes(color = Type))+
  scale_color_stata()+
  hori_x_theme

# Looking at all N-based transformations
# N-based transformations (with AAs)
N = trans[grep("N|ine|phan|urea|biotinyl|co-enzyme|uracil|amin|adenyl|Aspart", row.names(trans)),]
to.remove = grep("N/A|glucose-N-|Na_", row.names(N))
N = N[-to.remove,]
CHO.remove = row.names(N)

N = data.frame(Sample = colnames(N), N = colSums(N), Type = factors$Type, ID = factors$ID, Location = factors$Location)

N.plot = ggplot(data = N, aes(x = Type, y = N))+
  geom_boxplot(aes(color = Type))+
  xlab(NULL)+
  ylab("N-based Trans. (%)")+
  hori_x_theme

N.plot

# S-transformations
S = trans[grep("S|Cysteine|Cystine|glutathi|Methionine|co-enzyme|biotinyl|sulfate", row.names(trans)),]
to.remove = grep("Serine", row.names(S))
S = S[-to.remove,]
CHO.remove = c(CHO.remove, row.names(S))

S = data.frame(Sample = colnames(S), S = colSums(S), Type = factors$Type)

S.plot = ggplot(data = S, aes(x = Type, y = S))+
  geom_boxplot(aes(color = Type))+
  xlab(NULL)+
  ylab("S-based Trans. (%)")+
  hori_x_theme

S.plot

# P-transformations
P = trans[grep("P|co-enzyme|phosphate|Phosphate|adenylate", row.names(trans)),]
to.remove = grep("2ndIP|Phenylalanine|Proline", row.names(P))
P = P[-to.remove,]
CHO.remove = c(CHO.remove, row.names(P))

P = data.frame(Sample = colnames(P), P = colSums(P), Type = factors$Type)

P.plot = ggplot(data = P, aes(x = Type, y = P))+
  geom_boxplot(aes(color = Type))+
  xlab(NULL)+
  ylab("P-based Trans. (%)")+
  hori_x_theme

P.plot

# Non-NSP transformations
CHO = trans
to.remove = which(row.names(CHO) %in% CHO.remove)
to.remove = c(to.remove, grep("Cl|Na|N/A", row.names(CHO))) # Removing transformations with chloride, sodium, or that are ambiguous
CHO = CHO[-to.remove,]
CHO = CHO[grep("C|A|a|e", row.names(CHO)),] # Selecting those transformtions which are carbon containing

CHO = data.frame(Sample = colnames(CHO), CHO = colSums(CHO), Type = factors$Type)

CHO.plot = ggplot(data = CHO, aes(x = Type, y = CHO))+
  geom_boxplot(aes(color = Type))+
  xlab(NULL)+
  ylab("Non-Nutrient Trans. (%)")+
  hori_x_theme

CHO.plot

# Hydrogen transformations
OH = trans
OH = OH[-to.remove,]
OH = OH[-grep("C|A|a|e|P", row.names(OH)),] # Selecting those transformtions which are O/H, but not C containing

OH = data.frame(Sample = colnames(OH), OH = colSums(OH), Type = factors$Type)

OH.plot = ggplot(data = OH, aes(x = Type, y = OH))+
  geom_boxplot(aes(color = Type))+
  xlab(NULL)+
  ylab("O/H Trans. (%)")+
  hori_x_theme

OH.plot

# Combining plots
ggarrange(CHO.plot, N.plot, S.plot, P.plot)

### Calculating stats
stats = wilcox.test(CHO~Type, data = CHO)
stats = data.frame(Comparison = "CHO", MWU.stat = stats$statistic, p.value = stats$p.value)

temp = wilcox.test(N~Type, data = N)
stats = rbind(stats, 
              data.frame(Comparison = "N", MWU.stat = temp$statistic, p.value = temp$p.value))

temp = wilcox.test(S~Type, data = S)
stats = rbind(stats, 
              data.frame(Comparison = "S", MWU.stat = temp$statistic, p.value = temp$p.value))

temp = wilcox.test(P~Type, data = P)
stats = rbind(stats, 
              data.frame(Comparison = "P", MWU.stat = temp$statistic, p.value = temp$p.value))

stats$p.value = p.adjust(stats$p.value, method = "fdr")
stats$p.value = round(stats$p.value, digits = 4)
