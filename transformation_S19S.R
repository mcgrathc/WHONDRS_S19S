#set working directory
setwd("C:\\Users\\mcgr323\\OneDrive - PNNL\\Documents\\WHONDRS\\Transformation")

#set workspace path
workspace_path <- "C:\\Users\\mcgr323\\OneDrive - PNNL\\Documents\\WHONDRS\\Transformation"

#read in csv
trans_matrix <- read.csv(file.path(workspace_path,"S19S_Wat-Sed_8.12_Trans_Profiles.csv"))

#create an empty dataframe for the relative abundances 
rel_abund <- data.frame(matrix(NA, nrow = nrow(trans_matrix), ncol = (ncol(trans_matrix)-2)))

#get the relative abundance of the transformation
rel_abund <-as.data.frame(sapply(names(trans_matrix[3:ncol(trans_matrix)]), function(x) {
  trans_matrix[x] / sum(trans_matrix[x])
}))

#add two columns (transformation name, mass) to relative abundance 
trans_rel_abun <- cbind(trans_matrix[1:2],rel_abund)

#get rid of the duplicate column names from sapply function 
names(trans_rel_abun) <- names(trans_matrix)

#create dataframe of unique sample names 
unique_sample_names <- as.data.frame(unique(names(trans_matrix[3:ncol(trans_matrix)])))

#spilt sample names by "_"
temp <- strsplit(unique_sample_names[,1], '_')

#create dataframe of separated values  
max_length <- max(sapply(temp,length))
temp2 <- as.data.frame(t(sapply(temp, function(x){
  c(x, rep(NA, max_length - length(x)))
})))

#make new columns for spatial location, sample type, and site ID 
dat.extended <- temp2 %>% 
  mutate(Spatial_location = case_when(temp2[,4] == "ICR.1" ~ 1,
                          temp2[,4] == "ICR.2" ~ 2,
                          temp2[,4] == "ICR.3" ~ 3,
                          temp2[,6] == "ICR.U" ~ 4, #upstream
                          temp2[,6] == "ICR.M" ~ 5, #midstream
                          temp2[,6] == "ICR.D" ~ 6),#downstream
         Sample_type = ifelse(temp2[,4] == "Sed", "sediment", "surface water"),
         Site_ID = temp2[,3])

sample_desc <- cbind(unique_sample_names[,1], dat.extended[,8:10])


