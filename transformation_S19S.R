#set working directory
setwd("C:\\Users\\mcgr323\\OneDrive - PNNL\\Documents\\WHONDRS\\Transformation")

#set workspace path
workspace_path <- "C:\\Users\\mcgr323\\OneDrive - PNNL\\Documents\\WHONDRS\\Transformation"

#read in csv
trans_matrix <- read.csv(file.path(workspace_path,"S19S_Wat-Sed_8.12_Trans_Profiles.csv"))

#count the unique sample names in the trans matrix
unique_sample_names <- length(unique(names(trans_matrix[3:ncol(trans_matrix)])))

rel_abund <- data.frame(matrix(NA, nrow = nrow(trans_matrix), ncol = (ncol(trans_matrix)-2)))

#get the relative abundance of the transformation
rel_abund <-as.data.frame(sapply(names(trans_matrix[3:ncol(trans_matrix)]), function(x) {
  trans_matrix[paste0(x, "_rel_abud")] <<- trans_matrix[x] / sum(trans_matrix[x])
}))

#add two columns (transformation name, mass) to relative abunance 
trans_rel_abun <- cbind(trans_matrix[1:2],rel_abund)

