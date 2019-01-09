##################################################################################################################################
#Written for Proteomics Capstone Group Project
#University of Maryland, University College
#Author: Justin Hughes
#Date Modified: 11/25/2018

#################################################################################################################################
# Load Required Libraries
#################################################################################################################################
library(rpart)
library(rlist)
#################################################################################################################################

#################################################################################################################################

#Import File
peptideFile <- "Capstone Data - Proteomics.csv"
#proteinFile <- "Protein Data Table CSV.csv"


#Read the CSV file into a dataset, note that there are headers in the file, and set the row names to the first column (Study ID)
peptideSet <- read.csv(peptideFile, header=TRUE,  dec=".",  strip.white=TRUE, row.names = c(1))
#proteinSet <- read.csv(proteinFile, header=TRUE,  dec=".",  strip.white=TRUE, row.names = c(1))

#If Protein Data used, reads the Protein List file in and creates a subset of those variables
#proteins <- read.csv (listFile, header = FALSE, na.strings = "NA", dec = ".", strip.white = TRUE)
#proteinList <- unlist(proteins, use.names = TRUE)
#ProteinData <- proteinSet[proteinList]


#Create Subset with No Pseudo2
peptideNoP2 <- subset(peptideSet, Original.Group.Name != "PSEUDO2")

#Create Subset with Control, Asymptomatic, and Mild conditions
peptideCAM <- subset(peptideSet, Original.Group.Name == "CONTROL" | Original.Group.Name == "ASYMPTOMATIC" | Original.Group.Name == "MILD")

#Create Subset with just Pseudo and Severe
peptidePS <- subset(peptideSet, Original.Group.Name == "PSEUDO" | Original.Group.Name == "SEVERE")


##########################################################################################################################################################
#Create Trees and pull peptides from each subset
##########################################################################################################################################################
# For Peptide Dataset, Remove low-frequency spectral counts (anything lower than 5 is considered noise)) - assumes we are using non-normalized data
#peptideData[peptideData <= 1] <- 0 # removing low freq counts on the peptides dataset renders most columns nonusable

#Create the Classification Tree Using the "Class" method
rpart.tree <- rpart(Original.Group.Name ~ ., data = peptideSet, method = "class",  control = rpart.control (cp=0.001))

#Create list with important peptide names
output <- attributes (rpart.tree$variable.importance) 
write.table(output, "Peptide Test Output.csv", append = FALSE, row.names = FALSE, col.names = FALSE)

#Repeat setting any counts <= 1 as 0 to reduce noise
peptideSet[peptideSet <= 1] <- 0
rpart.tree <- rpart(Original.Group.Name ~ ., data = peptideSet, method = "class",  control = rpart.control (cp=0.001))
output <- attributes (rpart.tree$variable.importance) 
write.table(output, "Peptide Test Output.csv", append = TRUE, row.names = FALSE, col.names = FALSE)


#No Pseudo 2 (NOP2)
rpart.tree <- rpart(Original.Group.Name ~ ., data = peptideNoP2, method = "class",  control = rpart.control (cp=0.001))
output <- attributes (rpart.tree$variable.importance)
write.table(output, "Peptide Test Output.csv", append = TRUE, row.names = FALSE, col.names = FALSE)

#NOP2 Filtering
#peptideNoP2[peptideNoP2 <= 1] <- 0
#rpart.tree <- rpart(Original.Group.Name ~ ., data = peptideNoP2, method = "class",  control = rpart.control (cp=0.001))
#output <- attributes (rpart.tree$variable.importance) 
#write.table(output, "Peptide Test Output.csv", append = TRUE, row.names = FALSE, col.names = FALSE)


#Control, Asymptomatic, Mild (CAM)
rpart.tree <- rpart(Original.Group.Name ~ ., data = peptideCAM, method = "class",  control = rpart.control (cp=0.001))
output <- attributes (rpart.tree$variable.importance)
write.table(output, "Peptide Test Output.csv", append = TRUE, row.names = FALSE, col.names = FALSE)

#CAM Filtering
peptideCAM[peptideCAM <= 1] <- 0
rpart.tree <- rpart(Original.Group.Name ~ ., data = peptideCAM, method = "class",  control = rpart.control (cp=0.001))
summary(rpart.tree)
output <- attributes (rpart.tree$variable.importance) 
write.table(output, "Peptide Test Output.csv", append = TRUE, row.names = FALSE, col.names = FALSE)


#Pseudo and Severe (PS)
rpart.tree <- rpart(Original.Group.Name ~ ., data = peptidePS, method = "class",  control = rpart.control (cp=0.001))
output <- attributes (rpart.tree$variable.importance)
write.table(output, "Peptide Test Output.csv", append = TRUE, row.names = FALSE, col.names = FALSE)

#PS Filtering
#peptidePS[peptidePS <= 1] <- 0
#rpart.tree <- rpart(Original.Group.Name ~ ., data = peptidePS, method = "class",  control = rpart.control (cp=0.001))
#output <- attributes (rpart.tree$variable.importance) 
#write.table(output, "Peptide Test Output.csv", append = TRUE, row.names = FALSE, col.names = FALSE)


###########################################################################################################################################################
#Remove duplicates and write file for final peptide list
##########################################################################################################################################################
#If Peptide Data Used, reads the Peptide List file in and creates a subset of those variables
peptideListFile = "Peptide Test Output.csv"
peptides <- read.csv (peptideListFile, header = FALSE, na.strings = "NA", dec = ".", strip.white = TRUE)
peptideList <- unlist(peptides, use.names = TRUE)

#Removing any duplicate peptide entries
peptideList <- unique (peptideList, incomparables = FALSE)

#Writing the final Peptide List to a CSV file for use in LDA Script
write.table(peptideList, "Final List 2.csv", append = FALSE, row.names = FALSE, col.names = FALSE)

peptideData <- peptideSet[peptideList]
##################################################################################################################################
#Final Tree Pruning
#################################################################################################################################
#Create the tree using the final Peptide Dataset
rpart.tree <- rpart(Original.Group.Name ~ ., data = peptideData, method = "class",  control = rpart.control (cp=0.001))

#Plot unpruned tree and set the branch lengths and aesthetic parameters
plot(rpart.tree, uniform=TRUE, branch=0.2, margin=0.05)
text(rpart.tree, all=TRUE, use.n=TRUE)
title("Final Peptide Classification Tree Unpruned")


#Plot the cp values and print it out
plotcp(rpart.tree)
printcp(rpart.tree)


#Prune the tree based on the values above
rpart.tree2 = prune(rpart.tree, cp = 0.01)

plot(rpart.tree2, uniform = TRUE, margin = 0.05)
text (rpart.tree2, use.n = TRUE, cex = 0.75)
title("Final Peptide Classification Tree - Pruned")

summary (rpart.tree2)
printcp (rpart.tree2)
