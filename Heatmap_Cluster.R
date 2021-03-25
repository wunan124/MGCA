#########################################################
### A) Installing and loading required packages
#########################################################

if (!require("gplots")) {
   install.packages("gplots", dependencies = TRUE)
   library(gplots)
   }
if (!require("RColorBrewer")) {
   install.packages("RColorBrewer", dependencies = TRUE)
   library(RColorBrewer)
   }


#########################################################
### B) Reading in data and transform it into matrix format
#########################################################

data <- read.csv ("~/R-Scripts/datasets/Heatmap.csv", comment.char="#")
rnames <- data[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames                  # assign row names


#########################################################
### C) Customizing and plotting the heat map
#########################################################

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("blue", "black", "red"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(-1,-0.25,length=100),  # for red
  seq(-0.24,0.24,length=100),           # for purple
  seq(0.25,1,length=100))             # for blue

# creates a 5 x 5 inch image
png("~/R-Scripts/images/heatmaps_in_r.png",    # create PNG for the heat map        
  width = 5*300,        # 5 x 300 pixels
  height = 5*300,
  res = 300,            # 300 pixels per inch
  pointsize = 5)        # smaller font size

  heatmap.2(mat_data,
  #cellnote = "none",  # same data set for cell labels
  main = "Cluster", # heat map title
  #notecex=NULL,
  
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier
  breaks=col_breaks,    # enable color transition at specified limits
  notecol="black",      # change font color of cell labels to black
  dendrogram="row",     # only draw a row dendrogram
  #Rowv="NA",		# turn off row clustering
  Colv="NA")            # turn off column clustering

dev.off()               # close the PNG device
