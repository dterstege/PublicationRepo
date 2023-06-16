#Brain Analysis Pipeline - Coronal Sections
#
#Semi-Automated Pipeline for the Segmentation and Registration of Brain-Wide Labels
#
#Requires WholeBrain software package (http://www.wholebrainsoftware.org/cms/install/)
#Ilastik segmentation and further image preprocessing in ImageJ are recommended
#
#Read through commented documentation prior to running the script
#
#Requires initial setting of a registration filter, but should remain consistent across image set
#   Tips:
#   Can test a filter by running segment(testimage, filter=testfilter, display = TRUE)
#     Adjust sliders as needed to yield a nice outline of the brain section
#     Once pleased with outline, hit escape key
#     Many values will flood the Console. Identify any with the prefix "$filter$" (ex. $filter$brain.threshold)
#     Make note of these values and enter them into your registration filter
#
#DO NOT NEED TO ADJUST THE SEGMENTATION FILTER FROM IMAGE TO IMAGE
#Outputs named using registration image name
#
#WholeBrain - Daniel Furth. 2015. WholeBrain: a Neuroanatomical Information System, Beta Version 0.0.1
#
#Created: 01-03-2021 Dylan Terstege
#         Epp Lab, University of Calgary
#
#Contact: dylan.terstege@ucalgary.ca

#Enter Allen Atlas Coordinate
#Data entry format: Images named "Test01_00.tif" (animalNum_imNum.tif)
animalNum<-"Test01"   #animal ID
imNum<-"00"           #image identifier
coord<- 3.145         #Allen Atlas Coordinate

#Load Packages
library(wholebrain)

#Load Images
#adjust accordingly
#Data entry format: parent folder with the same name as the animal ID
#                   subfolders named after the types of images (dapi/registration, label to be registered, and point mask)
regim<-paste(sep="","/Volumes/sampledir/",animalNum,"/regi/",animalNum,"_",imNum,".tif")
segim<-paste(sep="","/Volumes/sampledir/",animalNum,"/ilastik_thresholded/",animalNum,"A_",imNum,".tif")
maskim<-paste(sep="","/Volumes/sampledir/",animalNum,"/mask/",animalNum,"_",imNum,"_mask.tif")
outputDir<-paste(sep="","/Volumes/sampledir/outputs/",animalNum,"/")
setwd(outputDir)

#Segmentation
segfilter<- structure(list(alim = c(0, 49), threshold.range = c(0, 225), eccentricity = 1000, Max = 4035, Min = 2252, brain.threshold = 1801, resize = 0.2, blur = 5, downsample = 0.25), .Names = c("alim", "threshold.range", "eccentricity", "Max", "Min", "brain.threshold", "resize", "blur", "downsample"))
labelseg<- segment(segim, filter = segfilter, display = FALSE)
maskseg<- segment(maskim, filter = segfilter, display = FALSE) #Will likely be a time consuming step

#Registration
#   Adjust "regfilter" as needed
quartz()
regfilter<-structure(list(alim = c(1, 1000), threshold.range = c(65535, 65536), eccentricity = 1000, Max = 100, Min = 0, brain.threshold = 10, resize = 0.14, blur =5, downsample = 0.1), .Names = c("alim", "threshold.range", "eccentricity", "Max", "Min", "brain.threshold", "resize", "blur", "downsample"))
regi<-registration(regim, coordinate=coord, filter = regfilter, display=TRUE)
#   Run to this line: Command+Option+B
#   The following 3 lines of code are to be used to add/change points in the registration
#   Run these individually and repeat as many times as needed
regi<-change.corrpoints(regi, 1:32)
regi<-add.corrpoints(regi, 100)
regi<-registration(regim, coordinate = coord, filter = regfilter, correspondance = regi)
#  Run to end: Command+Option+E
#Putting It All Together and Generating Output Files
label_dataset<-inspect.registration(regi, labelseg, draw.trans.grid = FALSE)
mask_dataset<-inspect.registration(regi, maskseg, draw.trans.grid = FALSE)
id<-sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(regim)) #using registration image to name outputs
gbfile=paste(sep="",id,".Rdata")
save(labelseg, regi, label_dataset, file = gbfile)
label<-table(label_dataset$acronym)
mask<-table(mask_dataset$acronym)
fullcells=paste(sep="",id,"_cells.csv")
write.csv(label, file = fullcells)
fullgrid=paste(sep="",id,"_grid.csv")
write.csv(mask, file = fullgrid)
end_time<-Sys.time()
runtime=end_time-start_time

graphics.off()
rm(list=ls())
cat("\014")