##
##
## ms-mice.R: A script to analyze Dorothy Madamba (Bucheli Lab) dead mouse VOC GC/MS data
##
##

#Note: Change working directory to your own working directory
working.directory = getwd() # or a specific directory like "C:/Users/DCH009/OneDrive - Sam Houston State University/docs/Research/Rats/AIAEXPRT.AIA"
if (getwd()!=working.directory)
{
  print("Changing working directory to script defined working directory of ")
  print(working.directory)
  setwd(working.directory)
  
}


if (!exists("ionlist.alkyl")) ionlist.alkyl=c(43,57,71,85,99,113,127,141,155,169)
if (!exists("colors.alkyl")) colors.alkyl=rainbow(length(ionlist.alkyl))
computer.host<-Sys.info()["nodename"]

# get user's My Documents (or other Windows home) directory
home.directory = path.expand('~')

#Note: You may need to edit the lines below to make sure the script can find these two files 
R.directory = working.directory # or something lik epaste0(home.directory, "/Git/R/")
if (!exists("dch_general_functions.loaded")) source(paste0(R.directory, "HainesLabFunctions/dch_general_functions.R"))
if (!exists("ms.functions.loaded")) source(paste0(R.directory,"/HainesLabMS/ms-functions.R"))


########################### Program Proper ##################################
#############################################################################

pdf(paste("Mouse Decomp Results -",format(Sys.time(), "%Y%m%d-%H%M%S"),".PDF",sep=""), width=8.5, height=11)


if (loadFiles || !exists("spectra")) spectra<-ms.load() #you can set loadFiles=FALSE to 
                                                        # skip loading the files if 
                                                        # spectra object exists already

tryCatch({decomp.ions<-read.csv("compounds.csv",header=TRUE,stringsAsFactors=FALSE)}, 
   warning=function(w){decomp.ions=NULL},
   error=function(e){decomp.ions=NULL})

# Define some useful groups of spectra
order1=c(3,10,5,4,8,9,1,6,7,2) #groups together
order2=c(10,5,8,9,6,7)  #phases together
order.last.samples=c(7,9,5)
order.p3a=c(6,8,10)
order.mask=c(FALSE,TRUE,TRUE,FALSE,TRUE,TRUE,FALSE,TRUE, TRUE,FALSE)

############################################################
## Make some overview plots
############################################################

par(mar = c(4, 4, 4, 4) + 0.2)

ionlist.temp<-sort(unique(c(stockionlist, decomp.ions[,"Ion.1"], decomp.ions[,"Ion.2"],decomp.ions[,"Ion.3"], 
                            decomp.ions[,"Ion.4"], decomp.ions[,"Ion.5"])))
ms.stackedplot(spectra, ionlist=ionlist.temp, which.spectra=order1, main="All Samples", showfilenames = TRUE)
ms.stackedplot(spectra, ionlist=ionlist.alkyl, which.spectr=order1, main="All Samples: Alkyl series", 
               eiccolors=colors.alkyl, ticcolor="grey", showfilenames = TRUE)


ms.stackedplot(spectra, ionlist=ionlist.temp, which.spectra=order2, main="All Samples", showfilenames = TRUE)
ms.stackedplot(spectra, ionlist=ionlist.alkyl, which.spectr=order2, main="All Samples: Alkyl series", 
               eiccolors=colors.alkyl, ticcolor="grey", showfilenames = TRUE)



#### plots added 9-19-17 after discussion with Dorothy Madamba
ms.stackedplot(spectra, ionlist=0, which.spectra=order.last.samples, main="", showfilenames = TRUE)
ms.stackedplot(spectra, ionlist=ionlist.temp, which.spectra=order.last.samples, main="", showfilenames = TRUE)

ms.stackedplot(spectra, ionlist=0, which.spectra=order.p3a, main="", showfilenames = TRUE)
ms.stackedplot(spectra, ionlist=ionlist.temp, which.spectra=order.p3a, main="", showfilenames = TRUE)


############################################################
## Now work on barplots shoing individiul compounds from compounds.csv
############################################################

x.days=c(0,6,12,20,31,45,68)
y.values=x.days*NA
y.matrix=y.values
par(mar = c(8, 6, 4, 4) + 0.2, mfrow=c(1,1)) #add room for the rotated labels

for (i in 1:length(decomp.ions[,1])) # for every compound
{ 
  y.values=x.days*NA
  print(paste0(decomp.ions[i,"Long.Name"],"  (", decomp.ions[i,"Ion.1"], "):"))
  if (is.numeric(decomp.ions[i,"Ion.1"])&&(is.numeric(decomp.ions[i,"Retention.Time"]))) 
            y.values=integrate.peaks(spectra, peakmz=decomp.ions[i,"Ion.1"], 
                                                         peakrt=decomp.ions[i,"Retention.Time"], 
                                                         peakRTTolerance=decomp.ions[i,"Retention.Delta"])
  print("Integration done.")
  print(par("mfrow"))
  print(par("mar"))
  print(par("oma"))
  if (i==1) y.matrix=y.values
  if (i>1) y.matrix=cbind(y.matrix, y.values)
  par(mar = c(12, 5, 4, 3) + 0.2, mfrow=c(1,1)) #add room for the rotated labels
  #print("Integration done.")
  print(par("mfrow"))
  print(par("mar"))
  print(par("oma"))
  if (max(y.values, na.rm=TRUE)>0) barplot(y.values[order1], beside=TRUE,names.arg=spectra$print.names[order1,"print.name"], 
                                           las=2,cex.names=1, cex.axis=1, main=decomp.ions[i,"Short.Name"]
          )#args.legend = list(x = "topright", bty = "n", inset=c(-0.15, 0)))
  if (max(y.values, na.rm=TRUE)>0) barplot(y.values[order2], beside=TRUE,names.arg=spectra$print.names[order2,"print.name"], 
                                           las=2,cex.names=1, cex.axis=1, main=decomp.ions[i,"Short.Name"]
  )#args.legend = list(x = "topright", bty = "n", inset=c(-0.15, 0)))
  print("First plot done.") 
    par(mar = c(4, 4, 4, 4) + 0.2)
    #plot the extracted ion chromatogram in vicinity of the compounds ion #1 in compounds.csv
    plot.eic(spectra, which.spectra=order1,ion.to.plot=decomp.ions[i,"Ion.1"], 
           xlimits=c(max(0,decomp.ions[i,"Retention.Time"]-decomp.ions[i,"Retention.Delta"]-0.2),
                     decomp.ions[i,"Retention.Time"]+decomp.ions[i,"Retention.Delta"]+0.2),
           main=paste0(decomp.ions[i,"Short.Name"], ": All Samples"),
           plot.reference.lines=TRUE, reference.lines=c(decomp.ions[i,"Retention.Time"],decomp.ions[i,"Retention.Time"]-decomp.ions[i,"Retention.Delta"],
            decomp.ions[i,"Retention.Time"]+decomp.ions[i,"Retention.Delta"])
            )
    plot.eic(spectra, which.spectra=order2,ion.to.plot=decomp.ions[i,"Ion.1"], 
             xlimits=c(max(0,decomp.ions[i,"Retention.Time"]-decomp.ions[i,"Retention.Delta"]-0.2),
                       decomp.ions[i,"Retention.Time"]+decomp.ions[i,"Retention.Delta"]+0.2),
             main=paste0(decomp.ions[i,"Short.Name"], ": All Samples"),
             plot.reference.lines=TRUE, reference.lines=c(decomp.ions[i,"Retention.Time"],decomp.ions[i,"Retention.Time"]-decomp.ions[i,"Retention.Delta"],
                                                          decomp.ions[i,"Retention.Time"]+decomp.ions[i,"Retention.Delta"])
             )
    print(y.values)
  
   
  
}

############################################################
## Now work on group barplots with relative bars
############################################################


colnames(y.matrix)=decomp.ions[,"Short.Name"]
rownames(y.matrix)=spectra$cdfFiles

y.matrix=y.matrix[order1,]

y.matrix.rel=sweep(y.matrix, 2, apply(y.matrix[,], 2, max), FUN="/")*100 
                                               ## for each ion. set highest intens for water sample to 100
#y.matrix.rel=y.matrix.rel[,c(1,5,2,6,3,4)]
par(mar = c(8, 6, 4, 4) + 0.2, xpd=TRUE)
print("Relative Group plots coming")
print("Individual plots done")

## Multiplanel barplot
par(mar = c(8, 8, 4, 4) + 0.2, mfrow=c(1,1))

#################################################
#### PCA Time
#################################################
par(xpd=FALSE) #turn off writing outside plotting area if it is on
scale.setting=FALSE
#mypca.all <-prcomp(y.matrix[,], center=TRUE, scale=scale.setting)
mypca.all <- pca.report(t(y.matrix), center=TRUE, scale=TRUE, 
                        column.factors=spectra$print.names[,"Phase.Group"],
                        column.factors2=spectra$print.names[,"Run.of.Day"]
                        )
#mypca.all.2 <- pca.report(t(y.matrix[order.mask,]), center=TRUE, scale=TRUE, 
  #                        column.factors=spectra$print.names[order.mask,"Phase.Group"],
   #                       column.factors2=spectra$print.names[order.mask,"Run.of.Day"]
#)

dev.off()
## later check in void: 31 (methanol/methoxy?), 42/55/69 at 8.5min, 119 at 16min