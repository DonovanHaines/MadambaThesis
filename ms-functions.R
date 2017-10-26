#######################################################################################################
#######################################################################################################
##                                                                                                   ##
##            ms reading and stack plotting routines by                                              ##
##                                         Dr. Donovan C. Haines                                     ##
##                                                                                                   ##
##                                                                                                   ##
#######################################################################################################
#######################################################################################################

library(ncdf)
library(Peaks)
library(qgraph)
if (!exists("ms.functions.loaded")) ms.functions.loaded=TRUE

if (!exists("loadFiles")) loadFiles=TRUE  #set to false if spectra list of matrices has already been loaded into memory
#if (!exists("stackednamesize")) stackednamesize=0.5  # for cex font size in placing filenames on plot
#if (!exists("dayconversion")) stackeddayconversion=TRUE     # if dayconversion, try to label with day instead of filename
if (!exists("stockionlist")) stockionlist=c(34,43,48,57,58,71, 85, 93, 97, 99, 106, 113, 127, 133, 151, 193, 207, 209, 223, 247, 281)  #define a stock ion list that ionlist can be set to (and is default)
if (!exists("stockeiccolors")) stockeiccolors=c("red", "green", "blue", "darkgoldenrod", "purple", "orange", "skyblue", "turquoise", "wheat", "navy", "peru", "rosybrown", "gray10", "gray50", "gray 85" )  #define some useful colors that can be assigned to eiccolors
if (!exists("stocklighteiccolors")) stocklighteiccolors=c("red", "green", "darkgoldenrod", "purple", "orange", "skyblue", "turquoise", "wheat", "peru", "rosybrown", "gray10", "gray50", "gray 85" )  #define some useful colors that can be assigned to eiccolors
#if (!exists("eiccolors")) eiccolors = append(stocklighteiccolors, rainbow(length(stockionlist)))  #add 30 random colors to incr number avail

if (!exists("stackedalignchrom")) stackedalignchrom=FALSE
if (!exists("stackedalignrt")) stackedalignrt=5.05
if (!exists("stackedalignrtwindow")) stackedalignrtwindow=0.1

if (!exists("stackedplotalignrt")) stackedplotalignrt=TRUE


#######################################################################################################
#######################################################################################################
##                                                                                                   ##
##            Function cdf.to.csv.matrix: Read CDFs to matrix then save as csv                      ##
##                                                                                                   ##
##                                                                                                   ##
#######################################################################################################
#######################################################################################################

cdf.to.csv.matrix <- function(filename, 
                              eiclowerbound=-0.3, # this is the Chemstation default for binning ions that are to tenths when doing extracted ion chromatograms
                                                  # so m/z 20.6 is rounded to 20 but 20.7 is 21
                              verbosity=0)  # set verbosity to higher numbers for more info printed
{
  # filename is name of file
  #eic lowerbound is a limit for binning ions as simple rounding isn't what Chemstation does 
  # this is the Chemstation default for binning ions that are to tenths when doing extracted ion chromatograms
  # so m/z 20.6 is rounded to 20 but 20.7 is 21
  # verbosity: if true, print extra information about what the function is doing (debug etc) 
  
  ncin <- open.ncdf(filename)
  
  total_intensity=get.var.ncdf(ncin, "total_intensity")
  scans=length(total_intensity)
  if (verbosity>0) print(paste0("  Scans: ", scans))
  min_mass=min(get.var.ncdf(ncin,"mass_values"))
  max_mass=max(get.var.ncdf(ncin,"mass_values"))
  
  ionlist=seq(round(min_mass-0.101,0),round(max_mass-0.101,0))
  
  if (verbosity>0) print(paste0("  Masses from ", round(min_mass,2)," to ", round(max_mass,2)))
  
  scan_acquisition_time=max(get.var.ncdf(ncin, "scan_acquisition_time"))
  intervals=get.var.ncdf(ncin, "scan_acquisition_time")[2]-get.var.ncdf(ncin, "scan_acquisition_time")[1]
  
  #convert rt to min instead of sec (later should update this to be an option, for now not optional)
  scan_acquisition_time=scan_acquisition_time/60
  print(paste0("  Acq Time: ", round(scan_acquisition_time,3), "    (Delta: ", intervals,")"))
  
  eic_intensities=matrix(0, nrow=scans, ncol=length(ionlist)+2)
  
  scan_index=get.var.ncdf(ncin, "scan_index")
  mass_values=get.var.ncdf(ncin, "mass_values")
  intensity_values=get.var.ncdf(ncin, "intensity_values")
  progress=1:10
  progress<-round(length(intensity_values)/10*progress,0)
  
  if (verbosity>0) print(paste0("Generating EIC matrix, please be patient. Processing ", length(intensity_values), " individual entries."))
  for (entrynum in 1:length(intensity_values)) # for each mass, intensity pair in the file
  {
    if (entrynum %in% progress) print(paste0(round(entrynum/length(intensity_values)*100,0)," percent   (",entrynum," of ",length(intensity_values),")"))
     {
      ionnum=findInterval(mass_values[entrynum]+(1+eiclowerbound),ionlist) 
       {
        scannumber=findInterval(entrynum-1, scan_index)   
        eic_intensities[scannumber, ionnum]=eic_intensities[scannumber, ionnum] + intensity_values[entrynum]
      }
    }
  }


  eic_intensities[,length(ionlist)+1]<-total_intensity

 eic_intensities[,length(ionlist)+2]<-rowSums(eic_intensities[,1:length(ionlist)])

  colnames(eic_intensities)<-c(paste0("mz",ionlist),"TIC","TIC.Recalc")
  rownames(eic_intensities)<-round(get.var.ncdf(ncin, "scan_acquisition_time")/60,4) #convert to min
  close.ncdf(ncin)

  if (verbosity>0) print(paste0("Individual file eic matrix for the file:", filename, " : "))
  if (verbosity>0) print(eic_intensities[1:10,1:20])    #report the first part of the matrix for troubleshooting

  write.csv(eic_intensities,file=sub("\\.[[:alnum:]]+$",".CSV",filename))
  return(TRUE)
}    #### end of function cdf.to.csv.matrix

#######################################################################################################
#######################################################################################################
##                                                                                                   ##
##            Function cdf.to.csv.matrix.all: Read all CDFs to matrices then save as csv                      ##
##                                                                                                   ##
##                                                                                                   ##
#######################################################################################################
#######################################################################################################

cdf.to.csv.matrix.all <- function(file.directory=getwd(), eiclowerbound=-0.3, verbosity=0)
{
  # Define some parameters but don't overwrite existing values if there are some (allows parameters to be adjusted and retained run to run)
  
  #if (!exists("eiclowerbound")) eiclowerbound=-0.3   # this is the Chemstation default for binning ions that are to tenths when doing extracted ion chromatograms
  # so m/z 20.6 is rounded to 20 but 20.7 is 21
  
  
  # Make a list of all CDF files in the path (currently non-recursive only)
  
  cdfFiles=dir(file.directory, pattern="\\.CDF$")   
  cdfFiles=sort(sub("\\.[[:alnum:]]+$","",cdfFiles))
  
  donecdfFiles=dir(file.directory, pattern="\\.CSV")
  donecdfFiles=sort(sub("\\.[[:alnum:]]+$","",donecdfFiles))
  
  filesToConvert=setdiff(cdfFiles, donecdfFiles)
  #The list of donecdfFiles with extensions stripped
  if (length(filesToConvert)>0)
  {
    for (filenum in 1:length(filesToConvert))
    {
      print(paste("Reading file #", filenum, ". ", filesToConvert[filenum], ".CDF"))
      cdf.to.csv.matrix(paste0(getwd(),"/",filesToConvert[filenum],".CDF"), eiclowerbound=eiclowerbound, verbosity=verbosity)  
      
    }
    
  }
  
  print("Script ms-cdf-rmatrix.R done")
  return(length(filesToConvert)) 
}#### end of function cdf.to.csv.matrix.all

#cdf.to.csv.matrix.all(getwd(),verbosity=TRUE);




#######################################################################################################
#######################################################################################################
##                                                                                                   ##
##            Function ms.load: Load csv matrices and return 'spectra' object(list)                  ##
##                                                                                                   ##
##                                                                                                   ##
#######################################################################################################
#######################################################################################################


ms.load <- function(ms.dir=getwd(), verbosity=0, stackeddayconversion=TRUE)
{
#############################################################
#List the files and generate CSV versions of CDFs if needed
#############################################################

 #get a list of CDF Files
 cdfFiles=dir(ms.dir, pattern="\\.CDF$")   #CSV exported from CDF by ms-cdf-to-rmatrix.R
 cdfFiles=sort(sub("\\.[[:alnum:]]+$","",cdfFiles))  #sort and remove extension
 if (verbosity>0) print(cdfFiles)
 donecdfFiles=dir(ms.dir, pattern="\\.CSV$")
 donecdfFiles=sort(sub("\\.[[:alnum:]]+$","",donecdfFiles))
 if (verbosity>0) print(cdfFiles)
 filesToConvert=setdiff(cdfFiles, donecdfFiles)
 #The list of donecdfFiles with extensions stripped
 if (length(filesToConvert)>0) # First need to generate some missing .CSV versions of .CDF files in current workign directory
 {
  print("Calling .CSV generating script.")
  cdf.to.csv.matrix.all(getwd(), verbosity=0)
  if (verbosity>0) print(cdfFiles)
 }
if (verbosity>0) print(cdfFiles)
 #############################################################
 #Read the files into the list named 'spectra', which has cdfFiles as element #1 and the file data as a second list starting at #2
 #############################################################

 spectra<-list(cdfFiles=cdfFiles, spectra=list(NULL), ionlist=list(NULL), retention.adjust=list(NULL))

 for (filenum in 1:length(spectra$cdfFiles))
 {
  tempobj<-read.csv(paste0(cdfFiles[filenum],".CSV"),header=TRUE)
  print(paste0("Reading csv file", sub("\\.[[:alnum:]]+$",".CSV",cdfFiles[filenum]), " (File ", filenum," of ", length(spectra$cdfFiles),")") )
  spectra$spectra[[filenum]]<-tempobj   #add spectral matrix
  #spectra$ionlist[[filenum]]<-          #add list of ions for columns (since loading put in X's etc)
  spectra$ionlist[[filenum]]<-as.numeric(gsub("mz","",(colnames(spectra$spectra[[filenum]]))[2:(dim(spectra$spectra[[filenum]])[2]-2)]))
  spectra$retention.adjust[[filenum]]=spectra$spectra[[filenum]][,1]*0
 }
 print("Files read into spectra object.")
 if (length(unique(sapply(spectra$ionlist, "length")))>1)
  print(paste0("Warning! Ionlists vary in length: ", sapply(spectra$ionlist, "length")))
 spectra$masterionlist=unique(sort((sapply(spectra$ionlist,"unique"))))
 #if (!exists("stackedrtadj")) stackedrtadj=rep(0, length(spectra$cdfFiles))

 #Determine days if day conversion is requested
 if (stackeddayconversion)
  {
   spectra$cdfFilesDays=as.Date(strtrim(spectra$cdfFiles,8), "%Y%m%d")
   spectra$cdfFilesDays=spectra$cdfFilesDays-min(spectra$cdfFilesDays)+1
 }
 
 #### load more nicely printing names if they exist

 tryCatch(
   {
    nick.name.temp<-read.csv("printnames.csv",header=TRUE, stringsAsFactors=FALSE)  ### if a file printnames.csv exists, load it (assumes file has two cols, one Filenames and the other Print.Names)
    spectra$print.names=nick.name.temp
    #print("Print.names loaded")
    #print(nick.name.temp)
   },
   error=function(e)
   {
    spectra$print.names=cbind(spectra$cdfFiles, spectra$cdfFiles, stringsAsFactors=FALSE) # if the file can't be loaded or doesnt' exist, 
                                                                  # make both columns the filename
    print("Error loading print.names")
   },
   warning=function(w)
   {
     spectra$print.names=cbind(spectra$cdfFiles, spectra$cdfFiles) # if the file can't be loaded or doesnt' exist, 
                                                                   # make both columns the filename
     print("Warning loading print.names")
   }
   )
 
 
 spectra$print.names<-spectra$print.names[order(spectra$print.names[,"filename"]),] # order the same as cdfFiles, based on filename
 if (verbosity>5) print(spectra$print.names)
 
 if (sum(!(spectra$cdfFiles%in%spectra$print.names[,"filename"]))>0)
 {# some cdfFiles don't have info in print.names table
   print(" **** Warning! Some spectra cdfFiles don't have entries in the print.names table")
   for (missing.entry in 1:sum(!(spectra$cdfFiles%in%spectra$print.names[,"filename"])))
   {
     spectra$print.names<-rbind(spectra$print.names, 
                                c(spectra$cdfFiles[which(!(spectra$cdfFiles%in%spectra$print.names[2]))[missing.entry]],
                                  spectra$cdfFiles[which(!(spectra$cdfFiles%in%spectra$print.names[2]))[missing.entry]],
                                  rep(NA,dim(spectra$print.names)[2]-2)))
     
   }
 }
 if (sum(!(spectra$print.names[,"filename"]%in%spectra$cdfFiles))>0)
 {# some cdfFiles don't have info in print.names table
   print(" **** Warning! Some spectra print.names files don't have corresponding cdfFiles!!!")
 }
 return(spectra)
}   ############################# end ms loading function





#######################################################################################################
#######################################################################################################
##                                                                                                   ##
##            Function ms.stackedplot: Create eic and tic stacked plots from 'spectra' list object   ##
##                                                                                                   ##
##                                                                                                   ##
#######################################################################################################
#######################################################################################################

ms.stackedplot <- function (spectra=NULL, #a list object made by the loading routines
                            ionlist=stockionlist, #list of ions for extraction ion chromatograms
                            eiccolors=rainbow(length(ionlist)),
                            which.spectra=0, #a vector of indices of the spectra in 'spectra' to plot
                            chromatogramxlim=c(-1,-1), 
                            stackedshoweic=TRUE, 
                            stackedshowtic=TRUE, 
                            verbosity=0, #set to higher number to be more verbose in printing info
                            relativeeic=FALSE, #plot EICs as relative to max value for the eic (usually want false)
                            stackzoom=1, # zoom factor, increase to zoom y-axis
                            retentiontimelist=NA, # list of retention time at which to draw vertical lines
                            showfilenames=TRUE,
                            stackeddayconversion=FALSE,
                            stackednamesize=0.5,
                            namesattop=FALSE,
                            namesatright=FALSE,
                            namesatcenter=FALSE,
                            namesoffset=0,
                            margins=c(4,4,4,4),
                            showaxes=TRUE,
                            tic.line.width=1, #width of TIC line
                            yaxisside=2,
                            showyaxis=TRUE,
                            xaxisside=1,# 2 is left, 4 is right
                            showxaxis=TRUE,
                            showbox=TRUE,
                            showtitle=TRUE,
                            ticcolor="black",
                            main="",
                            correct.split=FALSE,
                            groups=0,  #if replicates, each is a group as designated by this factor vector
                            ...)
{
  old.par=par(no.readonly=TRUE) #record graphic parameters so we can politely put them back that way when done
  if (correct.split==FALSE) split.factors<-rep(1.0,length(spectra$cdfFiles))
  if (correct.split==TRUE) split.factors<-as.double(spectra$print.names[,"split"])
  
 #find maximum TIC y-value overall
  maxtic=0
  for (i in 1:length(spectra$cdfFiles))  ###need tic max for each and every ion in ionlist
  {
    if((which.spectra[1]==0)||(i %in% which.spectra))
    {
      maxtic=#max(maxtic,max(spectra$spectra[[i]][,"TIC"]))
            max( c(maxtic,
                   split.factors[i]*max(spectra$spectra[[i]][(findInterval(chromatogramxlim[1],spectra$spectra[[i]][,1]):
                                    findInterval(chromatogramxlim[2],spectra$spectra[[i]][,1])),
                                    "TIC"]    )))
      if ((chromatogramxlim[1]==-1) && (chromatogramxlim[2]==-1))
        maxtic=max( c(maxtic, split.factors[i]*max(spectra$spectra[[i]][,"TIC"])))
    }
  }
 #maxtic=max(unlist(lapply(spectra$spectra,FUN=max)))   #get the max value of all the matrices, which should be the max TIC
 maxticplus=maxtic*1.1
 if (verbosity>3) print(paste0("Max TIC: ", maxtic))
 # determine number of panels required
 if (length(groups)<2) number.of.panels=length(which.spectra)
 if (length(groups)>1)
 {
   if (which.spectra[1]!=0) number.of.panels=nlevels(groups[which.spectra])
   if (which.spectra[1]==0) number.of.panels=nlevels(groups[])
 }
 if (which.spectra[1]==0) number.of.panels=length(spectra$cdfFiles)
 if (length(groups)>1) 
   {
    
    if ((which.spectra[1])!=0)
    {
      #trim the groups vector down to just those entries to be printed
      if(length(groups)==length(spectra$cdfFiles)) groups<-as.factor(as.character(groups[which.spectra]))
    }
   number.of.panels=length(levels(groups))
   }
 if (verbosity>5) print(paste0("Number of panels: ", number.of.panels))
 
 mintime=0 #min(spectra$spectra[[1]][,1])
 maxtime=-1 #max(spectra$spectra[[1]][,1])
 for (i in 1:length(spectra$cdfFiles))
 {
   if((which.spectra[1]==0)||(i %in% which.spectra))
   {
     mintime=min(mintime,min(spectra$spectra[[i]][,1]))
     maxtime=max(maxtime,max(spectra$spectra[[i]][,1]))
     
   }
 }
 #reset graph margins to a reasonable value
 par(mar = margins + 0.2)
 
 print(paste0("Min time: ", round(mintime,3),"  Max time: ", round(maxtime,3)))
 if (chromatogramxlim[2]==-1) chromatogramxlim=c(round(mintime,3),round(maxtime,3))
 if (stackedshoweic) 
 {
   
   if (verbosity>0)
   {
     print(spectra$maxeic[1:min(10, length(spectra$maxeic))])
     print(spectra$maxeic.trimmed[1:min(10, length(spectra$maxeic))])
   }  
   
  spectra$maxeic         = spectra$masterionlist*0
  spectra$maxeic.trimmed = spectra$maxeic
  
  if (verbosity>0)
  {
    print(spectra$maxeic[1:min(10, length(spectra$maxeic))])
    print(spectra$maxeic.trimmed[1:min(10, length(spectra$maxeic))])
    print(paste0("Chromatogramxlim is ", chromatogramxlim))
  }
  
  for (i in 1:length(spectra$cdfFiles))  ###Note:fix; need eic max for each and every ion in ionlist
   {
    if((which.spectra[1]==0)||(i %in% which.spectra))
    {
     for (j in 1:length(spectra$masterionlist))
     {
      spectra$maxeic[j]=
             max( c(
                  spectra$maxeic[j],
                   split.factors[i]*max(spectra$spectra[[i]][,findInterval(spectra$masterionlist[j],spectra$ionlist[[i]])+1]    )))
      spectra$maxeic.trimmed[j]=
             max( c(
                  spectra$maxeic.trimmed[j],
                   split.factors[i]*max(spectra$spectra[[i]][(findInterval(chromatogramxlim[1],spectra$spectra[[i]][,1]):findInterval(chromatogramxlim[2],spectra$spectra[[i]][,1])),findInterval(spectra$masterionlist[j],spectra$ionlist[[i]])+1]    )))
    
    
     }
    }
   }
 print("Maximum EIC Intensities")
 print("Ions   : ")
 print (ionlist)
 print("Max Int: (all then trimmed) ")
 print(spectra$maxeic[1:min(10, length(spectra$maxeic))])
 print(spectra$maxeic.trimmed[1:min(10, length(spectra$maxeic))])
 if (verbosity>20) 
  {
   print("5 sec pause...")
   Sys.sleep(5)
  }
 }



 scaleadj=1 #length(cdfFiles)

 stackedfontspacer=(maxtic*number.of.panels/stackzoom)*0.011 #was 0.13 before 1-31-17
                     #place labels this fraction of whole graph below traces

 if (stackedshoweic) maxeicall=max(spectra$maxeic.trimmed[findInterval(ionlist,spectra$masterionlist)])
 if (verbosity>0) print(paste0("Stackedshoweic is ", stackedshoweic," and Max eic overall is ", maxeicall))
 #plot(stacked_scan_acquisition_times[,1], 0, type="l", xlab="Time (min)", ylab="TIC Intensity", ylim=c(0,maxticplus*length(cdfFiles)/stackzoom), xlim=chromatogramxlim)
 if (verbosity>7) print(paste0("Groups as of right now"))
 if (verbosity>7) print(groups)
 if (verbosity>7) print(paste0("Which.spectra as of right now"))
 if (verbosity>7) print(which.spectra)
 if (verbosity>7) print(paste0("Groups[which.spectra] as of right now"))
 if (verbosity>7) print(groups[which.spectra])
 if (verbosity>9) print(paste0("spectra$print.names[,\"group.auto\"]"))
 if (verbosity>9) print(spectra$print.names[,"group.auto"])
 
 
 ############################ MAIN LOOP
 ###########################################
 
 printed.filenum=0
  for (filenum in 1:length(spectra$cdfFiles)) # for every file in the spectra list (note: not the which.spectra list)
  {
   
   if((which.spectra[1]==0)||(filenum %in% which.spectra)) #if it is a spectrum we are supposed to plot
   {
    ordered.file=which(which.spectra==filenum)  #which entry of which.spectrum are we working on
    printed.filenum=printed.filenum+1 #how many spectra have we printed
    if (verbosity>0) print(paste0("File ",filenum," counts, adding as printed.filenum: ", printed.filenum,"."))
    which.panel<-printed.filenum
    if (verbosity>9) print(paste0("Initial which.panel set to ", which.panel))
    
    if (length(groups)>1) which.panel<-which(as.integer(groups[ordered.file])==unique(as.integer(groups)))
      #length(levels(as.factor(as.integer(groups[1:printed.filenum])))) #how many levels up in groups is panel number
    
    if (verbosity>9) print(paste0("Which.panel set to ", which.panel, " after check for groups"))
    verticaloffset=(which.panel-1)*maxticplus/stackzoom #how high do we plot it based on which panel it goes to
    
    if(verbosity>5) print(paste0("Panel ", which.panel, " Vertical offset: ", verticaloffset))
    if (length(groups)==1) 
      {
       if (which.spectra[1]!=0) verticaloffset=(match(filenum, which.spectra)-1)*maxticplus/stackzoom #if no groups but which.spectra specifies subset to plot
      }
    
    graphtop=maxticplus/stackzoom*number.of.panels
    print(paste0("Initial plotting... stackzoom: ", stackzoom, "   verticaloffset: ", verticaloffset))
    #abline(h=verticaloffset) 
    if (printed.filenum==1)#(filenum==1)||(filenum==min(which.spectra))) # if this is the first spectrum we are plotting, then we need to generate the whole plot not just add line
      {
       if (stackedshowtic) # REMOVE AFTER DEBUG
         {
         plot(spectra$spectra[[filenum]][,1]#spectra$spectra[[filenum]][,"TIC"]<graphtop/number.of.panels,1]
              +spectra$retention.adjust[[filenum]],#[(spectra$spectra[[filenum]][,"TIC"]<graphtop/number.of.panels)], 
              #split.factors[filenum]
                sapply(spectra$spectra[[filenum]][,"TIC"], function(x) min(x*split.factors[filenum],maxtic/stackzoom))
                +verticaloffset,
              type="l", 
              xlab="Time (min)", ylab="TIC Intensity", ylim=c(0,graphtop), 
              xlim=chromatogramxlim, lwd=tic.line.width, col=ticcolor,
              axes=FALSE, ann=FALSE, ...)
         if (showxaxis) {axis(xaxisside); title(xlab="Time (min)")}
         if (showyaxis) {axis(yaxisside); title(ylab="Intensity")}
         if (showbox) box()
         if (showtitle) title(main)
         #abline(h=maxtic*0.99, col="red")
         #abline(h=maxtic/stackzoom*0.99, col="red", lty=2)
         
         }
    }
  # print(filenum)
  #  print(length(spectra$spectra[[filenum]][,"TIC"]))
  #  print(length(spectra$spectra[[filenum]][(spectra$spectra[[filenum]][,"TIC"]<maxticplus/stackzoom),1]
  #                +spectra$retention.adjust[[filenum]][(spectra$spectra[[filenum]][,"TIC"]<maxticplus/stackzoom)]))
  #  print(length(spectra$spectra[[filenum]][(spectra$spectra[[filenum]][,"TIC"]<maxticplus/stackzoom),1]))
  #  print(length(spectra$retention.adjust[[filenum]][(spectra$spectra[[filenum]][,"TIC"]<maxticplus/stackzoom)]))          
  #  
  #  print(length(split.factors[filenum]*spectra$spectra[[filenum]][(spectra$spectra[[filenum]][,"TIC"]<maxticplus/stackzoom),"TIC"]*stackzoom+verticaloffset))
  #  print(paste0("Split.factor[filenum] length is ", length(split.factors[filenum]), " and value is ", split.factors[filenum]))
  #  print(length(spectra$spectra[[filenum]][(spectra$spectra[[filenum]][,"TIC"]<maxticplus/stackzoom),"TIC"]))
  #  if (length(spectra$spectra[[filenum]]
  #              [(spectra$spectra[[filenum]][,"TIC"]<maxticplus/stackzoom),1]
  #              +spectra$retention.adjust[[filenum]])!=  
  #       length(split.factors[filenum]*spectra$spectra[[filenum]][(spectra$spectra[[filenum]][,"TIC"]<maxticplus/stackzoom),"TIC"]
  #              *stackzoom+verticaloffset))
  #     {
  #      print("Boom! Goes the dynamite.")
  #     }
    
    print(paste0("maxticplus is ", maxticplus, " and stackzoom is ", stackzoom))
    #if (filenum>1) {if (stackedshowtic) lines(spectra$spectra[[filenum]][(spectra$spectra[[filenum]][,"TIC"]<maxticplus/stackzoom),1]
    #                                          +spectra$retention.adjust[[filenum]][(spectra$spectra[[filenum]][,"TIC"]<maxticplus/stackzoom)], 
    #                   split.factors[filenum]*spectra$spectra[[filenum]][(spectra$spectra[[filenum]][,"TIC"]<maxticplus/stackzoom),"TIC"]*stackzoom+verticaloffset, 
    #                                          col=ticcolor, lwd=tic.line.width)}
     if (stackedshoweic)
     {
      for (k in 1:length(ionlist))
       {
        if (!stackedshowtic) 
         {
          if ((filenum==1)||(filenum==which.spectra[1])) 
           {
            if (k==1) 
              {
              plot(spectra$spectra[[filenum]][,1]
                   +spectra$retention.adjust[[filenum]][(spectra$spectra[[filenum]][,"TIC"]<maxticplus/stackzoom)], 
                          split.factors[filenum]*spectra$spectra[[filenum]][,1+findInterval(ionlist[k],spectra$ionlist[[filenum]])]*stackzoom*maxtic/maxeicall+verticaloffset, 
                             col=eiccolors[k], type="l", xlab="Time (min)", ylab="TIC Intensity", ylim=c(0,graphtop), xlim=chromatogramxlim, 
                             axes=FALSE, ann=FALSE,...)
              if (showxaxis) {axis(xaxisside); title(xlab="Time (min)")}
                if (showyaxis) {axis(yaxisside); title(ylab="Intensity")}
              if (showbox) box()
              if (showtitle) title(main)
            }
           }
         }
        lines(spectra$spectra[[filenum]][,1]+spectra$retention.adjust[[filenum]], 
             split.factors[filenum]*spectra$spectra[[filenum]][,1+findInterval(ionlist[k],spectra$ionlist[[filenum]])]*stackzoom*maxtic/maxeicall+verticaloffset, col=eiccolors[k])
       }

     }
    #abline(h=maxtic, col="red")
    #abline(h=maxticplus, col="green")
    
    
    
    
    #plot lines over the top
    if (stackedshowtic) lines(spectra$spectra[[filenum]][,1]#spectra$spectra[[filenum]][,"TIC"]<graphtop/number.of.panels,1]
                              +spectra$retention.adjust[[filenum]],#[(spectra$spectra[[filenum]][,"TIC"]<graphtop/number.of.panels)], 
                              #split.factors[filenum]
                              sapply(spectra$spectra[[filenum]][,"TIC"], function(x) min(x*split.factors[filenum],maxtic/stackzoom))
                              +verticaloffset, type="l", 
                              lwd=tic.line.width, col=ticcolor)
    #######if (stackedshowtic) lines(stacked_scan_acquisition_times[,filenum]+stackedrtadj[filenum], stacked_total_intensities[,filenum]*stackzoom+verticaloffset, lwd=tic.line.width*2)
  
    tracelabel=""
    
    if ((showfilenames)&&("Print.Name" %in% colnames(spectra$print.names))) tracelabel=paste0(tracelabel, spectra$print.names[filenum,"Print.Name"], " ")
    if ((showfilenames)&&("print.name" %in% colnames(spectra$print.names))) tracelabel=paste0(tracelabel, spectra$print.names[filenum,"print.name"], " ")
    
    if ((showfilenames)&&(length(groups)>1)) 
    {
      tracelabel=paste(spectra$print.names[which.spectra[ which(groups==groups[which(which.spectra==filenum)])],"print.name"], collapse=", ") 
    }
    if (stackeddayconversion) tracelabel=paste0(tracelabel, "Day ", spectra$cdfFilesDays[filenum], " ")
    print(paste0("Filenum: ", filenum, "   Tracelabel: ", tracelabel)) ### DEBUG 
    if ((showfilenames)&&(!namesattop)&&(!namesatright)&&(!namesatcenter)) 
           text(namesoffset*(chromatogramxlim[2]-chromatogramxlim[1])+chromatogramxlim[1], 
                verticaloffset-stackedfontspacer-stackednamesize*maxticplus/stackzoom*0.1, 
                tracelabel, cex=stackednamesize, adj=c(0,0))
    if ((showfilenames)&&(namesattop)&&(!namesatright)&&(!namesatcenter)) text(namesoffset*(chromatogramxlim[2]-chromatogramxlim[1])+chromatogramxlim[1], verticaloffset+maxticplus/stackzoom-stackedfontspacer-stackednamesize*maxticplus/stackzoom*0.1, tracelabel, cex=stackednamesize, adj=c(0,0))
    if ((showfilenames)&&(!namesattop)&&(namesatright)&&(!namesatcenter)) text(namesoffset*(chromatogramxlim[2]-chromatogramxlim[1])+chromatogramxlim[2]-(chromatogramxlim[2]-chromatogramxlim[1])*0.01, verticaloffset-stackedfontspacer-stackednamesize*maxticplus/stackzoom*0.1, tracelabel, cex=stackednamesize, adj=c(1,0))
    if ((showfilenames)&&(namesattop)&&(namesatright)&&(!namesatcenter)) text(namesoffset*(chromatogramxlim[2]-chromatogramxlim[1])+chromatogramxlim[2]-(chromatogramxlim[2]-chromatogramxlim[1])*0.01, verticaloffset+maxticplus/stackzoom-stackedfontspacer-stackednamesize*maxticplus/stackzoom*0.1, tracelabel, cex=stackednamesize, adj=c(1,0))
    if ((showfilenames)&&(!namesattop)&&(!namesatright)&&(namesatcenter)) 
      text(namesoffset*(chromatogramxlim[2]-chromatogramxlim[1])+mean(c(chromatogramxlim[2], chromatogramxlim[1])), 
           verticaloffset-stackedfontspacer-stackednamesize*maxticplus/stackzoom*0.1, 
           tracelabel, cex=stackednamesize)
    if ((showfilenames)&&(namesattop)&&(!namesatright)&&(namesatcenter)) text(namesoffset*(chromatogramxlim[2]-chromatogramxlim[1])+mean(c(chromatogramxlim[2], chromatogramxlim[1])), verticaloffset+maxticplus/stackzoom-stackedfontspacer-stackednamesize*maxticplus/stackzoom*0.1, tracelabel, cex=stackednamesize)
    
    
   abline(h=verticaloffset, lwd=1) 
     if (((filenum==1)&&(which.spectra[1]==0))||(filenum==which.spectra[1])) # first file so write labels
      {
       eiclabelsperline=50
       eiclabellinespace=stackedfontspacer
       eiclabelspacing=(chromatogramxlim[2]-chromatogramxlim[1])/(eiclabelsperline+2)
       
       #if (stackedshowtic) text(chromatogramxlim[1], graphtop*1.02+eiclabellinespace, "TIC", col="black", cex=0.5)
       if (stackedshoweic)
       {
       for (eiccount in 1:length(ionlist))
        {
          eiclabelrow=(eiccount-1)%/%eiclabelsperline
          eiclabelcol=(eiccount-1)%%eiclabelsperline
          if (ionlist[eiccount]!=0) text(chromatogramxlim[1]+(eiclabelcol)*eiclabelspacing, graphtop*1.02-eiclabellinespace*(eiclabelrow-1), ionlist[eiccount], col=eiccolors[eiccount], cex=0.5)
          #if (eiccount<eiclabelsperline) text(chromatogramxlim[1]+(eiccount)*eiclabelspacing, graphtop*1.02+eiclabellinespace, ionlist[eiccount], col=eiccolors[eiccount], cex=0.5)
          # else if (eiccount<eiclabelsperline*2) {text(chromatogramxlim[1]+(eiccount-eiclabelsperline-1)*eiclabelspacing, graphtop*1.02,                    ionlist[eiccount], col=eiccolors[eiccount], cex=0.5)}
          # else if (eiccount<eiclabelsperline*3){text(chromatogramxlim[1]+(eiccount-eiclabelsperline*2-1)*eiclabelspacing, graphtop*1.02-eiclabellinespace*1,               ionlist[eiccount], col=eiccolors[eiccount], cex=0.5) }
          # else if (eiccount<eiclabelsperline*4){text(chromatogramxlim[1]+(eiccount-eiclabelsperline*3-1)*eiclabelspacing, graphtop*1.02-eiclabellinespace*2,              ionlist[eiccount], col=eiccolors[eiccount], cex=0.5) }
          # else if (eiccount<eiclabelsperline*5){text(chromatogramxlim[1]+(eiccount-eiclabelsperline*4-1)*eiclabelspacing, graphtop*1.02-eiclabellinespace*3,               ionlist[eiccount], col=eiccolors[eiccount], cex=0.5) }
          # else if (eiccount<eiclabelsperline*6){text(chromatogramxlim[1]+(eiccount-eiclabelsperline*5-1)*eiclabelspacing, graphtop*1.02-eiclabellinespace*4,            ionlist[eiccount], col=eiccolors[eiccount], cex=0.5) }
        }
       if ((stackedshowtic)&&(ionlist[1]!=0)) text(chromatogramxlim[1]+(eiclabelcol+1)*eiclabelspacing, graphtop*1.02-eiclabellinespace*(eiclabelrow-1), "TIC", cex=0.5, adj=c(0,NA), col=ticcolor)
       if (ionlist[1]!=0) text(chromatogramxlim[1]+(eiclabelcol+2)*eiclabelspacing, graphtop*1.02-eiclabellinespace*(eiclabelrow-1), paste0("(EICs are x", round(maxtic/maxeicall,2), ")"), cex=0.5, adj=c(0,NA))
             
       }      
      }
     
   }
 }
 #title(paste0("CDF Files in ", stackpath))
 if (verbosity>1) print(!is.na(retentiontimelist[1]))
 if (!is.na(retentiontimelist[1]))
 {
  for (i in 1:length(retentiontimelist)) abline(v=retentiontimelist[i], lty="dotted")
 }

 if ((stackedalignchrom)**(stackedplotalignrt)) abline(v=stackedrtmax[1], lty="dotted", col="gray")
par(old.par)
}
#######################################################################################################
#######################################################################################################
##                                                                                                   ##
##            Function integrate.peaks: return integration for all peaks                             ##
##                                                                                                   ##
##                                                                                                   ##
#######################################################################################################
#######################################################################################################

################## Function integrate.peaks (to be improved later)
integrate.peaks<-function (spectra, peakmz=117, peakrt=25.5, which.spectra=0, 
                           peakRTTolerance=0.2, verbosity=0, plot.integration=FALSE, subtract.background=FALSE)
{
#inegrate defaults mz117 at 24.4 min as indole
 intsum=rep(0,length(spectra$cdfFiles))
if (which.spectra==0) which.spectra=1:length(spectra$cdfFiles)
#peakmz=117
#peakrt=24.4

for (filenum in 1:length(spectra$cdfFiles))
{ 
 if (verbosity>5) print(paste0("Int.peaks file number ", filenum, " of ", length(spectra$cdfFiles)))
 intsum[filenum]<- integrate.peak(spectra, which.spectrum=filenum, peakmz=peakmz, peakrt=peakrt,
                                  peakRTTolerance=peakRTTolerance, verbosity=verbosity,
                                  plot.integration=plot.integration, subtract.background=subtract.background)
 if (verbosity>5) print(intsum[filenum])
}

 if (FALSE)#plot.integration)
 {
   plot.eic(spectra, which.spectra=which.spectra, ion.to.plot=peakmz, 
            xlimits=c(peakrt-10*peakRTTolerance, peakrt+10*peakRTTolerance), verbosity=verbosity,
            plot.reference.lines=TRUE, 
            reference.lines=c(peakrt, peakrt+peakRTTolerance, peakrt-peakRTTolerance)
            )
   
 }
 if (verbosity>1) print(intsum)
return(intsum[])
}
#######################################################################################################
#######################################################################################################
##                                                                                                   ##
##            Function integrate.peak: return integration for one peak                               ##
##                                                                                                   ##
##                                                                                                   ##
#######################################################################################################
#######################################################################################################

integrate.peak<-function (spectra, which.spectrum=1, peakmz=117, peakrt=25.5, peakRTTolerance=0.2, 
                          verbosity=0, plot.integration=FALSE,method="simple",
                          subtract.background=FALSE)
{
  if (method=="simple") return(integrate.peak.simple(spectra, which.spectrum, peakmz, peakrt, peakRTTolerance, 
                                              verbosity, plot.integration, subtract.background=subtract.background))
  if (method=="gaussian") return(integrate.peak.gaussian(spectra, which.spectrum, peakmz, peakrt, peakRTTolerance, 
                                                      verbosity, plot.integration))
}


integrate.peak.simple<-function (spectra, which.spectrum=1, peakmz=117, peakrt=25.5, peakRTTolerance=0.2, 
                          verbosity=0, plot.integration=FALSE,method="simple", 
                          subtract.background=subtract.background)
{
  #inegrate defaults mz117 at 24.4 min as indole
  sum.total=0
  
  #peakmz=117
  #peakrt=24.4
  
  
    sum.total<- sum(spectra$spectra[[which.spectrum]][
      ((spectra$spectra[[which.spectrum]][,1]>=(peakrt-peakRTTolerance))&(spectra$spectra[[which.spectrum]][,1]<=(peakrt+peakRTTolerance))),
      1+findInterval(peakmz,spectra$ionlist[[which.spectrum]])])
    
    number.x.points<-sum((spectra$spectra[[which.spectrum]][,1]>=(peakrt-peakRTTolerance))&(spectra$spectra[[which.spectrum]][,1]<=(peakrt+peakRTTolerance)))
    
    background.min <- min(spectra$spectra[[which.spectrum]][
      ((spectra$spectra[[which.spectrum]][,1]>=(peakrt-peakRTTolerance))&(spectra$spectra[[which.spectrum]][,1]<=(peakrt+peakRTTolerance))),
      1+findInterval(peakmz,spectra$ionlist[[which.spectrum]])])
    background.sum <- number.x.points*background.min
    if (!subtract.background) background.sum=0
    
    if (verbosity>10) print(spectra$spectra[[which.spectrum]][
      ((spectra$spectra[[which.spectrum]][,1]>=(peakrt-peakRTTolerance))&(spectra$spectra[[which.spectrum]][,1]<=(peakrt+peakRTTolerance))),
      1+findInterval(peakmz,spectra$ionlist[[which.spectrum]])])
    
    if (plot.integration)
    {
      plot.eic(spectra, which.spectra=which.spectrum, ion.to.plot=peakmz, 
               xlimits=c(peakrt-2*peakRTTolerance, peakrt+2*peakRTTolerance), verbosity=verbosity,
               plot.reference.lines=TRUE, 
               reference.lines=c(peakrt, peakrt+peakRTTolerance, peakrt-peakRTTolerance),
               plot.segment.lines = TRUE,
               segment.lines = c(peakrt-peakRTTolerance, background.min,
                                 peakrt+peakRTTolerance, background.min))
      
    }
  
  
  
  if (verbosity>0) print(paste0("Returning ", sum.total, " from integrate.peak function."))
  return(sum.total-background.sum)
}
integrate.peak.gaussian<-function (spectra, which.spectrum=1, peakmz=117, peakrt=25.5, peakRTTolerance=0.2, 
                                 verbosity=0, plot.integration=FALSE,method="simple")
{
  #inegrate defaults mz117 at 24.4 min as indole
  sum.total=0
  
  #peakmz=117
  #peakrt=24.4
  
  
  sum.total<- sum(spectra$spectra[[which.spectrum]][
    ((spectra$spectra[[which.spectrum]][,1]>=(peakrt-peakRTTolerance))&(spectra$spectra[[which.spectrum]][,1]<=(peakrt+peakRTTolerance))),
    1+findInterval(peakmz,spectra$ionlist[[which.spectrum]])])
  
  if (verbosity>10) print(spectra$spectra[[which.spectrum]][
    ((spectra$spectra[[which.spectrum]][,1]>=(peakrt-peakRTTolerance))&(spectra$spectra[[which.spectrum]][,1]<=(peakrt+peakRTTolerance))),
    1+findInterval(peakmz,spectra$ionlist[[which.spectrum]])])
  
  if (plot.integration)
  {
    plot.eic(spectra, which.spectra=which.spectrum, ion.to.plot=peakmz, 
             xlimits=c(peakrt-10*peakRTTolerance, peakrt+10*peakRTTolerance), verbosity=verbosity,
             plot.reference.lines=TRUE, 
             reference.lines=c(peakrt, peakrt+pearkRTTolerance, peakrt-pearkRTTolerance))
    
  }
  
  
  
  if (verbosity>0) print(paste0("Returning ", sum.total, " from integrate.peak function."))
  return(sum.total)
}

#######################################################################################################
#######################################################################################################
##                                                                                                   ##
##            Function plot.eic: Plot local eic traces for selected spectra                          ##
##                                                                                                   ##
##                                                                                                   ##
#######################################################################################################
#######################################################################################################


plot.eic<-function(spectra, 
                   which.spectra=1, #which spectra numbers to plot (vector)
                   ion.to.plot=34,  #which ion to plot, 0 means TIC
                   xlimits=c(0,30), #time limits to plot
                   eic.colors=rainbow(length(which.spectra)), #colors to use
                   vertical.offset=0, # absolute vertical offset for each line
                   vertical.offset.relative=0,#vertical offset for each line, 
                                              #as a fraction of the y.max
                   horizontal.offset=0, #absolute horizontal offset for each line
                   show.legend=TRUE,
                   verbosity=0,  # set to higher integer to print more info for debugging
                   plot.reference.lines=FALSE,
                   reference.lines=c(0,0,0),
                   reference.line.colors=rep("black", length(reference.lines)),
                   reference.lty=rep(2, length(reference.lines)),
                   plot.segment.lines=FALSE,
                   segment.lines=c(0,0,0,0),
                   segment.lty=3,
                   segment.color="black",
                   use.par=FALSE,
                   par.to.use=NULL,
               
                   ...) #catch any other arguments and set to plot functions
{  
  #old.par=par(no.readonly=TRUE)
  if (is.na(ion.to.plot))
  {
    print("Error! Ion to plot is NA in plot.eic! Printing blank plot")
    plot.new()
    return()
    
  }
  if (ion.to.plot>0) ion.col=paste0("mz",ion.to.plot)
  if (ion.to.plot==0) ion.col="TIC"  #if 0 use TIC
  
  max.y=0
  min.y=0
  for (i in 1:length(which.spectra))
   {
    max.y<- max( c(max.y,
      max(((i-1)*vertical.offset)+(vertical.offset.relative*(i-1)*max.y)+spectra$spectra[[(which.spectra[i])]]
             [(findInterval(
                    xlimits[1],
                       spectra$spectra[[(which.spectra[i])]][,1]):findInterval(xlimits[2],spectra$spectra[[(which.spectra[i])]][,1])),
               ion.col]    )))
    min.y<-min( c(min.y,
                min(((i-1)*vertical.offset)+(vertical.offset.relative*(i-1)*max.y)+spectra$spectra[[(which.spectra[i])]]
                    [(findInterval(xlimits[1],spectra$spectra[[(which.spectra[i])]][,1]):findInterval(xlimits[2],spectra$spectra[[(which.spectra[i])]][,1])),
                      ion.col]    )))
  }
  new.xlimits=xlimits
  new.xlimits[1]=min(xlimits[1],xlimits[1]+horizontal.offset)
  new.xlimits[2]=max(xlimits[2],xlimits[2]+horizontal.offset)
  
  new.ylimits=c(min(0,min.y),max.y)
  
  
  if (verbosity>0)
  {
    print(paste0("Plotting eic/tic:"))
    print(paste0("  Spectra: ", which.spectra))
    print(paste0("  Ion: ", ion.to.plot))
    if (ion.to.plot==0) print ("      (0 means TIC is plotted)")
    print(paste0("xlimits: ", xlimits))
    print(paste0("offsets: vertical: ", vertical.offset, "  horizontal: ", horizontal.offset))
    print(paste0("         vertical.relative: ", vertical.offset.relative))
    
    print(paste0("Show legend? ", show.legend))
    print(paste0("Calculated min.y as ", min.y, " and max.y as ", max.y))
  }  
  
  if (show.legend) layout(rbind(1,2), heights=c(5,1))
  par(mar = c(4, 4, 4, 4) + 0.2)
  if (use.par) par(par.to.use)
  plot(spectra$spectra[[which.spectra[1]]][,1], spectra$spectra[[which.spectra[1]]][,ion.col], pch=NA, lty=1, type='l', 
       xlim=new.xlimits, 
       ylim=new.ylimits, 
       col=eic.colors[1],
       xlab="Time (min)", 
       ylab=paste0("Intensity (",ion.col,")"),...)
  for (spectrum.no in 2:length(which.spectra))
  {
    points(spectra$spectra[[which.spectra[spectrum.no]]][,1]+horizontal.offset*(spectrum.no-1), 
           spectra$spectra[[which.spectra[spectrum.no]]][,ion.col]
                   +vertical.offset*(spectrum.no-1)
                   +vertical.offset.relative*(spectrum.no-1)*max.y, 
           type="l", col=eic.colors[spectrum.no])
  }
  if (plot.reference.lines)
  {
    for (line.number in 1:length(reference.lines))
    {
      abline(v=reference.lines[line.number], lty=reference.lty[line.number],
             col=reference.line.colors[line.number])
    }
  }
  if (plot.segment.lines)
  {
    for (line.number in 1:((length(segment.lines)/2-1)))
    {
      segments(segment.lines[i], segment.lines[i+1], segment.lines[i+2], segment.lines[i+3],
               lty=segment.lty,
             col=segment.color)
    }
  }
  
  
  if (show.legend) 
    {
     par(mar=c(0.5,0.5,0.5,0.5))
     plot.new()
     if (length(which.spectra)<12) legend("bottom",legend=spectra$print.names[which.spectra,2], 
            lty=1, col=eic.colors, cex=0.8, bty="n", 
            ncol=max(1,floor(length(which.spectra)/3)))
       else legend("bottom",legend=spectra$print.names[which.spectra,2], 
                   lty=1, col=eic.colors, cex=0.4, bty="n", 
                   ncol=6,
                   inset=c(0,0.4))
     par(mar = c(4, 4, 4, 4) + 0.2)
  }
  #par(old.par)
}

#######################################################################################################
#######################################################################################################
##                                                                                                   ##
##            Function plot.tic: Plot local tic traces for selected spectra                          ##
##                                                                                                   ##
##                                                                                                   ##
#######################################################################################################
#######################################################################################################


plot.tic<-function(spectra, which.spectra=1, xlimits=c(0,30), tic.colors="black",
                   vertical.offset=0, vertical.offset.relative=0, horizontal.offset=0, show.legend=TRUE,...)
{ 
  plot.eic(spectra, which.spectra=which.spectra, ion.to.plot=0, xlimits=xlimits, eic.colors=tic.colors,
          vertical.offset=vertical.offset, vertical.offset.relative=vertical.offset.relative, horizontal.offset=horizontal.offset, show.legend=show.legend, ...)
  
  #note: old code archived below, was replaced when I realized that if I use ion #0 to mean TIC
  #I can keep from having to add options to both functions (plot.tic and plot.eic) when I add
  #functionality
  
  #max.y=0
  #for (i in 1:length(which.spectra))
  #{
  #  max.y<- max( c(max.y,
  #                 max(spectra$spectra[[(which.spectra[i])]]
  #                     [(findInterval(xlimits[1],
  #                                    spectra$spectra[[(which.spectra[i])]][,1]):
  #                                          findInterval(xlimits[2],spectra$spectra[[(which.spectra[i])]][,1])), "TIC"])))    
  #} 
  #  
  #
  #plot(spectra$spectra[[which.spectra[1]]][,1], spectra$spectra[[which.spectra[1]]][,"TIC"], pch=NA, lty=1, type='l', xlim=xlimits, ylim=c(0,max.y), col=tic.colors[1],
  #     xlab="Time (min)", ylab=paste0("Intensity (TIC)"),...)
  #for (spectrum.no in 2:length(which.spectra))
  #{
  #  points(spectra$spectra[[which.spectra[spectrum.no]]][,1], spectra$spectra[[which.spectra[spectrum.no]]][,"TIC"], type="l", col=tic.colors[spectrum.no])
  #}
  #legend("topright",legend=spectra$print.names[which.spectra,2], lty=1, col=tic.colors, cex=0.8)
}



#######################################################################################################
#######################################################################################################
##                                                                                                   ##
##            Function find.peak.eic: find retentiont time of local max for eic                      ##
##                                                                                                   ##
##                                                                                                   ##
#######################################################################################################
#######################################################################################################


find.peak.eic<-function(spectra, which.spectrum=1, ion=0, xlimits=c(0,30), finding.method="simple", verbosity=0, ploteic=FALSE)
{  
  #print(finding.method)
  #print(class(finding.method))
  #print(typeof(finding.method))
  if (ion>0) ion.col=paste0("mz",ion)
  if (ion==0) ion.col="TIC" #if ion==0 then use TIC instead of EIC
  max.y=0
  index.low=findInterval(xlimits[1],spectra$spectra[[(which.spectrum)]][,1])
  index.high=findInterval(xlimits[2],spectra$spectra[[(which.spectrum)]][,1])
  max.x=which.max(spectra$spectra[[(which.spectrum)]][index.low:index.high, ion.col])
  max.y <-    max(spectra$spectra[[(which.spectrum)]][index.low:index.high, ion.col])   
  min.y <-    min(spectra$spectra[[(which.spectrum)]][index.low:index.high, ion.col])     
    
  index.of.match=which(spectra$spectra[[(which.spectrum)]]
                         [(findInterval(xlimits[1],spectra$spectra[[(which.spectrum)]][,1]):
                           findInterval(xlimits[2],spectra$spectra[[(which.spectrum)]][,1])),
                           ion.col]==max.y)
  if (verbosity>0) 
  {
    print(paste0("find.peak.eic for spectrum #", which.spectrum, 
                 "(", spectra$print.names[which.spectrum,2],")",
                 " ion ", ion,
                 " between ", xlimits[1], " and ", xlimits[2], "."))
    print(paste0("findInterval(xlimits[1],spectra$spectra[[(which.spectrum)]][,1]) is ", findInterval(xlimits[1],spectra$spectra[[(which.spectrum)]][,1])))
    print(paste0("findInterval(xlimits[2],spectra$spectra[[(which.spectrum)]][,1]) is ", findInterval(xlimits[2],spectra$spectra[[(which.spectrum)]][,1])))
    print(paste0("Maximum intensity is ", max.y, "."))
    print(paste0("Index of match is ", max.x+index.low-1))
    print(paste0(" whose y-value is ",spectra$spectra[[(which.spectrum)]][max.x+index.low-1,ion.col] ))
    print(paste0(" and whose x-value is ", spectra$spectra[[(which.spectrum)]][max.x+index.low-1,1]))
   #if (verbosity>5) print(spectra$spectra[[which.spectrum]][max.x+index.low-1,c(1,ion.col)])
  }
  
  
  if (finding.method=="simple") 
    {
     if (ploteic)
      {
        plot.eic(spectra, ion=ion, which.spectra=which.spectrum, xlimits=xlimits, main="find.peak.eic")
        abline(v=spectra$spectra[[(which.spectrum)]][max.x+index.low-1,1])
      }
    return(spectra$spectra[[(which.spectrum)]][max.x+index.low-1,1])
    }
  if (finding.method=="gaussian") 
  {
    x=spectra$spectra[[which.spectrum]][index.low:index.high,1]
    r=spectra$spectra[[which.spectrum]][index.low:index.high,ion.col]
    if (verbosity>2) 
     {
      print(x) 
      print(r)
     }
    tab=data.frame(x,r)
    
    if (verbosity>2) print(paste0("Initial guesses:    sigma, ", 0.02, 
                                  "   k, ", (spectra$spectra[[(which.spectrum)]][max.x+index.low-1,ion.col]),
                                  "   mu, ", (spectra$spectra[[(which.spectrum)]][max.x+index.low-1,1])))
    if (ploteic)
     {
      plot.eic(spectra, ion=ion, which.spectra=which.spectrum, xlimits=xlimits, main="find.peak.eic pre-fit")
      lines(as.vector(x), min.y+(spectra$spectra[[(which.spectrum)]][max.x+index.low-1,ion.col])*
              exp(-1/2*(x-(spectra$spectra[[(which.spectrum)]][max.x+index.low-1,1]))^2/(0.02^2)),  col="green")
     }
    
    tryCatch(
      {
        res <- nls( r ~ bg+k*exp(-1/2*(x-mu)^2/sigma^2), 
                 start=c(bg=min.y,mu=(spectra$spectra[[(which.spectrum)]][max.x+index.low-1,1]),
                         sigma=0.02,k=(spectra$spectra[[(which.spectrum)]][max.x+index.low-1,ion.col])) , data = tab)
      },
      error=function(e)
      { 
        print("Failed to fit Gaussian, finding simple max instead.")
        return(find.peak.eic(spectra=spectra, 
                                       which.spectrum=which.spectrum, ion=ion, xlimits=xlimits, 
                                       finding.method="simple", verbosity=verbosity, ploteic=ploteic))
      },
      warning=function(w)
      {})
    
    if (verbosity>0) print (res)
    if (ploteic)
    {
     plot.eic(spectra, ion=ion, which.spectra=which.spectrum, xlimits=xlimits, main="find.peak.eic fit")
     lines(as.vector(x), predict(res), col="blue")
    }
    return(coef(res)[2])
    }
  }

#######################################################################################################
#######################################################################################################
##                                                                                                   ##
##            Function ms.barplot: Plot barplot of integrated peak for selected spectra              ##
##                                                                                                   ##
##                                                                                                   ##
#######################################################################################################
#######################################################################################################

ms.barplot<-function(spectra, filenum=1, peakrt=1, time.range=0, verbosity=0, mz.min=0, mz.max=1000)
{ 
  temp.mean<-spectra$ionlist[[filenum]]*0.0
  for (mz in spectra$ionlist[[filenum]][1]:spectra$ionlist[[filenum]][length(spectra$ionlist[[filenum]])])
   {
     #print(mz)
    
     temp.mean[(mz-(spectra$ionlist[[filenum]][1]))] <- 
         mean(spectra$spectra[[filenum]][
            ((spectra$spectra[[filenum]][,1]>=(peakrt-time.range))&(spectra$spectra[[filenum]][,1]<=(peakrt+time.range))),
             1+findInterval(mz,spectra$ionlist[[filenum]])], na.rm=TRUE)
   }
  if (verbosity>1) print(temp.mean)
  
  barplot(temp.mean, names=spectra$print.names[[filenum]])
  
  return(temp.mean)
}
  
#######################################################################################################
#######################################################################################################
##                                                                                                   ##
##            Function plot.mass.spectrum: Plot barplot type ms              ##
##                                                                                                   ##
##                                                                                                   ##
#######################################################################################################
#######################################################################################################

plot.mass.spectrum<-function(spectra, which.spectrum=1, peakrt=1, verbosity=0, mz.min=0, mz.max=10000, 
                             zoom=1, return.spectrum=FALSE, ...)
{ 
  rt.index=findInterval(peakrt,spectra$spectra[[(which.spectrum)]][,1])
  mz.index.low  = findInterval(mz.min,spectra$ionlist[[(which.spectrum)]])+1
  mz.index.high = findInterval(mz.max,spectra$ionlist[[(which.spectrum)]])
  
  y.values=t(spectra$spectra[[which.spectrum]][rt.index,((mz.index.low+1):(mz.index.high+1))])
  if (verbosity>10) print((spectra$spectra[[which.spectrum]][rt.index,((mz.index.low+1):(mz.index.high+1))])[,1:10])
  if (verbosity>10) print(t(spectra$spectra[[which.spectrum]][rt.index,((mz.index.low+1):(mz.index.high+1))])[1:10,])
  
  
  bar.names=spectra$ionlist[[which.spectrum]][mz.index.low:mz.index.high]
  
  if (verbosity>5) print(cbind(bar.names,y.values)[1:10,])
  if (verbosity>0) print(paste0("rt.index is ", rt.index))
  if (verbosity>0) print(paste0("mz.index.low is ", mz.index.low, " and .high is ", mz.index.high))
  if (verbosity>0) print(paste0("Dim of y.values is ", dim(y.values)))
  if (verbosity>0) print(paste0("Dim of bar.names is ", length(bar.names)))
  
  barplot(y.values[,1], names=bar.names, ylim=c(0,max(y.values[,1])/zoom, ...))
  
  if (return.spectrum) return(cbind(bar.names,y.values[,1]))
  
}



plot.with.labels<-function(x,y,labels,...)
{
  plot(x,y,... )
  textxy(x,y,labels)
  
}

plot.versus.day<-function(column="H2S", which.spectra=0)
{
  if (which.spectra==0) which.spectra=1:length(spectra$cdfFilesDays)
  plot.with.labels(spectra$print.names[which.spectra,"Day"],
                   y.matrix[which.spectra, column],
                   y.matrix[which.spectra, column], xlab="Days", ylab=column,
                   xlim=c(min(spectra$print.names[which.spectra,"Day"]), max(spectra$print.names[which.spectra,"Day"])*1.1),
                   ylim=c(min(0,min(y.matrix[which.spectra,column])),max(y.matrix[which.spectra,column])*1.1))
}


plot.peakgrid<-function(spectra=NULL, which.spectra=1:length(spectra$cdfFiles), ion.to.plot=0, retention.time=5, retention.delta=0.2, 
                        number.rows=3, number.columns=2, width.mult=10, ...)
{
  ppg=number.rows*number.columns # panels per graph page
  
  par(mfrow = c(number.rows, number.columns))
  par(cex = 0.5)
  par(mar = c(6, 4, 4, 4), oma = c(4, 1, 4, 0.5))
  par(tcl = -0.25)
  par(mgp = c(3, 0.3, 0))
  
  for (spectrum in 1:length(spectra$cdfFiles))
  {
    if (spectrum %in% which.spectra) plot.eic(spectra, 
              which.spectra=spectrum, #which spectra numbers to plot (vector)
              ion.to.plot=ion.to.plot,  #which ion to plot, 0 means TIC
              xlimits=c(retention.time-width.mult*retention.delta,
                        retention.time+width.mult*retention.delta), #time limits to plot
              eic.colors="black", #colors to use
              vertical.offset=0, # absolute vertical offset for each line
              vertical.offset.relative=0,#vertical offset for each line, 
              #as a fraction of the y.max
              horizontal.offset=0, #absolute horizontal offset for each line
              show.legend=FALSE,
              verbosity=0,  # set to higher integer to print more info for debugging
              plot.reference.lines=FALSE,
              reference.lines=c(0,0,0),
              ...)
  }
  
}
#######################################################################################################
#######################################################################################################
##                                                                                                   ##
##            Function generate.ionlist                                                              ##
##                                                                                                   ##
##                                                                                                   ##
#######################################################################################################
#######################################################################################################
generate.ionlist<-function(spectra=NULL, which.spectra=0, number=20, lower.limit=0, upper.limit=10000, verbosity=0, plot.it=FALSE)
{
  min.ion=max(lower.limit,min(spectra$masterionlist))
  max.ion=min(upper.limit,max(spectra$masterionlist))
  ion.max=rep(0,(max.ion-min.ion)+1)
  ion.min=rep(0,(max.ion-min.ion)+1)
  ion.diff=rep(0,(max.ion-min.ion)+1)
  ions=min.ion:max.ion
  for (counter in min.ion:max.ion)
  
  if (which.spectra==0) which.spectra=1:length(spectra$cdfFiles)
  for (ion in min.ion:max.ion)
  {
   for (i in 1:length(which.spectra))
   {
    ion.max[ion-min.ion+1]<- max( c(ion.max[ion-min.ion+1], max(spectra$spectra[[(which.spectra[i])]][,paste0("mz",ion)])))
    ion.min[ion-min.ion+1]<- min( c(ion.max[ion-min.ion+1], min(spectra$spectra[[(which.spectra[i])]][,paste0("mz",ion)])))
                                            
   }                   
  }
  ion.diff=ion.max-ion.min
  if (verbosity>8) print(ion.min)
  if (verbosity>8) print(ion.max)
  if (verbosity>8) print(ion.diff)
  if (verbosity>8) print(ions)
  if (plot.it==TRUE)
  {
    length(ions)
    length(ion.diff)
    plot(ions, ion.diff, pch=20)
  }
  if (verbosity>8) print("Sorted differences:")
  if (verbosity>8) print(sort(unique(ion.diff), decreasing=TRUE))
  
  if (verbosity>8) print("Indices of sorted differences:")
  if (verbosity>8) print(which(ion.diff == sort(unique(ion.diff), decreasing=TRUE)))
  
  if (verbosity>8) print("Returning ionlist")
  ion.list=ions[which(ion.diff == sort(unique(ion.diff),decreasing=TRUE))]
  return(ions[order(ion.diff, decreasing=TRUE)[1:number]])
  
}

#######################################################################################################
#######################################################################################################
##                                                                                                   ##
##            Function ms.overlayplot: Create tic overlay plots from 'spectra' list object   ##
##                                                                                                   ##
##                                                                                                   ##
#######################################################################################################
#######################################################################################################

ms.overlayplot <- function (spectra=NULL, #a list object made by the loading routines
                            tic.colors=-1,
                            stackzoom=1,
                            which.spectra=0, #a vector of indices of the spectra in 'spectra' to plot
                            chromatogramxlim=c(-1,-1), 
                            verbosity=0, #set to higher number to be more verbose in printing info
                            zoom=1, # zoom factor, increase to zoom y-axis
                            retentiontimelist=NA, # list of retention time at which to draw vertical lines
                            showlegend=TRUE,
                            margins=c(10,6,4,4),
                            showaxes=TRUE,
                            tic.line.width=1, #width of TIC line
                            yaxisside=2,
                            showyaxis=TRUE,
                            xaxisside=1,# 2 is left, 4 is right
                            showxaxis=TRUE,
                            showbox=TRUE,
                            showtitle=TRUE,
                            ticcolor="black",
                            main="",
                            correct.split=FALSE,
                            groups=0,  #if replicates, each is a group as designated by this factor vector
                            ...)
{
  if ((length(which.spectra)>1)&&(colors=1)) colors=rainbow(length(which.spectra))
  if ((length(which.spectra)<2)&&(colors=1)) colors=rainbow(length(spectra$cdfFiles))
  old.par=par(no.readonly=TRUE) #record graphic parameters so we can politely put them back that way when done
  if (correct.split==FALSE) split.factors<-rep(1.0,length(spectra$cdfFiles))
  if (correct.split==TRUE) split.factors<-as.double(spectra$print.names[,"split"])
  if (tic.colors==-1)
  {
    if (length(groups)>1) tic.colors=rainbow(length(unique(as.integer(groups[which.spectra]))))
     else tic.colors=rainbow(length(spectra$cdfFiles))
  }
  if (tic.colors==0)
  {
    if (length(groups)>1) tic.colors=rep("black",length(unique(as.integer(groups[which.spectra]))))
    else tic.colors=rep("black", length(spectra$cdfFiles))
  }
  #find maximum TIC y-value overall
  maxtic=0
  for (i in 1:length(spectra$cdfFiles))  ###need tic max for each and every ion in ionlist
  {
    if((which.spectra[1]==0)||(i %in% which.spectra))
    {
      maxtic=#max(maxtic,max(spectra$spectra[[i]][,"TIC"]))
        max( c(maxtic,
               split.factors[i]*max(spectra$spectra[[i]][(findInterval(chromatogramxlim[1],spectra$spectra[[i]][,1]):
                                                    findInterval(chromatogramxlim[2],spectra$spectra[[i]][,1])),
                                                 "TIC"]    )))
      if ((chromatogramxlim[1]==-1) && (chromatogramxlim[2]==-1))
        maxtic=max( c(maxtic, split.factors[i]*max(spectra$spectra[[i]][,"TIC"])))
    }
  }
  #maxtic=max(unlist(lapply(spectra$spectra,FUN=max)))   #get the max value of all the matrices, which should be the max TIC
  maxticplus=maxtic*1.1
  if (verbosity>3) print(paste0("Max TIC: ", maxtic))
  
  mintime=0 #min(spectra$spectra[[1]][,1])
  maxtime=-1 #max(spectra$spectra[[1]][,1])
  for (i in 1:length(spectra$cdfFiles))
  {
    if((which.spectra[1]==0)||(i %in% which.spectra))
    {
      mintime=min(mintime,min(spectra$spectra[[i]][,1]))
      maxtime=max(maxtime,max(spectra$spectra[[i]][,1]))
      
    }
  }
  
  # determine number of pools required
  if (length(groups)<2) number.of.pools=length(which.spectra)
  if (length(groups)>1)
  {
    if (which.spectra[1]!=0) number.of.pools=nlevels(groups[which.spectra])
    if (which.spectra[1]==0) number.of.pools=nlevels(groups[])
  }
  if (which.spectra[1]==0) number.of.pools=length(spectra$cdfFiles)
  if (length(groups)>1) 
  {
    
    if ((which.spectra[1])!=0)
    {
      #trim the groups vector down to just those entries to be printed
      if(length(groups)==length(spectra$cdfFiles)) groups<-as.factor(as.character(groups[which.spectra]))
    }
    number.of.pools=length(levels(groups))
  }
  if (verbosity>5) print(paste0("Number of pools: ", number.of.pools))
  
  #reset graph margins to a reasonable value
  par(mar = margins + 0.2, oma=c(1,1,1,1))
  
  print(paste0("Min time: ", round(mintime,3),"  Max time: ", round(maxtime,3)))
  if (chromatogramxlim[2]==-1) chromatogramxlim=c(round(mintime,3),round(maxtime,3))
  
  scaleadj=1 #length(cdfFiles)
  
  stackedfontspacer=(maxtic*number.of.pools/stackzoom)*0.011 #was 0.13 before 1-31-17
  #place labels this fraction of whole graph below traces
  
  if (verbosity>7) print(paste0("Groups as of right now"))
  if (verbosity>7) print(groups)
  if (verbosity>7) print(paste0("Which.spectra as of right now"))
  if (verbosity>7) print(which.spectra)
  if (verbosity>7) print(paste0("Groups[which.spectra] as of right now"))
  if (verbosity>7) print(groups[which.spectra])
  if (verbosity>9) print(paste0("spectra$print.names[,\"group.auto\"]"))
  if (verbosity>9) print(spectra$print.names[,"group.auto"])
  
  
  ############################ MAIN LOOP
  ###########################################
  
  printed.filenum=0
  for (filenum in 1:length(spectra$cdfFiles)) # for every file in the spectra list (note: not the which.spectra list)
  {
    
    if((which.spectra[1]==0)||(filenum %in% which.spectra)) #if it is a spectrum we are supposed to plot
    {
      ordered.file=which(which.spectra==filenum)  #which entry of which.spectrum are we working on
      printed.filenum=printed.filenum+1 #how many spectra have we printed
      if (verbosity>0) print(paste0("File ",filenum," counts, adding as printed.filenum: ", printed.filenum,"."))
      which.pool<-printed.filenum
      if (verbosity>9) print(paste0("Initial which.pool set to ", which.pool))
      
      if (length(groups)>1) which.pool<-which(as.integer(groups[ordered.file])==unique(as.integer(groups)))
      #length(levels(as.factor(as.integer(groups[1:printed.filenum])))) #how many levels up in groups is pool number
      
      if (verbosity>9) print(paste0("Which.pool set to ", which.pool, " after check for groups"))
      verticaloffset=0 #(which.pool-1)*maxticplus/stackzoom #how high do we plot it based on which pool it goes to
      
      if(verbosity>5) print(paste0("pool ", which.pool, " Vertical offset: ", verticaloffset))
      #if (length(groups)==1) 
      #{
      #  if (which.spectra[1]!=0) verticaloffset=(match(filenum, which.spectra)-1)*maxticplus/stackzoom #if no groups but which.spectra specifies subset to plot
      #}
      
      graphtop=maxticplus#*number.of.pools/stackzoom
      
      #abline(h=verticaloffset) 
      if (printed.filenum==1)#(filenum==1)||(filenum==min(which.spectra))) # if this is the first spectrum we are plotting, then we need to generate the whole plot not just add line
      {
        #if (stackedshowtic) 
        {
          plot(spectra$spectra[[filenum]][,1]+spectra$retention.adjust[[filenum]][1], 
               split.factors[filenum]*spectra$spectra[[filenum]][,"TIC"]*stackzoom+verticaloffset,
               type="l", 
               xlab="Time (min)", ylab="TIC Intensity", ylim=c(0,graphtop), 
               xlim=chromatogramxlim, lwd=tic.line.width, col=tic.colors[which.pool],
               axes=FALSE, ann=FALSE,...)
          if (showxaxis) {axis(xaxisside); title(xlab="Time (min)")}
          if (showyaxis) {axis(yaxisside); title(ylab="Intensity")}
          if (showbox) box()
          if (showtitle) title(main)
        }
      }
      
      if (filenum>1) {lines(spectra$spectra[[filenum]][,1]+spectra$retention.adjust[[filenum]], 
                                                split.factors[filenum]*spectra$spectra[[filenum]][,"TIC"]*stackzoom+verticaloffset, 
                                                col=tic.colors[which.pool], lwd=tic.line.width)}
      
      #if (stackedshowtic) 
        lines(spectra$spectra[[filenum]][,1]+spectra$retention.adjust[[filenum]], split.factors[filenum]*spectra$spectra[[filenum]][,"TIC"]*stackzoom+verticaloffset, type="l", 
                                lwd=tic.line.width, col=tic.colors[which.pool])
      #if (stackedshowtic) lines(stacked_scan_acquisition_times[,filenum]+stackedrtadj[filenum], stacked_total_intensities[,filenum]*stackzoom+verticaloffset, lwd=tic.line.width*2)
      
      tracelabel=""
      
      #if (showfilenames) 
        tracelabel=paste0(tracelabel, spectra$print.names[filenum,"print.name"], " ")
      if ((length(groups)>1)) 
      {
        tracelabel=paste(spectra$print.names[which.spectra[ which(groups==groups[which(which.spectra==filenum)])],"print.name"], collapse=", ") 
      }
      #if (stackeddayconversion) tracelabel=paste0(tracelabel, "Day ", spectra$cdfFilesDays[filenum], " ")
      
      
      #abline(h=verticaloffset, lwd=1) 
        
      
      
    }
  }
  #title(paste0("CDF Files in ", stackpath))
  if (verbosity>1) print(!is.na(retentiontimelist[1]))
  if (!is.na(retentiontimelist[1]))
  {
    for (i in 1:length(retentiontimelist)) abline(v=retentiontimelist[i], lty="dotted")
  }
  
  if ((stackedalignchrom)**(stackedplotalignrt)) abline(v=stackedrtmax[1], lty="dotted", col="gray")
  
  if (length(groups)<2)
  {
   if (length(which.spectra)<12) legend("bottom",legend=spectra$print.names[which.spectra,2], 
                                       lty=1, col=tic.colors, cex=0.8, bty="n", xpd=TRUE,
                                       inset=c(0,-0.15),
                                       ncol=max(1,floor(length(which.spectra)/3)))
   else legend("bottom",legend=spectra$print.names[which.spectra,2], 
              lty=1, col=tic.colors, cex=0.6, bty="n", xpd=TRUE,
              ncol=5,
              inset=c(0,-0.15))
  } else
  {
    if (length(which.spectra)<12) legend("bottom",legend=levels(groups)[unique(as.integer(groups))], 
                                         lty=1, col=tic.colors[1:length(unique(as.integer(groups)))], 
                                         cex=0.8, bty="n", xpd=TRUE, 
                                         inset=c(0,-0.2),
                                         ncol=3)
    else legend("bottom",legend=levels(groups)[unique(as.integer(groups))], 
                lty=1, col=tic.colors[1:length(unique(as.integer(groups)))], 
                cex=0.6, bty="n", xpd=TRUE,
                ncol=5, 
                inset=c(0,-0.2))
  }
  par(old.par)
  return("Overlay Plot Finished.")
}


#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
show.integration<-function(spectra, decomp.ions.show)
{
  for (cmpd.no in 1:dim(decomp.ions.show)[1])
  {
    plot.eic(spectra, which.spectra=1:length(spectra$cdfFiles), 
             ion.to.plot=decomp.ions.show[cmpd.no,"Ion.1"], #was decomp.ions.rev, not sure why
             xlimits=c(decomp.ions.show[cmpd.no,"Retention.Time"]-5*decomp.ions.show[cmpd.no,"Retention.Delta"],
                       decomp.ions.show[cmpd.no,"Retention.Time"]+5*decomp.ions.show[cmpd.no,"Retention.Delta"]),
             plot.reference.lines=TRUE,
             reference.lines=c(decomp.ions.show[,"Retention.Time"],
                               decomp.ions.show[cmpd.no,"Retention.Time"],
                               decomp.ions.show[cmpd.no,"Retention.Time"]-decomp.ions.show[cmpd.no,"Retention.Delta"],
                               decomp.ions.show[cmpd.no,"Retention.Time"]+decomp.ions.show[cmpd.no,"Retention.Delta"]
                               ),
             reference.line.colors=c(rep("black", length(decomp.ions.show[,"Retention.Time"])),"black", "red", "red"),
             reference.lty=c(rep(2, length(decomp.ions.show[,"Retention.Time"])),decomp.ions.show[,"Retention.Time"],1,3,3),
             main=decomp.ions.show[cmpd.no,"Long.Name"]
             
             )
    mtext(paste0("Ion ", decomp.ions.show[cmpd.no,"Ion.1"],
               " m/z, RT = ", decomp.ions.show[cmpd.no,"Retention.Time"],
               " min"))
  }
}


#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
test.integration<-function(spectra, decomp.ions.show)
{
  dch.create.pdf("Test Integrations ")
  old.par<-par(no.readonly=TRUE)
  
  #par(mfrow=c(3,4), cex=0.6, cex.axis=0.3, cex.lab=0.5, cex.main=0.6, cex.sub=0.4)
  
  #new.pars=par(no.readonly=TRUE)
  
  for (cmpd.no in 1:dim(decomp.ions.show)[1])
  {
    for (ion.no in 1:3)
    {
    plot.eic(spectra, which.spectra=1:length(spectra$cdfFiles), 
             ion.to.plot=decomp.ions.rev[cmpd.no, paste0("Ion.", ion.no)], 
             xlimits=c(decomp.ions.show[cmpd.no,"Retention.Time"]-5*decomp.ions.show[cmpd.no,"Retention.Delta"],
                       decomp.ions.show[cmpd.no,"Retention.Time"]+5*decomp.ions.show[cmpd.no,"Retention.Delta"]),
             plot.reference.lines=TRUE,
             reference.lines=c(decomp.ions.show[,"Retention.Time"],
                               decomp.ions.show[cmpd.no,"Retention.Time"],
                               decomp.ions.show[cmpd.no,"Retention.Time"]-decomp.ions.show[cmpd.no,"Retention.Delta"],
                               decomp.ions.show[cmpd.no,"Retention.Time"]+decomp.ions.show[cmpd.no,"Retention.Delta"]
             ),
             reference.line.colors=c(rep("black", length(decomp.ions.show[,"Retention.Time"])),"black", "red", "red"),
             reference.lty=c(rep(2, length(decomp.ions.show[,"Retention.Time"])),decomp.ions.show[,"Retention.Time"],1,3,3),
             main=decomp.ions.show[cmpd.no,"Long.Name"]#,
             #use.par=TRUE,
             #par.to.use=new.pars
             
    )
    #mtext(paste0("Ion ", decomp.ions.rev[cmpd.no,paste0("Ion.", ion.no)],
    #             " m/z, RT = ", decomp.ions.rev[cmpd.no,"Retention.Time"],
    #             " min"), side=3)
    }
  }
  par(old.par)
  dev.off()
}

#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################

fatty.acid.tms.analysis<-function(spectra, which.spectrum=1, use.pdf=TRUE)
{
  if (use.pdf) dch.create.pdf("Fatty Acid Analysis of GC-MS Data ")
  #par(mfrow=c(5,1), mar=c(4,4,2,2))
  key.ion=c(242+(0:4)*28)
  key.ion=c(key.ion, key.ion-2, key.ion-4, key.ion-6, key.ion-8, key.ion+16)
  ms.stackedplot(spectra, which.spectra=which.spectrum, ionlist=key.ion, main="Fatty Acid Analsyis", stackedshoweic=TRUE)
}


#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
integration.matrix<-function(spectra, decomp.ions.show, debug.level=10)
{
  par(mfrow=c(1,1))
  par(cex = 0.5)
  par(mar = c(3.5, 3.5, 0, 0), oma = c(1, 1, 1, 1))
  par(tcl = -0.25)
  par(mgp = c(3, 0.3, 0))
  for (cmpd.no in 1:dim(decomp.ions.show)[1])
  {
    for (ion in 1:5)
    {
     if (!is.na(decomp.ions.show[cmpd.no,paste0("Ion.",ion)]))
     {
      if (debug.level>5) print(paste0("Cmpd ", cmpd.no,"  Ion ", ion,
                               "  ", decomp.ions.show[cmpd.no,paste0("Ion.",ion)]))
      plot.eic(spectra, which.spectra=1:length(spectra$cdfFiles), 
             ion.to.plot=decomp.ions.show[cmpd.no,paste0("Ion.",ion)], 
             xlimits=c(decomp.ions.show[cmpd.no,"Retention.Time"]-5*decomp.ions.show[cmpd.no,"Retention.Delta"],
                       decomp.ions.show[cmpd.no,"Retention.Time"]+5*decomp.ions.show[cmpd.no,"Retention.Delta"]),
             plot.reference.lines=TRUE,
             reference.lines=c(decomp.ions.show[,"Retention.Time"],
                               decomp.ions.show[cmpd.no,"Retention.Time"],
                               decomp.ions.show[cmpd.no,"Retention.Time"]-decomp.ions.show[cmpd.no,"Retention.Delta"],
                               decomp.ions.show[cmpd.no,"Retention.Time"]+decomp.ions.show[cmpd.no,"Retention.Delta"]
             ),
             reference.line.colors=c(rep("grey30", length(decomp.ions.show[,"Retention.Time"])),"green", "red", "red"),
             reference.lty=c(rep(3, length(decomp.ions.show[,"Retention.Time"])),decomp.ions.show[,"Retention.Time"],1,3,3),
             main=decomp.ions.show[cmpd.no,"Long.Name"],
             sub=paste0("Ion ", decomp.ions.show[cmpd.no,paste0("Ion.",ion)],
                        " m/z, RT = ", decomp.ions.show[cmpd.no,"Retention.Time"],
                        " min")
             
                )
     #mtext(paste0("Ion ", decomp.ions.show[cmpd.no,paste0("Ion.",ion)],
    #             " m/z, RT = ", decomp.ions.show[cmpd.no,"Retention.Time"],
     #            " min"), side=3)
     }
     #else (plot.new())
     #plot.new() 
    }
  }
}

#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
image.spectrum<-function(spectra, spectrum.to.plot=1, debug.level=10)
{
  par(fig=c(0,1,0,0.3), mar=c(2,2,0,2), oma=c(0,0,0,0), new=TRUE)
  plot(spectra$spectra[[spectrum.to.plot]][,1],
       (spectra$spectra[[spectrum.to.plot]][,"TIC"]), 
       lty=1, type="l", xaxs="i", xlab="Retention Time (min)", 
       ylim=c(0,1.05*max(spectra$spectra[[spectrum.to.plot]][,"TIC"])))
  
  par(fig=c(0,1,0.3,1), mar=c(0,2,2,2), oma=c(0,0,0,0), new=TRUE)
  plot(1, type="n", xlab="", ylab="", 
       xlim=c(min(spectra$spectra[[spectrum.to.plot]][,1]), max(spectra$spectra[[spectrum.to.plot]][,1])),
       ylim=c(0, 600))
  max.int=max(spectra$spectra[[spectrum.to.plot]][,-c(1,dim(spectra$spectra[[spectrum.to.plot]])[2])])
  my.colors=rainbow(100)
  
  for (col in 2:(dim(spectra$spectra[[spectrum.to.plot]])[2]-1))
  {
   mz.value=gsub('mz([0-9]+).*','\\1',colnames(spectra$spectra[[spectrum.to.plot]])[col])
   for (row in 1:(dim(spectra$spectra[[spectrum.to.plot]])[1]))
   {
     points(spectra$spectra[[spectrum.to.plot]][row,1],
            mz.value,
            pch=1, cex=0.1,
            col=my.colors[floor(spectra$spectra[[spectrum.to.plot]][row, col]/max.int*100)])
   }
  }

  

  
}

#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
write.tic.csv <- function (spectra, filename="TIC matrix.csv", which.spectra=0, 
                           return.data=TRUE, verbosity=0)
{
  tic.matrix.temp<-tic.matrix(spectra, which.spectra, verbosity)
  write.csv(tic.matrix.temp, file=filename)
  if (return.data) return(tic.matrix.temp)
}

tic.matrix <- function (spectra, which.spectra=0, verbosity=0)
{
 if (which.spectra==0) which.spectra=seq(1, length(spectra$spectra))
 
 for (i in 1:length(which.spectra))
 {
   if (verbosity>15) print (paste0("write.tsv.csv iteration #", i))
                            
   if (i==1)
     {
      tic.matrix<-cbind(spectra$spectra[[i]][,1], spectra$spectra[[i]][,"TIC"])
      if (verbosity>50) print(paste0("Dim of initial tic.matrix = ", dim(tic.matrix)))
      colnames(tic.matrix)<-c("Temp.Time", "Temp.TIC") # to avoid colname setting below from freaking out
     } 
   else
     {
     tic.matrix<-cbind(tic.matrix, spectra$spectra[[i]][,1], spectra$spectra[[i]][,"TIC"])
     }
   if (verbosity>10) print (tic.matrix[1:5,])
   if (verbosity>20) print (colnames(tic.matrix))
   
   colnames(tic.matrix)[2*i-1]<-paste0(spectra$print.names$print.name[i], ".Time")
   colnames(tic.matrix)[2*i-0]<-paste0(spectra$print.names$print.name[i], ".TIC")
   
   
 }
 return(tic.matrix)
}