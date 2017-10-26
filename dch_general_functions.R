#######################################################################################################
#######################################################################################################
##                                                                                                   ##
##            dch_general_functions.R                                                                ##
##                 Some general function by Dr. Haines, SHSU                                         ##
##                                                                                                   ##
#######################################################################################################
#######################################################################################################

library(StatDA)
if (!exists("dch_general_functions.loaded")) dch_general_functions.loaded=TRUE
dch.summary <- function(x)
{
  print(paste0("Statistical Summary for ",deparse(substitute(x))))
  edaplot(x, H.freq=FALSE)
  print(paste0("Mean : ", mean(x)))
  print(paste0("Median : ", median(x)))
  print(paste0("Standard deviation : ", sd(x)))
  print(paste0("Quantiles : "))
  print(quantile(x))
  print(paste0("Interquartile Range (IQR) : ", IQR(x)))
  print(paste0("s(IQR) = 0.7413*IQR : ", 0.7413*IQR(x)))
  print(paste0("Median Absolute Deviation (MAD) : ", mad(x)))
  print(paste0("s(MAD) = 1.483 * MAD : ", 1.483*mad(x)))
  
}

dch.compare.distributions <- function(x,y, 
                                      x.name=deparse(substitute(x)),
                                      y.name=deparse(substitute(y)),
                                      verbosity=0)
{ 
  
  P.x.main<-print(paste0("Statistical Summary for ",x.name))
  P.y.main<-print(paste0("Statistical Summary for ",y.name))
  
  P.xlim.min<-min(min(x),min(y))
  P.xlim.max<-max(max(x),max(y))
  P.xlim.margin<-(P.xlim.max-P.xlim.min)*0.05
  
  P.xlim.min<-P.xlim.min-P.xlim.margin
  P.xlim.max<-P.xlim.max+P.xlim.margin
  par(mfrow=c(2,1))
  edaplot(x, H.freq=FALSE, P.xlim=c(P.xlim.min,P.xlim.max), P.main=P.x.main, P.xlab=x.name)
  edaplot(y, H.freq=FALSE, P.xlim=c(P.xlim.min,P.xlim.max), P.main=P.y.main, P.xlab=y.name)
  
  x.shapiro<-shapiro.test(x)
  y.shapiro<-shapiro.test(y)
  print(paste0("Shapiro-Wilk Normality Test on ", x.name))
  print(x.shapiro)
  if (x.shapiro$p.value<0.05) print("WARNING: Likely not a normal distribution!!!!!")
  print("")
  print(paste0("Shapiro-Wilk Normality Test on ", y.name))
  print(y.shapiro)
  if (y.shapiro$p.value<0.05) print("WARNING: Likely not a normal distribution!!!!!")
  print("")
  
  my.t<-t.test(x,y)
  
  print(my.t)
  par(mfrow=c(1,1))
  return(my.t$p.value)
}

temp.shapiro<-function()
{
  test.data<-rnorm(5000)
  edaplot(test.data, H.freq=FALSE)
  print(shapiro.test(test.data))
}
#######################################################################################################
#######################################################################################################
##                                                                                                   ##
##            Function pca.report                                  ##
##                                                                                                   ##
##                                                                                                   ##
#######################################################################################################
#######################################################################################################
pca.report<-function(data,
                     scale=TRUE,
                     center=TRUE,
                     plots=TRUE,
                     number.rows.strip=5,
                     row.factors=rep(1, dim(data)[1]),
                     column.factors=rep(1, dim(data)[2]),
                     column.factors2=rep(1, dim(data)[2]),
                     scores.to.plot=100,
                     order.by.rotation=FALSE,
                     pcs.to.plot=25,
                     rotation.strips.to.plot=100, #maximum strip to plot
                     verbosity=0,
                     barplot.margins=c(20, 4, 4, 4))
{
  my.pca<-prcomp(data,scale=scale,center=center)
  text.page(capture.output(summary(my.pca)))
  if (verbosity>0) print(summary(my.pca))
  if (plots)
  {
    par(mfrow = c(1,1))
    par(cex = 0.5)
    par(mar = c(8, 8, 8, 8), oma = c(4, 1, 4, 0.5))
    par(tcl = -0.25)
    par(mgp = c(3, 0.3, 0))
   plot(my.pca, main="Scree Plot")
    
   
    ## Score plots
   number.rows=1
   number.columns=1
   ppg=number.rows*number.columns # panels per graph page
   
   par(mfrow = c(number.rows, number.columns))
   par(cex = 0.8)
   par(mar = barplot.margins, oma = c(4, 1, 4, 0.5))
   if (verbosity>1) print(barplot.margins)
   par(tcl = -0.25)
   par(mgp = c(3, 0.3, 0))
   
    number.bars=min(dim(my.pca$x)[1],scores.to.plot)
    for (pc in 1:min(pcs.to.plot, dim(my.pca$x)[2]))
    {
      if (!order.by.rotation) 
         {
     
      barplot(my.pca$x[1:number.bars,pc], names.arg=rownames(my.pca$x)[1:number.bars],
              main=paste0("Scores for PC", pc), 
              ylab="Score", las=2, cex.names = 0.6)
      
         } else
         {
           #temp.rotation  <- my.pca$rotation[order(my.pca$rotation[,pc],decreasing=TRUE),pc]
           barplot(my.pca$x[order(my.pca$x[,pc],decreasing=TRUE),pc], names.arg=rownames(my.pca$x)[order(my.pca$x[,pc],decreasing=TRUE)],
                   main=paste0("Scores for PC", pc), 
                   ylab="Rotation", las=2, cex.names = 0.6)
         }
    }
    ##Rotation strip charts
   
    number.rows=5
    number.columns=2
    ppg=number.rows*number.columns # panels per graph page
    
    par(mfrow = c(number.rows, number.columns))
    par(cex = 0.8) #was 0.5
    par(mar = c(4, 6, 0.5, 4), oma = c(4, 1, 4, 0.5))
    par(tcl = -0.25)
    par(mgp = c(3, 0.3, 0))
   for (i in 1:min(dim(my.pca$rotation)[2], rotation.strips.to.plot))
   {
    
     plot(my.pca$rotation[,i], as.integer(as.factor(column.factors)), 
         col=as.integer(as.factor(column.factors)), 
         cex=0.7+(2.5*as.integer(as.factor(column.factors2)))/nlevels(as.factor(column.factors2)),
         ann=FALSE, ylim=c(0,1+nlevels(as.factor(column.factors))), yaxt='n',
         xlim=c(-(max(abs(my.pca$rotation[,i])))*1.1,(max(abs(my.pca$rotation[,i])))*1.1) #center x-axis at 0
          )
    title(xlab=paste0("PC", i), line=2) 
    boxplot(my.pca$rotation[,i]~as.integer(as.factor(column.factors)), 
            names=paste0(levels(as.factor(column.factors))," "),
            ylim=c(-(max(abs(my.pca$rotation[,i])))*1.1,(max(abs(my.pca$rotation[,i])))*1.1) ,#center x-axis at 0
            horizontal=TRUE, las=1)
    #}
    #}
    #}
   
    #mtext("Black is Phase I Chicken, Red is Phase I Soil, Green is Phase II Human, Blue is Phase II Soil. Time series small to large symbol.", outer=TRUE, side=1, cex=0.7)
  }
  
  par (mar=c(12,4,4,4)+0.5, xpd=TRUE, mfrow=c(1,1))
  ################## now biplots

  
  for (page in 1:floor((dim(my.pca$rotation)[2]-number.columns)/number.rows))
  {
    number.rows=6
    number.columns=5
    ppg=number.rows*number.columns # panels per graph page
    
    par(mfrow = c(number.rows, number.columns))
    par(cex = 0.5)
    par(mar = c(3.5, 3.5, 0, 0), oma = c(4, 2, 4, 4))
    par(tcl = -0.25)
    par(mgp = c(3, 0.3, 0))
    for (row in 1:number.rows)
    {
      for (column in 1:number.columns)
      { 
        if (verbosity>8) print(paste0("Row ", row, " Column ", column))
        if (
            ((((page-1)*number.rows)+row)<=dim(my.pca$rotation)[2])&&
            ((((page-1)*number.rows)+row+column)<=dim(my.pca$rotation)[2])
            )
        {
         plot(my.pca$rotation[,(page-1)*number.rows+row], my.pca$rotation[,(page-1)*number.rows+row+column], 
             col=as.integer(as.factor(column.factors)), 
             cex=0.5+(as.integer(as.factor(column.factors2))/nlevels(as.factor(column.factors2))),
             ann=FALSE,
             xlim=c(-(max(abs(my.pca$rotation[,(page-1)*number.rows+row])))*1.1,(max(abs(my.pca$rotation[,(page-1)*number.rows+row])))*1.1), #center x-axis at 0)
             ylim=c(-(max(abs(my.pca$rotation[,(page-1)*number.rows+row+column])))*1.1,(max(abs(my.pca$rotation[,(page-1)*number.rows+row+column])))*1.1)) #center x-axis at 0
         title(xlab=paste0("PC", (page-1)*number.rows+row), line=1.5) 
         title(ylab=paste0("PC", (page-1)*number.rows+row+column), line=1.5)
        }
        else plot.new()
      }
    }
    #Make a blank full page overaly plot so we can add legend
    par(fig=c(0, 1, 0, 1), oma=c(1, 1, 1, 1), mar=c(0, 0, 0, 0), new=TRUE)
    plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
    legend("bottom", pch=1, col=1:nlevels(as.factor(column.factors)),
           legend=levels(as.factor(column.factors)), ncol=min(5,nlevels(as.factor(column.factors))))
    mtext("Rotations for Top PCs", outer=TRUE, line=2)
    #mtext("Black is Phase I, Red is Phase II, Green is special. Time series small to large symbol.", outer=TRUE, side=1, cex=0.7)
  }
  }
  return(my.pca)
}


#######################################################################################################
#######################################################################################################
##                                                                                                   ##
##            Function pca.report.2                                                         ##
##                                                                                                   ##
##                                                                                                   ##
#######################################################################################################
#######################################################################################################
pca.report.2<-function(data,
                     scale=TRUE,
                     center=TRUE,
                     plots=TRUE,
                     number.rows.strip=5,
                     row.factors=rep(1, dim(data)[1]),
                     row.factors2=rep(1, dim(data)[1]),
                     column.factors=rep(1, dim(data)[2]),
                     column.factors2=rep(1, dim(data)[2]),
                     scores.to.plot=100,
                     pcs.to.plot=25,
                     order.by.rotation=FALSE,
                     rotation.strips.to.plot=0, #maximum strips to plot
                     score.strips.to.plot=100, #maximum strips to plot
                     verbosity=0,
                     barplot.margins=c(20, 4, 4, 4))
{
  my.pca<-prcomp(data,scale=scale,center=center)
  text.page(capture.output(summary(my.pca)))
  if (verbosity>0) print(summary(my.pca))
  if (plots)
  {
    par(mfrow = c(1,1))
    par(cex = 0.5)
    par(mar = c(8, 8, 8, 8), oma = c(4, 1, 4, 0.5))
    par(tcl = -0.25)
    par(mgp = c(3, 0.3, 0))
    
    ####################################################################
    #### Custom Scree code derived from http://www.sthda.com/english/wiki/principal-component-analysis-in-r-prcomp-vs-princomp-r-software-and-data-mining
    
    # Eigenvalues
    eig <- (my.pca$sdev)^2
    # Variances in percentage
    variance <- eig*100/sum(eig)
    # Cumulative variances
    cumvar <- cumsum(variance)
    eig.my.pca <- data.frame(eig = eig, variance = variance,
                                        cumvariance = cumvar)
    
    df.bar<-barplot(eig.my.pca[, 2], names.arg=1:nrow(eig.my.pca), 
            main = "Variance",
            xlab = "Principal Component",
            ylab = "Percentage of variance",
            col ="steelblue")
    # Add connected line segments to the plot
    
    lines(x = df.bar, 
          eig.my.pca[, 2], 
          type="b", pch=19, col = "red")
    
    ###plot(my.pca, main="Scree Plot")
    
    #####################################################################
    #####################################################################
    #####################################################################
    #####################################################################
    
    ## Rotation plots
    number.rows=1
    number.columns=1
    ppg=number.rows*number.columns # panels per graph page
    
    par(mfrow = c(number.rows, number.columns))
    par(cex = 0.8)
    par(mar = barplot.margins, oma = c(4, 1, 4, 0.5))
    if (verbosity>1) print(barplot.margins)
    par(tcl = -0.25)
    par(mgp = c(3, 0.3, 0))
    
    number.bars=min(dim(my.pca$rotation)[1],scores.to.plot)
    for (pc in 1:min(pcs.to.plot, dim(my.pca$rotation)[2]))
    {
      #if (!order.by.rotation) 
      #   {
          barplot(my.pca$rotation[1:number.bars,pc], names.arg=rownames(my.pca$rotation)[1:number.bars],
              main=paste0("Rotation for PC", pc), 
              ylab="Rotation", las=2, cex.names = 0.6)
      #   } else
      #   {
      #     temp.rotation  <- my.pca$rotation[order(my.pca$rotation[,pc],decreasing=TRUE),pc]
      #     barplot(my.pca$rotation[order(my.pca$rotation[,pc],decreasing=TRUE),pc], names.arg=rownames(my.pca$rotation)[order(my.pca$rotation[,pc],decreasing=TRUE)],
      #             main=paste0("Rotation for PC", pc), 
      #             ylab="Rotation", las=2, cex.names = 0.6)
      #   }
    }
    
    
    #####################################################################
    #####################################################################
    ##Score strip charts
    
    number.rows=5
    number.columns=2
    ppg=number.rows*number.columns # panels per graph page
    
    par(mfrow = c(number.rows, number.columns))
    par(cex = 0.8) #was 0.5
    par(mar = c(4, 6, 0.5, 4), oma = c(4, 1, 4, 0.5))
    par(tcl = -0.25)
    par(mgp = c(3, 0.3, 0))
    for (i in 1:min(dim(my.pca$x)[2], score.strips.to.plot))
    {
      
      plot(my.pca$x[,i], as.integer(as.factor(row.factors)), 
           col=as.integer(as.factor(row.factors)), 
           cex=0.7+(2.5*as.integer(as.factor(row.factors2)))/nlevels(as.factor(row.factors2)),
           ann=FALSE, ylim=c(0,1+nlevels(as.factor(row.factors))), yaxt='n',
           xlim=c(-(max(abs(my.pca$x[,i])))*1.1,(max(abs(my.pca$x[,i])))*1.1) #center x-axis at 0
      )
      title(xlab=paste0("PC", i), line=2) 
      boxplot(my.pca$x[,i]~as.integer(as.factor(row.factors)), 
              names=paste0(levels(as.factor(row.factors))," "),
              ylim=c(-(max(abs(my.pca$x[,i])))*1.1,(max(abs(my.pca$x[,i])))*1.1) ,#center x-axis at 0
              horizontal=TRUE, las=1)
      #}
      #}
      #}
      
      mtext("Score Strip Plots.", outer=TRUE, side=1, cex=0.7)
    }
    
    par (mar=c(12,4,4,4)+0.5, xpd=TRUE, mfrow=c(1,1))
    
    
    ################## now biplots (two pc's, scores)
        for (page in 1:floor((dim(my.pca$x)[2]-number.columns)/number.rows))
    {
      number.rows=6
      number.columns=5
      ppg=number.rows*number.columns # panels per graph page
      
      par(mfrow = c(number.rows, number.columns))
      par(cex = 0.5)
      par(mar = c(3.5, 3.5, 0, 0), oma = c(4, 2, 4, 4))
      par(tcl = -0.25)
      par(mgp = c(3, 0.3, 0))
      for (row in 1:number.rows)
      {
        for (column in 1:number.columns)
        { 
          if (verbosity>8) print(paste0("Row ", row, " Column ", column))
          if (
            ((((page-1)*number.rows)+row)<=dim(my.pca$x)[2])&&
            ((((page-1)*number.rows)+row+column)<=dim(my.pca$x)[2])
          )
          {
            plot(my.pca$x[,(page-1)*number.rows+row], my.pca$x[,(page-1)*number.rows+row+column], 
                 col=as.integer(as.factor(row.factors)), 
                 cex=0.5+(as.integer(as.factor(row.factors2))/nlevels(as.factor(row.factors2))),
                 ann=FALSE,
                 xlim=c(-(max(abs(my.pca$x[,(page-1)*number.rows+row])))*1.1,(max(abs(my.pca$x[,(page-1)*number.rows+row])))*1.1), #center x-axis at 0)
                 ylim=c(-(max(abs(my.pca$x[,(page-1)*number.rows+row+column])))*1.1,(max(abs(my.pca$x[,(page-1)*number.rows+row+column])))*1.1)) #center x-axis at 0
            title(xlab=paste0("PC", (page-1)*number.rows+row), line=1.5) 
            title(ylab=paste0("PC", (page-1)*number.rows+row+column), line=1.5)
          }
          else plot.new()
        }
      }
      
      
      #####################################################################
      #Make a blank full page overaly plot so we can add legend
      par(fig=c(0, 1, 0, 1), oma=c(1, 1, 1, 1), mar=c(0, 0, 0, 0), new=TRUE)
      plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
      legend("bottom", pch=1, col=1:nlevels(as.factor(row.factors)),
             legend=levels(as.factor(row.factors)), ncol=min(5,nlevels(as.factor(row.factors))))
      mtext("Scores/Loadings for Top PCs", outer=TRUE, line=2)
      #mtext("Black is Phase I, Red is Phase II, Green is special. Time series small to large symbol.", outer=TRUE, side=1, cex=0.7)
    }
    par(mfrow=c(1,1))
  }
  return(my.pca)
}


#######################################################################################################
#######################################################################################################
##                                                                                                   ##
##            Function title.page                                  ##
##                                                                                                   ##
##                                                                                                   ##
#######################################################################################################
#######################################################################################################
title.page<-function(title="", subtitle="", 
                     title.size=4, 
                     title.x=0, title.y=40, 
                     title.color="black",
                     title.width=25,
                     subtitle.width=50,
                     subtitle.size=2, 
                     subtitle.x=0, subtitle.y=-40,
                     subtitle.color="black",
                     alignment=c(0.5,0.5))
{
 par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0))
 plot(0, 0, xlim=c(-100,100), ylim=c(-100,100),type='n', bty='n', xaxt='n', yaxt='n')## set up a fake graph that is invisible with 0,0 at center of page
 title2=paste(strwrap(title, width=title.width), collapse="\n")
 subtitle2=paste(strwrap(subtitle, width=subtitle.width), collapse="\n")
 
 text(title.x,title.y+1,title2, cex=title.size, adj=alignment, col=title.color)
   
 text(subtitle.x,subtitle.y,subtitle2, cex=subtitle.size, adj=alignment, col=subtitle.color)
 
}


#######################################################################################################
#######################################################################################################
##                                                                                                   ##
##            Function text.page : Write text on graphics device page                                ##
##                                                                                                   ##
##                                                                                                   ##
#######################################################################################################
#######################################################################################################
text.page<-function(text="", title="", 
                     title.size=2, 
                     title.x=0, title.y=90, 
                     title.color="black", 
                     text.size=1, 
                     text.x=-90, text.y=75,
                     text.color="black",
                     text.alignment=c(0,0.5),
                     title.alignment=c(0.5,0.5),
                     lines.per.page=40,
                     auto.lines.per.page=TRUE,
                     verbosity=10)
{
   # a good ref for fonts is http://cran.stat.auckland.ac.nz/doc/Rnews/Rnews_2004-2.pdf
  if (verbosity>9)
  {
    print("Entering text.page function with text:")
    print(text)
  }
  row.spacer=text.size*5
  if (auto.lines.per.page) lines.per.page=floor(160/row.spacer)
  
  for (page in 1:ceiling(length(text)/lines.per.page))
    
  {
    if (verbosity>8) print(paste0("Printing text page ", page, " at ", lines.per.page," lines per page."))
    par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0))
    plot(0, 0, xlim=c(-100,100), ylim=c(-100,100),type='n', bty='n', xaxt='n', yaxt='n')## set up a fake graph that is invisible with 0,0 at center of page
    if (page>1) title2<-paste0(title," (Cont.)") else title2<-title
    text(title.x,title.y,title2, cex=title.size, adj=title.alignment, col=title.color)#, fontface="bold")
    
   for (row in 1:lines.per.page)
   {
    if (((page-1)*lines.per.page+row)<=length(text))
      {
       text(text.x,text.y-row.spacer*(row-1),
        text[row+(page-1)*lines.per.page], 
        cex=text.size, adj=text.alignment, 
        col=text.color, family="mono")
      }
   }
  }
}
#######################################################################################################
#######################################################################################################
##                                                                                                   ##
##            Function table.to.text.page : Write table on graphics device page                                ##
##                                                                                                   ##
##                                                                                                   ##
#######################################################################################################
#######################################################################################################

table.to.text.page<-function(table, title="", 
                             title.size=2, 
                             title.x=0, title.y=90, 
                             title.color="black", 
                             text.size=1, 
                             text.x=-90, text.y=75,
                             text.color="black",
                             text.alignment=c(0,0.5),
                             title.alignment=c(0.5,0.5),
                             lines.per.page=60,
                             auto.lines.per.page=TRUE,
                             max.print.temp=12000, #about 200 pages of printing
                             width.temp=110, #columns of letters wide each page
                             verbosity=0)
{
  old.max<-getOption("max.print")
  old.width<-getOption("width")
  options(max.print = max.print.temp)
  options(width=width.temp)
 text.page(capture.output(table),
           title=title,
           title.x=title.x, title.y=title.y,
           title.color=title.color,
           text.size=text.size,
           text.x=text.x, text.y=text.y,
           text.color=text.color,
           text.alignment=text.alignment,
           title.alignment=title.alignment,
           lines.per.page=lines.per.page,
           auto.lines.per.page=auto.lines.per.page,
           verbosity=verbosity) 
 options(max.print = old.max)
 options(width = old.width)
}
#######################################################################################################
#######################################################################################################
##                                                                                                   ##
##            Function dch.create.pdf : Open a PDF file with timestamp                                ##
##                                                                                                   ##
##                                                                                                   ##
#######################################################################################################
#######################################################################################################

dch.create.pdf<-function(filename="Results from", timestamp=TRUE)
{
  if (timestamp) filename<-paste0(filename, " ",format(Sys.time(), "%Y%m%d-%H%M%S"),".PDF")
     else filename<-paste0(filename,".PDF")
  pdf(filename, width=8.5, height=11)
  
}

#######################################################################################################
#######################################################################################################
##                                                                                                   ##
##            Function label.plot : Write table on graphics device page                                ##
##                                                                                                   ##
##                                                                                                   ##
#######################################################################################################
#######################################################################################################

dch.label.plot<-function(names="", x.pos=0, y.pos=0, textsize=0.5, 
                         tick.type="dotted.line",
                         tick.length=0.008, #length of arrow/tick in fraction of graph area
                         tick.color="black",
                         tick.margin=0.002,
                         tick.lty=1,
                         label.pos=0.15
                         )
{
  list.length<-min(length(names), length(x.pos))
  xl=par("usr")[1:2]
  yl=par("usr")[3:4]
  
  
  for (lab in 1:list.length)
  {
    if (length(y.pos)<list.length)
    {
      label.y<-yl[2]-(label.pos)*(yl[2]-yl[1])  #position to hang text labels below
      tick.y1<-yl[2]-(label.pos+tick.margin)*(yl[2]-yl[1])  #position to hang arrows below
      tick.y2<-yl[2]-(label.pos+tick.margin+tick.length)*(yl[2]-yl[1])
    } else
    {
      label.y<-y.pos[lab]+(tick.margin+tick.length)*(yl[2]-yl[1])  #position to hang text labels below
      tick.y1<-y.pos[lab]  #position to hang arrows below
      tick.y2<-y.pos[lab]-(tick.length)*(yl[2]-yl[1])
    }
    
    text(x.pos[lab], label.y, names[lab], adj=c(0,0.5), cex=textsize, srt=90)
    if (tick.type=="dotted.line")
    {
      segments(x.pos[lab], tick.y1, x.pos[lab], tick.y2, col=tick.color, lty=tick.lty)
    }
  }
  
  

}

#++++++++++++++++++++++++++++++++++++++++++++
#generate a plot of point shapes which R knows about.
# adapted from http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
# and http://www.cookbook-r.com/Graphs/Shapes_and_line_types/
#++++++++++++++++++++++++++++++++++++++++++++
dch.symbols<-function(){
  oldpar<-par(no.readonly=TRUE)
  par(font=2, mar=c(0.5,0,0,0))
  y=rev(c(rep(1,6),rep(2,5), rep(3,5), rep(4,5), rep(5,5)))
  x=c(rep(1:5,5),6)
  plot(x, y, pch = 0:25, cex=1.5, ylim=c(1,14), xlim=c(0,27), 
       axes=FALSE, xlab="", ylab="", bg="blue")
  text(x, y, labels=0:25, pos=3)
 
  
  
  # Draw the lines
  for (i in 0:6) {
    segments(4,i+7, 7,i+7, lty=i, lwd=2)
  }
  # Add labels
  text(0, 7, "0. 'blank'"   ,  adj=c(0,.5))
  text(0, 8, "1. 'solid'"   ,  adj=c(0,.5))
  text(0, 9, "2. 'dashed'"  ,  adj=c(0,.5))
  text(0, 10, "3. 'dotted'"  ,  adj=c(0,.5))
  text(0, 11, "4. 'dotdash'" ,  adj=c(0,.5))
  text(0, 12, "5. 'longdash'",  adj=c(0,.5))
  text(0, 13, "6. 'twodash'" ,  adj=c(0,.5))
  
  for (col.no in 1:24)
  {
    text(9, 14-col.no/2, paste0("Color #", col.no), col=col.no, cex=0.9)
  }
  
  rainbow.num=48
  for (col.no in 1:rainbow.num)
  {
    text(13, 14-col.no*(12/rainbow.num), paste0("Rainbow(", rainbow.num,") #", col.no), 
         col=rainbow(rainbow.num)[col.no], cex=0.5)
  }
  
  cm.num=48
  for (col.no in 1:rainbow.num)
  {
    text(16, 14-col.no*(12/cm.num), paste0("cm.colors(", cm.num,") #", col.no), 
         col=cm.colors(cm.num)[col.no], cex=0.5)
  }
  topo.num=48
  for (col.no in 1:rainbow.num)
  {
    text(19, 14-col.no*(12/topo.num), paste0("topo.colors(", topo.num,") #", col.no), 
         col=topo.colors(topo.num)[col.no], cex=0.5)
  }
  terrain.num=48
  for (col.no in 1:rainbow.num)
  {
    text(22, 14-col.no*(12/terrain.num), paste0("terrain.colors(", terrain.num,") #", col.no), 
         col=terrain.colors(terrain.num)[col.no], cex=0.5)
  }
  heat.num=48
  for (col.no in 1:rainbow.num)
  {
    text(25, 14-col.no*(12/heat.num), paste0("heat.colors(", heat.num,") #", col.no), 
         col=heat.colors(heat.num)[col.no], cex=0.5)
  }
  par(oldpar )
}


###################################################
###################################################
###################################################
######                                      #######
######    Add error bars                    #######
######                                      #######
###################################################
###################################################
###################################################

#based on http://stackoverflow.com/questions/13032777/scatter-plot-with-error-bars


add.error.bars<-function(x, y, error.size)
{
  arrows(x, y-error.size, x, y+error.size, length=0.05, angle=90, code=3)
}
add.error.bars.x<-function(x, y, error.size)
{
  arrows(x-error.size, y, x+error.size, y, length=0.05, angle=90, code=3)
}


######################################################################################################
######################################################################################################
######################################################################################################
######                                                                                         #######
######                             Qgraph Correlations                                         #######
######                Graphical network analysis of correlations of a matrix                   #######
######                                                                                         #######
######                                                                                         #######
######################################################################################################
######################################################################################################
######################################################################################################

library(qgraph)

#prints two side-by side graphs, one with text labels and one with column numbers, also prints column names
qgraph.corr<-function(datam, names=NULL, name.table=colnames(datam),name.table.cex=0.5,
                      layout="spring", #instead of default "circle",
                      print.name.graph=TRUE, print.number.graph=TRUE, 
                      print.name.table=TRUE, color="white")
{
  par(mfrow=c(1,1))
  if (is.null(names)&&(print.name.graph)) qgraph((cor(datam)), layout=layout, color=color)
  if (!is.null(names)&&(print.name.graph)) qgraph((cor(datam)), layout=layout, labels=names, color=color)
  if (print.number.graph) qgraph((cor(datam)), layout=layout, labels=c(1:dim(datam)[2]), color=color)
  if (print.name.table) table.to.text.page(name.table, text.size=0.5)
  
  par(mfrow=c(1,1))
}

qgraph.corr.pdf<-function(datam, filename="Correlation QGraph", date.stamp=TRUE,
                          layout="spring", #instead of default "circle",
                          names=NULL, name.table=colnames(datam), name.table.cex=0.5,
                          print.name.graph=TRUE, print.number.graph=TRUE, 
                          print.name.table=TRUE,color="white")
{ 
  if (date.stamp==TRUE) filename<-paste(filename, " ", format(Sys.time(), "%Y%m%d-%H%M%S"),".PDF")
   else filename<-paste(filename, ".PDF") 
  pdf(filename, width=8.5, height=11)  
  qgraph.corr(datam, names, name.table, name.table.cex, layout,
              print.name.graph, print.number.graph, 
              print.name.table, color=color)
  dev.off()
}

## Function to get users dropbox folder location
## from https://www.r-bloggers.com/get-a-path-to-your-dropbox-folder/
get.dropbox.folder <- function() {
  
  if (Sys.info()['sysname'] == 'Darwin') {
    info <- RJSONIO::fromJSON(
      file.path(path.expand("~"),'.dropbox','info.json'))
  }
  if (Sys.info()['sysname'] == 'Windows') {
    info <- RJSONIO::fromJSON(
      if (file.exists(file.path(Sys.getenv('APPDATA'), 'Dropbox','info.json'))) {
        file.path(Sys.getenv('APPDATA'), 'Dropbox', 'info.json')
      } else {
        file.path(Sys.getenv('LOCALAPPDATA'),'Dropbox','info.json')
      }
    )
  }
  
  dropbox_base <- info$personal$path
  return(dropbox_base)
}

######################################################################################################
######################################################################################################
######                                                                                         #######
######                                                                                                     #######
######                Open current working directory in windows explorer                   #######
######                                                                                         #######
######                                                                                         #######
######################################################################################################
######################################################################################################
######################################################################################################

dch.openwd<-function()
{
  
  shell(paste0("explorer ", getwd()), intern=TRUE, translate=TRUE)

}

######################################################################################################
######################################################################################################
######                                                                                         #######
######                                                                                                     #######
######                Check memory usage                                                 #######
######                                                                                         #######
######                                                                                         #######
######################################################################################################
######################################################################################################
######################################################################################################

dch.memory<-function() # based on https://heuristically.wordpress.com/2010/01/04/r-memory-usage-statistics-variable/
  
{
# print aggregate memory usage statistics
print(paste('R is using', memory.size(), 'MB out of limit', memory.limit(), 'MB'))

# create function to return matrix of memory consumption
object.sizes <- function()
{
  return(rev(sort(sapply(ls(envir=.GlobalEnv), function (object.name) 
    object.size(get(object.name))))))
}

# print to console in table format
print(object.sizes())

# draw bar plot
barplot(object.sizes()/1024/1024, 
        main="Memory usage by object", ylab="MegaBytes", xlab="Variable name", 
        col=heat.colors(length(object.sizes())),
        las=2)

# draw dot chart
dotchart(object.sizes()[1:min(length(object.sizes()),30)]/1024/1024,
          main="Memory usage by object", xlab="MegaBytes")

# draw pie chart
pie(object.sizes(), main="Memory usage by object")
}
