## RandomizePV4.R - 5/16/2018; SLGarman - BLM NOC, Denver.  Version 4.  Uses spsample and bb=bbox() to locate random samples in GenPts(). And, this
##                                                          version can use (buffered) NHD lines to evaluate aquatic AIM designs.

## Determines if the dispersion of points in a survey design is different from random.  Developed specifically to evaluate the spatial balance of merged (i.e., expanded) GRTS designs.
## For a specified Frame and n GRTS points, derives x number of random sets of n points, and compares
## the mean nearest neighbor (MNN) distance of the input points with that of each set of randomly selected points.

## Records the proportion of random-point sets with a MNN equal to or greater than the survey-design (GRTS) points.
## This proportion can be treated as a P value to test:
##	Ho: GRTS design is not different from random
##	Ha: GRTS design is different from random

##      E.g., a proportion of 0.03 means that the input GRTS points are different from random.
##            A proportion of 0.20 means that the input GRTS points are Likely Not different from random - that is,
##              randomly selected point sets were more dispersed than the GRTS points in 20% of the replications.

## Mean NN is based on Arithmetic mean and Geometric mean.  Arithmetic-mean results are most applicable when there is 1 highly dominant polygon per feature.
## That is, the frame or a strata is represented almost entirely by 1 polygon.  Where a frame or stratum is represented by multiple disjunct
## (and highly dispersed) polygons, the Geometric mean should be considered (better represents the central tendency of highly skewed data).  If random proportions of
##  both measures are > 0.1, then best to retain Ho:, else good bet to go with Ha: especially if the Geometric mean is <= 0.05.


## This version only works on Polygonal sample frames and strata.  Additionally, this version includes the option to perform randomization tests on a strata basis.


## I.  AIM Terrestrial analyses.  For both Frame and Strata, requires polygons to be dissolve - only 1 feature per attribute.  The strata of each point in the point file
## is derived by overlaying the points on the strata file.  Random() [see bottom of the script] is called to process AIM terrestrial designs.

## II.  AIM Aquatic analyses.  This was an add-on.  The strata and frame in this case are the same, and consists of buffered NHD lines (buffer by say 0.5 m to create a poly-line
## file - easier to process using this algorithm). Unlike AIM Terrestrial frames and strata, NHD poly-lines SHOULD NOT be dissolved; multiple polygons per feature (i.e., stream order segment)
## are EXPECTED.  The user specifies the name of the strata field in the NHD 'strata' file.  It is assumed that the Aquatic points file contains a field called STRM_ORDR that defines
## the strata (stream order), and STRM_ORDR values match the strata values in the NHD poly-line file.  STRM_ORDER is directly used to determine the strata of points - the script
## doesn't overlay points onto the NHD file to determine strata like it does with Terrestrial AIM points (aquatic point locations are field-derived coordinates and do not
## always overlap even the buffered NHD poly-line).  RandomAquatic() is called to process AIM Aquatic designs.



# How to run this script.
#1)  Copy all of section I into the R console

#2)  Modify the I/O arguments in section II, then cut and paste arguments and the call to Random() to the R console.
## Specify the working directory where the following shapefiles reside
## Specify the polygonal frame (shapefile)
## Specify the GRTS points (the input points shapefile)
## Specify the number of random replications (500 should do the trick, can try more depending on computer speed)
## Specify the output file name (summary of results - short and sweet)
## Specify the random number seed or leave as is
## Specify a strata file for an analysis based on strata, else set this to NA.  For terrestrial AIM analyses, strata file must be dissolved.
## Specify the strata field name in the strata file, else set this to NA
## Regardless if a strata file/analysis is performed, can specify an entire frame analysis that ignores strata

#3)  OR use Batch mode to run the script after you modify I/O arguments.  Your system path must include the bin subdirectory of the R version you want to use.
#     Do the following in a command-line console:   R CMD BATCH randomizepv4.r capture_console_output

#     capture_console_output is an ASCII file that records the output you'd see in the R console when running from the console.
#     You can batch multiple jobs at once (rename your modified versions of randomizepv4.r & use different file names to capture the results),
#     grab lunch, then come back and check out all of the results.......


#I.  Copy all of Section I. into the R console.
#####################################################################
library(raster)
library(tidyverse)
library(stringr)
library(rgdal)
library(rgeos)
library(maptools)
library(digest)
library(sp)

################################################################
#  Extract area of polygons & build a cumulative Prob. Distribution indexed by polygon number (spdf is the strata or frame)
## This only works if polygons were dissolved - designed for polygonal frames (e.g., terrestrial AIM)
ExtractPolyArea <- function(spdf) {
  if (class(spdf) != "SpatialPolygonsDataFrame") {
    stop("spdf must be a spatial polygons data frame")
  }

  # Get the areas of the polygons in the SPDF
  areas <- sapply(X = spdf@polygons[[1]]@Polygons,
                  FUN = function(X) {
                    X@area
                  })

  # Make a data frame with the areas and the within-polygon ID
  areas_df <- data.frame(area = areas,
                         id = 1:length(areas))

  # Sort from largest to smallest area
  # (This can help speed up selecting from the probability distribution)
  areas_df <- areas_df[order(-areas_df[["area"]]), ]

  # Get the total area
  total_area <- sum(areas_df[["area"]])

  # Add proportional area to the data frame
  areas_df[["area_prop"]] <- areas_df[["area"]] / total_area

  # Add cumulative frequency distribution, which can be treated as a probability distribution
  areas_df[["cum_freq"]] <- cumsum(areas_df[["area_prop"]])

  # Return this data frame!
  return(areas_df)
}

## End of ExtractPolyArea()
#########################################################################################################
#  Extract area of polygons & build a cumulative Prob. Distribution indexed by polygon number (spdf is the strata or frame)
## This only works if polygons were NOT dissolved - designed for aquatic lines buffered to form a polygonal frame
ExtractPolyAreaAquatic <- function(spdf) {
  if (class(spdf) != "SpatialPolygonsDataFrame") {
    stop("spdf must be a spatial polygons data frame")
  }

  # Get the areas of the polygons
  areas <- sapply(X = spdf@polygons,
                  FUN = function(X) {
                    sapply(X = X@Polygons,
                           FUN = function(X) {
                             X@area
                           })
                  })

  # Make a data frame with the areas and the within-polygon ID
  areas_df <- data.frame(area = areas,
                         id = 1:length(areas))

  # Sort from largest to smallest area
  # (This can help speed up selecting from the probability distribution)
  areas_df <- areas_df[order(-areas_df[["area"]]), ]

  # Get the total area
  total_area <- sum(areas_df[["area"]])

  # Add proportional area to the data frame
  areas_df[["area_prop"]] <- areas_df[["area"]] / total_area

  # Add cumulative frequency distribution, which can be treated as a probability distribution
  areas_df[["cum_freq"]] <- cumsum(areas_df[["area_prop"]])

  # Return this data frame!
  return(areas_df)
}

###############################################################################
##  Select polygon from prob distribution.  randn is a urv.
SelectFrDistr<-function(ProbDistr,randn)
{
  for(i in 1:nrow(ProbDistr)) {
    if(randn<ProbDistr$VAR1[i]) {
      return(ProbDistr$VAR2[i])
      break
    }
  }
  return(1)		## A default
}
###########################################################
## Geometric mean calcs
gm_mean = function(x, na.rm=TRUE){
  ret<- exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  return(ret)
}
###########################################################
NN<-function(apts)	# Derive mean Nearest Neighbor using a pts file containing XMETERS and YMETERS
{
  a<--1*apts$XMETERS		## Assumes X is always negative
  b<-apts$YMETERS
  a<-outer(a,a,'-')
  a<-a*a
  b<-outer(b,b,'-')
  b<-b*b
  a<-a+b
  a<-sqrt(a)			## Euclidean distance -> Sqrt of a**2 + b**2
  a[a==0] <-999999999		## This just gets rid of the difference between the same pt, enabling min() function to find the true NN of the point

  Geo<-NULL
  RMuNN<-0
  for(i in 1:nrow(apts)) {
    nn<-min(a[,i])
    RMuNN<-RMuNN+nn
    Geo<-rbind(Geo,nn)
  }
  RMuNN<-(RMuNN/nrow(apts))	## Derive mean NN - arithmetic
  GMuNN<-gm_mean(Geo,T)		## Derive mean NN - geometric
  ret<-rbind(RMuNN,GMuNN)
  return(ret)
}
###########################################################################GenPts(number,thepts,tempS,StrataNN,StrataBox)
## Generate sets of random points, derive Mean NN of each set, compare with the specified NN (MeanNN) of the GRTS (input) pts, record the proportion
##          of random replicates where the random NN >= NN of the input pts.
GenPts<-function(number,		## Number of reps
                 thepts.spdf,		## The pts file
                 aoi.spdf,		## The AOI (frame or a stratum file)
                 MeanNN,		## Vector of NN distances, output from NN()
                 type)			## type ==2 for aquatic analysis where poly-lines are not dissolved, else 1

{
  projection=CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")
  projectionAL <- CRS("+proj=aea")
  target<-1


  cnt<-0
  samplesize<-number		## The requested number of random simulations
  proportion1<-0		## The proportion of random sets with a mean NN => than MNN of the input pts
  proportion2<-0
  if(type==1)ProbDistr<-ExtractPolyArea(aoi.spdf)	## Generate the cumulative Prob Distribution
  if(type==2)ProbDistr<-ExtractPolyAreaAquatic(aoi.spdf)

  while(number>0) {
    store.spdf<-NULL
    counter<-nrow(thepts.spdf)			## Number of random points we need to generate
    while(counter>0) {
      temp<-NULL
      for(i in 1:(nrow(thepts.spdf)*1.25) ) {  	## Pick up counter*1.25 draws then check them out.  With bbox() below, may select points outside of polygon area.
        ## The 25% adjustment helps to account for non-overlapping points.  By adding perhaps more than we need, we at least cut
        ## down on repeat conversion from Spatial Points to SpatialPointsDataFrame and the use of the over() function [seems to
        ## be a bit of a bottle neck in terms of time!]. If we exceed what we need, we get rid of the extra points below....
        set.seed(1)
        urv<-runif(1)				##  uniform random variate (urv) for selecting a polygon
        set.seed(1)
        opt<-SelectFrDistr(ProbDistr,urv)	##  Using the urv, determine the polygon number (opt) from the cumulative freq distribution
        if(type==1)poly<-aoi.spdf@polygons[[1]]@Polygons[[opt]]	## If dissolved (as it should be), then always access polygons[[1]].  Polygons[[x]] are the multiple polygons.
        if(type==2)poly<-aoi.spdf@polygons[[opt]]@Polygons[[1]] ## If not dissolved
        set.seed(1)
        z<-spsample(poly,n=1,"random",bb=bbox(poly))	## Use the bounding box in spsample to select just 1 random point.
        if(is.null(temp)){temp<-z}else{temp<-rbind(temp,z)}	## Probably doing something wrong here, but can't seem to bind temp if it isn't already set???????
      }


      ##  Translate from SpatialPoints to SpatialPointsDataFrame, overlay onto the aoi (a stratum or the Frame), and only retain overlapping points
      proj4string(temp)=CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")
      temp.spdf <- SpatialPointsDataFrame(temp, data.frame(id=1:length(temp)))
      goodpts<-over(temp.spdf,aoi.spdf)		## goodpts$CODE will be =1 if a point overlaps an actual polygon in aoi.spdf (a stratum or the frame)
      temp.spdf$CODE<-goodpts$CODE
      temp.spdf<-temp.spdf[temp.spdf$CODE %in% target ,]	## Only want overlapping points

      ## Bind the points that overlap aoi.spdf, and adjust counter.
      if(nrow(temp.spdf)>0) {
        if(is.null(store.spdf)) {store.spdf<-temp.spdf}else {store.spdf<-rbind(store.spdf,temp.spdf)}
        counter<-counter-nrow(temp.spdf)
        #print(counter)
      }
    }
    ## For reference only...
    #coordinates(store.df)=~X+Y
    #proj4string(store.df)=CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")
    #rand.spdf <- SpatialPointsDataFrame(store.df, data.frame(id=1:length(store.df)))
    ####

    rand.spdf<-store.spdf
    rand.spdf<-rand.spdf[1:nrow(thepts.spdf),]	## In case we have more random points than needed, just pick the first nrow(thepts.spdf) pts.  We
    ## still have a fully random sample since we effectively store the random points by accession, and
    ## we elimiinate points from the bottom up.
    number<-number-1				## Decrement rep counter
    temp.spdf<- spTransform(rand.spdf,projection)
    rand.spdf@data <- cbind(rand.spdf@data,temp.spdf@coords)
    names(rand.spdf)[names(rand.spdf) == "X"] <- "LONG"		## Noticed this can be X or x; so lets account for both possibilities.
    names(rand.spdf)[names(rand.spdf) == "Y"] <- "LAT"
    names(rand.spdf)[names(rand.spdf) == "x"] <- "LONG"
    names(rand.spdf)[names(rand.spdf) == "y"] <- "LAT"
    temp.spdf<- spTransform(rand.spdf,projectionAL)
    rand.spdf@data <- cbind(rand.spdf@data,temp.spdf@coords)
    names(rand.spdf)[names(rand.spdf) == "X"] <- "XMETERS"
    names(rand.spdf)[names(rand.spdf) == "Y"] <- "YMETERS"
    names(rand.spdf)[names(rand.spdf) == "x"] <- "XMETERS"
    names(rand.spdf)[names(rand.spdf) == "y"] <- "YMETERS"

    RMuNN<-NN(rand.spdf)					## Returns Arithmetic and Geometric NN
    if(RMuNN[1]>= MeanNN[1])proportion1<-proportion1+1	## Arithmetic NN
    if(RMuNN[2]>= MeanNN[2])proportion2<-proportion2+1	## Geometric NN

    ## Just for show: Let's you see the evolution of the P value as number of reps increase.  Output to console every 100 reps.
    #print(number)
    cnt<-samplesize-number
    a<-cnt/100
    b<-round(a,digits=0)
    if(a*100 == b*100 & cnt!=0){
      print(paste("Working on rep # = ",cnt," Evolving Arithmetic P value= ",proportion1/cnt,sep=" "))
      print(paste("Working on rep # = ",cnt," Evolving Geometric P value= ",proportion2/cnt,sep=" "))
    }

    ## In case you wanna look at the random pts - quickie diagnostic for special occasions.
    #writeOGR(rand.spdf, ".",output,driver="ESRI Shapefile",overwrite_layer=T)
    #q()


  }  ## While number >0

  proportion1<-proportion1/samplesize		## Convert count to proportion
  proportion2<-proportion2/samplesize
  ret<-c(proportion1,proportion2)		## s.l. P value based on Arithmetic mean, based on Geometric mean
  return(ret)
}
############################################################################################
Ingest<-function(layername)	## layer name of shapefile to ingest
{
  projection = CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")

  Aspdf <- readOGR(dsn=getwd(),layer=layername,stringsAsFactors=FALSE)  ## The sample frame
  Aspdf <- spTransform(Aspdf, projection)
  names(Aspdf@data) <- str_to_upper(names(Aspdf@data))
  return(Aspdf)
}
###############################################################################
GetXY<-function(thepts)
{
  projectionAL <- CRS("+proj=aea")

  # access the x, y coordinates of the pts and derive nearest neighbor
  thepts@data <- cbind(thepts@data, thepts@coords)
  names(thepts)[names(thepts) == "coords.x1"] <- "LONGITUDE"
  names(thepts)[names(thepts) == "coords.x2"] <- "LATITUDE"

  temp.spdf<- spTransform(thepts,projectionAL)
  thepts@data <- cbind(thepts@data, temp.spdf@coords)
  names(thepts)[names(thepts) == "coords.x1"] <- "XMETERS"
  names(thepts)[names(thepts) == "coords.x2"] <- "YMETERS"

  thepts<-thepts[, c("XMETERS","YMETERS")]
  return(thepts)
}
######################################################################
Random<-function(WD,		##working directory
                 frame,		## the sample frame as a spdf (shapefile)
                 pts,		## the GRTS points as a spdf  (shapefile)
                 reps,		## Number of sets of random points
                 output,	## Name of file to store results
                 strata,	## Name of strata file (analyses will be conducted at the strata level) or NA
                 stratafield,	## Strata-field name in strata or NA
                 doFrame) 	## Set to T for analysis across the entire FRAME (ignores strata if specified), else F

{
  projection = CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")
  projectionAL <- CRS("+proj=aea")

  setwd(WD)
  numberofreps<-reps

  blank<-"Randomization Results"
  write.table(blank,output,row.names=F,col.names=F,quote=F,append=F)


  ## Ingest all of the points corresponding to the entire Frame.  Even if doFrame==F, we'll use the full suite of pts for the strata-based analysis
  Framepts.spdf<-Ingest(pts)
  Framepts.spdf<-GetXY(Framepts.spdf)	## Access X & Y

  # quickie for experimenting with N
  # Framepts.spdf$DELETE<-0
  #for(i in 201:nrow(Framepts.spdf)) {
  #   Framepts.spdf$DELETE[i]<-1
  #}
  #Framepts.spdf<-Framepts.spdf[Framepts.spdf$DELETE==0 ,]



  if(doFrame==T) {		## If requested, first analyze entire frame
    ## ingest the frame
    frame.spdf<-Ingest(frame)
    frame.spdf$CODE<-1
    frame.spdf<-frame.spdf[, c("CODE")]			## Render Frame to just CODE

    # Derive Mean NN of the Frame points
    FrameNN<-NN(Framepts.spdf)		## Mean NN; [1]=arithmetic, [2]=geometric
    print(paste("Arithmetic Mean nearest neighbor of the input points= ",FrameNN[1]," for N= ",nrow(Framepts.spdf)," points",sep=" "))
    print(paste("Geometric Mean nearest neighbor of the input points= ",FrameNN[2]," for N= ",nrow(Framepts.spdf)," points",sep=" "))

    ## Do randomization test for the entire frame
    number<-reps
    thepts<-Framepts.spdf
    proportion<-GenPts(number,thepts,frame.spdf,FrameNN,1)	## Number of requested sim reps, the pts file, the aoi file (frame or a stratum file),
    ##   the previously calculated vector (arithmetic and geometric) of MNN distances of the input pts

    print(paste("Proportion of random draws with Arithmetic NN >=  NN of the input points = ",proportion[1],sep=" "))
    a<-data.frame(COL="Arithmetic Mean NN(meters) of the ",COL2="",COL3="",N=nrow(Framepts.spdf),COL4="input points = ",VAL=FrameNN[1])
    b<-data.frame(COL="Proportion of",COL2=as.character(numberofreps),COL3=" random draws with an Arithmetic mean NN >= mean NN of the",N=nrow(Framepts.spdf),COL4=" input points = ",VAL=proportion[1])
    ans<-rbind(a,b)

    print(paste("Proportion of random draws with Geometric NN >=  NN of the input points = ",proportion[2],sep=" "))
    a<-data.frame(COL="Geometric Mean NN(meters) of the ",COL2="",COL3="",N=nrow(Framepts.spdf),COL4="input points = ",VAL=FrameNN[2])
    b<-data.frame(COL="Proportion of",COL2=as.character(numberofreps),COL3=" random draws with a Geometric mean NN >= mean NN of the",N=nrow(Framepts.spdf),COL4=" input points = ",VAL=proportion[2])
    ans<-rbind(ans,a)
    ans<-rbind(ans,b)
    write.table(ans,output,row.names=F,col.names=F,quote=F,append=T)

  }else {
    frame.spdf<-NULL
  }



  ######### Analyze by strata if requested
  strata.spdf<-NULL
  # Ingest the strata file if specified, and analyze on a stratum basis
  if(!is.na(strata)) {
    strata.spdf<-Ingest(strata)
  }else {
    q()
  }

  if(!is.null(strata.spdf)) {
    names(strata.spdf)[names(strata.spdf)==stratafield]<-"STRATA"
    strata.spdf<-strata.spdf[, c("STRATA")]
    slist<-unique(strata.spdf$STRATA)

    for(i in 1:length(slist)) {
      tempS<-strata.spdf[strata.spdf$STRATA %in% slist[i] ,]
      tempS$CODE<-1
      tempT<-over(Framepts.spdf,tempS)
      Framepts.spdf$USE<-tempT$STRATA
      tempPTS<-Framepts.spdf[!is.na(Framepts.spdf$USE) ,]		## Now we have just the points within strata i (from slist)
      Framepts.spdf$USE<-NULL			## Clear for re-use
      if(nrow(tempPTS)>1) {			## If we have >1 points, do a randomization
        # Derive Mean NN of the Stratum points
        StrataNN<-NN(tempPTS)		## Mean NN; [1]=arithmetic, [2]=geometric
        print(paste("For Stratum = ",slist[i],":",sep=""))
        print(paste("Arithmetic Mean nearest neighbor of the input points= ",StrataNN[1]," for N= ",nrow(tempPTS)," points",sep=" "))
        print(paste("Geometric Mean nearest neighbor of the input points= ",StrataNN[2]," for N= ",nrow(tempPTS)," points",sep=" "))

        ## Do randomization test for this stratum
        number<-reps
        thepts<-tempPTS
        proportion<-GenPts(number,thepts,tempS,StrataNN,1)	## Number of requested sim reps, the pts file, the aoi file (frame or a stratum file),
        ##   the previously calculated vector (arithmetic and geometric) of NN distances of the input pts
        print(paste("Proportion of random draws with Arithmetic NN >=  NN of the input points = ",proportion[1],sep=" "))
        a<-data.frame(COL="Arithmetic Mean NN(meters) of the ",COL2="",COL3="",N=nrow(tempPTS),COL4="input points = ",VAL=StrataNN[1])
        b<-data.frame(COL="Proportion of",COL2=as.character(numberofreps),COL3=" random draws with an Arithmetic mean NN >= mean NN of the",N=nrow(tempPTS),COL4=" input points = ",VAL=proportion[1])
        ans<-rbind(a,b)

        print(paste("Proportion of random draws with Geometric NN >=  NN of the input points = ",proportion[2],sep=" "))
        a<-data.frame(COL="Geometric Mean NN(meters) of the ",COL2="",COL3="",N=nrow(tempPTS),COL4="input points = ",VAL=StrataNN[2])
        b<-data.frame(COL="Proportion of",COL2=as.character(numberofreps),COL3=" random draws with a Geometric mean NN >= mean NN of the",N=nrow(tempPTS),COL4=" input points = ",VAL=proportion[2])
        ans<-rbind(ans,a)
        ans<-rbind(ans,b)
        blank<-c(" ")
        write.table(blank,output,row.names=F,col.names=F,quote=F,append=T)
        blank<-paste("For Stratum= ",slist[i],":",sep="")
        write.table(blank,output,row.names=F,col.names=F,quote=F,append=T)
        write.table(ans,output,row.names=F,col.names=F,quote=F,append=T)
      }else {
        print(paste("For Stratum = ",slist[i],":",sep=""))
        print(paste("WARNING, number of points only = ",nrow(tempPTS),sep=""))
        blank<-" "
        write.table(blank,output,row.names=F,col.names=F,quote=F,append=T)
        blank<-paste("For Stratum = ",slist[i],":",sep="")
        write.table(blank,output,row.names=F,col.names=F,quote=F,append=T)
        blank<-paste("WARNING, number of points only = ",nrow(tempPTS),sep="")
        write.table(blank,output,row.names=F,col.names=F,quote=F,append=T)
      }
    }
  }
} ## End of Random()



######################################################################
##  This more-or-less dups Random() but is customized to work on buffered NHD polygons for aquatic analyses.
RandomAquatic<-function(WD,		##working directory
                        frame,		## the sample frame as a spdf (shapefile)
                        pts,		## the GRTS points as a spdf  (shapefile)
                        reps,		## Number of sets of random points
                        output,	## Name of file to store results
                        strata,	## Name of strata file (analyses will be conducted at the strata level) or NA
                        stratafield,	## Strata-field name in strata or NA
                        doFrame) 	## Set to T for analysis across the entire FRAME (ignores strata if specified), else F

{
  projection = CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")
  projectionAL <- CRS("+proj=aea")

  setwd(WD)
  numberofreps<-reps

  blank<-"Randomization Results"
  write.table(blank,output,row.names=F,col.names=F,quote=F,append=F)


  ## Ingest all of the points corresponding to the entire Frame.  Even if doFrame==F, we'll use the full suite of pts for the strata-based analysis
  Framepts.spdf<-Ingest(pts)
  a<-GetXY(Framepts.spdf)	## Access X & Y
  a$STRATA<-Framepts.spdf$STRM_ORDR			## Assumes pts have field STRM_ORDR
  Framepts.spdf<-a

  # quickie for experimenting with N
  # Framepts.spdf$DELETE<-0
  #for(i in 201:nrow(Framepts.spdf)) {
  #   Framepts.spdf$DELETE[i]<-1
  #}
  #Framepts.spdf<-Framepts.spdf[Framepts.spdf$DELETE==0 ,]



  if(doFrame==T) {		## If requested, first analyze entire frame
    ## ingest the frame
    frame.spdf<-Ingest(frame)
    frame.spdf$CODE<-1
    frame.spdf<-frame.spdf[, c("CODE")]			## Render Frame to just CODE

    # Derive Mean NN of the Frame points
    FrameNN<-NN(Framepts.spdf)		## Mean NN; [1]=arithmetic, [2]=geometric
    print(paste("Arithmetic Mean nearest neighbor of the input points= ",FrameNN[1]," for N= ",nrow(Framepts.spdf)," points",sep=" "))
    print(paste("Geometric Mean nearest neighbor of the input points= ",FrameNN[2]," for N= ",nrow(Framepts.spdf)," points",sep=" "))

    ## Do randomization test for the entire frame
    number<-reps
    thepts<-Framepts.spdf
    proportion<-GenPts(number,thepts,frame.spdf,FrameNN,2)	## Number of requested sim reps, the pts file, the aoi file (frame or a stratum file),
    ##   the previously calculated vector (arithmetic and geometric) of MNN distances of the input pts

    print(paste("Proportion of random draws with Arithmetic NN >=  NN of the input points = ",proportion[1],sep=" "))
    a<-data.frame(COL="Arithmetic Mean NN(meters) of the ",COL2="",COL3="",N=nrow(Framepts.spdf),COL4="input points = ",VAL=FrameNN[1])
    b<-data.frame(COL="Proportion of",COL2=as.character(numberofreps),COL3=" random draws with an Arithmetic mean NN >= mean NN of the",N=nrow(Framepts.spdf),COL4=" input points = ",VAL=proportion[1])
    ans<-rbind(a,b)

    print(paste("Proportion of random draws with Geometric NN >=  NN of the input points = ",proportion[2],sep=" "))
    a<-data.frame(COL="Geometric Mean NN(meters) of the ",COL2="",COL3="",N=nrow(Framepts.spdf),COL4="input points = ",VAL=FrameNN[2])
    b<-data.frame(COL="Proportion of",COL2=as.character(numberofreps),COL3=" random draws with a Geometric mean NN >= mean NN of the",N=nrow(Framepts.spdf),COL4=" input points = ",VAL=proportion[2])
    ans<-rbind(ans,a)
    ans<-rbind(ans,b)
    write.table(ans,output,row.names=F,col.names=F,quote=F,append=T)

  }else {
    frame.spdf<-NULL
  }



  ######### Analyze by strata if requested
  strata.spdf<-NULL
  # Ingest the strata file if specified, and analyze on a stratum basis
  if(!is.na(strata)) {
    strata.spdf<-Ingest(strata)
  }else {
    q()
  }

  if(!is.null(strata.spdf)) {
    names(strata.spdf)[names(strata.spdf)==stratafield]<-"STRATA"
    strata.spdf<-strata.spdf[, c("STRATA")]
    slist<-unique(strata.spdf$STRATA)

    for(i in 1:length(slist)) {
      tempS<-strata.spdf[strata.spdf$STRATA %in% slist[i] ,]
      tempS$CODE<-1
      tempPTS<-Framepts.spdf[Framepts.spdf$STRATA==slist[i] ,]
      if(nrow(tempPTS)>1) {			## If we have >1 points, do a randomization
        # Derive Mean NN of the Stratum points
        StrataNN<-NN(tempPTS)		## Mean NN; [1]=arithmetic, [2]=geometric
        print(paste("For Stratum = ",slist[i],":",sep=""))
        print(paste("Arithmetic Mean nearest neighbor of the input points= ",StrataNN[1]," for N= ",nrow(tempPTS)," points",sep=" "))
        print(paste("Geometric Mean nearest neighbor of the input points= ",StrataNN[2]," for N= ",nrow(tempPTS)," points",sep=" "))

        ## Do randomization test for this stratum
        number<-reps
        thepts<-tempPTS
        proportion<-GenPts(number,thepts,tempS,StrataNN,2)	## Number of requested sim reps, the pts file, the aoi file (frame or a stratum file),
        ##   the previously calculated vector (arithmetic and geometric) of NN distances of the input pts
        print(paste("Proportion of random draws with Arithmetic NN >=  NN of the input points = ",proportion[1],sep=" "))
        a<-data.frame(COL="Arithmetic Mean NN(meters) of the ",COL2="",COL3="",N=nrow(tempPTS),COL4="input points = ",VAL=StrataNN[1])
        b<-data.frame(COL="Proportion of",COL2=as.character(numberofreps),COL3=" random draws with an Arithmetic mean NN >= mean NN of the",N=nrow(tempPTS),COL4=" input points = ",VAL=proportion[1])
        ans<-rbind(a,b)

        print(paste("Proportion of random draws with Geometric NN >=  NN of the input points = ",proportion[2],sep=" "))
        a<-data.frame(COL="Geometric Mean NN(meters) of the ",COL2="",COL3="",N=nrow(tempPTS),COL4="input points = ",VAL=StrataNN[2])
        b<-data.frame(COL="Proportion of",COL2=as.character(numberofreps),COL3=" random draws with a Geometric mean NN >= mean NN of the",N=nrow(tempPTS),COL4=" input points = ",VAL=proportion[2])
        ans<-rbind(ans,a)
        ans<-rbind(ans,b)
        blank<-c(" ")
        write.table(blank,output,row.names=F,col.names=F,quote=F,append=T)
        blank<-paste("For Stratum= ",slist[i],":",sep="")
        write.table(blank,output,row.names=F,col.names=F,quote=F,append=T)
        write.table(ans,output,row.names=F,col.names=F,quote=F,append=T)
      }else {
        print(paste("For Stratum = ",slist[i],":",sep=""))
        print(paste("WARNING, number of points only = ",nrow(tempPTS),sep=""))
        blank<-" "
        write.table(blank,output,row.names=F,col.names=F,quote=F,append=T)
        blank<-paste("For Stratum = ",slist[i],":",sep="")
        write.table(blank,output,row.names=F,col.names=F,quote=F,append=T)
        blank<-paste("WARNING, number of points only = ",nrow(tempPTS),sep="")
        write.table(blank,output,row.names=F,col.names=F,quote=F,append=T)
      }
    }
  }
} ## End of RandomAquatic()
##  Paste all of the above into the R console.

##################################################################################################################
## II. Modify the following arguments, then call Random for Terrestrial AIM analysis, else RandomAquatic for Aquatic AIM analysis.
##     For Terrestrial AIM analyses, all shapefiles must be dissolved (1 row per feature in the shapefile).
##     For Aquatic AIM analyses, the strata/frame file (1 in the same) should not be dissolved (pertains to the standard way strata are generated in the aquatic world).


WD<-c("c:/projects/mastergrts")	# Working directory - all input files must reside in this directory
frame<-c("Frame")			# Name of the frame file in the working directory - this will be the same as strata for AIM Aquatic analyses
pts<-c("grtssample")			# Name of the GRTS pts file in the working directory.  A field called USE is created in the terrestrial overlay function, so
##    the incoming pts file shouldn't have a field named USE.
reps<-500				# No. of randomly selected sets of points
output<-c("test")		        # Name of file to record the results
set.seed(1)				# Set the random number seed or not!

############## Strata info
strata<-NA
#strata<-c("eagleSTRDISS")		# Set to NA if you don't want to analyze by strata (e.g., strata<-NA) , else enter the name of the strata shapefile.
#stratafield<-c("GRIDCODE")		# If strata is set, then specify the strata field in the strata shapefile (upper case), else stratafield is ignored.
# stratafield is the strata attribute within file=strata.
## NOTE:  The strata field is renamed STRATA in prcessing, so the strata file can't already have this field name.
## For Terrestrial AIM processing, the over() function is used to derive the 'actual' strata of points; thus,
## any strata designation within the pts file is not used!  THUS, the strata file used to generate in input GRTS pts should be same as used here.

doFrame<-T				# Set to T if you want to analyze across the enter Frame without regard to strata, else set to F.  Even if you set strata, setting this to T
#     will generate results where strata are ignored; in essence you do 2 analyses at once - frame-based and strata-based.  If strata  == NA,
#     make sure this is set to T, otherwise no analyses are performed!

## Nothing is returned from Random() - the ASCII-formatted results are output during processing to the file specified by output (set above)
Random(WD,frame,pts,reps,output,strata,stratafield,doFrame)		## Call this if doing Terrestrial AIM
#RandomAquatic(WD,frame,pts,reps,output,strata,stratafield,doFrame)	## Call this if doing Aquatic AIM
q()									## Use this to terminate when running a batch file

