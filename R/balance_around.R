##  GenDistLists() takes the distance matrix and orders the New points from closest to farthest for each existing GRTS point.

##  Somewhat tricky schema, so the following is provided to help us remember what the heck is going on here...

##  dist lists the distances between each of the existing N GRTS pts and the new pts.  Each column is an existing GRTS pts for a total of N columns.
#   NRows will equal the number of new points minus the number of existing GRTS pts.  Each row in dist is the accession point-number of the new points starting at N+1.

## disto lists the accession point-number from the nearest to the farthest.  Each column is the same as dist, each row however is ordered from closest to
## farthest and the numeric value is the accession point number.

## E.g., say dist has 52 columns and 104 rows (52 existing GRTS pts and 104 new points).  Column 1 is the first existing GRTS point.
##       Row 1 is the first new point which is row #53 in apts, row 104 is the last new point.
##       disto has 52 columns and 104 rows.  Again, the first column is the first existing GRTS point.  The first row of this column say is 9, meaning that
##       row 9 in dist (col=1, row=9) will have the smallest distance.  This is the smallest for the first existing GRTS point.  Column 2, row 1 is say 70 which
##       means dist(col=2, row=70) will show the smallest distance for the second existing GRTS point.  In apts, these smallest points will be entry 9+52 (61) and 70+52 (122). ETC..

GenDistLists<-function(a,extant,apts)
{
  displ<-extant+1
  dist<-NULL
  for(c in 1:extant) {
    z<-data.frame(VAR=a[displ:nrow(apts),c])
    zo<-data.frame(order(z))

    if(c==1) {
      dist<-z
      disto<-zo
    }else {
      dist<-cbind(dist,z)
      disto<-cbind(disto,zo)
    }
  }
  ret<-cbind(dist,disto)	## odd way to combine and return results, but it works!
  return(ret)
}
##################################################
## GetClosestPts - select the existing pts that spatially best matches the dispersion of a new GRTS draw.
GetClosestPts<-function(apts,stratafield)
{

  rel<-NULL
  names(apts)[names(apts)==stratafield]<-"THESTRATA"
  slist<-unique(apts$THESTRATA)
  for(k in slist) {
    checkE<-apts[apts$CODE==1 & apts$THESTRATA==k ,]
    checkN<-apts[apts$CODE==2 & apts$THESTRATA==k ,]
    checkE$USED<-0

    for(i in 1:nrow(checkN) ) {
      mdist<-99999999999
      savej<-0
      for(j in 1:nrow(checkE) ) {
        if(checkE$USED[j]!=1) {
          dist<- (checkN$XMETERS[i] - checkE$XMETERS[j])*   (checkN$XMETERS[i] - checkE$XMETERS[j])
          dist<-dist+ ( (checkN$YMETERS[i] - checkE$YMETERS[j])*  (checkN$YMETERS[i] - checkE$YMETERS[j]) )
          if(dist<mdist) {
            mdist<-dist
            savej<-j
          }## dist<
        }## checkE$USED
      }## for j
      checkE$USED[savej]<-1
    }## for i
    checkE<-checkE[checkE$USED==1 ,]
    if(k==slist[1]){
      rel<-checkE
    }else {
      rel<-rbind(rel,checkE)
    }
  } ## for k
  rel$USED<-NULL
  return(rel)
}
#############################################################################################
## FindClosest() - the work-horse that determines the unique New points that are the closest to existing points.
##                 After identifying these points, eliminates them from the combined points file (apts).

FindClosest<-function(apts,dist,disto,extant)
{
  ## Best to loop.  If shared points, keep the closest and increment the row index for the !closest.  Keep looping until no more shared pts.

  if(extant>1) {
    z<-as.numeric(disto[1,])		## Start with the closest
  }else z<-as.numeric(disto[1])

  zindex<-NULL
  delta<-nrow(apts)-extant
  zindex[1:delta]<-1			## Row index for each column (each existing GRTS pts); we start with row index ==1 (the closest)

  resolve<-1			## while loop continues until we have identified N unique New points, where N is the no. of existing GRTS points
  while(resolve) {
    resolve<-0
    for(i in 1:length(z)){
      tr<-0
      for(j in 1:length(z)) {
        if(i!=j) {
          if(z[i]==z[j]){		## Shared point (i.e., the closest New point is the same for 2 of the existing GRTS pts)
            resolve<-1
            tr<-1
            dist1<-dist[zindex[i],i]
            dist2<-dist[zindex[j],j]
            #print(paste(i,dist1,j,dist2,sep=" "))
            if(dist1<dist2) {			## Find the existing GRTS pt !closest to this New point, increment row index, then set
              ## z$VAR to the point accession number (of the New pts) of the next closest
              zindex[j]<-zindex[j]+1
              z[j]=as.numeric(disto[zindex[j],j])		## Pick up the next closest point.  row index is in zindex[]
            }else {
              zindex[i]<-zindex[i]+1
              z[i]=as.numeric(disto[zindex[i],i])
            }
          }
        }
      }
      #if(tr)print(paste("Not Unique ",z[i],sep=" "))
    }
  }

  ## Check to ensure we have selected nrow(apts) points (New pts that are closest to existing GRTS pts) to delete
  y<-unique(z)
  if(length(y) != extant ) {
    print(paste("ERROR, we have not eliminated sufficient pts", y, nrow(apts)-extant,sep=" "))
    q()
  }

  apts$DELETE<-0
  for(i in y) {
    apts$DELETE[i+extant]<-1
  }

  apts<-apts[apts$DELETE==0 ,]
  apts$DELETE<-NULL

  return(apts)
}
##########################################################################################
NN<-function(apts,extant,stratafield)	# Derive New points that are closest to existing GRTS points, such that a New point is selected only once.
  # Then eliminates these New points to produce an expanded, balanced design.
{

  if(is.na(stratafield)) {						## If NA, then we are dealing with just 1 frame area - all points will be considered
    a<-dist_matrix(apts)		## Matrix of distances among every point.
    a[a == 0] <- Inf
    ret<-GenDistLists(a,extant,apts)	## Generate the distance and New-point accession order lists for further processing
    dist<-ret[,1:extant]			## The 1:extant columns are Euclidean distance
    disto<-ret[,(extant+1):(extant*2)]	## The extant+1:extant*2 columns are the accession numbers of New points, where the accession number
    ## in apts is actually the number plus extant (e.g., accession = 9 is actually 9 + 52 when
    ## there are 52 existing GRTS points).
    ## Find the unique list of New points that are closest to each existing GRTS point, then eliminate them from apts
    apts<-FindClosest(apts,dist,disto,extant)
  }else {								## else Balance by strata
    savepts<-NULL
    names(apts)[names(apts)==stratafield]<-"THESTRATA"
    slist<-unique(apts$THESTRATA)
    for(i in slist) {
      tr<-0
      checkE<-apts[apts$CODE==1 & apts$THESTRATA==i ,]		## Compare number of existing and new points by strata
      checkN<-apts[apts$CODE==2 & apts$THESTRATA==i ,]
      if(nrow(checkE)>0) {		## Any existing points?
        if(nrow(checkN)>nrow(checkE)) {			## The number of New points must be > number of existing points, else nothing to do here
          Spts<-rbind(checkE,checkN)		## Stratum pts file
          a<-dist_matrix(Spts)		## All of the remaining dups the processing described directly above, but here its by stratum
          a[a == 0] <- Inf
          extant<-nrow(checkE)
          ret=GenDistLists(a,extant,Spts)
          dist<-ret[,1:extant]
          disto<-ret[,(extant+1):(extant*2)]
          Spts<-FindClosest(Spts,dist,disto,extant)
          #Spts$THESTRATA<-NULL		## clean up
          tr<-1
        }
      } # if nrow(checkE)

      if(!tr) {			## Ensures we carry-over all the existing pts and/or all New points that were not eliminated
        if(nrow(checkE)>0)Spts<-checkE
        if(nrow(checkN)>0)Spts<-checkN
        if(nrow(checkE)>0 & nrow(checkN)>0)Spts<-rbind(checkE,checkN)
        #Spts$THESTRATA<-NULL
      }
      if(is.null(savepts)) {
        savepts<-Spts
      }else {
        savepts<-rbind(savepts,Spts)
      }


    }## for i in slist
    apts<-savepts
  } #if else

  ## Clean up
  apts$XMETERS<-NULL
  apts$YMETERS<-NULL
  apts$LONGITUDE<-NULL
  apts$LATITUDE<-NULL

  return(apts)
}
#######################################################################################
##  Ingest inputs, call functions to access XY coords and to identify and eliminate New points, and output
##          expanded, balanced design.

BalancePTS<-function(layerE,		## Name of existing points shapefile
                     layerN,		## Name of New points shapefile
                     stratafield,  	## Name of the stratum field in the ingested point files.  If set, then points will be balanced on a stratum by stratum basis.
                     ## If stratafield=NA, then ingested point files are assumed to represent an entire frame and spatial balance is based on the entire
                     ## collection of existing and New points.

                     output,		## This is the output shapefile.
                     option		## 1 or 2- see top of script for explanation
)

{
  pts<-Ingest(layerE)			## Ingest the existing points
  pts2<-Ingest(layerN)			## Ingest the New points
  pts$CODE<-1				## CODE assignment
  pts2$CODE<-2
  pts$PREVDATE<-pts$DATEVISITE
  pts<-pts[ ,c("CODE",stratafield,"PLOTID","PLOTKEY","PRIMARYKEY","PROJECTNAM","PREVDATE")]
  pts2$PLOTID<-NA
  pts2$PLOTKEY<-NA
  pts2$PRIMARYKEY<-NA
  pts2$PROJECTNAM<-NA
  pts2$PREVDATE<-NA

  pts2<-pts2[ ,c("CODE",stratafield,"PLOTID","PLOTKEY","PRIMARYKEY","PROJECTNAM","PREVDATE")]

  pts<-rbind(pts,pts2)			## Bind the 2 point files together
  extant<-pts[pts$CODE==1 ,]
  extant<-nrow(extant)			## Determine the number of existing points

  pts<-GetXY(pts)			## Access point location in meters (UTM)
  if(option==1) pts<-NN(pts,extant,stratafield)	## This determines the number of New points to eliminate, and eliminates the points


  ## Option 2.
  #########TO derive the best spatial balance, skip the above call to NN and do the following.
  ##          This was designed for/is most useful whenever you are selecting x revisit points per stratum from a e.g. 5-yr design.
  ##          See selectpts.R which derives the number of points per strata and extracts the points. However, sometimes the original
  ##          design is not balanced, so selectpts.r output is very unbalanced.  Here we skip the extraction portion of selectpts.r and
  ##          do what we can to ID the best set of existing points (most spatially balanced) given a 'template' GRTS example (LayerN), where this
  ##          template has the exact number of points we want by strata.
  if(option==2) pts<-GetClosestPts(pts,stratafield)
  }

  ###########################################  WHERE WE FORMAT THE DATA.  This generally needs to be customized.
  new_indices_remaining <- pts@data[["TYPE"]] == "NEW"

  # Get the existing plot ids and renumber them
  plotids <- pts@data[new_indices_remaining, "PLOTID"]
  plotids <- gsub(plotids,
                  pattern = "\d*$",
                  replacement = "")
  pts@data[new_indices_remaining, "PLOTID"] <- paste0(plotids, 1:length(new_indices_remaining))

  # TODO: Rename within strata.
  # Does this mean using over() with stratification polygons to determine new stratification assignments for old points?


  # None of the points are considered visited now!
  pts@data[["DATEVIS"]] <- ""
  pts@data[["EVALSTA"]] <- "NotEval"
  pts@data[["FINAL_DESI"]] <- ""

  # Add coordinates!
  pts <- get_coords(pts,
                    x_var = "LAT",
                    y_var = "LONG",
                    projection = projection)

  pts$ORDERDes<-pts$COUNTER
  pts$COUNTER<-NULL
  writeOGR(pts, ".",output,driver="ESRI Shapefile",overwrite_layer=T)	# output the final shapefile
  return(pts)
}
