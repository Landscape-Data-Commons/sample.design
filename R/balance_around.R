    }
  }
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
keep_farthest <- function(existing_points_spdf,
                          new_points_spdf){
  # TODO: Sanitize
  if (!(class(existing_points_spdf) %in% "SpatialPointsDataFrame")) {
    stop("existing_points_spdf must be a spatial points data frame")
  }
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

  # Now that we have our indices to remove, let's do it as we combine points
  output <- rbind(existing_points_spdf[, common_varnames],
                  new_points_spdf[-removal_indices, common_varnames])

  return(output)
}


##  Ingest inputs, call functions to access XY coords and to identify and eliminate New points, and output
##          expanded, balanced design.
#' @param existing_points_spdf Spatial ponits data frame. The existing points that will be balanced around.
#' @param new_points_spdf Spatial points data frame. The points that will be compared against the existing points and selected from to create a balanced design.
#' @param stratafield Character string. The name of the variable in common between \code{existing_points_spdf} and \code{new_points_spdf} that contains stratum identities. This is used to balance by stratum. If \code{NULL} then balancing will not take strata into account. Defaults to \code{NULL}.
#' @param option Somethingorother
#' @param projection CRS object. The projection to force on the spatial objects. Defaults to \code{sp::CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")}.
#' @return A spatial points data frame containing all the points from \code{existing_points_spdf} and the selected points from \code{new_points_spdf}. The projection will match \code{projection}.
#' @export
balance_around <- function(existing_points_spdf,		## Name of existing points shapefile
                           new_points_spdf,		## Name of New points shapefile
                           stratafield = NULL,  	## Name of the stratum field in the ingested point files.  If set, then points will be balanced on a stratum by stratum basis.
                           ## If stratafield=NA, then ingested point files are assumed to represent an entire frame and spatial balance is based on the entire
                           ## collection of existing and New points.
                           projection = sp::CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")){
  # TODO: Sanitization (including reprojection)

  # Assign the codes that indicate if they're existing plots or freshly-drawn ones
  existing_pts_spdf@data[["TYPE"]] <- "EXISTING"
  new_points_spdf@data[["TYPE"]] <- "NEW"

  # Add in the fields that the points don't have for easy combination
  missing_vars_new <- names(existing_points_spdf@data)[!(names(existing_points_spdf@data) %in% names(new_points_spdf@data))]
  new_points_spdf@data[, missing_vars_new] <- NA
  missing_vars_existing <- names(new_points_spdf@data)[!(names(new_points_spdf@data) %in% names(existing_points_spdf@data))]
  existing_points_spdf@data[, missing_vars_existing] <- NA

  # Bind the 2 point files together
  pts <- rbind(existing_points_spdf, new_points_spdf)

  # What are the existing points' indices?
  extant_indices <- 1:nrow(existing_points_spdf@data)
  # Which are the new points' indices?
  new_indices <- (nrow(existing_points_spdf@data) + 1):(nrow(existing_points_spdf@data) + nrow(new_points_spdf@data))

  # Rename the date that the plot was sampled to "PREVDATE"
  pts@data[["PREVDATE"]] <- pts@data[["DATEVISITE"]]

  # Restrict to only relevant fields
  pts <- pts[ , c("TYPE",
                  stratafield,
                  "PLOTID",
                  "PLOTKEY",
                  "PRIMARYKEY",
                  "PANEL",
                  "PROJECTNAM",
                  "PREVDATE")]

  # Make these all NA for the new points
  new_points_spdf@data[new_indices, "PLOTID"] <- NA
  new_points_spdf@data[new_indices, "PLOTKEY"] <- NA
  new_points_spdf@data[new_indices, "PRIMARYKEY"] <- NA
  new_points_spdf@data[new_indices, "PROJECTNAM"] <- NA
  new_points_spdf@data[new_indices, "PREVDATE"] <- NA

  # Determine the number of existing points
  extant <- nrow(existing_points_spdf@data)

  # Add coordinates to the combined points for distance calculations
  pts <- get_coords(pts,
                    x_var = "XMETERS",
                    y_var = "YMETERS",
                    projection = sp::CRS("+proj=aea"))

  # This determines the number of New points to eliminate, and eliminates the points
    pts <- keep_farthest(existing_points_spdf = existing_points_spdf,
                         new_points_spdf = new_points_spdf)

  # Time to tweak the new points that were kept
  new_indices_remaining <- pts@data[["TYPE"]] == "NEW"

  # Get the existing plot ids and renumber them
  if (any(new_indices_remaining)) {
    plotids <- pts@data[new_indices_remaining, "PLOTID"]
    plotids <- gsub(plotids,
                    pattern = "\d*$",
                    replacement = "")
    pts@data[new_indices_remaining, "PLOTID"] <- paste0(plotids, 1:length(new_indices_remaining))
  }

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

  return(pts)
}
