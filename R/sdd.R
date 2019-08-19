terra_sample_frame <- function(spdf,
                               project_name,
                               type,
                               description,
                               group = NULL,
                               projection = sp::CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")){
  valid_types <- c("ESR", "LUP", "INTS", "LR", "TRT", "FUELS")
  if (!(type %in% valid_types)) {
    stop("The valid values for type are: '", paste(valid_types, collapse = "' '"), "'")
  }

  # Make sure that project_name is a character string
  project_name <- as.character(project_name)
  # Check for illegal characters in project_name
  if (grepl(project_name, pattern = "\\W")) {
    stop("Project name must be a character string containing only alphanumeric characters and _")
  }

  # Reproject if necessary
  if (!(class(polygons_spdf) %in% "SpatialPolygonsDataFrame")) {
    stop("polygons_spdf must be a spatial polygons data frame")
  }
  if (!identical(projection, polygons_spdf@proj4string)) {
    polygons_spdf <- sp::spTransform(polygons_spdf,
                                     projection)
  }

  # Make and add the the frame ID
  frame_id <- paste(project_name, type, sep = "_")
  spdf@data[["TERRA_SAMPLE_FRAME_ID"]] <- frame_id

  # Dissolve if necessary
  if (nrow(spdf@data) > 1) {
    spdf <- dissolve_spdf(spdf = spdf,
                          dissolve_field = "TERRA_SAMPLE_FRAME_ID")
  }


  spdf@data[["SAMPLE_FRAME_DESC"]] <- description
  spdf@data[["SAMPLE_FRAME_GROUP"]] <- group
  spdf@data[["SAMPLE_FRAME_AREA_SQKM"]] <- add_area(polygons = spdf,
                                                    area_ha = FALSE,
                                                    area_sqkm = TRUE,
                                                    area_acres = FALSE,
                                                    byid = TRUE)@data[["AREA.SQKM"]]
  spdf@data[["SAMPLE_FRAME_TYP"]] <- type

  return(spdf)
}

terra_sample_points <- function(spdf,
                                strata_spdf,
                                project_spdf,
                                plotid_field){

  input_spdf <- spdf
  # TERRA_PLOT_ID
  spdf@data[["TERRA_PLOT_ID"]] <- NULL
  # PLOT_KEY
  spdf@data[["PLOT_KEY"]] <- NULL
  # PLOT_NM
  spdf@data[["PLOT_NM"]] <- input_spdf@data[[plotid_field]]
  # TERRA_PRJCT_AREA_ID
  # TERRA_CLSTR_ID
  # TERRA_STRTM_ID
  # TERRA_STRTM_VAL_ID
  # TERRA_SAMPLE_FRAME_ID
  # TERRA_TERRADAT_ID
  # ADMIN_ST
  # ADM_UNIT_CD
  # DSGN_STRTM_NM
  # ACTL_STRTM_NM
  # MDCATY
  # SAMPLE_DSGN_ALGRTHM
  # PT_DRAW
  over_vector <- grepl(input_spdf@data[["PANEL"]], pattern = "oversample", ignore.case = TRUE)
  base_vector <- !over_vector
  spdf@data[["PT_DRAW"]][over_vector] <- "OVER"
  spdf@data[["PT_DRAW"]][base_vector] <- "BASE"
  # PANEL

  # INTL_PT_WGT

  # FINAL_DESIG
  spdf@data[["FINAL_DESIG"]] <- NULL
  # DT_VST
  spdf@data[["DT_VST"]] <- NULL
  # COMMENT
  # REVISIT
  # REVISIT_ORIG_PLOTID
  # REVISIT_TERRA_PK
  # TARGET
  spdf@data[["TARGET"]] <- "NO"
  # TARGETDESCRIPTION
  spdf@data[["TARGETDESCRIPTION"]] <- NULL
}

terra_strtfctn <- function(spdf){

}
