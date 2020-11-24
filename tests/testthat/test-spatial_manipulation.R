# Create a spatial points data frame to use
set.seed(420)
xcoords <- round(rnorm(n = 20, mean = 100, sd = 20))
set.seed(69)
ycoords <- round(rnorm(n = 20, mean = 100, sd = 20))

proj_wgs84 <- sp::CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")
proj_nad83 <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

points_spdf <- sp::SpatialPointsDataFrame(coords = cbind(xcoords, ycoords),
                                          data = data.frame("x" = xcoords,
                                                            "y" = ycoords,
                                                            stringsAsFactors = FALSE),
                                          proj4string = proj_nad83)

test_that("multiplication works", {
  expect_equal(2 * val, 4)
})
