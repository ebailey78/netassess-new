library(sp)
library(deldir)
library(maptools)
library(rgeos)

x <- read.csv("area.served/o3network.csv", stringsAsFactors = FALSE)
us <- readShapePoly("area.served/cb_2013_contus_nation_20m.shp")

ids <- x$MONITOR_ID
lats <- x$LATITUDE
longs <- x$LONGITUDE

tracts <- read.csv("area.served/2010_centroids.txt", stringsAsFactors = FALSE)
tracts <- SpatialPointsDataFrame(list(tracts$Longitude, tracts$Latitude), tracts)

voronoi <- function(ids, lats, longs, boundary) {
  
  df <- cbind(ids, lats, longs)
  points <- SpatialPointsDataFrame(list(longs, lats), data.frame(ids))
  
  crds <- points@coords
  
  if(!missing(boundary)) {
    bb <- bbox(boundary)
    rw <- as.numeric(t(bbox(boundary)))
    z <- deldir(crds[,1], crds[,2], rw=rw)  
  } else {
    z <- deldir(crds[,1], crds[,2])
  }
  
  w <- tile.list(z)
  polys <- vector(mode="list", length=length(w))
  
  for(i in seq(along = polys)) {
    pcrds = cbind(w[[i]]$x, w[[i]]$y)
    pcrds = rbind(pcrds, pcrds[1, ])
    polys[[i]] = Polygons(list(Polygon(pcrds)), ID=ids[i])
  }
  
  SP = SpatialPolygons(polys)
  voronoi = SpatialPolygonsDataFrame(SP, data=data.frame(x=crds[,1], y = crds[,2], row.names=sapply(slot(SP, "polygons"), function(x) slot(x, "ID"))))
  
  return(voronoi)
  
}

v <- voronoi(ids, lats, longs, us)
gg <- gIntersection(v, us, byid=TRUE)
plot(gg)

ov <- over(tracts, gg, returnList=TRUE)

