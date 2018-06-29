### The "grog" project
#
# What I did. What data has been collected. What analysis can be
# done. What still needs to be done.
#
# The 'grog' project is active and should be opened whenever 
# working on the project.

# Otherwise: load the subset of the EEBO-TCP index that includes
# the 28 files used in 'grog'
load("index.rda")

# Now load the semantic model
load("geo.rda")

# Load the network graph
load("g.rda")

# Load the corresponding communities
load("wt.rda")

# Load the g-geo key
load("geo_nodes.rda")

# Load compose_nodes function
load("compose_nodes.rda")

# Load word-communities matrix
load("geo_c.rda")

# Load place context matrix
# Rows: 1668 toponyms
# Columns: 2831 divs from geography corpus
load("pc.rda")

# Load eebo place context matrix
# Rows: 1822 toponyms
# Columns: 2002 keywords (same as eebo)
load("eebo_pc.rda")

# Load list of coordinates
load("coordinates.rda")

# Load data from A37760
el = read.csv("A37760 edges.csv")
g1 = graph.data.frame(el, directed = F)  ## Creates a tree graph
polygons = read.csv("A37760 polygons.txt")

# Load placefulness records
load("./distance matrices/placefulness.rda")

### Now that data is loaded, what to do

## Easy to do word maps

# So, this vec gets very strong correlations (> .70)

# Now we need to calculate point distances using A43514
dim(coordinates$A43514)


### How I built geo_c and re-wrote wt to match
# Build communities semantic model
comms = as.integer(names(sort(table(wt$membership), decreasing = T)))[1:200]
geo_c = matrix(0, nrow(geo), length(comms))
pnames = c()
for (i in 1:length(comms)) {
  pname = names(sort(degree(g)[wt$membership == comms[i]], decreasing = T))[1]
  pnames = c(pnames,pname)
}
pnames[1] = "Europe"
pnames[8] = "Greece"
colnames(geo_c) = pnames
rownames(geo_c) = rownames(geo)
for (i in 1:length(comms)) {
  print(i)
  vec = rep(0, nrow(geo_c))
  vec = try(compose_nodes(nodes = which(wt$membership == comms[i])))
  vec = as.numeric(vec)
  if (is.na(sum(vec))) vec = rep(0, nrow(geo_c))
  geo_c[,i] = vec
}
totals = apply(geo_c, 2, sum)
ord = which(totals > 200)
geo_c = geo_c[,ord]

# Re-write wt membership with comms
for (i in 1:length(comms)) {
  wt$membership[which(wt$membership == comms[ord][i])] = pnames[i]
}


for (i in 1:ncol(mat)) {
  print(i)
  pname = colnames(mat)[i]
  pnames = V(g)$name[neighbors(g, V(g)[pname])]
  mat[pnames,pname] = 1
}


## Working with A41559 polygons
centx = c()
centy = c()
for (i in 1:nrow(coordinates$A41559)) {
  row = as.numeric(coordinates$A41559[i,2:9])
  centx = c(centx, mean(c(row[5], row[7])))
  centy = c(centy, mean(c(row[1], row[3])))
}

plot(c(0,360), c(-90,90), type = "n")
for (i in 1:nrow(coordinates$A41559)) {
  row = as.numeric(coordinates$A41559[i,2:9])
  rect(xleft = row[5], xright = row[7], ybottom = row[1], ytop = row[3])
}
text(centx, centy, labels = coordinates$A41559$clean_heads)
points(coordinates$A37551$LON, coordinates$A37551$LAT)
points(coordinates$A43514$LON, coordinates$A43514$LAT)

## Plot A37760 polygons
png("texture.png")
plot(c(0,360), c(-50,75), type = "n")
points(coordinates$A37551$LON, coordinates$A37551$LAT, pch = 20, col = "blue")
points(coordinates$A43514$LON, coordinates$A43514$LAT, pch = 20, col = "blue")
for (i in 1:nrow(polygons)) {
  row = as.numeric(polygons[i,2:5])
  rect(xleft = row[1], xright = row[2], ybottom = row[3], ytop = row[4])
}
centx = c()
centy = c()
for (i in 1:nrow(polygons)) {
   row = as.numeric(polygons[i,2:5])
   centx = c(centx, mean(c(row[1], row[2])))
   centy = c(centy, mean(c(row[3], row[4])))
}
# text(centx, centy, labels = polygons$Name)
dev.off()



# Collecting A43514 data
fp = "A43514.xml"
parsedText = htmlTreeParse(fp, useInternalNodes = TRUE)
nodes = getNodeSet(parsedText, "//div[@type='table_of_locations']")
loc_tables = nodes

rows = getNodeSet(loc_tables[[1]], "//row")


library(rgl)
earthRadius <- 3600
spheres3d(0,0,0, radius = earthRadius, lit=T,color="lightblue")
spheres3d(0,0,0,radius= earthRadius + .01,lit=FALSE,color="black",front="lines")
cosLat = cos(lat * pi/180)
cosLon = cos(lon * pi/180)
sinLon = sin(lon * pi/180)
sinLat = sin(lat * pi/180)
x = (earthRadius * cosLon * cosLat)
y = (earthRadius * cosLat * sinLon)
z = (earthRadius * sinLat)
spheres3d(x,y,z,col="red",radius=25)

earthRadius <- 3600
spheres3d(0,0,0, radius = earthRadius, lit=T,color="lightblue")
spheres3d(0,0,0,radius= earthRadius + .01,lit=FALSE,color="black",front="lines")
for (i in 1:nrow(polygons)) {
  lat = as.numeric(polygons[i,3:4])
  lon = as.numeric(polygons[i,1:2])
  cosLat = cos(lat * pi/180)
  cosLon = cos(lon * pi/180)
  sinLon = sin(lon * pi/180)
  sinLat = sin(lat * pi/180)
  x = (earthRadius * cosLon * cosLat)
  y = (earthRadius * cosLat * sinLon)
  z = (earthRadius * sinLat)
  polygon3d(x, y, z, col = i)
}

var cosLat = Math.cos(lat * Math.PI / 180.0);
var sinLat = Math.sin(lat * Math.PI / 180.0);
var cosLon = Math.cos(lon * Math.PI / 180.0);
var sinLon = Math.sin(lon * Math.PI / 180.0);
var rad = 500.0;
marker_mesh.position.x = rad * cosLat * cosLon;
marker_mesh.position.y = rad * cosLat * sinLon;
marker_mesh.position.z = rad * sinLat;

