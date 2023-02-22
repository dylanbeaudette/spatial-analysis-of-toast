library(terra)
library(sf)
library(colorspace)
library(aqp)
library(farver)
library(MASS)
library(magick)
library(grid)
library(cluster)
library(ape)
library(gstat)
library(rasterVis)
library(lattice)


# extent set from pixels
x <- rast('toast-grid.jpg')
ext(x)
plot(x)

# # one-time only, digitize toast centers, kind of tedious
# xy <- terra::click(x, n = 90, xy = TRUE, col = 'green')
# write.csv(xy, file = 'points.csv', row.names = FALSE)
# 

# convert to CIELAB
x.lab <- x
values(x.lab) <- grDevices::convertColor(values(x) / 255, from = 'sRGB', to = 'Lab')
names(x.lab) <- c('L', 'A', 'B')

# check: ok
plot(x.lab)

# read digitized toast-centers
xy <- read.csv('points.csv')

# setup IDs
names(xy) <- c('x', 'y', 'r', 'g', 'b')
ids <- expand.grid(1:9, LETTERS[1:10])
xy$id <- apply(ids, 1, function(i) {paste(i, collapse = '')})


# init spatvect
p <- vect(xy, geom = c('x', 'y'))

# buffer for extracting values
b <- buffer(p, width  = 15)

# graphical check: ok
plot(x)
points(p)
lines(b, col = 'green')
text(crds(p), p$id, pos = 3, cex = 0.66, col = 'red', offset = 1.5)


# mean LAB values within buffers
e <- extract(x.lab, b, fun = mean, na.rm = TRUE)

# fix names, IDs
names(e) <- c('id', 'L', 'A', 'B')
e$toast.id <- xy$id

# check
head(e)

# back to sRGB
.srgb <- grDevices::convertColor(e[, c('L', 'A', 'B')], from = 'Lab', to = 'sRGB', from.ref.white = 'D65', to.ref.white = 'D65')

e$r <- .srgb[, 1]
e$g <- .srgb[, 2]
e$b <- .srgb[, 3]

# Munsell notation
m <- rgb2munsell(e[, c('r', 'g', 'b')])
m$m <- sprintf('%s %s/%s', m$hue, m$value, m$chroma)
m$col <- rgb(e[, c('r', 'g', 'b')], maxColorValue = 1)

# annotate with Munsell colors
plot(x, mar = c(0, 0, 0, 0))
text(xy$x, xy$y, labels = m$m, col = invertLabelColor(m$col), cex = 0.7, font = 2)

# most frequent colors?
sort(table(m$m), decreasing = TRUE)

# mean toast colors on the L ~ A plane
plot(e$A, e$L, type = 'n')
points(e$A, e$L, pch = 15, col = m$col, cex = 4)
text(e$A, e$L, labels = e$toast.id, col = invertLabelColor(m$col), cex = 0.66)





d <- farver::compare_colour(
  e[, c('L', 'A', 'B')], 
  e[, c('L', 'A', 'B')], 
  from_space='lab', 
  white_from = 'D65', 
  method='cie2000'
)

# copy over SPC ids
dimnames(d) <- list(e$toast.id, e$toast.id)

# convert to dist object
d <- as.dist(d)

# NMDS
mds <- sammon(d)

# mean colors arranged via NMDS
plot(mds$points[, 1], mds$points[, 2], type = 'n')
points(mds$points[, 1], mds$points[, 2], pch = 15, col = m$col, cex = 4)
text(mds$points[, 1], mds$points[, 2], labels = e$toast.id, col = invertLabelColor(m$col), cex = 0.66)


par(mar = c(0, 0, 1, 0), mfcol = c(1, 2))
plot(x, mar = c(0, 0, 1, 0), main = 'Original')
plot(x, mar = c(0, 0, 1, 0), alpha = 1, main = 'Weighted Mean CIELAB')
plot(b, col = m$col, add = TRUE)


dev.off()



## crop original image into slices 
unlink('slices', recursive = TRUE)
dir.create('slices', recursive = TRUE)

for(i in 1:nrow(b)) {
  f <- file.path('slices', sprintf("%03d.png", i))
  crop(x, b[i, ], filename = f, overwrite = TRUE) 
}

# remove side-car files
unlink(x = 'slices/*.xml')

par(mar = c(0.1, 0.1, 1.5, 0.1))
plot(mds$points[, 1], mds$points[, 2], type = 'n', asp  = 1, axes = FALSE, xlab = '', ylab = '', main = 'Non-metric Multidimensional Scaling of Toast')
box()
abline(h = 0, v = 0, lty = 2)

for(i in 1:90) {
  f <- file.path('slices', sprintf("%03d.png", i))
  im <- image_read(f)
  
  .off <- 1.5
  rasterImage(
    im, 
    xleft = mds$points[i, 1] - .off, 
    ybottom = mds$points[i, 2] - .off, 
    xright = mds$points[i, 1] + .off, 
    ytop = mds$points[i, 2] + .off
  )
  
}




## load previously cropped slices
loaf <- sprc(lapply(list.files('slices/', full.names = TRUE), rast))

# single slice of toast
.n <- 45
.sl <- loaf[.n]
plot(.sl)


## variogram models

# SpatRast -> sf -> sp -> gstat
# sample
s.pts <- spatSample(.sl, size = 1000, method = 'random', as.points = TRUE, lonlat = FALSE)
names(s.pts) <- c('r', 'g', 'b')

.lab <- grDevices::convertColor(as.data.frame(s.pts)[, c('r', 'g', 'b')] / 255, from = 'sRGB', to = 'Lab', from.ref.white = 'D65', to.ref.white = 'D65')

.lab <- data.frame(.lab)
names(.lab) <- c('L', 'A', 'B')

s.pts <- cbind(s.pts, .lab)

# remove NA
idx <- which(!is.na(s.pts$L))
s.pts <- s.pts[idx, ]

# plot(s.pts)

s.pts <- as_Spatial(st_as_sf(s.pts))

v <- automap::autofitVariogram(L ~ 1, s.pts)

plot(v)

# plot(v$var_model, cutoff = max(v$var_model$psill))




## conditional simulation
# requires reasonable variogram model
# requires random sampling

# make target grid as SpatialGridDataFrame
# terra -> raster
.r <- raster::raster(.sl)[[1]]
values(.r) <- 1
.r <- aggregate(.r, fact = 2)


g <- as(.r, 'SpatialGridDataFrame')

# there is no CRS, set to NA
proj4string(s.pts) <- ''

# conditional simulation
sim <- krige(formula = L ~ 1, locations = s.pts, newdata = g, model = v$var_model, nsim = 90, nmax = 10, beta = mean(s.pts$L))

# sp::spplot(sim, as.table = TRUE)

# convert to raster stack for plotting
sim <- raster::stack(sim)


# plot
.title <- sprintf('Conditional Simulation of Slice %s\nCIELAB L-coordinate', xy$id[.n])
levelplot(sim, scales = list(draw = FALSE), col.regions = viridisLite::viridis(100), layout = c(10, 9), names.attr = rep('', times = length(names(sim))), colorkey = TRUE, main = .title)

## 








# 
# # divisive hierarchical clustering
# dd <- diana(d)
# hc <- as.hclust(dd)
# dd.p <- as.phylo(hc)
# 
# par(mar=c(0,0,0,0), bg='black', fg='white')
# 
# plot(dd.p, show.tip.label = TRUE, cex = 0.8)






