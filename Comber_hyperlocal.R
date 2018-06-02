
### Load data and libraries
### set up data.sp etc
library(GISTools)
library(raster)
library(OpenStreetMap)
library(rgdal)
library(spdep)
library(GWmodel)
library(tidyverse)
library(SpatialTools)
library(repmis)
library(AICcmodavg)


# 1. load data from GitHub 
source_data("https://github.com/lexcomber/CAS_GW_Training/blob/master/Soils.RData?raw=True")
source_data("https://github.com/lexcomber/CAS_GW_Training/blob/master/boundary.RData?raw=True")
# or locally
proj. <- CRS("+proj=tmerc +lat_0=0 +lon_0=108 +k=1 +x_0=500000 +y_0=0 +ellps=krass +units=m +no_defs ")
data.sp <- SpatialPointsDataFrame(coords, data = data.frame(data), proj4string = proj.)

# transform the data
data.sp@data$TNPC <- log(data.sp@data$TNPC+0.0001)
data.sp@data$TPPC <- (data.sp@data$TPPC)^0.5
data.sp@data$SOCgkg <- log(data.sp@data$SOCgkg)
data.sp@data$ClayPC <- (data.sp@data$ClayPC)^0.5
data.sp@data$NO3Ngkg <- log(abs(data.sp@data$NO3Ngkg))
data.sp@data$NH4Ngkg <- log(data.sp@data$NH4Ngkg)

## Figure 1
crs.val <- CRS("+proj=longlat ")
# Map the data
fac = 0.001 # in case a bit is needed
tmp <- spTransform(data.sp, crs.val)
ul <- as.vector(cbind(bbox(tmp)[2,2]+fac, bbox(tmp)[1,1]-fac))
lr <- as.vector(cbind(bbox(tmp)[2,1]-fac, bbox(tmp)[1,2]+fac))
MyMap <- openmap(ul,lr,zoom = NULL, type = "bing")
save("MyMap", file = "MyMap.RData")
setEPS()
postscript("F1.eps", width = 5, height = 6) 
plot(MyMap, removeMargin=T) 
plot(spTransform(data.sp, osm()),add = T, col="white", pch = 19, cex =1)
scalebar(1000, type = "bar", xy = c(12286385, 4690886), divs = 4, col = "white", below = "m")
dev.off()

### 2. Linear Regression
terms <- names(data)[c(4:11)]
## OLS N percentage
regmod <- paste(terms[1], "~")
for ( i in 3: length(terms) ) {
	if ( i != length(terms) ) regmod <- paste(regmod, terms[i],  "+")
	if ( i == length(terms) ) regmod <- paste(regmod, terms[i])
}
regmod.n <- as.formula(regmod)
## OLS P percentage
regmod <- paste(terms[2], "~")
for ( i in 3: length(terms) ) {
	if ( i != length(terms) ) regmod <- paste(regmod, terms[i],  "+")
	if ( i == length(terms) ) regmod <- paste(regmod, terms[i])
}
regmod.p <- as.formula(regmod)

# OLS regression
mN <- lm(regmod.n, data = data.sp@data)
mP <- lm(regmod.p, data = data.sp@data)
# set up for stepwaise later
step.N <- as.formula(mN)
step.P <- as.formula(mP)
step.N1 <- (as.character(step.N)[[3]])
step.P1 <- (as.character(step.P)[[3]])
# OLS regression using GWR
EUDM <-  gw.dist(coordinates(data.sp))
gwr.lm <- gwr.basic(regmod.n, data = data.sp, 
  kernel = "boxcar", bw = max(EUDM+1000),
  adaptive = F, dMat = EUDM)
gwr.lm$GW.diagnostic$gw.R2
gwr.lm2 <- gwr.basic(regmod.p, data = data.sp, 
  kernel = "boxcar", bw = max(EUDM+1000),
  adaptive = F, dMat = EUDM)
# Tables 1 and 2
tab1 <- round(summary(mN)$coefficients, 4)
tab2 <- round(summary(mP)$coefficients, 4)
# do stepwise 
step.mod.N <- stepAIC(lm(regmod.n, data), trace = F)
tab1a <- round(summary(step.mod.N)$coefficients, 4)
index <- match(rownames(tab1a), rownames(tab1))
m.tmp <- matrix(ncol = 4, nrow = nrow(tab1))
m.tmp[index, ] <- tab1a 
colnames(m.tmp) <- colnames(tab1)
tab1 <- data.frame(tab1, m.tmp)

step.mod.P <- stepAIC(lm(regmod.p, data), trace = F)
tab2a <- round(summary(step.mod.P)$coefficients, 4)
index <- match(rownames(tab2a), rownames(tab2))
m.tmp <- matrix(ncol = 4, nrow = nrow(tab1))
m.tmp[index, ] <- tab2a 
colnames(m.tmp) <- colnames(tab2)
tab2 <- data.frame(tab2, m.tmp)

### R2 adjusted R2 and AIC
t1.i <- paste("R2:", round(unlist(summary(mN)[8]), 3), 
	", adj R2:", round(unlist(summary(mN)[9]), 3),
	", AIC:", round(AIC(mN), 1), sep = "")
t1.ii <- paste("R2:", round(unlist(summary(step.mod.N)[8]), 3), 
	", adj R2:", round(unlist(summary(step.mod.N)[9]), 3),
	", AIC:", round(AIC(step.mod.N), 1), sep = "")

t2.i <- paste("R2:", round(unlist(summary(mP)[8]), 3), 
	", adj R2:", round(unlist(summary(mP)[9]), 3),
	", AIC:", round(AIC(mP), 1), sep = "")
t2.ii <- paste("R2:", round(unlist(summary(step.mod.P)[8]), 3), 
	", adj R2:", round(unlist(summary(step.mod.P)[9]), 3),
	", AIC:", round(AIC(step.mod.P), 1), sep = "")
## Summarise OLS
tab2 <- tab2[,-c(3,7)]
tab1 <- tab1[,-c(3,7)]
write.csv(tab1, file = "Tab1.csv")
write.csv(tab2, file = "Tab2.csv")

### 3. GWR Regression
# 3.1 N
bw.n <- bw.gwr(regmod.n, data = data.sp, 
  kernel = "bisquare", adaptive = F, approach = "AIC") 
gwr.n <- gwr.basic(regmod.n, data = data.sp, bw = bw.n, 
  kernel = "bisquare", adaptive = F)
# test for collinearity
gw.col.n <- gwr.collin.diagno(regmod.n, data = data.sp, bw = bw.n, 
  kernel = "bisquare", adaptive = F)
  
# 3.2 P
bw.p <- bw.gwr(regmod.p, data = data.sp, 
  kernel = "bisquare", adaptive = F, approach = "AIC") 
gwr.p <- gwr.basic(regmod.p, data = data.sp, bw = bw.p, 
  kernel = "bisquare", adaptive = F)
# test for collinearity
gw.col.p <- gwr.collin.diagno(regmod.p, data = data.sp, bw = bw.p, 
  kernel = "bisquare", adaptive = F)

# sumamry tables
tab3 <- (gwr.n$SDF@data[, 1:(nrow(tab1))])
tab3 <- (apply(tab3, 2, summary))
tab3 <- t(round(tab3, 4))
tab3 <- cbind(tab3, IQR = tab3[,5]-tab3[,2], Global = tab1[,1])

tab4 <- (gwr.p$SDF@data[, 1:(nrow(tab2))])
tab4 <- (apply(tab4, 2, summary))
tab4 <- t(round(tab4,4))
tab4 <- cbind(tab4, IQR = tab4[,5]-tab4[,2], Global = tab2[,1])

write.csv(tab3, file = "Tab3.csv")
write.csv(tab4, file = "Tab4.csv")

#### Figures 2 and 3: GWR coefficient surfaces
# plot function
## set up plot
boundary@data$id = rownames(boundary@data)
boundary.points = fortify(boundary, region="id")
boundary.df = left_join(boundary.points, boundary@data, by="id")

plot.function.gg <- function(i = 2, sh = "GnBu" ) {
vals <- gwr.sp@data[,i]
tit <- names(gwr.sp@data)[i]
index <- gwr.sp@data[,i+19] > 1.96 | gwr.sp@data[,i+19] < -1.96 
gwr.sp.index <- gwr.sp[index,]
p <- ggplot(boundary.df) + 
      geom_polygon(aes(x=long, y=lat), colour="grey", fill="grey") +
      coord_equal() +
      geom_point(data = gwr.sp@data, aes(x = X, y = Y, colour= vals), size = 1) +
      scale_colour_distiller(type="seq", direction = 1, 
      	palette = sh) +
      geom_point(data = gwr.sp.index@data, aes(x = X, y = Y), size = 0.3) +
      labs(subtitle = tit) +
      theme(axis.title.x=element_blank(),
      	axis.text.x=element_blank(),
    	axis.ticks.x=element_blank(), 
    	axis.title.y=element_blank(),
    	axis.text.y=element_blank(),
    	axis.ticks.y=element_blank(),
    	panel.background = element_blank(), 
    	legend.direction = "horizontal", 
    	legend.position = "bottom", 
    	legend.key.width = unit(0.9, "cm"),
    	legend.title=element_blank())
return(p)}

## multiplot function
multiplot2 <- function(plot.list, file, cols=3, layout=NULL) {
  library(grid)
  numPlots = length(plot.list)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }
 if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plot.list[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

sh.list <- c("YlOrBr", "YlGnBu","PuBuGn", "YlOrRd", "OrRd", "BuGn") 
### Figure 2
gwr.sp <- gwr.n$SDF
X <- coordinates(gwr.sp)[,1]
Y <- coordinates(gwr.sp)[,2]
gwr.sp <- spCbind(gwr.sp, data.frame(X, Y))
setEPS()
postscript("F2.eps", width = 7, height = 7) 
par(mfrow = c(2,3))
par(mar = c(0,0,0,0))
p1 <- plot.function.gg(i = 2, sh = sh.list[1] ) 
p2 <- plot.function.gg(i = 3, sh = sh.list[2] ) 
p3 <- plot.function.gg(i = 4, sh = sh.list[3] ) 
p4 <- plot.function.gg(i = 5, sh = sh.list[4] ) 
p5 <- plot.function.gg(i = 6, sh = sh.list[5] ) 
p6 <- plot.function.gg(i = 7, sh = sh.list[6] ) 
multiplot2(list(p1,p2,p3,p4,p5,p6), cols = 3)
dev.off()

### Figure 3
gwr.sp <- gwr.p$SDF
X <- coordinates(gwr.sp)[,1]
Y <- coordinates(gwr.sp)[,2]
gwr.sp <- spCbind(gwr.sp, data.frame(X, Y))
setEPS()
postscript("F3.eps", width = 7, height = 7) 
par(mfrow = c(2,3))
par(mar = c(0,0,0,0))
p1 <- plot.function.gg(i = 2, sh = sh.list[1] ) 
p2 <- plot.function.gg(i = 3, sh = sh.list[2] ) 
p3 <- plot.function.gg(i = 4, sh = sh.list[3] ) 
p4 <- plot.function.gg(i = 5, sh = sh.list[4] ) 
p5 <- plot.function.gg(i = 6, sh = sh.list[5] ) 
p6 <- plot.function.gg(i = 7, sh = sh.list[6] ) 
multiplot2(list(p1,p2,p3,p4,p5,p6), cols = 3)
dev.off()

### 4. Hyper local GWR Regression
# This will take about 40 minutes to run

##### START Hyper local GWR Regression
#---DISTANCES
# do local model colibration of GWR
# set up inputs 
dMat <- dist2(coordinates(data.sp), coordinates(data.sp))
bwd.range <- c(seq(200,max(dMat),by=50))

# names(data)
terms <- names(data)[c(4:11)]
## OLS N percentage
regmod <- paste(terms[1], "~")
for ( i in 3: length(terms) ) {
  if ( i != length(terms) ) regmod <- paste(regmod, "+", terms[i],  "+")
  if ( i == length(terms) ) regmod <- paste(regmod, terms[i])
}
regmod.n <- as.formula(regmod)
#summary(lm (regmod.n, data.sp@data))
## OLS P percentage
regmod <- paste(terms[2], "~")
for ( i in 3: length(terms) ) {
  if ( i != length(terms) ) regmod <- paste(regmod, "+", terms[i],  "+")
  if ( i == length(terms) ) regmod <- paste(regmod, terms[i])
}
regmod.p <- as.formula(regmod)
#summary(lm (regmod.p, data.sp@data))

regmod.all.N <- as.formula(regmod.n)
regmod.all.N.orig <- regmod.all.N

regmod.all.P <- as.formula(regmod.p)
regmod.all.P.orig <- regmod.all.P

######## 1. N

# and outputs
aic.score.nd <- matrix(nrow=length(data.sp), ncol=length(bwd.range), data = 0)
model.mat.nd <- matrix(nrow=length(data.sp), ncol=length(bwd.range), data = 0)

## local Model and AIC
for ( i in 1:length(data.sp)) {
  # for each pt
  pt.i <- data.sp[i,]
  dMat.i <- dMat[i,]
  
  # for each bw
  # st <- Sys.time()
  for ( j in 1:length(bwd.range) ) {
    ## sort out the clostest
    index.j <- dMat.i < bwd.range[j]
      #j.integer.val <- round(j*(nrow(data.sp)/100), 0)
    #index.j <- order(dMat.i)[ 1:j.integer.val]  
    data.j <- data.sp[index.j, c(4,6:11)]
    #plot(data.sp)
    #plot(pt.i, pch = 19, add = T)
    #plot(data.j, col = "red", add = T)
    # get rid of any singletons
    regmod.all.N = regmod.all.N.orig
    tmp.test <- lapply(data.j@data, unique)
    tmp.test <- lapply(tmp.test, length)
    tmp.test <- as.vector(unlist(tmp.test))
    
    ## Get rid of any singularities 
    # droped from regmod 
    if (any(tmp.test == 1)) {
      index.jj <- which(tmp.test == 1)
      tit <- names(data.j)[index.jj]
      terms.j <- terms[2:length(terms)]
      index.jj <- match(tit, terms.j)
      terms.j <- terms.j[-index.jj]
      ltj <- length(terms.j)     
      regmod <- paste(terms[1], "~")
      for ( k in c(1:ltj) ) {
        if ( k != ltj ) regmod <- paste(regmod, terms.j[k],  "+")
        if ( k == ltj ) regmod <- paste(regmod, terms.j[k])
      }
      regmod.all.N <- regmod
    } 
    
    ## GWeight the data 
    dists.j <- gw.dist(coordinates(pt.i), coordinates(data.j))
    weights.j <- as.vector(gw.weight(dists.j,bwd.range[j], kernel = "bisquare", adaptive = F))
    #plot(dists.i, weights.i)
    #plot(data.j, pch = 19, cex = weights.j)
    data.j@data <- (data.j@data) * as.matrix(weights.j, ncol = 1)
 
 	#plot(data.j, pch = 19, cex = weights.j, add = T)
        
    if ( class(tryCatch(stepAIC(lm(regmod.all.N, data.j@data), 
    		trace = F), error = function(e) e))[1] != "simpleError" ) {
    		step.j <- stepAIC(lm(regmod.all.N, data.j@data), trace = F) 
    		regmod.j <- as.formula(step.j)
    		tmp <- lm(regmod.j, data.j@data)
		aic.score.nd[i,j] <-	AIC(tmp)
    		model.mat.nd[i,j] <- as.vector(as.character(regmod.j)[3]) } else {
		  
		step.j <- regmod.all.N
		regmod.j <- as.formula(step.j)
  		tmp <- lm(regmod.j, data.j@data)
		aic.score.nd[i,j] <-	AIC(tmp)
		model.mat.nd[i,j] <- as.vector(as.character(regmod.all.N)[3]) }
  }
  #Sys.time() - st
  if(i %% 5 == 0) cat(i)
}

######## 2. P

# and outputs
aic.score.pd <- matrix(nrow=length(data.sp), ncol=length(bwd.range), data = 0)
model.mat.pd <- matrix(nrow=length(data.sp), ncol=length(bwd.range), data = 0)

## local Model and AIC
for ( i in 1:length(data.sp)) {
  # for each pt
  pt.i <- data.sp[i,]
  dMat.i <- dMat[i,]
  
  # for each bw
  # st <- Sys.time()
  for ( j in 1:length(bwd.range) ) {
    ## sort out the clostest
    index.j <- dMat.i < bwd.range[j]
    
    #j.integer.val <- round(j*(nrow(data.sp)/100), 0)
    #index.j <- order(dMat.i)[ 1:j.integer.val]  
    data.j <- data.sp[index.j, c(5,6:11)]
    #plot(data.sp)
    #plot(pt.i, pch = 19, add = T)
    #plot(data.j, col = "red", add = T)
    # get rid of any singletons
    regmod.all.P = regmod.all.P.orig
    tmp.test <- lapply(data.j@data, unique)
    tmp.test <- lapply(tmp.test, length)
    tmp.test <- as.vector(unlist(tmp.test))
    
    ## Get rid of any singularities 
    # droped from regmod 
    if (any(tmp.test == 1)) {
      index.jj <- which(tmp.test == 1)
      tit <- names(data.j)[index.jj]
      terms.j <- terms[2:length(terms)]
      index.jj <- match(tit, terms.j)
      terms.j <- terms.j[-index.jj]
      ltj <- length(terms.j)     
      regmod <- paste(terms[1], "~")
      for ( k in c(1:ltj) ) {
        if ( k != ltj ) regmod <- paste(regmod, terms.j[k],  "+")
        if ( k == ltj ) regmod <- paste(regmod, terms.j[k])
      }
      regmod.all.P<- regmod
    } 
    
    ## GWeight the data 
     
    dists.j <- gw.dist(coordinates(pt.i), coordinates(data.j))
    weights.j <- as.vector(gw.weight(dists.j,bwd.range[j], 
    		kernel = "bisquare", adaptive = F))
    #plot(dists.i, weights.i)
    #plot(data.j, pch = 19, cex = weights.j)
    data.j@data <- (data.j@data) * as.matrix(weights.j, ncol = 1)
 
    if ( class(tryCatch(stepAIC(lm(regmod.all.P, data.j@data), 
    		trace = F), error = function(e) e))[1] != "simpleError" ) {
    		step.j <- stepAIC(lm(regmod.all.P, data.j@data), trace = F) 
    		regmod.j <- as.formula(step.j)
    		tmp <- lm(regmod.j, data.j@data)
		aic.score.pd[i,j] <-	AIC(tmp)
    		model.mat.pd[i,j] <- as.vector(as.character(regmod.j)[3]) } else {
		  
		step.j <- regmod.all.P
		regmod.j <- as.formula(step.j)
  		tmp <- lm(regmod.j, data.j@data)
		aic.score.pd[i,j] <-	AIC(tmp)
		model.mat.pd[i,j] <- as.vector(as.character(regmod.all.P)[3]) }
  }
  #Sys.time() - st
  if(i %% 5 == 0) cat(i)
}

colnames(aic.score.nd) <- bwd.range
rownames(aic.score.nd) <- rownames(data.sp@data)
colnames(aic.score.pd) <- bwd.range
rownames(aic.score.pd) <- rownames(data.sp@data)

colnames(model.mat.nd) <- bwd.range
rownames(model.mat.nd) <- rownames(data.sp@data)
colnames(model.mat.pd) <- bwd.range
rownames(model.mat.pd) <- rownames(data.sp@data)

save(list = c("model.mat.pd","model.mat.nd"), file = "fixed.model.mat_v4.RData")
save(list = c("aic.score.pd","aic.score.nd"), file = "fixed.aic.score_v4.RData")

##### END Hyper local GWR Regression
# load("fixed.model.mat_v4.RData")
# load("fixed.aic.score_v4.RData")

## evaluate - Figure 4
# helper function
my.which.min <- function(x) {
  y <- min(x, na.rm = T)
  x <- match(y, x)
  x
}
### STN
index <- is.infinite(aic.score.nd)
aic.score.nd[index] <- NA
bw.vals <- apply((aic.score.nd), 1, my.which.min)
bw.vals <- as.numeric(colnames(aic.score.nd)[bw.vals])
tab5 <-(data.frame(table(bw.vals)))
names(tab5)[1] <- "Bandwidth (m)"
#write.csv(tab5, file = "Tab5.csv")
bw.vals.N <- bw.vals
gwr.sp <- data.sp
X <- coordinates(gwr.sp)[,1]
Y <- coordinates(gwr.sp)[,2]
gwr.sp <- spCbind(gwr.sp, data.frame(X, Y, bw = bw.vals))
Bandwidth = bw.vals.N
setEPS()
postscript("F4a.eps", width = 6, height = 8) 
#png(filename = "F4a_rev.png", w = 6, h = 8, units = "in", res = 600)
ggplot(boundary.df) + 
  geom_polygon(aes(x=long, y=lat), colour=NA, fill="lightgrey") +
  coord_equal() +
  geom_point(data = gwr.sp@data, aes(x = X, y = Y, size = Bandwidth), 
  	show.legend = T) +
  labs(title = "STN") +
  theme(axis.title.x=element_blank(),
      	axis.text.x=element_blank(),
    	axis.ticks.x=element_blank(), 
    	axis.title.y=element_blank(),
    	axis.text.y=element_blank(),
    	axis.ticks.y=element_blank(),
    	panel.background = element_blank(), 
    	legend.direction = "horizontal", 
    	legend.position = "bottom", 
    	legend.key.width = unit(0.9, "cm"),    	
    	plot.title = element_text(size=22),
    	legend.title=element_text(size=14))  
dev.off()
### STP
index <- is.infinite(aic.score.pd)
aic.score.pd[index] <- NA
bw.vals <- apply((aic.score.pd), 1, my.which.min)
bw.vals <- as.numeric(colnames(aic.score.pd)[bw.vals])
tab6 <-(data.frame(table(bw.vals)))
names(tab6)[1] <- "Bandwidth (m)"
#write.csv(tab6, file = "Tab6.csv")
bw.vals.P <- bw.vals
gwr.sp <- data.sp
X <- coordinates(gwr.sp)[,1]
Y <- coordinates(gwr.sp)[,2]
gwr.sp <- spCbind(gwr.sp, data.frame(X, Y, bw = bw.vals))
Bandwidth = bw.vals.P
setEPS()
postscript("F4b.eps", width = 6, height = 8) 
#png(filename = "F4b_rev.png", w = 6, h = 8, units = "in", res = 600)
ggplot(boundary.df) + 
  geom_polygon(aes(x=long, y=lat), colour=NA, fill="lightgrey") +
  coord_equal() +
  geom_point(data = gwr.sp@data, aes(x = X, y = Y, size = Bandwidth), 
  	show.legend = T) +
  labs(title = "STP") +
      theme(axis.title.x=element_blank(),
      	axis.text.x=element_blank(),
    	axis.ticks.x=element_blank(), 
    	axis.title.y=element_blank(),
    	axis.text.y=element_blank(),
    	axis.ticks.y=element_blank(),
    	panel.background = element_blank(), 
    	legend.direction = "horizontal", 
    	legend.position = "bottom", 
    	legend.key.width = unit(0.9, "cm"), 
    	plot.title = element_text(size=22),
    	legend.title=element_text(size=14))
dev.off()

### Local variable selection - Table 5
### STN
bw.vals <- apply((aic.score.nd), 1, my.which.min)
bw.vals <- as.numeric(colnames(aic.score.nd)[bw.vals])

index <- match(bw.vals, colnames(model.mat.nd))
best.mods <- vector()
for (i in 1:689) {
  best.mods <- append(best.mods, model.mat.nd[i, index[i]])
}
#length(terms)
plot.mat.n <- matrix(ncol = (length(terms)-2), nrow = nrow(data.sp), data = NA)
terms.n <- terms[3:8]
colnames(plot.mat.n) <- terms.n
for (i in 1:length(terms.n)) {
  grep.i <- grep(terms.n[i], best.mods)
  plot.mat.n[grep.i,i] <- 1
}
STN <- colSums(plot.mat.n, na.rm = T)
#index.n <- (which(STN < 689 & STN > 0))
index.n <- which(STN > 0)

### STP
bw.vals <- apply((aic.score.pd), 1, my.which.min)
bw.vals <- as.numeric(colnames(aic.score.pd)[bw.vals])

index <- match(bw.vals, colnames(model.mat.pd))
best.mods <- vector()
for (i in 1:689) {
  best.mods <- append(best.mods, model.mat.pd[i, index[i]])
}
plot.mat.p <- matrix(ncol = (length(terms)-2), nrow = nrow(data.sp), data = NA)
terms.p <- terms[3:8]
colnames(plot.mat.p) <- terms.p
for (i in 1:length(terms.n)) {
  grep.i <- grep(terms.n[i], best.mods)
  plot.mat.p[grep.i,i] <- 1
}
STP <- colSums(plot.mat.p, na.rm = T)
#index.p <- (which(STP < 689 & STP > 0))
index.p <- which(STP > 0)

tab5 <- cbind(STN, STP)
write.csv(tab5, file = "Tab5.csv")

### Figures 5 and 6 
# Maps with covariate selection and t-values
# significance of covariates
dMat <- dist2(coordinates(data.sp), coordinates(data.sp))
#X <- as.matrix(cbind(1,data.sp@data[2:10]))
BKWcn <- function(X) {
  p <- dim(X)[2]
  Xscale <- sweep(X, 2, sqrt(colSums(X^2)), "/")
  Xsvd <- svd(Xscale)$d
  Xsvd[1] / Xsvd[p]
}
#BKWcn(X)

### Determine the t-values and CN from the models 
# N
bw.vals <- apply((aic.score.nd), 1, my.which.min)
bw.vals <- as.numeric(colnames(aic.score.nd)[bw.vals])
index <- match(bw.vals, colnames(model.mat.nd))
best.mods <- vector()
for (i in 1:689) {
  best.mods <- append(best.mods, model.mat.nd[i, index[i]])
}
t.mat.n <- matrix(data = NA, nrow = 698, ncol = 7)
colnames(t.mat.n) <- c(names(data.sp@data[6:11]), "CN")
for (i in 1:689) {
	
	regmod.i <- (best.mods[i])
	regmod.i <- gsub("[:+;]", "", regmod.i)
	regmod.i <- gsub("  ", " ", regmod.i)
	regmod.i <- gsub(" ", "+", regmod.i)	
	regmod.i <- as.formula(paste0("TNPC~", regmod.i))
	bw.i <- bw.vals[i]
	pt.i <- data.sp[i,]
  	dMat.i <- dMat[i,]
  
    index.i <- dMat.i < bw.i
    data.i <- data.sp[index.i, c(4,6:11)]
    ## GWeight the data 
    dists.i <- gw.dist(coordinates(pt.i), coordinates(data.i))
    weights.i <- as.vector(gw.weight(dists.i,bw.i, kernel = "bisquare", adaptive = F))
    data.i@data <- (data.i@data) * as.matrix(weights.i, ncol = 1)
    
    tmp <- lm(regmod.i, data.i@data)
    tval.i <- summary(tmp)$coefficients[-1,3]
    index.i <- match(names(tval.i),colnames(t.mat.n))
    t.mat.n[i,index.i] <- tval.i
    
    X <- as.matrix(cbind(1,data.i@data[,names(tval.i)]))
    t.mat.n[i, 7] <- BKWcn(X) 
    
}
# check
n = 114
model.mat.nd[n, index[n]]
t.mat.n[n,]
best.mods[n]

# P
bw.vals <- apply((aic.score.pd), 1, my.which.min)
bw.vals <- as.numeric(colnames(aic.score.pd)[bw.vals])
index <- match(bw.vals, colnames(model.mat.pd))
best.mods <- vector()
for (i in 1:689) {
  best.mods <- append(best.mods, model.mat.pd[i, index[i]])
}

t.mat.p <- matrix(data = NA, nrow = 698, ncol = 7)
colnames(t.mat.p) <- c(names(data.sp@data[6:11]), "CN")
for (i in 1:689) {
	
	regmod.i <- (best.mods[i])
	regmod.i <- gsub("[:+;]", "", regmod.i)
	regmod.i <- gsub("  ", " ", regmod.i)
	regmod.i <- gsub(" ", "+", regmod.i)	
	regmod.i <- as.formula(paste0("TPPC~", regmod.i))
	bw.i <- bw.vals[i]
	pt.i <- data.sp[i,]
  	dMat.i <- dMat[i,]
  
    index.i <- dMat.i < bw.i
    data.i <- data.sp[index.i, c(5,6:11)]
    ## GWeight the data 
    dists.i <- gw.dist(coordinates(pt.i), coordinates(data.i))
    weights.i <- as.vector(gw.weight(dists.i,bw.i, kernel = "bisquare", adaptive = F))
    data.i@data <- (data.i@data) * as.matrix(weights.i, ncol = 1)
    
    tmp <- lm(regmod.i, data.i@data)
    tval.i <- summary(tmp)$coefficients[-1,3]
    index.i <- match(names(tval.i),colnames(t.mat.p))
    t.mat.p[i,index.i] <- tval.i
  
    X <- as.matrix(cbind(1,data.i@data[,names(tval.i)]))
    t.mat.p[i, 7] <- BKWcn(X)  
}
# check
n = 10
model.mat.pd[n, index[n]]
t.mat.p[n,]
best.mods[n]

### Now set up the plot
ggplot.func <- function(data.i, tit, index, sh = "GnBu" ) {
  gwr.sp <- data.i
  X <- coordinates(gwr.sp)[,1]
  Y <- coordinates(gwr.sp)[,2]
  gwr.sp <- spCbind(gwr.sp, data.frame(X, Y))
  gwr.sp.index <- gwr.sp[index,] 

  vals <-  gwr.sp@data[,1]
  q <- quantile(vals, c(0.25, 0.75))
  vals[vals < q[1]] <- q[1]
  vals[vals > q[2]] <- q[2]

  p <- ggplot(boundary.df) +       
      geom_polygon(aes(x=long, y=lat), colour="grey", fill="grey") +
      coord_equal() +
      geom_point(data = gwr.sp@data, aes(x = X, y = Y, colour= vals), size = 1) +
      scale_colour_distiller(type="seq", direction = 1, 
      	palette = sh) +
      geom_point(data = gwr.sp.index@data, aes(x = X, y = Y), size = 0.3) +
      labs(subtitle = tit) +
      theme(axis.title.x=element_blank(),
      	axis.text.x=element_blank(),
    	axis.ticks.x=element_blank(), 
    	axis.title.y=element_blank(),
    	axis.text.y=element_blank(),
    	axis.ticks.y=element_blank(),
    	panel.background = element_blank(), 
    	legend.direction = "horizontal", 
    	legend.position = "bottom", 
    	legend.key.width = unit(0.9, "cm"),
    	legend.title=element_blank())

  return(p)
}

### Figure 5 STN 
sh.list <- c("YlOrBr", "YlGnBu","PuBuGn", "YlOrRd", "OrRd", "BuGn") 
hyperlocal.tab.n <- matrix(ncol = 6, nrow = 6)
rownames(hyperlocal.tab.n) <- names(index.n)
colnames(hyperlocal.tab.n) <- colnames(tab3)[1:6]
for (i in 1:length(index.n)) {
  col.i <- plot.mat.n[,index.n[i]]
  row.i <- which(!is.na(col.i))
  data.i <- data.sp[row.i,]
  tit <- names(index.n)[i]
  t.mat.i <- t.mat.n[row.i,index.n[i]]
  index.pos <- which(t.mat.i > 1.96)
  index.neg <- which(t.mat.i < -1.96)
  index <- append(index.pos, index.neg)
  data.i <- SpatialPointsDataFrame(data.i, data.frame(t.mat.i))
  p.i <- ggplot.func(data.i, tit, index, sh = sh.list[i])
  tit <- paste("p",i, sep="")
  assign(tit, p.i)
  hyperlocal.tab.n[i, ] <- summary(t.mat.i)
}
setEPS()
postscript("F5.eps", width = 7, height = 7) 
#png(filename = "F5_new.png", w = 7, h = 6, units = "in", res = 300)
par(mfrow = c(2,3))
par(mar = c(0,0,0,0))
multiplot2(list(p1,p2,p3,p4,p5,p6), cols = 3)
dev.off()

### Figure 6 STP 
hyperlocal.tab.p <- matrix(ncol = 6, nrow = 6)
rownames(hyperlocal.tab.p) <- names(index.p)
colnames(hyperlocal.tab.p) <- colnames(tab3)[1:6]

for (i in 1:length(index.p)) {
  col.i <- plot.mat.p[,index.p[i]]
  row.i <- which(!is.na(col.i))
  data.i <- data.sp[row.i,]
  tit <- names(index.p)[i]
  t.mat.i <- t.mat.p[row.i,index.p[i]]
  index.pos <- which(t.mat.i > 1.96)
  index.neg <- which(t.mat.i < -1.96)
  index <- append(index.pos, index.neg)
  data.i <- SpatialPointsDataFrame(data.i, data.frame(t.mat.i))
  p.i <- ggplot.func(data.i, tit, index, sh = sh.list[i])
  tit <- paste("pp",i, sep="")
  assign(tit, p.i)
  hyperlocal.tab.p[i, ] <- summary(t.mat.i)
}
setEPS()
postscript("F6.eps", width = 7, height = 7) 
par(mfrow = c(2,3))
par(mar = c(0,0,0,0))
#png(filename = "F6_new.png", w = 7, h = 6, units = "in", res = 300)
multiplot2(list(pp1,pp2,pp3,pp4,pp5,pp6), cols = 3)
dev.off()

## Tables 6 and 7
write.csv(round(hyperlocal.tab.n, 3), file = "Tab6.csv")
write.csv(round(hyperlocal.tab.p, 3), file = "Tab7.csv")


## Figure 7 - sctterplots of fit
### STN
bw.vals <- apply((aic.score.nd), 1, my.which.min)
bw.vals <- as.numeric(colnames(aic.score.nd)[bw.vals])

index <- match(bw.vals, colnames(model.mat.nd))
best.mods <- vector()
for (i in 1:689) {
  best.mods <- append(best.mods, model.mat.nd[i, index[i]])
}
#length(terms)
plot.mat.n <- matrix(ncol = (length(terms)-2), nrow = nrow(data.sp), data = NA)
terms.n <- terms[3:8]
colnames(plot.mat.n) <- terms.n
for (i in 1:length(terms.n)) {
  grep.i <- grep(terms.n[i], best.mods)
  plot.mat.n[grep.i,i] <- 1
}
### use bw.vals to set bw
### use plot.mat.n to get the covariates
dMat <- dist2(coordinates(data.sp), coordinates(data.sp))
pred.obs.n <- matrix(ncol = 2, nrow = nrow(data.sp))
for ( i in 1: nrow(data.sp)) {
	pt.i <- data.sp[i,]
	bw <- bw.vals[i]
	dMat.i <- dMat[i,]
	terms.i <- plot.mat.n[i,]
	terms.i <- colnames(plot.mat.n)[which(terms.i == 1)]
	regmod.i <- paste(terms[1], "~")
	for ( j in 1:length(terms.i) ) {
		if ( j != length(terms.i) ) regmod.i <- paste(regmod.i, terms.i[j],  "+")
		if ( j == length(terms.i) ) regmod.i <- paste(regmod.i, terms.i[j])
	}
	regmod.i <- as.formula(regmod.i)
	terms.i <- append(terms[1], terms.i)
    ## sort out the clostest
    index.i <- which(dMat.i < bw)
    data.i <- data.sp[index.i,terms.i]

  	## Get rid of any singularities from regmod 
    # regmod.all = regmod.i
    tmp.test <- apply(data.i@data, 2, unique)
    
    if (is.matrix(tmp.test)) tmp.test <- apply(tmp.test,2, length)
    if (is.list(tmp.test)) tmp.test <- lapply(tmp.test,length)
    
    tmp.test <- tmp.test[-1]
    #tmp.test <- (unlist(tmp.test))   
    #tmp.test[2] = 1
    #regmod.i
    if (any(tmp.test == 1)) {
      index.jj <- names(tmp.test)[which(tmp.test == 1)]
      terms.j <- terms.i
      index.jj <- match(index.jj, terms.j)
      terms.j <- terms.j[-index.jj]
      ltj <- length(terms.j)     
      regmod.i <- paste(terms.j[1], "~")
      for ( k in 2:ltj ) {
        if ( k != ltj ) regmod.i <- paste(regmod.i, terms.j[k],  "+")
        if ( k == ltj ) regmod.i <- paste(regmod.i, terms.j[k])
      }
    } 
    #regmod.i
    ## GWeight the data 
    dists.i <- gw.dist(coordinates(pt.i), coordinates(data.i))
    weights.i <- gw.weight(dists.i,bw, kernel = "bisquare", adaptive = F)
    #plot(dists.i, weights.i)
    #plot(data.j, pch = 19, cex = weights.j)
    data.i@data <- (data.i@data) * as.matrix(weights.i, ncol = 1)
  	# do GW regression
    regmod.i <- as.formula(regmod.i)
	tmp <- lm(regmod.i, data.i@data)
	# extract the pred & obs	
	index.i <- match(rownames(pt.i@data), rownames(data.i@data))
	pred.obs.pair.i <- c(tmp$fitted.values[index.i], data.i@data[index.i,1])
	pred.obs.n[i,] <- pred.obs.pair.i
	if (i %% 100 == 0) cat(i, "\n")
}

### STP
bw.vals <- apply((aic.score.pd), 1, my.which.min)
bw.vals <- as.numeric(colnames(aic.score.pd)[bw.vals])

index <- match(bw.vals, colnames(model.mat.pd))
best.mods <- vector()
for (i in 1:689) {
  best.mods <- append(best.mods, model.mat.pd[i, index[i]])
}
#length(terms)
plot.mat.p <- matrix(ncol = (length(terms)-2), nrow = nrow(data.sp), data = NA)
terms.p <- terms[3:8]
colnames(plot.mat.p) <- terms.p
for (i in 1:length(terms.n)) {
  grep.i <- grep(terms.n[i], best.mods)
  plot.mat.p[grep.i,i] <- 1
}

### use bw.vals to set bw
### use plot.mat.n to get the covariates
dMat <- dist2(coordinates(data.sp), coordinates(data.sp))
pred.obs.p <- matrix(ncol = 2, nrow = nrow(data.sp))
for ( i in 1: nrow(data.sp)) {
	pt.i <- data.sp[i,]
	bw <- bw.vals[i]
	dMat.i <- dMat[i,]
	terms.i <- plot.mat.p[i,]
	terms.i <- colnames(plot.mat.p)[which(terms.i == 1)]
	regmod.i <- paste(terms[2], "~")
	for ( j in 1:length(terms.i) ) {
		if ( j != length(terms.i) ) regmod.i <- paste(regmod.i, terms.i[j],  "+")
		if ( j == length(terms.i) ) regmod.i <- paste(regmod.i, terms.i[j])
	}
	regmod.i <- as.formula(regmod.i)
	terms.i <- append(terms[2], terms.i)
    ## sort out the clostest
    index.i <- which(dMat.i < bw)
    data.i <- data.sp[index.i,terms.i]

  	## Get rid of any singularities from regmod 
    # regmod.all = regmod.i
    tmp.test <- apply(data.i@data, 2, unique)
    
    if (is.matrix(tmp.test)) tmp.test <- apply(tmp.test,2, length)
    if (is.list(tmp.test)) tmp.test <- lapply(tmp.test,length)
    
    tmp.test <- tmp.test[-1]
    #tmp.test <- (unlist(tmp.test))   
    #tmp.test[2] = 1
    #regmod.i
    if (any(tmp.test == 1)) {
      index.jj <- names(tmp.test)[which(tmp.test == 1)]
      terms.j <- terms.i
      index.jj <- match(index.jj, terms.j)
      terms.j <- terms.j[-index.jj]
      ltj <- length(terms.j)     
      regmod.i <- paste(terms.j[1], "~")
      for ( k in 2:ltj ) {
        if ( k != ltj ) regmod.i <- paste(regmod.i, terms.j[k],  "+")
        if ( k == ltj ) regmod.i <- paste(regmod.i, terms.j[k])
      }
    } 
    #regmod.i
    ## GWeight the data 
    dists.i <- gw.dist(coordinates(pt.i), coordinates(data.i))
    weights.i <- gw.weight(dists.i,bw, kernel = "bisquare", adaptive = F)
    #plot(dists.i, weights.i)
    #plot(data.j, pch = 19, cex = weights.j)
    data.i@data <- (data.i@data) * as.matrix(weights.i, ncol = 1)
  	# do GW regression
    regmod.i <- as.formula(regmod.i)
	tmp <- lm(regmod.i, data.i@data)
	# extract the pred & obs	
	index.i <- match(rownames(pt.i@data), rownames(data.i@data))
	pred.obs.pair.i <- c(tmp$fitted.values[index.i], data.i@data[index.i,1])
	pred.obs.p[i,] <- pred.obs.pair.i
	if (i %% 100 == 0) cat(i, "\n")
}	
colnames(pred.obs.p) <- c("Fitted", "Observed")
colnames(pred.obs.n) <- c("Fitted", "Observed")
## Other Obs v Pred
## STN
pred.obs.n.lm <- cbind(mN$fitted, data.sp$TNPC)
colnames(pred.obs.n.lm) <- c("Fitted", "Observed")
pred.obs.n.gwr <- matrix(ncol = 2, nrow = nrow(data.sp))
for (i in 1:nrow(data.sp)){
	x.int.i <- gwr.n$SDF@data[i,1]
	x.B.i <- gwr.n$SDF@data[i,2:7]
	x.data.i <- data.sp@data[i,c(6:11)]
	y.obs <- data.sp@data[i,"TNPC"]
	y.pred <- sum(x.B.i * x.data.i) + x.int.i
	pred.obs.n.gwr[i, ] <- c(y.pred, y.obs)
}
colnames(pred.obs.n.gwr) <- c("Fitted", "Observed")
## STP
pred.obs.p.lm <- cbind(mP$fitted, data.sp$TPPC)
colnames(pred.obs.p.lm) <- c("Fitted", "Observed")
pred.obs.p.gwr <- matrix(ncol = 2, nrow = nrow(data.sp))
for (i in 1:nrow(data.sp)){
	x.int.i <- gwr.p$SDF@data[i,1]
	x.B.i <- gwr.p$SDF@data[i,2:7]
	x.data.i <- data.sp@data[i,c(6:11)]
	y.obs <- data.sp@data[i,"TPPC"]
	y.pred <- sum(x.B.i * x.data.i) + x.int.i
	pred.obs.p.gwr[i, ] <- c(y.pred, y.obs)
}
colnames(pred.obs.p.gwr) <- c("Fitted", "Observed")

lm_eqn <- function(df){
    y <- df[,1]
    x <- df[,2]
    m <- lm(y ~ x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
         list(a = format(coef(m)[1], digits = 3), 
              b = format(coef(m)[2], digits = 3), 
             r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));                 
}

plot.pred.obs <- function(pred.obs.n, tit = "Hyper Local", top = T) {
	tmp <- data.frame(pred.obs.n)
	y1 <- min(pred.obs.n[,1])
	y2 <- max(pred.obs.n[,1])
	if (top == T) y3 <- max(pred.obs.n[,1])
	if (top == F) y3 <- min(pred.obs.n[,1])
	x1 <- (max(pred.obs.n[,2]) - ((max(pred.obs.n[,2]) - min(pred.obs.n[,2]))/2))*1.1
	textsize = 3
	textsize2 = 15
  # additional for EPS
	y <- tmp[,1]
	x <- tmp[,2]
	m <- lm(y ~ x, tmp)
	r2 = format(summary(m)$r.squared, digits = 3)
	
	r2 <- substitute(tit~~italic(r)^2~"="~r2,  
	                 list(tit = tit,
	                   r2 = format(summary(m)$r.squared, digits = 3)))
	p <- ggplot(data = tmp, aes(x = Observed, y = Fitted)) +
		geom_point() +
		scale_x_continuous(name = "Observed") + 
		scale_y_continuous(name = "Fitted", limits = c(y1, y2)) +
		# geom_text(x = x1, y = y3, label = lm_eqn(tmp[, c(2, 1)]), 
		  #	parse = TRUE, size = textsize ) +
		geom_smooth(method = "lm", se = T) +
		theme_minimal() +
		theme_bw() +
		labs(subtitle = r2) +
		theme(plot.title = element_text(hjust = 0.5, size = textsize2))
	return(p)
}
p1 <- plot.pred.obs(pred.obs.n, tit = "STN Hyper Local GWR")
p2 <- plot.pred.obs(pred.obs.n.gwr, tit = "STN GWR")
p3 <- plot.pred.obs(pred.obs.n.lm, tit = "STN Linear Regression")

p4 <- plot.pred.obs(pred.obs.p, tit = "STP Hyper Local GWR", top = F)
p5 <- plot.pred.obs(pred.obs.p.gwr, tit = "STP GWR", top = F)
p6 <- plot.pred.obs(pred.obs.p.lm, tit = "STP Linear Regression", top = F)


summary(lm(Observed~Fitted, data = data.frame(pred.obs.p.lm)))
summary(lm(Observed~Fitted, data = data.frame(pred.obs.p.gwr)))
summary(lm(Observed~Fitted, data = data.frame(pred.obs.p)))

summary(lm(Observed~Fitted, data = data.frame(pred.obs.n.lm)))
summary(lm(Observed~Fitted, data = data.frame(pred.obs.n.gwr)))
summary(lm(Observed~Fitted, data = data.frame(pred.obs.n)))

setEPS()
postscript("F7a.eps", width = 4, height = 8) 
#png(filename = "F7a_new.png", w = 4, h = 8, units = "in", res = 300)
multiplot2(list(p3,p2,p1),cols =1)
dev.off()
setEPS()
postscript("F7b.eps", width = 4, height = 8) 
#png(filename = "F7b_new.png", w = 4, h = 8, units = "in", res = 300)
multiplot2(list(p6,p5,p4),cols =1)
dev.off()

#### Figure 8
### maps comparing GWR and hyperlocal AIC values
### this is addressing the question of which model GWR or the HL approach
### results in the best fitting (most parsimonious) model  
# N
bw.vals <- apply((aic.score.nd), 1, my.which.min)
bw.vals <- as.numeric(colnames(aic.score.nd)[bw.vals])
index <- match(bw.vals, colnames(model.mat.nd))
best.mods <- vector()
for (i in 1:689) {
  best.mods <- append(best.mods, model.mat.nd[i, index[i]])
}
R2.mat.n <- matrix(data = NA, nrow = 689, ncol = 8)
colnames(R2.mat.n) <- c("R2GWR", "R2HL", "R2aGWR", "R2aHL", "AICGWR", "AICHL", "AICcGWR", "AICcHL")
  
for (i in 1:689) {
	
	regmod.i <- (best.mods[i])
	regmod.i <- gsub("[:+;]", "", regmod.i)
	regmod.i <- gsub("  ", " ", regmod.i)
	regmod.i <- gsub(" ", "+", regmod.i)	
	regmod.i <- as.formula(paste0("TNPC~", regmod.i))
	bw.i <- bw.vals[i]
	
	pt.i <- data.sp[i,]
  	dMat.i <- dMat[i,]
  
    index.i <- dMat.i < bw.i
    data.i <- data.sp[index.i, c(4,6:11)]
    ## GWeight the data 
    dists.i <- gw.dist(coordinates(pt.i), coordinates(data.i))
    weights.i <- as.vector(gw.weight(dists.i,bw.i, kernel = "bisquare", adaptive = F))
    data.i@data <- (data.i@data) * as.matrix(weights.i, ncol = 1)    
    tmp <- lm(regmod.i, data.i@data)
    #AIC.mat.n[i, 2] <- AIC(tmp)
    R2.mat.n[i, 2] <- unlist(summary(tmp)[8]) 
    R2.mat.n[i, 4] <- unlist(summary(tmp)[9])
    R2.mat.n[i, 6] <- AIC(tmp) 
    R2.mat.n[i, 8] <- AICc(tmp) 
 
    bw.i <- bw.n
    regmod.i <- regmod.n
    index.i <- dMat.i < bw.i
    data.i <- data.sp[index.i, c(4,6:11)]
    ## GWeight the data 
    dists.i <- gw.dist(coordinates(pt.i), coordinates(data.i))
    weights.i <- as.vector(gw.weight(dists.i,bw.i, kernel = "bisquare", adaptive = F))
    data.i@data <- (data.i@data) * as.matrix(weights.i, ncol = 1)
    tmp <- lm(regmod.i, data.i@data)
    #AIC.mat.n[i, 1] <- AIC(tmp)
    R2.mat.n[i, 1] <- unlist(summary(tmp)[8]) 
    R2.mat.n[i, 3] <- unlist(summary(tmp)[9])
    R2.mat.n[i, 5] <- AIC(tmp) 
    R2.mat.n[i, 7] <- AICc(tmp) 
 
}
head(R2.mat.n)
tail(R2.mat.n)
# does the HL have better R2
sum(R2.mat.n[, 2] > R2.mat.n[, 1], na.rm = T)

# P
bw.vals <- apply((aic.score.pd), 1, my.which.min)
bw.vals <- as.numeric(colnames(aic.score.pd)[bw.vals])
index <- match(bw.vals, colnames(model.mat.pd))
best.mods <- vector()
for (i in 1:689) {
  best.mods <- append(best.mods, model.mat.pd[i, index[i]])
}
R2.mat.p <- matrix(data = NA, nrow = 689, ncol = 8)
colnames(R2.mat.p) <- c("R2GWR", "R2HL", "R2aGWR", "R2aHL", "AICGWR", "AICHL", "AICcGWR", "AICcHL")

for (i in 1:689) {
	
	regmod.i <- (best.mods[i])
	regmod.i <- gsub("[:+;]", "", regmod.i)
	regmod.i <- gsub("  ", " ", regmod.i)
	regmod.i <- gsub(" ", "+", regmod.i)	
	regmod.i <- as.formula(paste0("TPPC~", regmod.i))
	bw.i <- bw.vals[i]
	
	pt.i <- data.sp[i,]
  	dMat.i <- dMat[i,]
  
    index.i <- dMat.i < bw.i
    data.i <- data.sp[index.i, c(5,6:11)]
    ## GWeight the data 
    dists.i <- gw.dist(coordinates(pt.i), coordinates(data.i))
    weights.i <- as.vector(gw.weight(dists.i,bw.i, kernel = "bisquare", adaptive = F))
    data.i@data <- (data.i@data) * as.matrix(weights.i, ncol = 1)    
    tmp <- lm(regmod.i, data.i@data)
    #AIC.mat.p[i, 2] <- AIC(tmp)
    R2.mat.p[i, 2] <- unlist(summary(tmp)[8]) 
   	R2.mat.p[i, 4] <- unlist(summary(tmp)[9])
    R2.mat.p[i, 6] <- AIC(tmp) 
    R2.mat.p[i, 8] <- AICc(tmp) 
 
    bw.i <- bw.p
    regmod.i <- regmod.p
    index.i <- dMat.i < bw.i
    data.i <- data.sp[index.i, c(5,6:11)]
    ## GWeight the data 
    dists.i <- gw.dist(coordinates(pt.i), coordinates(data.i))
    weights.i <- as.vector(gw.weight(dists.i,bw.i, kernel = "bisquare", adaptive = F))
    data.i@data <- (data.i@data) * as.matrix(weights.i, ncol = 1)
    tmp <- lm(regmod.i, data.i@data)
    #AIC.mat.p[i, 1] <- AIC(tmp)
    R2.mat.p[i, 1] <- unlist(summary(tmp)[8]) 
    R2.mat.p[i, 3] <- unlist(summary(tmp)[9])
    R2.mat.p[i, 5] <- AIC(tmp) 
    R2.mat.p[i, 7] <- AICc(tmp) 

}
head(R2.mat.p)
tail(R2.mat.p)
# does the HL have better R2
sum(R2.mat.p[, 2] > R2.mat.p[, 1], na.rm = T)
dif <- sqrt((R2.mat.p[,1]-R2.mat.p[,2])^2)

#### Plots of where HL is greater than GWR
plot.winner.func <- function(winner, tit = "") {	
	p <- ggplot(boundary.df) + 
      geom_polygon(aes(x=long, y=lat), colour=NA, fill="grey90") +
      coord_equal() +
      geom_point(data = gwr.sp@data, aes(x = X, y = Y), size = winner*1.5, colour = winner) +
      labs(title = tit) +
      theme(axis.title.x=element_blank(),
      	axis.text.x=element_blank(),
    	axis.ticks.x=element_blank(), 
    	axis.title.y=element_blank(),
    	axis.text.y=element_blank(),
    	axis.ticks.y=element_blank(),
    	panel.background = element_blank(), 
    	plot.title = element_text(size=22))
    return(p)
}
# R2 differences

diff.plot.func <- function(winner, diff, tit = "") {
	p <- ggplot(boundary.df) + 
      geom_polygon(aes(x=long, y=lat), colour=NA, fill="grey90") +
      coord_equal() +
      geom_point(data = gwr.sp@data, aes(x = X, y = Y), size = diff*15, colour = winner) +
      labs(title = tit) +
      theme(axis.title.x=element_blank(),
      	axis.text.x=element_blank(),
    	axis.ticks.x=element_blank(), 
    	axis.title.y=element_blank(),
    	axis.text.y=element_blank(),
    	axis.ticks.y=element_blank(),
    	panel.background = element_blank(), 
    	plot.title = element_text(size=22))
    return(p)
}

winner <- (R2.mat.n[, 2] > R2.mat.n[, 1] ) + 1
diff <- abs((R2.mat.n[, 2] - R2.mat.n[, 1] ))
png(filename = "F8a_rev.png", w = 6, h = 8, units = "in", res = 600)
diff.plot.func(winner, diff, tit = "STN") 
dev.off()

winner <- (R2.mat.p[, 2] > R2.mat.p[, 1] ) + 1
diff <- abs((R2.mat.p[, 2] - R2.mat.p[, 1] ))
png(filename = "F8b_rev.png", w = 6, h = 8, units = "in", res = 600)
diff.plot.func(winner, diff, tit = "STP") 
dev.off()


winner <- (R2.mat.n[, 2] > R2.mat.n[, 1] ) + 1
Difference <- abs((R2.mat.n[, 2] - R2.mat.n[, 1] ))
Model = c("GWR", "Hyper-Local")
Model= Model[winner]
setEPS()
postscript("F8a_rev.eps", width = 6, height = 8) 
#png(filename = "F8a_rev.png", w = 6, h = 8, units = "in", res = 600)
ggplot(boundary.df) + 
      geom_polygon(aes(x=long, y=lat), colour=NA, fill="grey90") +
      coord_equal() +
      geom_point(data = gwr.sp@data, aes(x = X, y = Y, size = Difference, 
      	colour = Model), show.legend = T) +    
      scale_color_manual(breaks = c("GWR", "Hyper-Local"),
                        values=c("black", "red")) + 
      scale_size(range = c(1,5)) +  
      labs(title = "STN") +
      theme(axis.title.x=element_blank(),
      	axis.text.x=element_blank(),
    	axis.ticks.x=element_blank(), 
    	axis.title.y=element_blank(),
    	axis.text.y=element_blank(),
    	axis.ticks.y=element_blank(),
    	panel.background = element_blank(), 
    	legend.direction = "horizontal", 
    	legend.position = "bottom", 
    	legend.key.width = unit(0.9, "cm"), 
    	plot.title = element_text(size=14),
    	legend.title=element_text(size=14))
dev.off()

winner <- (R2.mat.p[, 2] > R2.mat.p[, 1] ) + 1
Difference <- abs((R2.mat.p[, 2] - R2.mat.p[, 1] ))
Model = c("GWR", "Hyper-Local")
Model= Model[winner]
setEPS()
postscript("F8b_rev.eps", width = 6, height = 8) 
#png(filename = "F8b_rev.png", w = 6, h = 8, units = "in", res = 600)
ggplot(boundary.df) + 
      geom_polygon(aes(x=long, y=lat), colour=NA, fill="grey90") +
      coord_equal() +
      scale_color_manual(breaks = c("GWR", "Hyper-Local"),
                        values=c("black", "red")) +  
      scale_size(range = c(1,5), breaks = c(0.005, 0.01, 0.02)) +  
      geom_point(data = gwr.sp@data, aes(x = X, y = Y, size = Difference, 
      	colour = Model), show.legend = T) +      
      labs(title = "STP") +
      theme(axis.title.x=element_blank(),
      	axis.text.x=element_blank(),
    	axis.ticks.x=element_blank(), 
    	axis.title.y=element_blank(),
    	axis.text.y=element_blank(),
    	axis.ticks.y=element_blank(),
    	panel.background = element_blank(), 
    	legend.direction = "horizontal", 
    	legend.position = "bottom", 
    	legend.key.width = unit(0.9, "cm"), 
    	plot.title = element_text(size=14),
    	legend.title=element_text(size=14))
dev.off()

##### END

