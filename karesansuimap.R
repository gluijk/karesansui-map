# Karesansui (Japanese dry garden) style maps
# www.overfitting.net
# https://www.overfitting.net/2024/04/mapas-estilo-karesansui-con-r.html


library(terra)  # read GeoTIFF, reprojection, crop and resample
library(tiff)  # save 16-bit TIFF's
library(png)  # save 8-bit PNG's


hillshademap=function(DEM, dx=25, dlight=c(0, 2, 3), gamma=1) {
    # hillshademap() inputs DEM data and outputs a hillshade matrix
    #
    # DEM: digital elevation map matrix
    # dx: DEM resolution/cell size (same units as elevation values)
    # dlight: lighting direction (3D vector defined from observer to light source):
    #   dlight=c(0, 2, 3)  # sunrise
    #   dlight=c(0, 0, 1)  # midday
    #   dlight=c(0,-2, 3)  # sunset
    # gamma: optional output gamma lift
    
    DIMY=nrow(DEM)
    DIMX=ncol(DEM)
    # If array turn its first dimension into a genuine matrix
    if (!is.matrix(DEM)) DEM=matrix(DEM[,,1], nrow=DIMY, ncol=DIMX)
    
    dlightM=sum(dlight^2)^0.5
    
    # Vectorial product to calculate n (orthogonal vector)
    nx = 2*dx*(DEM[1:(DIMY-2), 2:(DIMX-1)] - DEM[3:DIMY,     2:(DIMX-1)])
    ny = 2*dx*(DEM[2:(DIMY-1), 1:(DIMX-2)] - DEM[2:(DIMY-1), 3:DIMX])
    nz = 4*dx^2
    nM = (nx^2 + ny^2 + nz^2)^0.5
    
    # Dot product to calculate cos(theta)
    dn = dlight[1]*nx + dlight[2]*ny + dlight[3]*nz  # (DIMY-2)x(DIMX-2) matrix
    
    # Reflectance (=cos(theta))
    hillshadepre=dn/(dlightM*nM)
    hillshadepre[hillshadepre<0]=0  # clip negative values
    
    # Add 1-pix 'lost' borders
    hillshademap=matrix(0, nrow=DIMY, ncol=DIMX)
    hillshademap[2:(DIMY-1), 2:(DIMX-1)]=hillshadepre
    rm(hillshadepre)
    hillshademap[c(1,DIMY),]=hillshademap[c(2,DIMY-1),]
    hillshademap[,c(1,DIMX)]=hillshademap[,c(2,DIMX-1)]
    
    return(hillshademap^(1/gamma))
}


# Blur
# https://stackoverflow.com/questions/70429190/how-can-i-perform-neighborhood-analysis-in-terra-or-raster-and-keep-the-same-na
arrayblur=function(img, radius=11) {
    # radius: radius of the circular averaging window
    require(terra)
    
    # Build circular kernel
    D=2*radius+1  # D will always be an odd number as required by focal()
    w=matrix(1, nrow=D, ncol=D)
    w[(row(w)-(radius+1))^2 + (col(w)-(radius+1))^2 >= (radius+1)^2]=NA
    writePNG(w, "blurkernel.png")
    
    raster=rast(img)  # process as raster
    rasterblur=focal(raster, w=w, fun='mean', na.rm=TRUE, na.policy='omit')
    
    if (is.matrix(img)) return (matrix(as.array(rasterblur), nrow=nrow(rasterblur)))
    else return (as.array(rasterblur))  # convert back to array
}


###########################################################

# 1. PROCESS GEOTIFF DATA TO GET THE DEM AS A MATRIX

baleares=rast("PNOA_MDT200_ETRS89_HU31_Baleares.tif")  # (985x1544 px)
baleares
plot(baleares)
RESOLUTION=res(baleares)[1]  # 200m grid resolution


# Convert to matrix and save as TIFF
DEM=matrix(as.array(baleares), nrow=nrow(baleares))
hist(DEM, breaks=1000)
DEM[is.nan(DEM)]=0
DEM[DEM<0]=0
HEIGHT=max(DEM)
DEM=DEM/max(DEM)
writeTIFF(DEM, "baleares.tif", compression='LZW', bits.per.sample=16)

DEM=readTIFF("baleares2.tif")  # added borders to Full HD (1080x1920 px)
DIMY=nrow(DEM)
DIMX=ncol(DEM)


###########################################################

# 2. PROCESS MATRIX TO OBTAIN MAP CONTOURS AND HILLSHADE

# Calculate solid map contour
solid=DEM
solid[solid>0]=1  # set areas to 1 (land)
writePNG(solid, "mapsolid.png")
         
         
# Calculate outline map from solid map
outline=solid*0
# 1 pixel thickness outline
outline[2:(DIMY-1), 2:(DIMX-1)]=
 abs(solid[1:(DIMY-2), 2:(DIMX-1)] -
     solid[2:(DIMY-1), 2:(DIMX-1)]) +
 abs(solid[2:(DIMY-1), 1:(DIMX-2)] -
     solid[2:(DIMY-1), 2:(DIMX-1)])
# increase to 2 pixel thickness outline
# outline[2:(DIMY-1), 2:(DIMX-1)]=outline[2:(DIMY-1), 2:(DIMX-1)]+
#     outline[1:(DIMY-2), 2:(DIMX-1)]+outline[2:(DIMY-1), 3:(DIMX-0)]
# increase to 3 pixel thickness outline
# outline[2:(DIMY-1), 2:(DIMX-1)]=outline[2:(DIMY-1), 2:(DIMX-1)]+
#     outline[1:(DIMY-2), 2:(DIMX-1)]+outline[3:(DIMY-0), 2:(DIMX-1)]+
#     outline[2:(DIMY-1), 1:(DIMX-2)]+outline[2:(DIMY-1), 3:(DIMX-0)]
outline[outline!=0]=1

writePNG(outline, "mapoutline.png")


# Calculate grayscale hillshade
MIX=0.7  # two light sources are mixed to fill dark areas a bit
hill=hillshademap(DEM, dx=RESOLUTION/HEIGHT, dlight=c(1, 2, 3))
hillfill=hillshademap(DEM, dx=RESOLUTION/HEIGHT, dlight=c(1, 3, 2))
hill=hill*MIX+hillfill*(1-MIX)
gamma=1
hill=(hill/max(hill))^(1/gamma)  # darken hillshade a bit

# Save hillshade
writeTIFF(hill, "hillshade.tif",
          bits.per.sample=16, compression="LZW")

# Display hillshade
image(t(hill[nrow(hill):1,]), useRaster=TRUE,
      col=c(gray.colors(256, start=0, end=1, gamma=2.2)),
      asp=nrow(hill)/ncol(hill), axes=FALSE)


###########################################################

# 3. OBTAIN KARESANSUI GROOVES (2 ALGORITHMS)

# Islands
GAPISLANDS=20  # islands grooves separation
NHORIZ=7
mode='accurate'  # mode=c('fast', 'accurate')  # Time difference of 4.127238 mins
mode='fast'  # mode=c('fast', 'accurate')  # Time difference of 13.39945 secs

if (mode=='accurate') {  # circles calculated from inner contour (R increases)
    horizon=outline*0
    indices=unname(which(outline==1, arr.ind=TRUE))  # contour pixels
    for (j in 1:NHORIZ) {
        print(paste0("Calculating groove ", j, "/", NHORIZ, "..."))
        R=j*GAPISLANDS  # radius in pixels from contour
        horizontmp=outline*0
        for (i in 1:nrow(indices)) {
            x0=indices[i,2]
            y0=indices[i,1]
            if (x0>R & y0>R & x0<DIMX-R & y0<DIMY-R)
                for (x in round(x0-R):round(x0+R)) {
                    for (y in round(y0-R):round(y0+R)) {
                        if ( ((x-x0)^2 + (y-y0)^2 )^0.5 < R) horizontmp[y,x]=1
                    }
                }
        }
        horizontmp[solid==1]=1  # fill inner area of horizontmp to 1
        horizon=horizon+horizontmp  # accumulate horizon ranges
    }
    rm(horizontmp)
} else {  # circles calculated from current contour (R=GAPISLANDS is fixed)
    outlinetmp=outline  # we already calculated solid and outline
    solidtmp=solid
    horizon=outlinetmp*0
    R=GAPISLANDS  # radius in pixels from current contour
    for (j in 1:NHORIZ) {
        print(paste0("Calculating groove ", j, "/", NHORIZ, "..."))
        
        # Calculate new contour
        if (j>1) {
            outlinetmp=solidtmp*0
            # 1 pixel thickness outline
            outlinetmp[2:(DIMY-1), 2:(DIMX-1)]=
                abs(solidtmp[1:(DIMY-2), 2:(DIMX-1)] -
                    solidtmp[2:(DIMY-1), 2:(DIMX-1)]) +
                abs(solidtmp[2:(DIMY-1), 1:(DIMX-2)] -
                    solidtmp[2:(DIMY-1), 2:(DIMX-1)])
            outlinetmp[outlinetmp!=0]=1
        }

        # Loop contour pixels
        indices=unname(which(outlinetmp==1, arr.ind=TRUE))  # current contour pixels
        horizontmp=outlinetmp*0
        for (i in 1:nrow(indices)) {
            x0=indices[i,2]
            y0=indices[i,1]
            if (x0>R & y0>R & x0<DIMX-R & y0<DIMY-R)
                for (x in round(x0-R):round(x0+R)) {
                    for (y in round(y0-R):round(y0+R)) {
                        if ( ((x-x0)^2 + (y-y0)^2 )^0.5 < R) horizontmp[y,x]=1
                    }
                }
        }
        horizontmp[solidtmp==1]=1  # fill inner area of horizontmp to 1
        horizon=horizon+horizontmp  # accumulate horizon ranges
        solidtmp=horizontmp
    }
    rm(outlinetmp, solidtmp, horizontmp)
}

writePNG(horizon/max(horizon), "horizon.png")

# Now paint in white all bands of the same odds/evens as band max(horizon) (land)
horizonbw=horizon*0
if (NHORIZ%%2) whites=which(horizon%%2==1) else whites=which(horizon%%2==0)
horizonbw[whites]=1
writePNG(horizonbw, "horizonbw.png")


# Sea
GAPSEA=14  # sea grooves separation

backgnd=horizon*0
backgnd[row(backgnd)%%(GAPSEA*2) < GAPSEA]=1
writePNG(backgnd, "backgnd.png")


#######################################################

# 4. BLUR OBTAINED GROOVES AND GENERATE A HILLSHADE

# It is better to calculate separate hillshades from each DEM and then mix them
# than mixing the input DEM and calculate a global hillshade because it can
# easily produce discontinuities

# Calculate grayscale hillshade of islands
# Blur radius should not be < GAPISLANDS/2
horizonblur=arrayblur(horizonbw, radius=GAPISLANDS/2)

MIX=0.7  # two light sources are mixed to fill dark areas a bit
hill=hillshademap(horizonblur, dx=0.1, dlight=c(1, 2, 3))
hillfill=hillshademap(horizonblur, dx=0.1, dlight=c(1, 3, 2))
hill=hill*MIX+hillfill*(1-MIX)
gamma=1/4
hillhorizon=(hill/max(hill))^(1/gamma)  # darken hillshade a bit


# Calculate grayscale hillshade of sea
# Blur radius should not be < GAPSEA/2
backgndblur=arrayblur(backgnd, radius=GAPSEA/2)

MIX=0.7  # two light sources are mixed to fill dark areas a bit
hill=hillshademap(backgndblur, dx=0.1, dlight=c(1, 2, 3))
hillfill=hillshademap(backgndblur, dx=0.1, dlight=c(1, 3, 2))
hill=hill*MIX+hillfill*(1-MIX)
gamma=1/4
hillsea=(hill/max(hill))^(1/gamma)  # darken hillshade a bit


# Mix and save hillshade
hillhorizon[horizon==0]=hillsea[horizon==0]
writeTIFF(hillhorizon, "hillshadeshogunNEW.tif",
          bits.per.sample=16, compression="LZW")
