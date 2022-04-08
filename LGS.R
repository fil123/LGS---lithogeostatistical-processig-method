#install.packages("GISTools")
#install.packages("fossil")
#install.packages("fields")
#install.packages("SpatialTools")
#install.packages("agricolae")
library(GISTools)
library(rgdal)
library(sp)
library(maptools)
library(rgeos)
library(raster)
library(SpatialTools)
library(readxl)
library(plyr)
library(agricolae)

MatrixGeol<-function( rnmeasI,centroid,mx,radius){
  for (i in 1:nrow(rnmeasI)){
    m1<-matrix(c(rnmeasI$X[i],rnmeasI$Y[i]),ncol=2)
    #cambio sistema di coordinate
    coord1<-m1
    coord1<-as.data.frame(coord1)
    colnames(coord1)<-c("X","Y")
    coord1<-as.data.frame(coord1)
    coordinates(coord1) <- ~X+Y
    coord2<-coord1
    proj4string(coord2)<-CRS("+init=epsg:32632")
    coord2<-spTransform(coord2,"+init=epsg:32632")
    m1<-matrix(c(coord2$X,coord2$Y),ncol=2)
    m2<-matrix(c(centroid$x,centroid$y),ncol=2)
    if (dist2(m1,m2)<=radius){
      mx[i,] <-c(rnmeasI$X[i],rnmeasI$Y[i],log(rnmeasI$CONCENTRAZIONE[i]),rnmeasI$LIV_2[i])
    }
  }  
  colnames(mx)<-c("X","Y","CONC","LIV_2")  
  mx<-as.data.frame(mx)
  print(mx)
  my<-ddply(mx,~LIV_2,summarise,
            mean=mean( as.numeric(as.character(CONC)) ),
            sd=sd(as.numeric(as.character(CONC)) ))
  return(my)
}


matX<-function(m1,m2,m3,m4,mn,mx){
  radiusC<-function(mc,mi){
    for (i in 1:nrow(mi)) {
      #print(nrow(mi))
      litx<-mi[i,1]
      for (j in 1:nrow(mc)) {
        lit1<-mc[j,1]
        if(as.character(lit1)==as.character(litx)&(is.na(lit1)==FALSE&is.na(litx)==FALSE)){
          #print(lit1==litx)
          #print(as.character(lit1))
          #print(as.character(litx))
          mi[i,2]<-(mc[j,2])
          mi[i,3]<-(mc[j,3])
        } 
        
      }
      
    }
    return(mi)  
  }
  radiusC1<-function(mc,mi){
    for (i in 1:nrow(mi)) {
      litx<-mi[i,1]
      #print(is.na(mi[i,2]))
      #print(is.na(mi[i,3]))
      #readline()
      if(is.na(mi[i,2])|is.na(mi[i,3])){
        #readline()  
        for (j in 1:nrow(mc)) {
          lit1<-mc[j,1]
          if(as.character(lit1)==as.character(litx)&(is.na(lit1)==FALSE&is.na(litx)==FALSE)){
            #print(lit1==litx)
            #print(as.character(lit1))
            #print(as.character(litx))
            #readline()
            mi[i,2]<-(mc[j,2])
            mi[i,3]<-(mc[j,3])
            #print(mi[i,2])
            #print(mi[i,3])
          } 
          
        }  
      }  
    }
    
    return(mi)  
  }
  
  c1<-m1
  c2<-m2
  c3<-m3
  c4<-m4
  cn<-mn
  cx<-mx
  #print(mx)
  #---------------------------------------------CONTROLLO SU RAGGIO 1------------------------  
  cx<-radiusC(c1,cx)
  cx<-as.data.frame(cx)
  #---------------------------------------------CONTROLLO SU RAGGIO 2------------------------  
  cx<-radiusC(c2,cx)
  #---------------------------------------------CONTROLLO SU RAGGIO 3------------------------ 
  cx<-radiusC(c3,cx)
  #---------------------------------------------CONTROLLO SU RAGGIO 4------------------------  
  cx<-radiusC(c4,cx)
  #---------------------------------------------CONTROLLO SU TABELLA MEDIE GEOLOGIE------------------------ 
  cx<-radiusC(cn,cx)
  
  #print(cx)
  return(cx)
}

RNCONCLITO <- read_excel("/Users/utente/Desktop/Radon TUTTO/Radon Data Analisys/RNCONCLITO.xlsx")
RNGEO <- read_excel("/Users/utente/Desktop/Radon TUTTO/Radon Data Analisys/RN_GEO_STAT.xlsx")
rnmeas<-RNCONCLITO
head(rnmeas)
coordinates(rnmeas) <- ~X+Y
files <- dir(pattern = "*.shp")

m<-1
#length(files)
for (m in 1:length(files)){
  #0 aprire shp file dall'IGM-Geologia
  name<-gsub(".shp", "", files[m])
  CA <-readOGR(".",layer=name)
  #1 ricavo il centroide
  a<-gCentroid(CA)
  plot(CA)
  points(a)
  #matrice delle distanze tra il centride e i punti di misura
  m3<-matrix(data=NA,ncol=4,nrow=nrow(rnmeas),byrow = TRUE)
  #MatrixGeol<-function( rnmeasI,centroid,mx)
  m4<-MatrixGeol(rnmeas,a,m3,28000)
  m5<-MatrixGeol(rnmeas,a,m3,56000)
  m6<-MatrixGeol(rnmeas,a,m3,84000)
  m7<-MatrixGeol(rnmeas,a,m3,112000)
  #mn<-head(CA[,1],n=-1)
  mn<-matrix(CA$LIV_2,nrow = length(CA$LIV_2) ,ncol=1,byrow = TRUE)
  mn<-as.data.frame(mn)
  print(mn)
  #readline(prompt="Press [enter] to continue")
  mn<-cbind(mn,"mean"=NA,"sd"=NA)
  colnames(mn)[1]<-"LIV_2"
  print(mn)
  #readline(prompt="Press [enter] to continue")
  A<-matX(m4,m5,m6,m7,RNGEO,mn)
  print(A)
  merged <- merge(CA, A, by='LIV_2')
  a1<-toString(A[,1])
  a2<-toString(head(merged$LIV_2))
  print(a1)
  print(a2)
  #readline(prompt="Press [enter] to continue")
  filenameshp <- paste0("MERGED", name, ".shp")
  shapefile(merged, filenameshp)
  #readline(prompt="Press [enter] to continue")
  }
#UNISCI GLI SHAPE FILE ELABORATI

libs <- c("rgdal", "maptools", "gridExtra","sp")
lapply(libs, require, character.only = TRUE)

# Import IGM data
filesM <- dir(pattern = "MERGE*.shp|.shp$")
filesM<-grep(pattern="MERGE",x=filesM,value = TRUE)
length(filesM)

#length(filesM)
i<-1
for (i in 1:length(filesM)){
                    f<-shapefile(filesM[i])
                    # Create a generic raster, set the extent to the same as wetlands
                    r.raster <- raster()  
                    #Use extent (from raster package) to read bounds of vector and assign to the raster:
                    extent(r.raster) <- extent(f)
                    res(r.raster) <- 1000 # set cell size to 1000 metres
                    # Make a raster of the wetlands:
                    rnconcMiles <- rasterize(f, r.raster,"mean")
                    st<-gregexpr(pattern ='_',filesM[i])
                    en<-gregexpr(pattern ='.s',filesM[i])
                    nr<-substring(filesM[i],st,en)
                    filenameras <- paste0("RAS", nr, "tif")
                    plot(rnconcMiles)
                    x<-writeRaster(rnconcMiles,filenameras,overwrite=TRUE)

}

