## setup R for the analyses in this book

packs <- c("aTSA", "boot", "brms", "changepoint", "cowplot", "DHARMa", "dplR", "dplyr", "effects", "effectsize","faraway", "future.apply", "ggbiplot", "ggfortify", "ggplot2", "gridExtra", "Hmisc", "Kendall", "latex2exp", "lattice", "lme4", "lmerTest", "maps", "MASS", "mgcv", "mice", "mvtnorm", "nortest", "parameters", "purrr", "pwr", "randomForest", "ranger", "raster", "readxl", "reshape2", "retrodesign", "rgbif", "simstudy", "sp", "TidyDensity", "tidyr", "vegan")

istherealready <- installed.packages()
toinstall <- setdiff(packs, rownames(istherealready))
install.packages(toinstall, depend=T)

istherealready[rownames(istherealready) %in% packs, c(3), drop=F]


#### For information: setting during simulations ####
R version 4.3.1 (2023-06-16)

             Version     
aTSA         "3.1.2"     
boot         "1.3-28.1"  
brms         "2.20.4"    
changepoint  "2.2.4"     
cowplot      "1.1.2"     
DHARMa       "0.4.6"     
dplR         "1.7.6"     
dplyr        "1.1.4"     
effects      "4.2-2"     
effectsize   "0.8.6"     
faraway      "1.0.8"     
future.apply "1.11.1"    
ggbiplot     "0.6.2"     
ggfortify    "0.4.16"    
ggplot2      "3.4.4"     
gridExtra    "2.3"       
Hmisc        "5.1-1"     
Kendall      "2.2.1"     
latex2exp    "0.9.6"     
lattice      "0.22-5"    
lme4         "1.1-35.1"  
lmerTest     "3.1-3"     
maps         "3.4.2"     
MASS         "7.3-60.0.1"
mgcv         "1.9-1"     
mice         "3.16.0"    
mvtnorm      "1.2-4"     
nortest      "1.0-4"     
parameters   "0.21.3"    
purrr        "1.0.2"     
pwr          "1.3-0"     
randomForest "4.7-1.1"   
ranger       "0.16.0"    
raster       "3.6-26"    
readxl       "1.4.3"     
reshape2     "1.4.4"     
retrodesign  "0.2.1"     
rgbif        "3.7.9"     
simstudy     "0.7.1"     
sp           "2.1-2"     
TidyDensity  "1.3.0"     
tidyr        "1.3.0"     
vegan        "2.6-4"


#### Code per chapter ####
The chapter-wise R code is a bit richer than what is presented in the book. Typically we use base plotting to get a first impression and ggplot2 only for production-level plots at the end.
Out-commented code is typically aimed to check that everything worked smoothly, has the right dimensions etc. You can ignore it for reproducing the examples, but it may serve additional guide on what we did during simulation.

Occassionally we save objects and load them again. This is because the for/replicate that yielded this object took a long time to compute and we do not want to re-do the same analysis again and again. This object is not part of the supplementary code but will be produced with the code.


#### Data used in book ####
The data used in the book are in the folder "Data", each with their own subfolder, apart from those without further information or invented for this book or other (e.g. dogs.csv from Dormann 2020). Sources are given in the book.