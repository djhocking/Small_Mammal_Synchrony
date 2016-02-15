######################################################
# Small Mammal Abundance
# White Mountain National Forest
# Daniel J. Hocking
# with Ryan Stephens, Becca Rowe, Mariko Yamasaki
# 2013
######################################################

#----------import data and load packages--------
getwd() # check working directory

NAIN <- read.table('Data/NAIN_abundance_All_Years.csv', header = TRUE, sep=',')
MIPI <- read.table('Data/MIPI_abundance_All_Years.csv', header = TRUE, sep=',')
MYGA <- read.table('Data/MYGA_abundance_All_Years.csv', header = TRUE, sep=',')
PELE <- read.table('Data/PELE_abundance_All_Years.csv', header = TRUE, sep=',')
PEMA <- read.table('Data/PEMA_abundance_All_Years.csv', header = TRUE, sep=',')
SOCI <- read.table('Data/SOCI_abundance_All_Years.csv', header = TRUE, sep=',')
SODI <- read.table('Data/SODI_abundance_All_Years.csv', header = TRUE, sep=',')
ZAHU <- read.table('Data/ZAHU_abundance_All_Years.csv', header = TRUE, sep=',')
SOFU <- read.table('Data/SOFU_abundance_All_Years.csv', header = TRUE, sep=',')
SOHO <- read.table('Data/SOHO_abundance_All_Years.csv', header = TRUE, sep=',')
SOPA <- read.table('Data/SOPA_abundance_All_Years.csv', header = TRUE, sep=',')
SYCO <- read.table('Data/SYCO_abundance_All_Years.csv', header = TRUE, sep=',')
BLBR <- read.table('Data/BLBR_abundance_All_Years.csv', header = TRUE, sep=',')

Habitat <- read.table('Data/Local_Level_All_Years.csv', header = TRUE, sep=',')
#Day <- read.table('Data/Day_Detection_All_Years.csv', header = TRUE, sep=',')
Precip <- read.table('Data/Daily_Precip_Detection_All_Years.csv', header = TRUE, sep=',')
#Snap <- read.table('Data/Snap_Detection_All_Years.csv', header = TRUE, sep=',')
#Season <- read.table('Data/Spring_Fall_Detection_All_Years.csv', header = TRUE, sep=',')
Trap <- read.table('Data/TrapType_All_Years.csv', header = TRUE, sep=',')
Landscape <- read.table('Data/Landscape_Level_Variables_All_Years.csv', header = TRUE, sep=',')
Temp <- read.table('Data/Daily_Temp_Detection_All_Years.csv', header = TRUE, sep=',')

#install.packages(c('car', 'rjags'))
library(ggplot2)
library(car)
library(rjags)
library(parallel)

#---------------------- Summarize data--------------
summary(NAIN)
summary(Habitat)
summary(Day)
summary(Precip)
summary(Snap) # number of snap trap nights. Zeros for pitfall traps
summary(Season)
summary(Landscape)

str(NAIN)
str(Habitat)
str(Day)
str(Precip)
str(Landscape)



head(NAIN, 10)

# Naive Abundance
NAIN.naive <- apply(NAIN[,-c(1)], 1, sum, na.rm = T)
data.frame(NAIN$Plot, NAIN.naive)

# Scatterplot matrix and correlations
panel.hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="gray", ...)
}

panel.cor <- function(x, y, method = "pearson", use = "pairwise.complete.obs", digits = 2, prefix="", cex.cor) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  x = x
  y = y
  Cor <- switch(method,
                spearman = cor(x, y, method = "spearman", use = use),
                pearson = cor(x, y, method = "pearson", use = use),
                kendell = cor(x, y, method = "kendall", use = use),
                stop("\nThe type of correlation must be spearman, pearson, or kendall\n"))
  txt <- format(c(Cor, 0.123456789), digits = digits)[1]
  txt <- paste(prefix, txt, sep="")
  #if(missing(cex.cor)) cex <- 0.9/strwidth(txt)
  text(0.5, 0.5, txt, cex = 1.5)
  #text(0.9, 0.9, Cor.Type, cex = 1)
}

Pairs <- data.frame(NAIN.naive, Habitat[,-c(1,2,3,4,5,8, 12,13, 11)])

pairs(Pairs, upper.panel=panel.smooth, diag.panel=panel.hist, lower.panel=panel.cor)

Pairs.100 <- data.frame(Habitat[ ,"Elevation"], Landscape[ , -1])
pairs(Pairs.100, upper.panel=panel.smooth, diag.panel=panel.hist, lower.panel=panel.cor)

Pairs.500 <- data.frame(Habitat[ ,"Elevation"], Landscape[ , -1])
pairs(Pairs.500, upper.panel=panel.smooth, diag.panel=panel.hist, lower.panel=panel.cor)

Pairs.Landscape <- data.frame(Habitat[ ,"Elevation"], Landscape[ , -c(1)])
pairs(Pairs.Landscape, upper.panel=panel.smooth, diag.panel=panel.hist, lower.panel=panel.cor)

Pairs.local.100 <- data.frame( Habitat[,c("Hardwood_Ha", "Softwood_Ha", "Stream_Dis")],  Landscape[ , -1])
pairs(Pairs.local.100, upper.panel=panel.smooth, diag.panel=panel.hist, lower.panel=panel.cor)

#------Check potential covariates for normality---------
# should have rough normality before standardizing??? - no
par(mfrow = c(1, 2))
boxplot(Habitat$Elevation) ; hist(Habitat$Elevation)  # good
boxplot(Habitat$Slope) ; hist(Habitat$Slope) 
boxplot(log(Habitat$Slope + 0.5)) ; hist(log(Habitat$Slope + 0.5))
boxplot(asin(Habitat$Slope/90)) ; hist(asin(Habitat$Slope/90))
boxplot(Habitat$DRAINAGE) ; hist(Habitat$DRAINAGE) 
boxplot(log(Habitat$DRAINAGE + 0.5)) ; hist(log(Habitat$DRAINAGE + 0.5))
boxplot(Habitat$Herb) ; hist(Habitat$Herb)
boxplot(asin(Habitat$Herb/100)) ; hist(asin(Habitat$Herb/100))
boxplot(log(Habitat$Herb + 0.5)) ; hist(log(Habitat$Herb + 0.5))

par(mfrow = c(1, 1))

#------------Standardize covariates------------
elev <- Habitat$Elevation
elev.s <- (elev - mean(elev))/sd(elev)
slope <- Habitat$Slope
slope.s <- (slope - mean(slope))/sd(slope)
drainage <- Habitat$DRAINAGE
drainage.s <- (drainage - mean(drainage))/sd(drainage)
hardwood <- Habitat$Hardwood_Ha
hardwood.s <- (hardwood - mean(hardwood))/sd(hardwood)
softwood <- Habitat$Softwood_Ha
softwood.s <- (softwood - mean(softwood))/sd(softwood)
herb <- Habitat$Herb
herb.s <- (herb - mean(herb))/sd(herb)
cwd <- Habitat$CWD
cwd.s <- (cwd - mean(cwd))/sd(cwd)
litter <- Habitat$Litter
litter.s <- (litter - mean(litter))/sd(litter)
stems <- Habitat$Stem_Ha
stems.s <- (stems - mean(stems))/sd(stems)
age <- Habitat$Age_local
age.s <- (age - mean(age))/sd(age)

Landscape.s <- apply(as.matrix(Landscape[ , -1]), MARGIN = 2, FUN = function(x)(x - mean(x)) / sd(x))

Temp.s <- apply(as.matrix(Temp[ , -1]), MARGIN = 2, FUN = function(x)(x - mean(x, na.rm = T)) / sd(x, na.rm = T))

Temp.s[is.na(Temp.s)] <- 0

# Prepare and standardize data
precip.m <- mean(as.matrix(Precip[,2:17]), na.rm = TRUE)
precip.sd <- sd(unlist(as.vector(Precip[ ,2:17])), na.rm = TRUE)

precip.s <- as.matrix((Precip[,2:17] - precip.m)/precip.sd) # Scale over all precip
dimnames(precip.s) <- NULL

precip.s[is.na(precip.s)] <- 0

nplots <- dim(NAIN)[1]

#day.m <- mean(as.matrix(Day[,2:17]), na.rm = TRUE)
#day.sd <- sd(unlist(as.vector(Day[ ,2:17])), na.rm = TRUE)

day.s <- as.matrix((Day[,2:17] - day.m)/day.sd) # Scale over all day
mean(as.matrix(day.s), na.rm = TRUE) # check that = 0
sd(unlist(as.vector(day.s)), na.rm = TRUE) # check that = 1
dimnames(day.s) <- NULL
str(day.s)
dim(day.s)

year2 <- seq(0, length.out=nplots)
year3 <- seq(0, length.out=nplots)
for(i in 1:nplots){
  if(Habitat$Year[i] == 1995){
    year2[i] <- 0
    year3[i] <- 0
  }
  if(Habitat$Year[i] == 1996){
    year2[i] <- 1
    year3[i] <- 0
  }
  if(Habitat$Year[i] == 1997) {
    year2[i] <- 0
    year3[i] <- 1
  }
}

#jday <- matrix(NA, 108, 16)
#for(i in 1:nplots){
 # for(j in 1:16){
  #  jday[i, j] <- as.factor(paste(Habitat$Year[i],Day[i, j], sep = "-"))
  #}
#}
#length(unique(levels(fday))) # 156 of 1728 plot-observations

cbind(Habitat$Year, year2, year3)

Pairs <- data.frame(NAIN.naive, Habitat[, c("Site", "Slope", "Elevation", "Hardwood_Ha", "Softwood_Ha", "CWD")], year2, year3)

pairs(Pairs, upper.panel=panel.smooth, diag.panel=panel.hist, lower.panel=panel.cor)

plot(jitter(as.numeric(unlist(Precip[ , 2:17]))), jitter(as.numeric(unlist(NAIN[,2:17]))))

plot(jitter(as.numeric(unlist(Temp[ , 2:17]))), jitter(as.numeric(unlist(NAIN[,2:17]))))
lines(smooth.spline(as.numeric(unlist(Temp[ , 2:17])), as.numeric(unlist(NAIN[,2:17]))), col = "red")

m1 <- lm(as.numeric(unlist(NAIN[,2:17])) ~ as.numeric(unlist(Temp[ , 2:17])))
summary(m1)

t1 <- as.numeric(unlist(Temp[ , 2:17]))

m2 <- lm(as.numeric(unlist(NAIN[,2:17])) ~ t1 + I(t1^2))
summary(m2)
