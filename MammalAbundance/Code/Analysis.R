

source('Code/Poisson_100.R')

dm100p.NAIN <- dm100p(NAIN, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/NAIN_P_Diagnostic.pdf", outfile2 = "Output/NAIN_Table.csv") # 
dm100p.MIPI <- dm100p(MIPI, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/MIPI_P_Diagnostic.pdf", outfile2 = "Output/MIPI_Table.csv") # 
dm100p.MYGA <- dm100p(MYGA, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/MYGA_P_Diagnostic.pdf", outfile2 = "Output/MYGA_Table.csv") # 
dm100p.PELE <- dm100p(PELE, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/PELE_P_Diagnostic.pdf", outfile2 = "Output/PELE_Table.csv") # 
dm100p.PEMA <- dm100p(PEMA, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/PEMA_P_Diagnostic.pdf", outfile2 = "Output/PEMA_Table.csv") # 
dm100p.SOCI <- dm100p(SOCI, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/SOCI_P_Diagnostic.pdf", outfile2 = "Output/SOCI_Table.csv") #
dm100p.SODI <- dm100p(SODI, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/SODI_P_Diagnostic.pdf", outfile2 = "Output/SODI_Table.csv") #
dm100p.ZAHU <- dm100p(ZAHU, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/ZAHU_P_Diagnostic.pdf", outfile2 = "Output/ZAHU_Table.csv") #
dm100p.SOFU <- dm100p(SOFU, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/SOFU_P_Diagnostic.pdf", outfile2 = "Output/SOFU_Table.csv") #
dm100p.SOHO <- dm100p(SOHO, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/SOHO_P_Diagnostic.pdf", outfile2 = "Output/SOHO_Table.csv") # 
dm100p.SOPA <- dm100p(SOPA, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/SOPA_P_Diagnostic.pdf", outfile2 = "Output/SOPA_Table.csv") # 
dm100p.SYCO <- dm100p(SYCO, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/SYCO_P_Diagnostic.pdf", outfile2 = "Output/SYCO_Table.csv") #
dm100p.BLBR <- dm100p(BLBR, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/BLBR_P_Diagnostic.pdf", outfile2 = "Output/BLBR_Table.csv") # 

source('Code/helper_functions.R')
abund.NAIN <- abund(model = dm100p.NAIN, output = "Output/N_NAIN.csv", sep = ",")
abund.MIPI <- abund(model = dm100p.MIPI, output = "Output/N_MIPI.csv", sep = ",")
abund.MYGA <- abund(model = dm100p.MYGA, output = "Output/N_MYGA.csv", sep = ",")
abund.PELE <- abund(model = dm100p.PELE, output = "Output/N_PELE.csv", sep = ",")
abund.PEMA <- abund(model = dm100p.PEMA, output = "Output/N_PEMA.csv", sep = ",")
abund.SOCI <- abund(model = dm100p.SOCI, output = "Output/N_SOCI.csv", sep = ",")
abund.SODI <- abund(model = dm100p.SODI, output = "Output/N_SODI.csv", sep = ",")
abund.ZAHU <- abund(model = dm100p.ZAHU, output = "Output/N_ZAHU.csv", sep = ",")
abund.SOFU <- abund(model = dm100p.SOFU, output = "Output/N_SOFU.csv", sep = ",")
abund.SOHO <- abund(model = dm100p.SOHO, output = "Output/N_SOHO.csv", sep = ",")
abund.SOPA <- abund(model = dm100p.SOPA, output = "Output/N_SOPA.csv", sep = ",")
abund.SYCO <- abund(model = dm100p.SYCO, output = "Output/N_SYCO.csv", sep = ",")
abund.BLBR <- abund(model = dm100p.BLBR, output = "Output/N_BLBR.csv", sep = ",")


############# 500m ##############

dm500p.NAIN <- dm500p(NAIN, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/NAIN_P_Diagnostic_500p.pdf", outfile2 = "Output/NAIN_Table_500p.csv") # 
dm500p.MIPI <- dm500p(MIPI, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/MIPI_P_Diagnostic_500p.pdf", outfile2 = "Output/MIPI_Table_500p.csv") # 
dm500p.MYGA <- dm500p(MYGA, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/MYGA_P_Diagnostic_500p.pdf", outfile2 = "Output/MYGA_Table_500p.csv") # 
dm500p.PELE <- dm500p(PELE, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/PELE_P_Diagnostic_500p.pdf", outfile2 = "Output/PELE_Table_500p.csv") # 
dm500p.PEMA <- dm500p(PEMA, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/PEMA_P_Diagnostic_500p.pdf", outfile2 = "Output/PEMA_Table_500p.csv") # 
dm500p.SOCI <- dm500p(SOCI, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/SOCI_P_Diagnostic_500p.pdf", outfile2 = "Output/SOCI_Table_500p.csv") #
dm500p.SODI <- dm500p(SODI, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/SODI_P_Diagnostic_500p.pdf", outfile2 = "Output/SODI_Table_500p.csv") #
dm500p.ZAHU <- dm500p(ZAHU, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/ZAHU_P_Diagnostic_500p.pdf", outfile2 = "Output/ZAHU_Table_500p.csv") #
dm500p.SOFU <- dm500p(SOFU, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/SOFU_P_Diagnostic_500p.pdf", outfile2 = "Output/SOFU_Table_500p.csv") #
dm500p.SOHO <- dm500p(SOHO, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/SOHO_P_Diagnostic_500p.pdf", outfile2 = "Output/SOHO_Table_500p.csv") # 
dm500p.SOPA <- dm500p(SOPA, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/SOPA_P_Diagnostic_500p.pdf", outfile2 = "Output/SOPA_Table_500p.csv") # 
dm500p.SYCO <- dm500p(SYCO, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/SYCO_P_Diagnostic_500p.pdf", outfile2 = "Output/SYCO_Table_500p.csv") #
dm500p.BLBR <- dm500p(BLBR, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/BLBR_P_Diagnostic_500p.pdf", outfile2 = "Output/BLBR_Table_500p.csv") #

############# Local ##############
source('Code/Poisson_local.R')

dmp.NAIN <- dmp(NAIN, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/Local/NAIN_P_Diagnostic_p.pdf", outfile2 = "Output/Local/NAIN_Table_p.csv") # 
dmp.MIPI <- dmp(MIPI, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/Local/MIPI_P_Diagnostic_p.pdf", outfile2 = "Output/Local/MIPI_Table_p.csv") # 
dmp.MYGA <- dmp(MYGA, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/Local/MYGA_P_Diagnostic_p.pdf", outfile2 = "Output/Local/MYGA_Table_p.csv") # 
dmp.PELE <- dmp(PELE, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/Local/PELE_P_Diagnostic_p.pdf", outfile2 = "Output/Local/PELE_Table_p.csv") # 
dmp.PEMA <- dmp(PEMA, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/Local/PEMA_P_Diagnostic_p.pdf", outfile2 = "Output/Local/PEMA_Table_p.csv") # 
dmp.SOCI <- dmp(SOCI, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/Local/SOCI_P_Diagnostic_p.pdf", outfile2 = "Output/Local/SOCI_Table_p.csv") #
dmp.SODI <- dmp(SODI, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/Local/SODI_P_Diagnostic_p.pdf", outfile2 = "Output/Local/SODI_Table_p.csv") #
dmp.ZAHU <- dmp(ZAHU, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/Local/ZAHU_P_Diagnostic_p.pdf", outfile2 = "Output/Local/ZAHU_Table_p.csv") #
dmp.SOFU <- dmp(SOFU, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/Local/SOFU_P_Diagnostic_p.pdf", outfile2 = "Output/Local/SOFU_Table_p.csv") #
dmp.SOHO <- dmp(SOHO, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/Local/SOHO_P_Diagnostic_p.pdf", outfile2 = "Output/Local/SOHO_Table_p.csv") # 
dmp.SOPA <- dmp(SOPA, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/Local/SOPA_P_Diagnostic_p.pdf", outfile2 = "Output/Local/SOPA_Table_p.csv") # 
dmp.SYCO <- dmp(SYCO, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/Local/SYCO_P_Diagnostic_p.pdf", outfile2 = "Output/Local/SYCO_Table_p.csv") #
dmp.BLBR <- dmp(BLBR, n.burn = 50000, n.it = 50000, n.thin = 50, outfile = "Output/Local/BLBR_P_Diagnostic_p.pdf", outfile2 = "Output/Local/BLBR_Table_p.csv") #





##################### Extra ###############



#-----------NAIN ---------okay - some autocorr------
# prepare data for this model


# Add random year effect?
# Take out things that correlate with year

nmixp.nain1 <- nmixp(NAIN, n.burn = 1000, n.it = 1000, n.thin = 1, K = 1, outfile = "Output/NAIN_diagnostics_nmixp1.pdf")
nmixp.nain2 <- nmixp(NAIN, n.burn = 1000, n.it = 1000, n.thin = 1, K = 2, outfile = "Output/NAIN_diagnostics_nmixp2.pdf")

dmp.nain <- dmp(NAIN, n.burn = 1000, n.it = 1000, n.thin = 1, outfile = "Output/NAIN_diagnostics.pdf")
dmp.myga <- dmp(MYGA, n.burn = 1000, n.it = 1000, n.thin = 1, outfile = "Output/MYGA_diagnostics.pdf")

dm100p.nain <- dm100p(NAIN, n.burn = 1000, n.it = 1000, n.thin = 1, outfile = "Output/NAIN_diagnostics.pdf")
dm100pod.nain <- dm100pod(NAIN, n.burn = 1000, n.it = 3000, n.thin = 1)
dm100zip.nain <- dm100zip(NAIN, n.burn = 3000, n.it = 5000, n.thin = 3)
dm100zipod.nain <- dm100zipod(NAIN, n.burn = 3000, n.it = 5000, n.thin = 3)

dm100pgama.nain <- dm100pgama(NAIN, n.burn = 3000, n.it = 3000, n.thin = 1, outfile = "NAIN_diagnostics.pdf")

dm100pgama.myga <- dm100pgama(MYGA, n.burn = 1000, n.it = 3000, n.thin = 3, outfile = "MYGA_diagnostics.pdf")

dm100p.year.myga <- dm100p.year(MYGA, n.burn = 1000, n.it = 3000, n.thin = 3, outfile = "Output/MYGA_diagnostics.pdf")

library(ggmcmc)
ggs_traceplot(ggs(dm100p.nain[ , c("p0", "p.precip","p.trap","a.N", "b.soft100", "b.age100", "b.stream100", "b.trap", "b.year2", "b.year3", "b.elev", "sigma.gam0", "gam1", "sigma.site")]))

ggmcmc(ggs(dm100p.nain[ , c("p0", "p.precip","p.trap","a.N", "b.soft100", "b.age100", "b.stream100", "b.trap", "b.year2", "b.year3", "b.elev", "sigma.gam0", "gam1", "sigma.site")]), file = "dm100p_NAIN.pdf")

plot(dm100p.nain[ , c("p0", "p.precip","p.trap","a.N", "b.soft100", "b.age100", "b.stream100", "b.trap", "b.year2", "b.year3", "b.elev", "sigma.gam0", "gam1", "sigma.site")]) # "gam.a",

plot(dm100pod.nain[ , c("p0", "p.precip","p.trap","a.N", "b.soft100", "b.age100", "b.stream100", "b.trap", "b.year2", "b.year3", "b.elev", "sigma.gam0", "gam1", "sigma.site", "sigma.delta")])

plot(dm100zip.nain[ , c("p0", "p.precip","p.trap","a.N", "b.soft100", "b.age100", "b.stream100", "b.trap", "b.year2", "b.year3", "b.elev", "sigma.gam0",  "gam.a", "gam1", "sigma.site", "omega")])

plot(dm100zipod.nain[ , c("p0", "p.precip","p.trap","a.N", "b.soft100", "b.age100", "b.stream100", "b.trap", "b.year2", "b.year3", "b.elev", "sigma.gam0", "gam.a", "gam1", "sigma.site", "sigma.delta", "omega")])



par(mfrow=c(1,1))
summary(dm100.list.nain[ , c("p0", "p.precip","p.trap","a.N", "b.hard100", "b.soft100", "b.age100", "b.stream100", "b.trap", "b.year2", "b.year3", "b.elev", "sigma.gam0", "gam1", "sigma.site")])

autocorr.plot(dm100.list.nain[ , c("p0", "p.precip","p.trap","a.N", "b.hard100", "b.soft100", "b.age100", "b.stream100", "b.trap", "b.year2", "b.year3", "b.elev", "sigma.gam0", "gam1", "sigma.site")])

acfplot(dm100.list.nain[ , c("p0", "p.precip","p.trap","a.N", "b.hard100", "b.soft100", "b.age100", "b.stream100", "b.trap", "b.year2", "b.year3", "b.elev", "sigma.gam0", "gam1", "sigma.site")])



# Check fit
for(i in 1:3) bayesP.nain <- mean(dm3.list.nain[, "fit.new",][[i]] > dm3.list.nain[, "fit",][[i]]) # 
print(bayesP.nain, dig = 3) # 

par(mfrow=c(1,1))
plot(as.matrix(dm3.list.nain[, "fit",]), as.matrix(dm3.list.nain[, "fit.new",])) # 
abline(0, 1, col = 'red')

print(gelman.diag(x=dm3.list.nain[ , c("p0", "p.precip","p.trap","a.N", "b.drainage", "b.stems", "b.softwood", "b.hardwood", "b.herb","b.cwd", "b.litter", "b.trap", "b.year2", "b.year3", "b.elev", "sigma.gam0", "gam1", "sigma.site")]), dig=3) #

gelman.diag(dm3.list.nain) # multivariate psrf = can't run

geweke.diag(dm3.list.nain[ , c("p0", "p.precip","p.trap","a.N", "b.drainage", "b.stems", "b.softwood", "b.hardwood", "b.herb","b.cwd", "b.litter", "b.trap", "b.year2", "b.year3", "b.elev", "sigma.gam0", "gam1", "sigma.site")])

heidel.diag(dm3.list.nain)

rejectionRate(dm3.list.nain[ , c("p0", "p.precip","p.trap","a.N", "b.drainage", "b.stems", "b.softwood", "b.hardwood", "b.herb","b.cwd", "b.litter", "b.trap", "b.year2", "b.year3", "b.elev", "sigma.gam0", "gam1", "sigma.site")])

rejectionRate(dm3.list.nain)

accept.nain <- 1 - rejectionRate(dm3.list.nain)

HPDinterval(dm3.list.nain)

# Calculate and organize summaries for publication tables--------------
Quants.nain <- apply(as.matrix(dm3.list.nain[ , c("p0", "p.precip","p.trap","a.N", "b.drainage", "b.stems", "b.softwood", "b.hardwood", "b.herb","b.cwd", "b.litter", "b.trap", "b.year2", "b.year3", "b.elev", "sigma.gam0", "gam1", "sigma.site")]), 2, FUN = quantile, probs = c(0.025, 0.5, 0.975))

Means.nain <- apply(as.matrix(dm3.list.nain[ , c("p0", "p.precip","p.trap","a.N", "b.drainage", "b.stems", "b.softwood", "b.hardwood", "b.herb","b.cwd", "b.litter", "b.trap", "b.year2", "b.year3", "b.elev", "sigma.gam0", "gam1", "sigma.site")]), 2, FUN = mean)

SDs.nain <- apply(as.matrix(dm3.list.nain[ , c("p0", "p.precip","p.trap","a.N", "b.drainage", "b.stems", "b.softwood", "b.hardwood", "b.herb","b.cwd", "b.litter", "b.trap", "b.year2", "b.year3", "b.elev", "sigma.gam0", "gam1", "sigma.site")]), 2, FUN = sd)

nain.variables <- c("p-intercept", "Precip", "Trap Type", "N-intercept", "Drainage", "Number Stems", "Softwood basal", "Hardwood basal", "Herbaceous", "CWD", "Litter", "Trap Type", "Year (1996)", "Year (1997)", "Elevation", "Random gamma SD", "Autoreg N", "Site SD")

nain.summary <- data.frame(nain.variables, Means.nain, SDs.nain, Quants.nain["2.5%", ], Quants.nain["50%", ], Quants.nain["97.5%", ])

colnames(nain.summary) <- c("Variable", "Mean", "SD", "2.5%", "Median", "97.5%")

write.table(nain.summary, file = "nain_Summary.csv", sep = ",", col.names = NA, row.names = TRUE)

N.nain <- matrix(NA, 108, 6)
for(i in 1:108){
  for(j in 1:2){
    foo <- apply(as.matrix(dm100p.nain[, c(paste("N[", i,",", j, "]", sep = ""))]), 2, FUN= quantile, probs = c(0.5, 0.025, 0.975))
    foo <- as.integer(foo)
    for(k in 1:3){
      if(j == 1) {
        N.nain[i,k] <- foo[k]
      }
      if(j == 2){
        N.nain[i,k+3] <- foo[k]
      }
    }
  }
}

ymax1 <- apply(as.matrix(NAIN[,2:9, 1]),1,sum)
ymax2 <- apply(NAIN[,10:17, 2],1,sum)
colnames(N.nain) <- c("Median1", "CI_2.5", "CI_97.5", "Median2", "CI_2.5", "CI_97.5")
cbind(ymax1, N.nain[ ,"Median1"], ymax2, N.nain[, "Median2"])

write.table(N.nain, file = "N_nain.csv", col.names = NA, row.names = TRUE, sep = ",")

range(N.nain[, "Median"]) # 
mean(N.nain[, "Median"]) # 

#---------------Figures------------------
# Abundance

# Elevation effects
plot(Data$elev, N.nain[ , 2])

elev.range <- seq(from = min(Data$elev), to = 2025, length.out = 1000)
elev.sd <- sd(Data$elev)
elev.mean <- mean(Data$elev)
elev.range.std <- (elev.range - elev.mean)/elev.sd

nain.elev.nbar <- exp(Means.nain["alpha.lam"] + Means.nain["beta1.lam"]*elev.range.std)
nain.elev.LCI <- exp(Quants.nain[1, "alpha.lam"] + Quants.nain[1, "beta1.lam"]*elev.range.std)
nain.elev.UCI <- exp(Quants.nain[3, "alpha.lam"] + Quants.nain[3, "beta1.lam"]*elev.range.std)

bitmap("Plot_nain_Elev.tiff", height = 12, width = 17, units = 'cm', type = "tifflzw", res = 300)
tiff("Plot_nain_Elev.tiff", height = 12, width = 17, units = 'cm', compression = "lzw", res = 300)
postscript("Plot_nain_Elev-courier.eps", width = 8, height = 8, horizontal = FALSE, onefile = FALSE, paper = "special", colormodel = "cmyk", family = "Courier")
par(mar=c(3.5,3.5,1,1), mgp=c(2,0.7,0), tck=-0.01)
plot(elev.range, nain.elev.nbar, type = 'n', xlab = 'Elevation (m)', ylab = expression(paste(italic(D.)," ",italic(wrighti), " abundance")))
polygon(c(elev.range, rev(elev.range)), c(nain.elev.UCI, rev(nain.elev.LCI)), col = 'light gray', border = NA)
lines(elev.range, nain.elev.nbar)
dev.off()

# Slope
slope.range <- seq(from = min(Data$Slope), to = max(Data$Slope), length.out = 1000)
slope.sd <- sd(Data$Slope)
slope.mean <- mean(Data$Slope)
slope.range.std <- (slope.range - slope.mean)/slope.sd

nain.slope.nbar <- exp(Means.nain["alpha.lam"] + Means.nain["beta3.lam"]*slope.range.std)
nain.slope.LCI <- exp(Quants.nain[1, "alpha.lam"] + Quants.nain[1, "beta3.lam"]*slope.range.std)
nain.slope.UCI <- exp(Quants.nain[3, "alpha.lam"] + Quants.nain[3, "beta3.lam"]*slope.range.std)

bitmap("Plot_nain_slope.tiff", height = 12, width = 17, units = 'cm', type = "tifflzw", res = 300)
tiff("Plot_nain_slope.tiff", height = 12, width = 17, units = 'cm', compression = "lzw", res = 300)
postscript("Plot_nain_slope.eps", width = 8, height = 8, horizontal = FALSE, onefile = FALSE, paper = "special", colormodel = "cmyk", family = "Times")
par(mar=c(3.5,3.5,1,1), mgp=c(2,0.7,0), tck=-0.01)
plot(slope.range, nain.slope.nbar, type = 'n', xlab = 'Slope', ylab = expression(paste(italic(D.)," ",italic(wrighti), " abundance")))
polygon(c(slope.range, rev(slope.range)), c(nain.slope.UCI, rev(nain.slope.LCI)), col = 'light gray', border = NA)
lines(slope.range, nain.slope.nbar)
dev.off()

# Combined
bitmap("Plot_nain_abund.tiff", height = 12, width = 17, units = 'cm', type = "tifflzw", res = 300)
tiff("Plot_nain_abund.tiff", height = 12, width = 17, units = 'cm', compression = "lzw", res = 300)
postscript("Plot_nain_abund.eps", width = 8, height = 8, horizontal = FALSE, onefile = FALSE, paper = "special", colormodel = "cmyk", family = "Times")
par(mfrow = c(1,2), mar=c(3.5,3.5,1,1), mgp=c(2,0.7,0), tck=-0.01)
plot(slope.range, nain.slope.nbar, type = 'n', xlab = 'Slope', ylab = expression(paste(italic(D.)," ",italic(wrighti), " abundance")))
polygon(c(slope.range, rev(slope.range)), c(nain.slope.UCI, rev(nain.slope.LCI)), col = 'light gray', border = NA)
lines(slope.range, nain.slope.nbar)

plot(elev.range, nain.elev.nbar, type = 'n', ylim = c(0, 150), xlab = 'Elevation (m)', ylab = '')
polygon(c(elev.range, rev(elev.range)), c(nain.elev.UCI, rev(nain.elev.LCI)), col = 'light gray', border = NA)
lines(elev.range, nain.elev.nbar)
dev.off()
par(mfrow = c(1,1))

# Detection
# Precip effects
precip.range <- seq(from = min(Precip, na.rm = T), to = max(Precip, na.rm = T), length.out = 1000)
precip.sd <- sd(as.matrix(Precip), na.rm = TRUE)
precip.mean <- mean(as.matrix(Precip), na.rm = TRUE)
precip.range.std <- (precip.range - precip.mean)/precip.sd

nain.precip.nbar <- 1 / (1 + exp(-1*(Means.nain["alpha.p"] + Means.nain["beta3.p"]*precip.range.std)))
nain.precip.LCI <- 1 / (1 + exp(-1 * (Quants.nain[1, "alpha.p"] + Quants.nain[1, "beta3.p"]*precip.range.std)))
nain.precip.UCI <- 1 / (1 + exp(-1 * (Quants.nain[3, "alpha.p"] + Quants.nain[3, "beta3.p"]*precip.range.std)))

# Humidity effects
RH.range <- seq(from = min(RH, na.rm = T), to = max(RH, na.rm = T), length.out = 1000)
RH.sd <- sd(as.matrix(RH), na.rm = TRUE)
RH.mean <- mean(as.matrix(RH), na.rm = TRUE)
RH.range.std <- (RH.range - RH.mean)/RH.sd

nain.RH.nbar <- 1 / (1 + exp(-1 * (Means.nain["alpha.p"] + Means.nain["beta10.p"]*RH.range.std)))
nain.RH.LCI <- 1 / (1 + exp(-1 * (Quants.nain[1, "alpha.p"] + Quants.nain[1, "beta10.p"]*RH.range.std)))
nain.RH.UCI <- 1 / (1 + exp(-1* (Quants.nain[3, "alpha.p"] + Quants.nain[3, "beta10.p"]*RH.range.std)))

# Temperature effects
Temp.range <- seq(from = min(Temp, na.rm = T), to = max(Temp, na.rm = T), length.out = 1000)
Temp.sd <- sd(as.matrix(Temp), na.rm = TRUE)
Temp.mean <- mean(as.matrix(Temp), na.rm = TRUE)
Temp.range.std <- (Temp.range - Temp.mean)/Temp.sd

nain.Temp.nbar <- 1 / (1 + exp(-1 * (Means.nain["alpha.p"] + Means.nain["beta1.p"]*Temp.range.std + Means.nain["beta2.p"]*Temp.range.std^2)))
nain.Temp.LCI <- 1 / (1 + exp(-1 * (Quants.nain[1, "alpha.p"] + Quants.nain[1, "beta1.p"]*Temp.range.std)))
nain.Temp.UCI <- 1 / (1 + exp(-1* (Quants.nain[3, "alpha.p"] + Quants.nain[3, "beta1.p"]*Temp.range.std)))

# Combined
bitmap("Plot_nain_detection.tiff", height = 12, width = 17, units = 'cm', type = "tifflzw", res = 300)
tiff("Plot_nain_detection.tiff", height = 12, width = 17, units = 'cm', compression = "lzw", res = 300)
postscript("Plot_nain_detection.eps", width = 8, height = 8, horizontal = FALSE, onefile = FALSE, paper = "special", colormodel = "cmyk", family = "Times")
par(mfrow = c(1,2), mar=c(3.5,3.5,1,1), mgp=c(2,0.7,0), tck=-0.01)
plot(precip.range, nain.precip.nbar, type = 'n', ylim = c(0, 0.6), xlab = '24-hour precipitation', ylab = expression(paste(italic(D.)," ",italic(wrighti), " detection probability")))
polygon(c(precip.range, rev(precip.range)), c(nain.precip.UCI, rev(nain.precip.LCI)), col = 'light gray', border = NA)
lines(precip.range, nain.precip.nbar)

plot(RH.range, nain.RH.nbar, type = 'n', ylim = c(0, 0.6), xlab = 'Relative humidity', ylab = '')
polygon(c(RH.range, rev(RH.range)), c(nain.RH.UCI, rev(nain.RH.LCI)), col = 'light gray', border = NA)
lines(RH.range, nain.RH.nbar)
dev.off()
par(mfrow = c(1,1))
