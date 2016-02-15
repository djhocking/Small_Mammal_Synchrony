rm(list=ls())# clears workspace

##load necessary packages
library(dplyr)
library(ggplot2)
library(grid)


## Set working directory and read in data files
setwd("C:/Users/Ryan/Documents/WMNF/Species_Tables")

##loads file into R
species_data <- read.csv("All_Species.csv",header=T) #Variables from model by species
head(species_data, n = 18)

##only selects the needed variables
var<-filter(species_data,Variable == "b.age100"| Variable == "b.elev"|Variable == "b.hard100"|
  Variable == "b.mix100"|Variable == "b.stream100"|Variable == "b.trap"|Variable == "b.year2"|
  Variable == "b.year3"|Variable == "p.precip"|Variable == "p.temp"|Variable == "p.trap"|Variable =="p.temp2") 
         
#renames Detection and abundance variables
var$Variable <-as.character(var$Variable)#Change to character to re-name
  var$Variable[var$Variable == 'p.trap'] <- 'Trap type (p)'
  var$Variable[var$Variable == 'b.age100'] <- 'Forest age'            
  var$Variable[var$Variable == 'b.elev'] <- 'Elevation'
  var$Variable[var$Variable == 'b.hard100'] <- 'Hardwood forest (%)'
  var$Variable[var$Variable == 'b.mix100'] <- 'Mixed forest (%)'
  var$Variable[var$Variable == 'b.stream100'] <- 'Stream/wetland (%)'
  var$Variable[var$Variable == 'b.trap'] <- 'Trap type'
  var$Variable[var$Variable == 'b.year2'] <- '1996'
  var$Variable[var$Variable == 'b.year3'] <- '1997'
  var$Variable[var$Variable == 'p.precip'] <- 'Precipitation (p)'
  var$Variable[var$Variable == 'p.temp'] <- 'Temperature (p)'
  var$Variable[var$Variable == 'p.temp2'] <- 'Temperature ^2 (p)'
head(var)

#reorders the levels of the variable for arrangement in facet_wrap
var$Variable1 <- factor(var$Variable, levels=c("Trap type (p)", "Precipitation (p)",
              "Temperature (p)", "Temperature ^2 (p)", "Trap type","Elevation",
              "Stream/wetland (%)","Hardwood forest (%)","Mixed forest (%)",
              "Forest age","1996","1997"))

#renames species
var$Species <-as.character(var$Species)#Change to character to re-name
  var$Species[var$Species == 'SOCI'] <- 'Sorex cinereus'
  var$Species[var$Species == 'MYGA'] <- 'Myodes gapperi'
  var$Species[var$Species == 'PEMA'] <- 'Peromyscus maniculatus'
  var$Species[var$Species == 'BLBR_S'] <- 'Blarina brevicauda(S)'
  var$Species[var$Species == 'BLBR_P'] <- 'Blarina brevicauda(P)'
  var$Species[var$Species == 'PELE'] <- 'Peromyscus leucopus'
  var$Species[var$Species == 'NAIN'] <- 'Napaeozapus insignis'
  var$Species[var$Species == 'SOFU'] <- 'Sorex fumeus'
  var$Species[var$Species == 'SYCO'] <- 'Synaptomys cooperi'
  var$Species[var$Species == 'SOHO'] <- 'Sorex hoyi'
  var$Species[var$Species == 'MIPE'] <- 'Microtus pennsylvanicus'
  var$Species[var$Species == 'ZAHU'] <- 'Zapus hudsonius'
  var$Species[var$Species == 'SODI'] <- 'Sorex dispar'
  var$Species[var$Species == 'SOPA'] <- 'Sorex palustris'
  var$Species[var$Species == 'MICH'] <- 'Microtus chrotorrhinus'
  var$Species[var$Species == 'TAST'] <- 'Tamias striatus'

#reorders the levels of the variable for arrangement in facet_wrap
var$Species1 <- factor(var$Species, levels=c("Tamias striatus","Microtus chrotorrhinus","Sorex palustris",
            "Sorex dispar","Zapus hudsonius","Microtus pennsylvanicus","Sorex hoyi","Synaptomys cooperi",
            "Sorex fumeus","Napaeozapus insignis","Peromyscus leucopus", "Blarina brevicauda(P)",
            "Blarina brevicauda(S)", "Peromyscus maniculatus", "Myodes gapperi","Sorex cinereus"))
  
#subsets the variables that are significant 
sig <- subset(var, Mean < 0 & X97.50.< 0 | Mean > 0 & X2.50.> 0)
head(sig)

head(var, n=18)

###############################################################################################
#Facets to display all species and variables#http://www.cookbook-r.com/Graphs/Facets_(ggplot2)/
###############################################################################################
Species<-ggplot(var, aes(x=Mean, y=Species1)) + geom_point(shape=1, size=4)+#plots data
  facet_wrap( ~ Variable1, scales="free_x")+#orders variables and scales each one indepenently 
  geom_errorbarh(height=.1, aes(xmin=X2.50., xmax=X97.50.))+#error bars 
  geom_point(data=var, color="white", size=4)+#makes points white
  geom_point(data=var, colour="black",pch=21,size=4)+#puts black border around the cirlce to make "open cicles"
  geom_point(data=sig, size=4)+#plots points that are significant and gives them a solid fill
  geom_vline(xintercept = 0)+#adds a verticle line at zero 

#removes background but keeps the boarder lines
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+#removes gridlines, but keeps border
  theme(panel.grid.major = element_line(colour = "gray95",size = 2.5))+#makes a very light gray line through each point 
  theme(panel.grid.major.x = element_blank())+#removes verticle line

#Formats the text
  theme(legend.key = element_blank())+#removes box around legend
  theme(strip.background = element_blank())+#removes gray background from titles
  theme(strip.text.x = element_text(size = 18, colour = "black"))+#Font size of variable lables
  theme(axis.text.y = element_text(face = "italic", size = 12))+#Font size for y axis
  
#makes shrews bold and italic and rodents just italic - basically renames them too
  scale_y_discrete(breaks=levels(factor(var$Species1)),
                   labels = c(expression(italic("Tamias striatus")),
                              expression(italic("Microtus chrotorrhinus")),
                              expression(bolditalic("Sorex palustris")),
                              expression(bolditalic("Sorex dispar")),
                              expression(italic("Zapus hudsonius")),
                              expression(italic("Microtus pennsylvanicus")),
                              expression(bolditalic("Sorex hoyi")),
                              expression(italic("Synaptomys cooperi")),
                              expression(bolditalic("Sorex fumeus")),
                              expression(italic("Napaeozapus insignis")),
                              expression(italic("Peromyscus leucopus")),
                              expression(bolditalic("Blarina brevicauda(P)*")),
                              expression(bolditalic("Blarina brevicauda(S)*")),
                              expression(italic("Peromyscus maniculatus")),
                              expression(italic("Myodes gapperi")),
                              expression(bolditalic("Sorex cinereus"))))+  
  theme(axis.text.x = element_text(size = 12.5))+#font size for x axis
  theme(panel.margin = unit(1, "lines"))+#adds space inbetween graphs (need package "grid" for this)
  theme(axis.title = element_blank())#Removes the x and y legends

#This changes the labels so that "p" for detection variables are in italics 
#http://stackoverflow.com/questions/19282897/how-to-add-expressions-to-labels-in-facet-wrap
facet_wrap_labeller <- function(gg.plot,labels=NULL) {
  #works with R 3.0.1 and ggplot2 0.9.3.1
  require(gridExtra)
  
  g <- ggplotGrob(gg.plot)
  gg <- g$grobs      
  strips <- grep("strip_t", names(gg))
  
  for(ii in seq_along(labels))  {
    modgrob <- getGrob(gg[[strips[ii]]], "strip.text", 
                       grep=TRUE, global=TRUE)
    gg[[strips[ii]]]$children[[modgrob$name]] <- editGrob(modgrob,label=labels[ii])
  }
  
  g$grobs <- gg
  class(g) = c("arrange", "ggplot",class(g)) 
  g
}

Species1<-facet_wrap_labeller(Species, labels = c(expression(paste("Trap type (", italic("p"),")")),#renames the first 4 panes so they have special characters
                                                  expression(paste("Precipitation (", italic("p"),")")),
                                                  expression(paste("Temperature (", italic("p"),")")),
                                                  expression(paste(Temperature^2," ","(", italic("p"),")"))))
####################################################################################################

pdf("Species_Response.pdf",height=11,width=16,useDingbats=F)#makes pdf of the graph
Species1
dev.off()




