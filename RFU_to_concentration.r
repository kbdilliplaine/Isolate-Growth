#For conversion of RFU to cell concentration based using count data.
library("ggplot2")
library("readxl")
library("dplyr")
library("tidyr")
library("ggthemes")
library(tidyverse)
library(data.table)
library(growthrates)
library(broom)

#Scan files and xml

#set function file as xlsx
f <- list.files(pattern="xlsx$")
#Read in all files and combine them into a single one
TOTAL <- map_df(f, read_excel)
TOTAL$Red=lapply(strsplit(as.character(TOTAL$RawData), " "), as.numeric)
TOTAL$Red=lapply(TOTAL$Red, mean)
TOTAL$Red=as.numeric(TOTAL$Red)
#change name2 to well
names(TOTAL)[names(TOTAL) == "Name2"] <- "Well"
#Split plate name
TOTAL=separate(TOTAL,Name, into=c("Species","Light"), sep="_")
#Set factors
TOTAL[, c(1,2,27,28)] = lapply(TOTAL[, c(1,2,27,28)], as.factor)

#convert to date
TOTAL$ReadTime=as.Date(TOTAL$ReadTime, format=c("%I:%M %p %m/%d/%Y"))
#Set the 0 Date
startdate <- as.Date("18/01/2021","%d/%m/%Y")
#Make a new column for the Day of Experiment
TOTAL$Days=as.numeric(difftime(TOTAL$ReadTime,startdate ,units="days"))
Order <- factor(c("Blk","3", "10", "20", "50", "100","125"))
TOTAL$Light <- factor(TOTAL$Light, levels=Order)
#make new DF based on mean of each days blank. Subtract blank from values.
blank=TOTAL %>%
    group_by(Days, Col) %>%
    summarize(blank=mean(Red[Species=="blk"]))
TOTAL=merge(TOTAL, blank)
TOTAL$norm=TOTAL$Red-TOTAL$blank

#subset where norm is blanked RFU
myvars = c("Species", "Light","Col", "Row", "Well", "Days", "norm")
sub = TOTAL[myvars]
sub=sub %>% filter(!is.na(Light) & !is.na(Days))
sub=droplevels(sub)
sub$Row=as.numeric(sub$Row)
Mean=sub %>% group_by(Species, Light, Days, Col) %>% summarize( AD=mean(norm[Row < 5]), EH=mean(norm[Row >4]))

Mean2=gather(Mean, key="Rows", value="RFU", c("AD", "EH"))
Mean2$Species=recode_factor(Mean2$Species, exp03="E3", exp10="E10", exp18="E18", exp22="E22" )
#import cell counts
setwd("./Cell Counts")
Hand=read_excel("Cell_Count_Hand.xlsx", sheet=1)
Harvest=read_excel("Cell_Harvest.xlsx", sheet=1)
H2=merge(Hand, Harvest, all = TRUE)

D1=read.csv("Processed_Summary2.csv")
D1$File=gsub('.TIF', '', D1$Slice)
D2= D1 %>% separate(File, into=c("Species", "Days", "Light", "Rep", "Tile"), sep="_")
D3=D2 %>% group_by(Species, Light, Days, Rep) %>% summarize(Count=sum(Count), Tiles=n())
D3$Count=D3$Count/0.89
D3=as.data.frame(D3)
D3[, c(2:6)] = lapply(D3[, c(2:6)], as.numeric)
H3=H2 %>%
        full_join(D3, by = c("Species", "Light", "Days", "Rep")) %>%
        mutate(Count = coalesce(Count.x, Count.y), Tiles = coalesce(Tiles.x, Tiles.y))
        ###Merging in RFU
H4=merge(H3, Mean2, all.x=T)
H4$avgtile=(H4$Count/H4$Tiles)/(1-(H4$Lugols_vol/(H4$Lugols_vol+H4$Sample_vol)))
H4$concml=round(H4$avgtile*1000, 0)
H4$rfu_cell=H4$RFU/(H4$avgtile*250)
Initial=H4 %>% filter(Days=="0")
Initialmean=Initial %>% group_by(Species, Light, Days) %>% summarise(conc=mean(concml), sd=sd(concml))
Initialmean$conc=round(Initialmean$conc, 0)
Initialmean$sd=round(Initialmean$sd, 0)

InitialmeanS=Initial %>% group_by(Species, Days) %>% summarise(conc=mean(concml), sd=sd(concml))
InitialmeanS$conc=round(InitialmeanS$conc, 0)
InitialmeanS$sd=round(InitialmeanS$sd, 0)

InitialmeanS=Initial %>% group_by(Species, Days) %>% summarise(median=median(concml)/10, mean=mean(concml)/10)
#dividing by 10 gives you the concentration in the actual well


#linear regression with r2
library(ggpubr)
H4$Species=recode_factor(H4$Species, "E3"= "A. septentrionalis", "E10"="F. cylindrus", "E18"="S. hyperborea str. 1", "E22" ="S. hyperborea str. 2" )

H4$Days=as.factor(H4$Days)
RFU_CellLinreg=
    ggplot(H4, aes(concml, RFU))+geom_smooth(method = "lm", se=F, na.rm=T)+
    geom_point(stat="identity", aes(colour=Days))+
    labs(y=expression(Relative~fluorescence),
         x=expression(Concentration~"("*cells~"mL"^"-1"*")" ),
         colour="Day of \ncollection")+
    theme(strip.text.y = element_text( face="bold.italic", size=15),
          strip.text.x = element_text(face="bold", size=15),
          panel.border = element_rect(fill=NA, colour = "black", size=1),
          axis.text = element_text( colour="black", size=10),
          legend.title=element_text(face="bold", size=20),
          legend.text =element_text(size=15),
          panel.spacing=unit(1.5,"lines"),
          axis.title = element_text(face="bold", size=25))+
    facet_grid(Species~Light, scales="free")+ggpubr::stat_cor(aes(label = ..rr.label..), color = "red", geom = "label")

#Force through 0 intercept
  coef_cts=H4 %>%
    group_by(Species, Light) %>%
    do({mod = lm(concml ~ RFU-1, data = .)
    tidy(mod)})
  coef_cts$term=recode_factor(coef_cts$term, "(Intercept)"="b", RFU="m")
  coef_cts=pivot_wider(  coef_cts, id_cols=c(Species, Light), names_from = term, values_from = c(estimate, std.error, statistic, p.value))
  coef_cts$Species=recode_factor(coef_cts$Species, E3="Attheya sp.", E10="Fragilariopsis sp.", E22="Synedropsis sp. 1", E18 = "Synedropsis sp. 2")




























setwd("C:/Users/Dilliplaine/Desktop/CMI 2020/WAF_Experiment/Master Data Sheet")
Full=read_excel("Master_uK_well.xlsx", sheet=1)
Full$Species=recode_factor(Full$Species, exp03="Attheya sp.", exp10="Fragilariopsis sp.", exp22="Synedropsis sp. 1", exp18 = "Synedropsis sp. 2")
slope$Species=recode_factor(slope$Species, E3="Attheya sp.", E10="Fragilariopsis sp.", E22="Synedropsis sp. 1", E18 = "Synedropsis sp. 2")
Full2=merge(Full, slope, all = TRUE)
Full2=rename(Full2, growthrate=mumax)

Full2$Carry=(Full2$K*Full2$Slope)+Full2$Intercept
ggplot(Full2, aes(Light, Slope, color=as.factor(Light))) +
    geom_point(stat="identity")+
    facet_wrap(~Species)

meank=Full2 %>%
    group_by(Species, Light, WAF) %>%
    summarize(meank=mean(Carry), sdk=sd(Carry))

    ggplot(meank, aes(WAF, meank)) +
    geom_point(stat="identity") +geom_line(stat="identity")+
    facet_wrap(~Species+Light, scales="free")

    ggplot(meank, aes(WAF, meank, color=as.factor(Light))) +
    geom_point(stat="identity") +geom_line(stat="identity")+
    geom_errorbar(aes(ymin=meank-sdk, ymax=meank+sdk), width=.2) +
    facet_wrap(~Species)

    ggplot(meank, aes(WAF, meank, color=as.factor(Light))) +
                      geom_point(stat="identity") +geom_smooth(span=0.8, se=F)+
                      geom_errorbar(aes(ymin=meank-sdk, ymax=meank+sdk), width=.2) +
                      labs(y=expression(Cells~ml^"-1"),
                           x=expression(Water~Accommodated~Fraction),
                           color=expression(Irradiance~"("*paste(mu, "mols")~"photons"~m^"-2"~s^"-1"*")" ))+
                      scale_x_continuous(labels = function(x) paste0(x, "%"))+
                      scale_color_viridis_d(option="C")+
                      facet_wrap(~Species)+theme_few()+theme(axis.title = element_text(size=20, face="bold"),
                                                             plot.title = element_text(hjust = 0.5, size=20, face="bold"),
                                                             strip.text.x = element_text(size = 15, face="bold.italic"),
                                                             axis.text = element_text(face="bold", size=10),
                                                             legend.title=element_text(face="bold", size=15),
                                                             legend.title.align=0.5,
                                                             legend.text=element_text(face="bold", size=13),
                                                             legend.position = "bottom")+
                      guides(color = guide_legend(override.aes = list(size=2)), size = guide_legend(title.position="top", title.hjust = 0.5))
