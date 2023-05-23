library("ggplot2")
library("readxl")
library("dplyr")
library(tidyverse)
library(data.table)
library(growthrates)
library(ggthemes)
library(drc)
library(bmd)
library(sandwich)
library(egg)

Oil_conc=read_excel("./Oil_Quant_Summarized.xlsx")
Oil_conc %>%
group_by(Experiment, Species, Test) %>%
summarize(conc=mean(Cal_Conc), sd=sd(Cal_Conc))

#Set to measured values of petroleum hydrocarbons
PTPH=9.66
PRRO=0.554
PDRO=2.99

setwd("./Master Data Sheet")
Wells=read_csv("./u_k_properror.csv")
Wells$Species=as.factor(Wells$Species)
Wells$Species=recode_factor(Wells$Species, "Attheya sp."= "A. septentrionalis", "Fragilariopsis sp."="F. cylindrus", "Synedropsis sp. 1"="S. hyperborea str. 1", "Synedropsis sp. 2" ="S. hyperborea str. 2" )
#Calling it wells because this contains growth and carrying capacity on a well by well basis.
Wells=rename(Wells, u=mumax, k=K)
Wells=Wells%>%
group_by(Species, Light, WAF, Row)%>%
mutate(uwellerror=sqrt(sumabserr2u), kwellerror=sqrt(sumabserr2k))

Wells=Wells%>%
group_by(Species, Light, WAF, Row)%>%
mutate(uwellerror=sqrt(sumabserr2u), kwellerror=sqrt(sumabserr2k))%>%
arrange(Light) #must arrange otherwise the order for the plot is out of wack.
cust6=c("#0D0887FF", "#6300A7FF", "#A62098FF", "#E76F5AFF", "#FDAD32FF") #Cust 6 wins for now.
par(mar=c(5,5.5,1,1))
Fcol=c("#0D0887FF", "#6300A7FF", "#A62098FF", "#E76F5AFF", "#FDAD32FF")
Fpch=c(1,2,3,4,5)
Flt=c(1,2,3,4,5)

Wells$WAF=round(Wells$WAF/100*PTPH, 1)

#GROWTH determine fit by scores and by interpretation of what may be appropriate. Drop light levels when necessary.
Fragmod=drm(formula = u ~ WAF, Light,
            data = Wells %>% filter(Species=="F. cylindrus"),
            fct = CRS.4b())

Syn2mod= drm(formula = u ~ WAF, Light,
data = Wells %>% filter(Species=="S. hyperborea str. 2") %>% filter( Light %in% c("3","10", "50", "125")),
fct = CRS.4c(), robust="median")
#Must drop L20

S2lin20=lm(uavg~WAF, Wells%>% filter(Species=="S. hyperborea str. 2" & Light %in% c("20") & Row=="1"))

Syn1mod= drm(formula = u ~ WAF, Light,
        data = Wells %>% filter(Species=="S. hyperborea str. 1")%>% filter( Light %in% c("3","20", "50", "125")),
        fct = CRS.4b(), robust="median")
S1lin10=lm(uavg~WAF, Wells%>% filter(Species=="S. hyperborea str. 1" & Light %in% c("10")& Row=="1"))

Atmod= drm(formula = u ~ WAF, Light,
          data = Wells %>% filter(Species=="A. septentrionalis")%>% filter( Light %in% c("3","10","20", "125")),
          fct = W1.3())
#Need to drop 50 and use a linear regression
Atlin50=lm(uavg~WAF, Wells%>% filter(Species=="A. septentrionalis" & Light %in% c("50") & Row=="1"))

#Repeat for maximum cell concentration
Fragmodk=drm(formula = k ~ WAF, Light,
             data = Wells %>% filter(Species=="F. cylindrus"),
             fct = cedergreen(fixed = c(NA, 9458, NA, NA, NA), names = c("b", "c", "d", "e", "f"), alpha=0.25))

Atmodk=drm(formula = k ~ WAF, Light,
                      data = Wells %>% filter(Species=="A. septentrionalis"),
                      fct = LL.4(fixed = c(NA, 4440, NA, NA), names = c("b", "c", "d", "e")))

s1modk=drm(formula = k ~ WAF, Light,
                      data = Wells %>% filter(Species=="S. hyperborea str. 1"),
                      fct = W1.4(fixed = c(NA, 2123, NA, NA), names = c("b", "c", "d", "e")))
s2modk=drm(formula = k ~ WAF, Light,
                      data = Wells %>% filter(Species=="S. hyperborea str. 2"),
                      fct = W1.4(fixed = c(NA, 11470, NA, NA), names = c("b", "c", "d", "e")))

Fr=Wells%>% filter(Species=="F. cylindrus")
At=Wells%>% filter(Species=="A. septentrionalis")
s1=Wells%>% filter(Species=="S. hyperborea str. 1")
s2=Wells%>% filter(Species=="S. hyperborea str. 2")

s2kfits <- expand.grid(WAF=exp(seq(log(max(Wells$WAF/100)), log(max(Wells$WAF)), length=1000)),Light=c("3","10","20","50","125"))
        pm <- predict(s2modk, newdata=s2kfits, interval="confidence")
        s2kfits$p <- pm[,1]
        s2kfits$pmin <- pm[,2]
        s2kfits$pmax <- pm[,3]
        s2kfits$Species <- "S. hyperborea str. 2"
s2ufits <- expand.grid(WAF=exp(seq(log(max(Wells$WAF/100)), log(max(Wells$WAF)), length=1000)),Light=c("3","10","50","125"))
        pm <- predict(Syn2mod, newdata=s2ufits, interval="confidence")
        s2ufits$p <- pm[,1]
        s2ufits$pmin <- pm[,2]
        s2ufits$pmax <- pm[,3]
        s2ufits$Species <- "S. hyperborea str. 2"
s2ulinfits=expand.grid(WAF=exp(seq(log(max(Wells$WAF/100)), log(max(Wells$WAF)), length=1000)),Light=c("20"))
  pm <- predict(S2lin20, newdata=s2ulinfits, interval="confidence")
  s2ulinfits$p <- pm[,1]
  s2ulinfits$pmin <- pm[,2]
  s2ulinfits$pmax <- pm[,3]
  s2ulinfits$Species <- "S. hyperborea str. 2"

s1kfits <- expand.grid(WAF=exp(seq(log(max(Wells$WAF/100)), log(max(Wells$WAF)), length=1000)),Light=c("3","10","20","50","125"))
        pm <- predict(s1modk, newdata=s1kfits, interval="confidence")
        s1kfits$p <- pm[,1]
        s1kfits$pmin <- pm[,2]
        s1kfits$pmax <- pm[,3]
        s1kfits$Species <- "S. hyperborea str. 1"
s1ufits <- expand.grid(WAF=exp(seq(log(max(Wells$WAF/100)), log(max(Wells$WAF)), length=1000)),Light=c("3","20","50","125"))
        pm <- predict(Syn1mod, newdata=s1ufits, interval="confidence")
        s1ufits$p <- pm[,1]
        s1ufits$pmin <- pm[,2]
        s1ufits$pmax <- pm[,3]
        s1ufits$Species <- "S. hyperborea str. 1"
s1ulinfits=expand.grid(WAF=exp(seq(log(max(Wells$WAF/100)), log(max(Wells$WAF)), length=1000)),Light=c("10"))
  pm <- predict(S1lin10, newdata=s1ulinfits, interval="confidence")
  s1ulinfits$p <- pm[,1]
  s1ulinfits$pmin <- pm[,2]
  s1ulinfits$pmax <- pm[,3]
  s1ulinfits$Species <- "S. hyperborea str. 1"


Atkfits <- expand.grid(WAF=exp(seq(log(max(Wells$WAF/100)), log(max(Wells$WAF)), length=1000)),Light=c("3","10","20","50","125"))
        pm <- predict(Atmodk, newdata=Atkfits, interval="confidence")
        Atkfits$p <- pm[,1]
        Atkfits$pmin <- pm[,2]
        Atkfits$pmax <- pm[,3]
        Atkfits$Species <- "A. septentrionalis"
Atufits <- expand.grid(WAF=exp(seq(log(max(Wells$WAF/100)), log(max(Wells$WAF)), length=1000)),Light=c("3","10","20","125"))
        pm <- predict(Atmod, newdata=Atufits, interval="confidence")
        Atufits$p <- pm[,1]
        Atufits$pmin <- pm[,2]
        Atufits$pmax <- pm[,3]
        Atufits$Species <- "A. septentrionalis"
Atulinfits=expand.grid(WAF=exp(seq(log(max(Wells$WAF/100)), log(max(Wells$WAF)), length=1000)),Light=c("50"))
  pm <- predict(Atlin50, newdata=Atulinfits, interval="confidence")
  Atulinfits$p <- pm[,1]
  Atulinfits$pmin <- pm[,2]
  Atulinfits$pmax <- pm[,3]
  Atulinfits$Species <- "A. septentrionalis"

Frufits <- expand.grid(WAF=exp(seq(log(max(Wells$WAF/100)), log(max(Wells$WAF)), length=1000)),Light=c("3","10","20","50","125"))
        pm <- predict(Fragmod, newdata=Frufits, interval="confidence")
        Frufits$p <- pm[,1]
        Frufits$pmin <- pm[,2]
        Frufits$pmax <- pm[,3]
        Frufits$Species <- "F. cylindrus"
Frkfits <- expand.grid(WAF=exp(seq(log(max(Wells$WAF/100)), log(max(Wells$WAF)), length=1000)),Light=c("3","10","20","50","125"))
        pm <- predict(Fragmodk, newdata=Frkfits, interval="confidence")
        Frkfits$p <- pm[,1]
        Frkfits$pmin <- pm[,2]
        Frkfits$pmax <- pm[,3]
        Frkfits$Species <- "F. cylindrus"

linfits=list(Atulinfits, s1ulinfits, s2ulinfits)
linufits= Reduce(function(x, y) merge(x, y, all=TRUE), linfits)
linufits$measure="u"

ufitslist <- list(Frufits, s2ufits, s1ufits, Atufits)
ufits= Reduce(function(x, y) merge(x, y, all=TRUE), ufitslist)
ufits$measure="u"
kfitslist<- list(s2kfits, s1kfits, Frkfits, Atkfits)
kfits= Reduce(function(x, y) merge(x, y, all=TRUE), kfitslist)
kfits$measure="k"

ukfits=merge(ufits, kfits, all=TRUE)
Wells$WAF[Wells$WAF==0]=max(Wells$WAF/100)

legend.title=expression(Irradiance~"("*paste(mu, "mols")~"photons"~m^"-2"~s^"-1"*")" )

#Facet K
kplot=ggplot(ukfits %>% filter(measure=="k"), aes(WAF, p/100000, color=Light))+
    theme_few() +
    labs(y=expression(Concentration~"(cells x 10"^"5"*ml^"-1"*")"),
x=expression(Total~Petroleum~Hydrocarbons~"("*mg~l^"-1"*")"))+
    scale_color_manual( legend.title, values=cust6)+scale_fill_manual(legend.title, values=cust6)+
    geom_line(aes(WAF, p/100000), size=1.2)+
    #  scale_x_log10()+
    # scale_x_continuous(labels = function(x) paste0(x, "%"))+
    ylim(c(0,NA))+
    geom_ribbon(aes(ymin=pmin/100000, ymax=pmax/100000, fill=Light), alpha=0.1, colour = NA)+
    theme(axis.title = element_text(size=20, face="bold"),
          plot.title = element_text(hjust = 0.5, size=20, face="bold"),
          # strip.text.x = element_text(size = 15, face="italic"),
          strip.background = element_blank(),
         strip.text.x = element_blank(),
          axis.text = element_text(face="bold", size=12, color="black"),
          legend.title=element_text(face="bold", size=15),
          legend.title.align=0.5,
          legend.text=element_text(face="bold", size=13),
          legend.position = "bottom",
          legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'),
        panel.border = element_rect(fill=NA, colour = "black", size=1))+
    guides(color = guide_legend(override.aes = list(size=2)), size = guide_legend(title.position="top", title.hjust = 0.5))+
    geom_point(data = Wells, mapping = aes( x=WAF, y = kavg/100000, color = as.factor(Light)))+
    geom_errorbar(data=Wells, mapping=aes(x=WAF, y=kavg/100000, ymin=kavg/100000-kerror/100000, ymax=kavg/100000+kerror/100000, color= as.factor(Light)), width=0)+
    facet_grid(~Species)
    k2=tag_facet(kplot, tag_pool = c("A", "B","C","D"), open="", close="", size=5,  y=4.5, vjust=0.01)

uplot=ggplot(ukfits %>% filter(measure=="u"), aes(WAF, p, color=Light))+
            theme_few() +
            labs(y=expression(Growth~"(d"^"-1"*")"),
x=expression(Total~Petroleum~Hydrocarbons~"("*mg~l^"-1"*")"))+
            scale_color_manual(legend.title, values=cust6)+scale_fill_manual(legend.title, values=cust6)+
            geom_line(data=linufits, aes(WAF,p), size=1, linetype="dashed", color="black")+
            geom_line(aes(WAF, p), size=1.2)+
            #scale_x_log10()+
            ylim(c(0,NA))+
            geom_ribbon(aes(ymin=pmin, ymax=pmax, fill=Light), alpha=0.1, colour = NA)+
            theme(axis.title = element_text(size=20, face="bold"),
                  plot.title = element_text(hjust = 0.5, size=20, face="bold"),
                  # strip.text.x = element_text(size = 15, face="italic"),
                  strip.background = element_blank(),
                 strip.text.x = element_blank(),
                  axis.text = element_text(face="bold", size=12, color="black"),
                  legend.title=element_text(face="bold", size=15),
                  legend.title.align=0.5,
                  legend.text=element_text(face="bold", size=13),
                  legend.position = "bottom",
                  legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'),
                  panel.border = element_rect(fill=NA, colour = "black", size=1))+
            guides(color = guide_legend(override.aes = list(size=2)), size = guide_legend(title.position="top", title.hjust = 0.5))+
geom_point(data = Wells, mapping = aes( x=WAF, y = uavg, color = as.factor(Light)))+
geom_errorbar(data=Wells, mapping=aes(x=WAF, y=uavg, ymin=uavg-uerror, ymax=uavg+uerror, color= as.factor(Light)), width=0)+
facet_grid(~Species)

up2=tag_facet(uplot, tag_pool = c("A", "B","C","D"), open="", close="", size=5, y=0.52, vjust=0.01)

#############Panel Combo
library(ggpubr)
library(gridExtra)
library(grid)

upan=ggplot(ukfits %>% filter(measure=="u"), aes(WAF, p, color=Light))+
    theme_few() +
    labs(y=expression(Growth~"(day"^"-1"*")"),
    x=expression(TPH~"("*mg~L^"-1"*")"))+
    scale_color_manual(legend.title, values=cust6)+scale_fill_manual(legend.title, values=cust6)+
    geom_line(aes(WAF, p), size=1.2)+
    ylim(c(0,NA))+
    geom_ribbon(aes(ymin=pmin, ymax=pmax, fill=Light), alpha=0.1, colour = NA)+
    theme(axis.title = element_text(size=15, face="bold"),
          axis.title.x=element_blank(),
          plot.title = element_text(hjust = 0.5, size=20, face="bold"),
          strip.text.x = element_text(size = 15, face="bold.italic"),
          axis.text.x=element_blank(),
          strip.text.y = element_text(size = 15, face="bold.italic"),
          axis.text = element_text(face="bold", size=12),
          legend.position="none" )+
  geom_errorbar(data=Wells, mapping=aes(x=WAF, y=uavg, ymin=uavg-uerror, ymax=uavg+uerror, color= as.factor(Light)), width=0)+
      geom_point(data = Wells, mapping = aes( x=WAF, y = uavg, color = as.factor(Light)))+
    facet_grid(~Species, scales ="free")

kpanwo=ggplot(ukfits %>% filter(measure=="k"), aes(WAF, p/100000, color=Light))+
                        theme_few() +
                        labs(y=expression(Growth~"(day"^"-1"*")"),
                             x=expression(Irradiance~"("*paste(mu, "mols")~"photons"~m^"-2"~s^"-1"*")" ))+
                        scale_color_manual( legend.title, values=cust6)+
                        scale_fill_manual(legend.title, values=cust6)+
                        geom_line(aes(WAF, p/100000), size=1.2)+
                        labs(y=expression(Cells~"(x 10"^"5"~ml^"-1"*")"),
                             x=expression(TPH~"("*mg~L^"-1"*")"))+
                        scale_y_continuous(labels = scales::number_format(accuracy = 0.1,
                                 decimal.mark = '.'))+
                        geom_ribbon(aes(ymin=pmin/100000, ymax=pmax/100000, fill=Light), alpha=0.1, colour = NA)+
                        theme(axis.title = element_text(size=15, face="bold"),
                              plot.title = element_text(hjust = 0.5, size=20, face="bold"),
                              strip.text.x = element_blank(),
                              strip.text.y = element_text(size = 15, face="bold.italic"),
                              axis.text = element_text(face="bold", size=12),
                              legend.position="none")+
                          guides(color = guide_legend(override.aes = list(size=2)), size = guide_legend(title.position="top", title.hjust = 0.5))+
                        geom_errorbar(data=Wells, mapping=aes(x=WAF, y=kavg/100000, ymin=kavg/100000-kerror/100000, ymax=kavg/100000+kerror/100000, color= as.factor(Light)), width=0)+
                        geom_point(data = Wells, mapping = aes( x=WAF, y = kavg/100000, color = as.factor(Light)))+
                        facet_grid(~Species, scales ="free")

kpan=ggplot(ukfits %>% filter(measure=="k"), aes(WAF, p/100000, color=Light))+
                        theme_few() +
                        labs(y=expression(Growth~"(day"^"-1"*")"),
                             x=expression(Irradiance~"("*paste(mu, "mols")~"photons"~m^"-2"~s^"-1"*")" ))+
                        scale_color_manual( legend.title, values=cust6)+
                        scale_fill_manual(legend.title, values=cust6)+
                        geom_line(aes(WAF, p/100000), size=1.2)+
                        labs(y=expression(Cells~"(x 10"^"5"~ml^"-1"*")"),
                             x=expression(TPH~"("*mg~L^"-1"*")"))+
                        scale_y_continuous(labels = scales::number_format(accuracy = 0.1,
                                 decimal.mark = '.'))+
                        geom_ribbon(aes(ymin=pmin/100000, ymax=pmax/100000, fill=Light), alpha=0.1, colour = NA)+
                        theme(axis.title = element_text(size=15, face="bold"),
                              plot.title = element_text(hjust = 0.5, size=20, face="bold"),
                              strip.text.x = element_blank(),
                              strip.text.y = element_text(size = 15, face="bold.italic"),
                              axis.text = element_text(face="bold", size=12),

                         legend.title=element_text(face="bold", size=15),
                         legend.title.align=0.5,
                         legend.text=element_text(face="bold", size=13),
                         legend.position = "bottom")+

                          guides(color = guide_legend(override.aes = list(size=2)), size = guide_legend(title.position="top", title.hjust = 0.5))+
                        geom_errorbar(data=Wells, mapping=aes(x=WAF, y=kavg/100000, ymin=kavg/100000-kerror/100000, ymax=kavg/100000+kerror/100000, color= as.factor(Light)), width=0)+
                        geom_point(data = Wells, mapping = aes( x=WAF, y = kavg/100000, color = as.factor(Light)))+
                        facet_grid(~Species)



get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

kpan_legend <- get_legend(kpan)

Panel1=ggarrange(upan,kpanwo, labels=c("A","B"), ncol=1, nrow=2)

Panel2=grid.arrange(Panel1,kpan_legend,
          ncol = 1, nrow = 2, heights=c(11, 1))

#IC10 and IC50
EDFr10 <- as.data.frame(ED(Fragmod, c(10,50), interval = "delta"))
EDFr10=cbind(EDFr10, read.table(text=row.names(EDFr10), sep=":",
header=FALSE, col.names = paste0("Light", 1:3), stringsAsFactors=FALSE))

EDAT <- as.data.frame(ED(Atmod, c(10,50),  interval = "delta"))
EDAT=cbind(EDAT, read.table(text=row.names(EDAT), sep=":",
                     header=FALSE, col.names = paste0("Light", 1:3), stringsAsFactors=FALSE))

EDS1 <- as.data.frame(ED(Syn1mod, c(10, 50), interval = "delta"))
EDS1=cbind(EDS1, read.table(text=row.names(EDS1), sep=":", header=FALSE,
col.names = paste0("Light", 1:3), stringsAsFactors=FALSE))

EDS2 <- as.data.frame(ED(Syn2mod, c(10, 50), interval = "delta"))
EDS2=cbind(EDS2, read.table(text=row.names(EDS2), sep=":",
header=FALSE, col.names = paste0("Light", 1:3), stringsAsFactors=FALSE))

EDAT$Light1[EDAT$Light1=="e"]="A. septentrionalis"
EDFr10$Light1[EDFr10$Light1=="e"]="F. cylindrus"
EDS1$Light1[EDS1$Light1=="e"]="S. hyperborea str. 1"
EDS2$Light1[EDS2$Light1=="e"]="S. hyperborea str. 2"


IC=merge(EDAT,EDFr10, all=T)
IC=merge(IC,EDS1, all=T)
IC=merge(IC,EDS2, all=T)
colnames(IC)[1] <- 'WAF'
colnames(IC)[5] <- 'Species'
colnames(IC)[6] <- 'Light'
colnames(IC)[7] <- 'Threshold'
IC$Species=as.factor(IC$Species)
IC$BV=IC$Species

IC = IC %>% filter(WAF<max(Wells$WAF))

IC_conc=IC %>%   dplyr::select(WAF, Species, Light, Threshold) %>% pivot_wider(names_from = Light, values_from = WAF)
IC_conc$conc="conc"
IC_SE=IC %>%   dplyr::select(`Std. Error`, Species, Light, Threshold) %>% pivot_wider(names_from = Light, values_from = `Std. Error`)
IC_SE$SE="SE"
IC_sum=merge(IC_conc,IC_SE, all=T)
IC_sum[,3:7]=round(IC_sum[,3:7],2)
library("writexl")
IC$Threshold=as.factor(IC$Threshold)
Supp_IC_Plot_u=
ggplot(IC, aes(Light, WAF, color=Threshold))+
    geom_errorbar(aes(ymin=ifelse(Lower < 0, 0, Lower), ymax=Upper), color="black", width=10)+
    geom_point(stat="identity", size=2.5)+
    ylim(0,10)+
    xlim(-5,130)+
    labs(y=expression(TPH~"("*mg~L^"-1"*")"),
         x=expression(Irradiance~"("*paste(mu, "mols")~"photons"~m^"-2"~s^"-1"*")" ))+
    facet_grid(~Species, scales="free")+theme_few()+theme(axis.title = element_text(size=20, face="bold"),
                                                          plot.title = element_text(hjust = 0.5, size=20, face="bold"),
                                                          strip.text.x = element_text(size = 15, face="bold"),
                                                          axis.text = element_text(face="bold", size=12),
                                                          legend.title=element_text(face="bold", size=15),
                                                          legend.title.align=0.5,
                                                          legend.text=element_text(face="bold", size=13),
                                                          legend.position = "bottom")+
    guides(color = guide_legend(override.aes = list(size=2), title="Inhibitory Concentration (%)"), size = guide_legend(title.position="top", title.hjust = 0.5))

###########################Maximum Concentration
EDFr10k <- as.data.frame(ED(Fragmodk, c(10,50), interval = "delta"))
EDFr10k=cbind(EDFr10k, read.table(text=row.names(EDFr10k), sep=":",
header=FALSE, col.names = paste0("Light", 1:3), stringsAsFactors=FALSE))

 EDATk <- as.data.frame(ED(Atmodk, c(10,50),  interval = "delta", clevel=c(50,125) ))
 EDATk=cbind(EDATk, read.table(text=row.names(EDATk), sep=":",
                      header=FALSE, col.names = paste0("Light", 1:3), stringsAsFactors=FALSE))

EDS1k <- as.data.frame(ED(s1modk, c(10, 50), interval = "delta"))
EDS1k=cbind(EDS1k, read.table(text=row.names(EDS1k), sep=":", header=FALSE,
col.names = paste0("Light", 1:3), stringsAsFactors=FALSE))

EDS2k <- as.data.frame(ED(s2modk, c(10, 50), interval = "delta"))
EDS2k=cbind(EDS2k, read.table(text=row.names(EDS2k), sep=":",
header=FALSE, col.names = paste0("Light", 1:3), stringsAsFactors=FALSE))

 EDATk$Light1[EDATk$Light1=="e"]="A. septentrionalis"
EDFr10k$Light1[EDFr10k$Light1=="e"]="F. cylindrus"
 EDS1k$Light1[EDS1k$Light1=="e"]="S. hyperborea str. 1"
EDS2k$Light1[EDS2k$Light1=="e"]="S. hyperborea str. 2"


IC=merge(EDFr10k,EDS2k, all=T)
 IC=merge(IC,EDS1k, all=T)
 IC=merge(IC,EDATk, all=T)
colnames(IC)[1] <- 'WAF'
colnames(IC)[5] <- 'Species'
colnames(IC)[6] <- 'Light'
colnames(IC)[7] <- 'Threshold'
IC$Species=as.factor(IC$Species)
IC$BV=IC$Species

ICk = IC %>% filter(WAF<max(Wells$WAF))

ICk_conc=ICk %>%   dplyr::select(WAF, Species, Light, Threshold) %>% pivot_wider(names_from = Light, values_from = WAF)
ICk_conc$conc="conc"
ICk_SE=ICk %>%   dplyr::select(`Std. Error`, Species, Light, Threshold) %>% pivot_wider(names_from = Light, values_from = `Std. Error`)
ICk_SE$SE="SE"
ICk_sum=merge(ICk_conc,ICk_SE, all=T)
ICk_sum[,3:7]=round(ICk_sum[,3:7],2)
library("writexl")

ICk$Threshold=as.factor(ICk$Threshold)
Supp_IC_Plot_k=
ggplot(ICk, aes(Light, WAF, color=Threshold))+
geom_errorbar(aes(ymin=ifelse(Lower < 0, 0, Lower), ymax=Upper), color="black", width=10)+
geom_point(stat="identity", size=2.5)+
ylim(0,10)+
xlim(-5,130)+
labs(y=expression(TPH~"("*mg~L^"-1"*")"),
     x=expression(Irradiance~"("*paste(mu, "mols")~"photons"~m^"-2"~s^"-1"*")" ))+
facet_grid(~Species, scales="free")+theme_few()+theme(axis.title = element_text(size=20, face="bold"),
                                                      plot.title = element_text(hjust = 0.5, size=20, face="bold"),
                                                      strip.text.x = element_text(size = 15, face="bold"),
                                                      axis.text = element_text(face="bold", size=12),
                                                      legend.title=element_text(face="bold", size=15),
                                                      legend.title.align=0.5,
                                                      legend.text=element_text(face="bold", size=13),
                                                      legend.position = "bottom")+
    guides(color = guide_legend(override.aes = list(size=2), title="Inhibitory Concentration (%)"), size = guide_legend(title.position="top", title.hjust = 0.5))
