library("ggplot2")
library("readxl")
library("dplyr")
library(tidyverse)
library(data.table)
library(growthrates)
library(ggthemes)
library(phytotools)
library(egg)
Wells=read_csv("./u_k_properror.csv")
Wells$Species=as.factor(Wells$Species)
Wells$Species=recode_factor(Wells$Species, "Attheya sp."= "A. septentrionalis", "Fragilariopsis sp."="F. cylindrus", "Synedropsis sp. 1"="S. hyperborea str. 1", "Synedropsis sp. 2" ="S. hyperborea str. 2" )
Full=Wells %>% filter(WAF=="0")

Full2=Wells %>% filter(WAF=="0" | WAF=="90")
library(paletteer)
Disc=paletteer_d("ggsci::default_jama", n=4)

At=Full%>%filter(Species=="A. septentrionalis")
Frag=Full%>%filter(Species=="F. cylindrus")
Syn1=Full%>%filter(Species=="S. hyperborea str. 1")
Syn2=Full%>%filter(Species=="S. hyperborea str. 2") #If we drop light 20, we get the photoinhibition.

#Assess fit using AIC and ssr
FragWebb=fitWebb(Frag$Light,Frag$mumax)
AtWebb=fitWebb(At$Light,At$mumax)
Syn1EP=fitEP(Syn1$Light,Syn1$mumax, fitmethod="Marq")
Syn2EP=fitEP(Syn2$Light,Syn2$mumax, fitmethod="Marq") #If including L20, I cannot use it because it errors out.

rm(PI)
PI <- data.frame(Light = seq(3,125, by = 1))
PI$"S. hyperborea str. 2"=with(Syn2EP, {
  P <- PI$Light/((1/(alpha[1]*eopt[1]^2))*PI$Light^2+(1/ps[1]-2/(alpha[1]*eopt[1]))*PI$Light+(1/alpha[1]))
})

PI$"S. hyperborea str. 1"=with(Syn1EP,{
  P <- PI$Light/((1/(alpha[1]*eopt[1]^2))*PI$Light^2+(1/ps[1]-2/(alpha[1]*eopt[1]))*PI$Light+(1/alpha[1]))
})

PI$"F. cylindrus"=  with(FragWebb,{
P <- alpha[1] * ek[1] * (1 - exp (-PI$Light / ek[1]))})

PI$"A. septentrionalis"=with(AtWebb, {
P <- alpha[1] * ek[1] * (1 - exp (-PI$Light / ek[1]))})

 PI=PI %>% gather(Species, value, -Light)

GvI=ggplot()+
    theme_few() +
    labs(y=expression(Growth~"(day"^"-1"*")"),
         x=expression(Irradiance~"("*paste(mu, "mols")~"photons"~m^"-2"~s^"-1"*")" ))+
    geom_line(data=PI, aes(Light, value, color=Species), size=1.2)+
    theme(axis.title = element_text(size=30, face="bold"),
          plot.title = element_text(hjust = 0.5, size=25, face="bold"),
          strip.text.x = element_text(size = 20, face="bold.italic"),
          axis.text = element_text(face="bold", size=25),
          legend.title=element_text(face="bold", size=25),
          legend.title.align=0.5,
          legend.text=element_text(face="italic", size=20),
          plot.margin=grid::unit(c(1,0,0,0), "cm"),
          legend.position = c(0.8, 0.25))+
    guides(color = guide_legend(override.aes = list(size=2)), size = guide_legend(title.position="top", title.hjust = 0.5))


PVI_facet=
ggplot()+
    theme_few() +
    labs(y=expression(Growth~"(d"^"-1"*")"),
         x=expression(Irradiance~"("*paste(mu, "mols")~"photons"~m^"-2"~s^"-1"*")" ))+
    geom_line(data=PI, aes(Light, value, color=Species), size=1.2)+
    scale_color_manual( values=Disc)+
    scale_fill_manual( values=Disc)+
    geom_point(data=Wells%>% filter(WAF=="0" ), aes(Light, uavg))+
    geom_errorbar(data=Wells%>% filter(WAF=="0" ),aes(x=Light, ymin=uavg-uerror, ymax=uavg+uerror), width=.4, color="black") +

    facet_grid(~Species)+theme_few()+theme(axis.title = element_text(size=20, face="bold"),
                                           plot.title = element_text(hjust = 0.5, size=25, face="bold"),
                                           axis.text = element_text(face="bold", size=12, colour="black"),
                                           strip.background = element_blank(),
                                          strip.text.x = element_blank(),
                                          legend.title=element_text(face="bold", size=15),
                                           legend.title.align=0.5,
                                           legend.text=element_text(face="italic", size=15),
                                           legend.position = "none",
                                         panel.border = element_rect(fill=NA, colour = "black", size=1))+
    guides(color = guide_legend(override.aes = list(size=2)), size = guide_legend(title.position="top", title.hjust = 0.5)) #+
    PVI2=tag_facet(PVI_facet, tag_pool = c("A", "B","C","D"), open="", close="", size=5)
