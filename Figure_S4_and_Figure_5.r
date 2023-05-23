#Adjusting the UTC so that it is centered near 12 so that I can average over 12 hours within an actual day.
library(neonUtilities)
library(raster)
library(lubridate)
library(ggplot2)
library("dplyr")
library(tidyverse)
library(readxl)
library(readr)
library(neonUtilities)
library(raster)

# Set global option to NOT convert all character variables to factors
options(stringsAsFactors=F)

#Read in the stacked files
par30 <- readTableNEON(
  dataFile="C:./PARPAR_30min.csv",
  varFile="C:./variables_00024.csv")
hrs <- hours(12)
par30$adjTime=par30$startDateTime-hrs
par30$DOY=lubridate::yday(par30$adjTime)-1
par30$ddate=decimal_date(par30$startDateTime)
par30$adjtime2=format(par30$adjTime,"%m-%d %H:%M")


#For March to June, solar noon =11:00-11:30 adjtime. Subset the data from 0600 to 17:30 adjTime
par30.2020 <- subset(par30,
                      startDateTime >=('2020-03-01 00:00') &
                          startDateTime <=('2020-05-15 23:30'))
par30.2022 <- subset(par30,
                      startDateTime >=('2022-03-01 00:00') &
                          startDateTime <=('2022-05-15 23:30'))
# combine these two years subsets.
par30_comb=rbind(par30.2020, par30.2022)

Trimmedlight=par30_comb%>%
    mutate(hour = format(ymd_hms(as.character(par30_comb$adjTime, tz = "UTC")),'%H:%M'))%>%
    mutate(year = format(ymd_hms(as.character(par30_comb$adjTime, tz = "UTC")),'%Y'))%>%
    mutate(yearless = format(ymd_hms(as.character(par30_comb$adjTime, tz = "UTC")),'%m-%d %H:%M'))%>%

    filter(verticalPosition=="040")%>%
    group_by(DOY)%>%
    filter(hour >=('06:00') &
               hour <=('17:30'))

Trimmedlight=na.omit(Trimmedlight)

PARDAvg=Trimmedlight%>%
    dplyr::select(DOY, startDateTime, PARMean, year)%>%
    group_by(DOY)%>%
    summarize( DayMean=mean(PARMean))

PARDAvg$Date=as.Date(PARDAvg$DOY, origin = "2020-01-01")

PARDAvg$Date <- PARDAvg$Date +
  hours(11) +
  minutes(30) +
  seconds(00)


  library(ggpubr)
  library(gridExtra)
  library(grid)
  library(ggthemes)

a2020=ggplot(data=Trimmedlight %>% filter(year=="2020"), aes(adjTime, PARMean))+
theme_few() +
  labs(x=expression(Day~of~Year))+
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1),
        axis.text = element_text( colour="black", size=15),
        axis.ticks.x=element_blank())+
        ylim(0,1500)+
  geom_line( size=0.4)

b2022=ggplot(data=Trimmedlight %>% filter(year=="2022"), aes(adjTime, PARMean))+
theme_few() +
  labs(x=expression(Day~of~Year))+
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1),
        axis.text = element_text( colour="black", size=15),
        axis.ticks.x=element_blank())+
        ylim(0,1500)+
  geom_line( size=0.4)

  cavg=ggplot(data=PARDAvg, aes(DOY, DayMean))+
        theme_few() +
        ylim(0,1500)+
        labs(x=expression(Day~of~Year))+
theme(axis.title.y=element_blank(),
panel.border = element_rect(fill=NA, colour = "black", size=1),
axis.text = element_text( colour="black", size=15),
axis.title.x=element_blank(),
axis.text.x = element_blank(),
axis.title= element_text(size=15),
axis.ticks.x=element_blank())+
        geom_line( size=0.4)

#Generate top portion of multipanel irradiance data.
panel1=
ggarrange(a2020, b2022, cavg, labels=c("A","B", "C"), hjust=-5.3, vjust= 2.1, ncol=1, nrow=3)
panel2=grid.arrange(panel1,   left = textGrob(expression(Irradiance~"("*paste(mu, "mols")~"photons"~m^"-2"~s^"-1"*")" ),
rot=90,  gp = gpar(col = "black", fontsize = 20)))


  PARDAvg$Date=as.Date(PARDAvg$DOY, origin = "2022-01-01")

  PARDAvg$Date <- PARDAvg$Date +
    hours(11) +
    minutes(30) +
    seconds(00)
#Need to get rid of year in the date format so I can plot the PAR data by year.
  #Add snow depths
PARDAvg$s1=0
PARDAvg$s2=.01
PARDAvg$s3=.02
PARDAvg$s4=.03
PARDAvg$s5=.04
PARDAvg$s6=.05
PARDAvg$s7=.06
PARDAvg$s8=.07
PARDAvg$s9=.08
PARDAvg$s10=.09
PARDAvg$s11=.1
PARDAvg$s12=.11
PARDAvg$s13=.12
PARDAvg$s14=.13
PARDAvg$s15=.14
PARDAvg$s16=.15
PARDAvg$s17=.16
PARDAvg$s18=.17
PARDAvg$s19=.18
PARDAvg$s20=.19
PARDAvg$s21=.20

#IMPORT Ice thickness   SIMB Data
#These data were derived from thermistor string
SIMB=read.csv("C:./SIMB_2014.csv", header = T, sep = ",", dec = ".")
SIMB=na.omit(SIMB)
Simsub=filter(SIMB, between(DOY, 59, 113))%>% #113 is the last day at the max ice thickness.
group_by(DOY)%>%
summarize( Hi=mean(Hi))

PAR_Merge=merge(PARDAvg,Simsub,by=c("DOY"))
plot(data=PAR_Merge, Hi~DOY, type="l")
Ice_Thick=ggplot(PAR_Merge)+
    geom_line(aes(DOY, Hi), size=1.2)+scale_y_reverse()

PAR_Merge=PAR_Merge %>% gather( snow, Hs, -c(DOY, DayMean, Date, Hi))
PAR_Merge=PAR_Merge[-c(5)]

#Set attenuation coefficients
ui=-3.5*.7 #Matthes 2019
us=-11
PAR_Merge$UI_PAR=round(PAR_Merge$DayMean*exp((PAR_Merge$Hi*ui)+(PAR_Merge$Hs*us)),2)

UI2=PAR_Merge %>% dplyr::filter( Hs==(0.00))
UI3=PAR_Merge %>% dplyr::filter( Hs==(0.10))
UI4=PAR_Merge %>% dplyr::filter( Hs==(0.20))

UI5= merge(UI2, UI3, all=T)
UI6=merge(UI5, UI4, all=T)

UI_PAR=
ggplot(UI6)+theme_few() +
    geom_line(aes(DOY, UI_PAR, linetype=as.factor(Hs)), size=0.4)+
    labs(x=expression(Day~of~Year), linetype=expression(Snow~depth~"(m)"))+
    theme(axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x = element_blank(),
          legend.position = c(0.2, 0.8),
          panel.border = element_rect(fill=NA, colour = "black", size=1),
          axis.text = element_text( colour="black", size=15),
          axis.ticks.x=element_blank())
ice_thick=
ggplot(PAR_Merge)+
    geom_line(aes(DOY, Hi), size=0.4)+scale_y_reverse()+theme_few()+
    labs(x=expression(Day~of~Year),
  )+
     theme(axis.text.x = element_text(size=15),
           axis.title= element_blank(),
          panel.border = element_rect(fill=NA, colour = "black", size=1),
          axis.text = element_text( colour="black", size=15),
          )

gA <- ggplotGrob(a2020)
gB <- ggplotGrob(b2022)
gc=ggplotGrob(cavg)
gd=ggplotGrob(UI_PAR)
ge=ggplotGrob(ice_thick)

#Full 5 panel plot of irradiance, ice thickness, and snow attenuation
library(cowplot)
final=
ggdraw() +
    draw_plot(plot_grid(gA,gB,gc,gd,ge, ncol = 1, align = 'hv', labels = "AUTO", hjust=-5, vjust=1.8))+draw_label("Day of Year", x=0.47, y=  0, vjust=-0.1, angle= 0)


library(readxl)
library(phytotools)
PAR_Merge=PAR_Merge %>% filter(between(DOY, 59, 113) )

setwd("C:./Master Data Sheet")
Wells=read_csv("C:/Users/Dilliplaine/Desktop/CMI 2020/WAF_Experiment/Master Data Sheet/u_k_properror.csv")
Wells$Species=as.factor(Wells$Species)
Wells$Species=recode_factor(Wells$Species, "Attheya sp."= "A. septentrionalis", "Fragilariopsis sp."="F. cylindrus", "Synedropsis sp. 1"="S. hyperborea str. 1", "Synedropsis sp. 2" ="S. hyperborea str. 2" )
Full=Wells %>% filter(WAF=="0")

At=Full%>%filter(Species=="A. septentrionalis")
Frag=Full%>%filter(Species=="F. cylindrus")
Syn1=Full%>%filter(Species=="S. hyperborea str. 1")
Syn2=Full%>%filter(Species=="S. hyperborea str. 2")


FragWebb=fitWebb(Frag$Light,Frag$mumax)
AtWebb=fitWebb(At$Light,At$mumax)
Syn1EP=fitEP(Syn1$Light,Syn1$mumax, fitmethod="Marq")
Syn2EP=fitEP(Syn2$Light,Syn2$mumax, fitmethod="Marq")

#Growthrate Calcs
PAR_Merge$"A. septentrionalis"=with(AtWebb, {
    P <- alpha[1] * ek[1] * (1 - exp (-PAR_Merge$UI_PAR / ek[1]))})

    PAR_Merge$"S. hyperborea str. 2"=with(Syn2EP, {
      P <- PAR_Merge$UI_PAR/((1/(alpha[1]*eopt[1]^2))*PAR_Merge$UI_PAR^2+(1/ps[1]-2/(alpha[1]*eopt[1]))*PAR_Merge$UI_PAR+(1/alpha[1]))
    })

    PAR_Merge$"S. hyperborea str. 1"=with(Syn1EP,{
      P <- PAR_Merge$UI_PAR/((1/(alpha[1]*eopt[1]^2))*PAR_Merge$UI_PAR^2+(1/ps[1]-2/(alpha[1]*eopt[1]))*PAR_Merge$UI_PAR+(1/alpha[1]))
    })

    PAR_Merge$"F. cylindrus"=  with(FragWebb,{
    P <- alpha[1] * ek[1] * (1 - exp (-PAR_Merge$UI_PAR / ek[1]))})

#G-I curves with oil.
Full_oil=Wells %>% filter(WAF=="90")
AtO=Full_oil%>%filter(Species=="A. septentrionalis")
FragO=Full_oil%>%filter(Species=="F. cylindrus")
Syn1O=Full_oil%>%filter(Species=="S. hyperborea str. 1")
Syn2O=Full_oil%>%filter(Species=="S. hyperborea str. 2")

FragEPO=fitEP(FragO$Light,FragO$mumax)
AtEPO=fitEP(AtO$Light,AtO$mumax)
Syn1JPO=fitJP(Syn1O$Light,Syn1O$mumax)
Syn2JPO=fitJP(Syn2O$Light,Syn2O$mumax)

PAR_Merge$"A. septentrionalis_Oil"=with(AtEPO, {
  P <- PAR_Merge$UI_PAR/((1/(alpha[1]*eopt[1]^2))*PAR_Merge$UI_PAR^2+(1/ps[1]-2/(alpha[1]*eopt[1]))*PAR_Merge$UI_PAR+(1/alpha[1]))
})

PAR_Merge$"F. cylindrus_Oil"=with(FragEPO, {
  P <- PAR_Merge$UI_PAR/((1/(alpha[1]*eopt[1]^2))*PAR_Merge$UI_PAR^2+(1/ps[1]-2/(alpha[1]*eopt[1]))*PAR_Merge$UI_PAR+(1/alpha[1]))
})

PAR_Merge$"S. hyperborea str. 2_Oil"=   with(Syn2JPO,{
  P <- alpha[1]*ek[1]*tanh(PAR_Merge$UI_PAR/ek[1])})

PAR_Merge$"S. hyperborea str. 1_Oil"=   with(Syn1JPO,{
    P <- alpha[1]*ek[1]*tanh(PAR_Merge$UI_PAR/ek[1])})


Growthrates=PAR_Merge[,c(1,5,7:14)]

#Growth calculations
Growthrates=read_excel("C:/Users/Dilliplaine/Desktop/CMI 2020/Manuscript 1/Growth_Rates_predicted.xlsx")
library(tidyr)
 Cells <- Growthrates %>% group_by(Hs) %>% mutate(across(-DOY, ~1000*cumprod(exp(as.numeric(.x)*.9))))
Cells=Cells%>%pivot_longer(cols=c("A. septentrionalis_Oil", "F. cylindrus_Oil", "S. hyperborea str. 2_Oil", "S. hyperborea str. 1_Oil", "A. septentrionalis", "F. cylindrus", "S. hyperborea str. 2", "S. hyperborea str. 1"), names_to="Species", values_to="Cells")
Cells=Cells %>% separate(col=Species, sep = "_", into = c("Species", "WAF"))
Cells$WAF[is.na(Cells$WAF)] <- 0.0
Cells$WAF[Cells$WAF == "Oil"] <- 8.7
library(paletteer)
Disc=paletteer_d("ggsci::default_jama", n=4)

library(egg)
new_labels <- c("0"="A", "8.7"="B")

#Area plot of modeled relative abundance
NEON_Area_D113=
ggplot(Cells %>% filter(DOY==113 & Hs <0.21), aes(fill=Species, y=Cells, x=Hs)) +
    geom_area(position="fill", stat="identity",linetype = 1, size =0.5 ,colour="black" )+
    scale_color_manual( values=Disc)+
    scale_fill_manual( values=Disc)+
    labs(y=expression(Relative~abundance),
         x=expression(Snow~depth~"(m)"))+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    #  ylim(c(0,NA))+
    theme(axis.title = element_text(size=20, face="bold"),
          plot.title = element_text(hjust = 0.5, size=20, face="bold"),
          strip.text.x = element_text(size = 20, face="bold", hjust=0),
          #strip.text.x = element_blank(),
          strip.text.y = element_text(size = 15, face="bold.italic"),
          axis.text = element_text(face="bold", size=15, colour = "black"),
          legend.text = element_text( face = "italic", size = 20),
          legend.title = element_blank(),
          axis.text.x = element_text(vjust = -1),
          legend.position="bottom",
          panel.spacing = unit(2.5, "lines"),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.title.x = element_text(vjust=-1),
          plot.margin = margin(8,18,0,2))+
    facet_grid(~WAF,labeller = labeller(WAF = new_labels))

#################################
#Cells to Biovolume
BV=data.frame(Species=c("A. septentrionalis", "F. cylindrus", "S. hyperborea str. 1", "S. hyperborea str. 2"),
BV=c(186, 142, 117, 80))
BV=BV %>% group_by(Species) %>% mutate(C=10^(log(BV, 10)*0.58+0.77))

Cells=merge(Cells,BV,by=c("Species"))
Cells$Carbonmg=(Cells$C*Cells$Cells)/1e+9

Cellsx=Cells %>% dplyr::filter(DOY==113 & Hs==c(0.00))
Cellsy=Cells %>% dplyr::filter(DOY==113 & Hs==c(0.10))
Cellsz= Cells %>% dplyr::filter(DOY==113 & Hs==c(0.2))
Cells2=  merge(Cellsx, Cellsy, all=T)
Cells2=merge(Cells2, Cellsz, all=T)
Cells2$Hs[Cells2$Hs == "0"] <- "C"
Cells2$Hs[Cells2$Hs == "0.1"] <- "D"
Cells2$Hs[Cells2$Hs == "0.2"] <- "E"

Cells4=Cells2 %>% group_by(WAF, Hs) %>% summarize(Carbonmg=sum(Carbonmg))

#Bar plot of carbon
Bar_Carbon=
ggplot(Cells4, aes( y=Carbonmg, x=WAF)) +
    geom_bar(stat="identity")+
    labs(y=expression(Biomass~"(mg"~"C)"),
        x=expression(TPH~"("*mg~L^"-1"*")"))+
    theme(axis.title = element_text(face="bold", size=20),
          plot.title = element_text(hjust = 0.5, size=20),
axis.text = element_text(face="bold", size=15, colour = "black"),
          legend.title = element_blank(),
          strip.text = element_text(hjust = -0),
          strip.text.x = element_text(size = 20, face="bold", hjust=0),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line.x.bottom= element_line(color="black"),
          axis.line.x.top= element_line(color="black"),
          axis.line.y.left= element_line(color="black"),
          axis.line.y.right= element_blank(),
          axis.text.x = element_text(vjust = -1),
          axis.title.x = element_text(vjust=-1),
          strip.background = element_blank(),
          plot.margin = margin(2,2,2,2))+
    facet_wrap(~Hs, ncol=1, scales="free_y")


#Combined figure
Figure_5=ggarrange(NEON_Area_D113, Bar_Carbon,
          ncol = 2, widths = c(4, 1))
