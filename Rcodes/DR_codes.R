# -----------------------------------------------------------------------------------------------------#
# Codes for replicating the paper Who has been saved more? Trends in gender differences in mortality
# through the revivorship approach.
# author: Vanessa di Lego
# Wittgenstein Centre for Demography and Global Human Capital(IIASA/OeAW/UniWien)
# Vienna Institute of Demography at the Austrian Academy of Sciences
# Vordere Zollamtstraße 3, 1030 Vienna, Austria
# -----------------------------------------------------------------------------------------------------#


# loading necessary packages

library(MortalityLaws)
library(dplyr)
library(tidyverse)
library(lubridate)
library(data.table)
library(purrr)
library(styler)
library(forcats)
library(broom)
library(here)
library(ggpubr)
library(ggExtra)
library(ggthemes)
library(ggrepel)
library(Hmisc)
library(suffrager)
library(devtools)
#install_github("alburezg/suffrager")   # this is for using the suffragette palette : )
library(ggthemes)

# Loading useful functions into environment
source(here("Rcodes","WhoWasSavedMore_DR","Functions_DR.R"))
options(scipen=999)

# Reading in the life tables in single ages, by single-year periods.
# Choosing countries with the longest and varied series available

cntr <- c("SWE", "DNK","FRATNP","JPN", "USA")

# females
LT_fem<- ReadHMD(what = "LT_f",
               countries = cntr,
               interval = "1x1",
               username = "add here your username",
               password = "add here your password",save = FALSE)$data

# males
LT_men<- ReadHMD(what = "LT_m",
               countries = cntr,
               interval = "1x1",
               username = "add here your username",
               password = "add here your password",save = FALSE)$data


# Life table closing at 90+ for Women

lt_f<-LT_fem%>%
  group_by(country,Year) %>%
  mutate(dx=case_when(Age>=90~sum(dx[Age>=90], na.rm = T),TRUE~dx),
         Lx=case_when(Age>=90~sum(Lx[Age>=90], na.rm = T),TRUE~Lx),
         mx=dx/Lx)%>%
  filter(Age<=90) %>%
  ungroup() %>%
  arrange(country,Year,Age) %>%
  group_by(country,Year) %>%
  group_modify(~life.table(.x$mx, sex="f"), .keep=T) %>%
  mutate(Age=0:90, Sex="f")%>%
  relocate(Age, .after = Year) %>%
  ungroup()

# Life table closing at 90+ for Men

lt_m<-LT_men%>%
  group_by(country,Year) %>%
  mutate(dx=case_when(Age>=90~sum(dx[Age>=90], na.rm = T),TRUE~dx),
         Lx=case_when(Age>=90~sum(Lx[Age>=90], na.rm = T),TRUE~Lx),
         mx=dx/Lx)%>%
  filter(Age<=90) %>%
  ungroup() %>%
  arrange(country,Year,Age) %>%
  group_by(country,Year) %>%
  group_modify(~life.table(.x$mx, sex="m"), .keep=T) %>%
  mutate(Age=0:90,
         Sex="m")%>%
  relocate(Age, .after = Year) %>%
  ungroup()

# just for Fig.2 with survival curves combine both sexes as retrieved from the HMD

lt_f$Sex<-"f"
lt_m$Sex<-"m"

lt_t<-rbind(lt_f,lt_m)
lt_t$country<- factor(lt_t$country,levels=c("SWE", "DNK","FRATNP","JPN", "USA"),
                      labels = c("Sweden", "Denmark","France","Japan", "USA"))
lt_t$Sex<- factor(lt_t$Sex, labels = c("Women", "Men"))


LT_fem$Sex<-"f"
LT_men$Sex<-"m"

LT_t<-rbind(LT_fem,LT_men)
LT_t$country<- factor(LT_t$country,levels=c("SWE", "DNK","FRATNP","JPN", "USA"),
                      labels = c("Sweden", "Denmark","France","Japan", "USA"))
LT_t$Sex<- factor(LT_t$Sex, labels = c("Women", "Men"))


LT_t_la<-LT_t%>%
  arrange(Year) %>%
  group_by(country,Age,Sex) %>%
  slice(c(1,n()))
View(LT_t_la)
LT_t_la<-LT_t_la %>%
             filter(Age==30)

# FIG 2:

lx<-ggplot(LT_t, aes(Age,lx, color=Year, group=Year))+
  geom_line(size=1.2, alpha=0.8)+
  scale_color_viridis_c()+
  facet_grid(Sex~country)+
  theme_pander(base_size = 24)+
  geom_point(data=LT_t_la %>%
               filter(Age==30),aes(x=Age, y=lx,group=Year),
             size=3, color="black", shape=17)+
  geom_text_repel(data=LT_t_la %>%
                     filter(Age==30 & Sex=="Women"),aes(x=Age, y=lx,group=Year, label=Year),
                   fontface = 'bold',
                   hjust = "right",
                   size = 6,
                  color="black",
                  box.padding = 1.5)+
  geom_vline(xintercept = as.numeric(LT_t_la$Age),
             color = "dark green", size = 0.4)

# saving as .tiff
tiff("lx.tiff", units="in", width=22, height=10, res=300)
lx
dev.off()
# saving as .png
png("lx.png", units="in", width=22, height=10, res=300)
lx
dev.off()

# ----------------------------------------------------------------------------------------------------------#
# first performing analysis by sex separately Countries have different years of availability.
# Denmark: 1835-2020; France: 1816-2018; Sweden: 1751-2019; Japan: 1947-2019; US: 1933-2019
# selecting years according to the literature on mortality transition
# ----------------------------------------------------------------------------------------------------------#

# These were the lx used for the stylized text on the main paper to build Fig 1., panel A.
# ggplot(LT_t %>% filter(country=="Denmark" & Year%in%c(1835,2019) & Sex=="Women")
#       , aes(Age,lx, group=as.factor(Year)))+
#        geom_line(size=1)+
# theme_classic()

# ----------------------------------------------------------------------------------------------------------#
#  I will cut all into 4 groups for each country using literature-driven criteria.
#  First selecting years and then combining into one
# ----------------------------------------------------------------------------------------------------------#

# ----------------------------------------------------------------------------------------------------------#
# estimating lambdas for all countries: using the formulas in Vaupel and Yashin (1987,a): ln(lx_new/lx_old)
# ----------------------------------------------------------------------------------------------------------#


# Sweden
lt_swe<-lt_t %>%
  filter(country=="Sweden")

lt_t_swe<-lt_swe %>%
  filter(Year %in% c(1751,1900,1975,2019)) %>%
  group_by(Age,Sex) %>%
  mutate(lambda=log(lx/lag(lx)),
         lambda=coalesce(lambda,0))%>%
  ungroup() %>%
  group_by(Age,Sex) %>%
  mutate(Year=as.numeric(Year),
         Year_lab=as.character(paste(shift(Year),Year,sep="-")))


# Denmark
lt_den<-lt_t %>%
  filter(country=="Denmark")

lt_t_den<-lt_den %>%
  filter(Year %in% c(1835,1920,1975,2020)) %>%
           group_by(Age,Sex) %>%
    mutate(lambda=log(lx/lag(lx)),
           lambda=coalesce(lambda,0))%>%
    ungroup() %>%
  group_by(Age,Sex) %>%
  mutate(Year=as.numeric(Year),
         Year_lab=as.character(paste(shift(Year),Year,sep="-")))

# France
lt_fra<-lt_t %>%
  filter(country=="France")

lt_t_fra<-lt_fra %>%
  filter(Year %in% c(1816,1900,1975,2018)) %>%
  group_by(Age,Sex) %>%
  mutate(lambda=log(lx/lag(lx)),
         lambda=coalesce(lambda,0))%>%
  ungroup() %>%
  group_by(Age,Sex) %>%
  mutate(Year=as.numeric(Year),
         Year_lab=as.character(paste(shift(Year),Year,sep="-")))

# Japan
lt_jap<-lt_t %>%
  filter(country=="Japan")

lt_t_jap<-lt_jap %>%
  filter(Year %in% c(1947,1975,2002,2019)) %>%
  group_by(Age,Sex) %>%
  mutate(lambda=log(lx/lag(lx)),
         lambda=coalesce(lambda,0))%>%
  ungroup() %>%
  group_by(Age,Sex) %>%
  mutate(Year=as.numeric(Year),
         Year_lab=as.character(paste(shift(Year),Year,sep="-")))


# USA
lt_usa<-lt_t %>%
  filter(country=="USA")

lt_t_usa<-lt_usa %>%
  filter(Year %in% c(1933,1955,1975,2019)) %>%
  group_by(Age,Sex) %>%
  mutate(lambda=log(lx/lag(lx)),
         lambda=coalesce(lambda,0))%>%
  ungroup() %>%
  group_by(Age,Sex) %>%
  mutate(Year=as.numeric(Year),
         Year_lab=as.character(paste(shift(Year),Year,sep="-")))


# join all countries

lt_cnt_all<-rbind(lt_t_swe,lt_t_den,lt_t_fra,lt_t_jap,lt_t_usa)
View(lt_cnt_all)


# building labels
label_year<-lt_cnt_all %>%
  #arrange(country,Year,Sex) %>%
  filter(lambda!=0) %>%
  group_by(country,Sex) %>%
  top_n(1, Age)

# ploting everything: Fig 3.

X11()
l_all<-ggplot(lt_cnt_all %>% filter(lambda!=0),
              aes(Age,lambda, color=Sex, group=Sex))+
  geom_line(aes(linetype=Sex), size=1.8)+
  scale_color_manual(values = rev(suf_palette("hanwell", n = 2,type = "discrete")))+
  facet_wrap(country~Year_lab, scales = "free")+
  labs(y=expression(paste("Intensity of lifesaving ", ~ (lambda[i]))), x="Age")+
  theme_bw()+
  geom_hline(yintercept = 0, size=1)+
  theme_pander(base_size = 18)+
  theme(legend.position = "bottom", axis.title.x = element_blank())
l_all

lt_t_all<-lt_cnt_all

# -----------------------------------------------------------------------------------------#
# estimating the i times revivorship function for averting deaths up to 10 times
# -----------------------------------------------------------------------------------------#

setDT(lt_t_all)

setkey(lt_t_all, country,Year,Age,Sex)

res_all<-function (x,y,z) {
  res<-(lag(x)*(y^z))/factorial(z)
  return(res)
}

# applying formula for i states of revivorship (replicating the paper, until 10th time)

lt_t_all2<-lt_t_all[ , lx_0:= res_all(lx,lambda,0), by=.(country,Sex,Age)]
lt_t_all2<-lt_t_all[ , lx_1:= res_all(lx,lambda,1), by=.(country,Sex,Age)]
lt_t_all2<-lt_t_all[ , lx_2:= res_all(lx,lambda,2), by=.(country,Sex,Age)]
lt_t_all2<-lt_t_all[ , lx_3:= res_all(lx,lambda,3), by=.(country,Sex,Age)]
lt_t_all2<-lt_t_all[ , lx_4:= res_all(lx,lambda,4), by=.(country,Sex,Age)]
lt_t_all2<-lt_t_all[ , lx_5:= res_all(lx,lambda,5), by=.(country,Sex,Age)]
lt_t_all2<-lt_t_all[ , lx_6:= res_all(lx,lambda,6), by=.(country,Sex,Age)]
lt_t_all2<-lt_t_all[ , lx_7:= res_all(lx,lambda,7), by=.(country,Sex,Age)]
lt_t_all2<-lt_t_all[ , lx_8:= res_all(lx,lambda,8), by=.(country,Sex,Age)]
lt_t_all2<-lt_t_all[ , lx_9:= res_all(lx,lambda,9), by=.(country,Sex,Age)]
lt_t_all2<-lt_t_all[ , lx_10:= res_all(lx,lambda,10), by=.(country,Sex,Age)]


lt_t_all2<-lt_t_all2 %>%
  mutate(lx_New=rowSums(lt_t_all2[, c(15:25)], na.rm=TRUE))


lt_t_all3<-lt_t_all2 %>%
   mutate(lx_3=rowSums(lt_t_all2[, c(18:25)], na.rm=TRUE))%>%
select(-c(19:25)) %>%
  mutate(lx_0=coalesce(lx_0,0),
         lx_1=coalesce(lx_1,0),
         lx_2=coalesce(lx_2,0),
         lx_3=coalesce(lx_3,0))

# putting everything in long format
library(magrittr)
library("dplyr")

res_long<-lt_t_all3 %>%
  group_by(country,Year,Sex) %>%
  pivot_longer(cols = c(lx_0:lx_3,lx_New),
               names_to = "lx_res",
               values_to = "res")%>%
  separate(lx_res, into=c("lx_res","i"), sep="_")


# ----------------------------------------------------------------------------------------------------#
# plot revivorship and survivorship differences for each country, then combining similar-transitions
# ----------------------------------------------------------------------------------------------------#

# Sweden

# survivorship
surv_swe<-ggplot()+
  geom_line(data=res_long %>% filter(lambda!=0 & i=="New" & country=="Sweden" & Age>0),
            aes(Age,res, color=factor(i), group=factor(i)), color="black",linetype="dashed", size=1.6)+
  facet_grid(Sex~Year_lab)+
  geom_line(data=res_long%>% filter(lambda!=0 & i==c("0")& country=="Sweden"),
            aes(Age,res,
                color=factor(i), group=factor(i)),size=1.6)+
  scale_color_viridis_d(drop = TRUE, guide="none")+
  geom_area(data=res_long%>% filter(lambda!=0 & i==c("0")& country=="Sweden"),
            aes(Age,res,
                color=factor(i), group=factor(i), fill=factor(i)),position = "identity") +
  scale_fill_manual(values = alpha(rev(suf_palette("hanwell", n = 2,type = "discrete")), 0.5), guide="none")+
  facet_grid(Sex~Year_lab)+
  theme_pander(base_size = 14)+
  theme(legend.position = "none",
        plot.title = element_text(size = 12))+
  #strip.text.x = element_blank())+
  labs(y=expression(paste("Revivorship", ~ italic(l)[i]), x="Age"))+
  ggtitle(label="Survivorship, Sweden")

#revivorship
res_swe<-ggplot()+
  geom_line(data=res_long%>% filter(lambda!=0 & i==c("1","2","3")& country=="Sweden"),
            aes(Age,res,
                color=factor(i), group=factor(i)),size=1.6)+
  scale_color_viridis_d(name="Number of times being saved")+
  geom_area(data=res_long%>% filter(lambda!=0 & i==c("1","2","3")& country=="Sweden"),
            aes(Age,res,
                color=factor(i), group=factor(i), fill=factor(i)),position = "identity") +
  scale_fill_viridis_d(name="Number of times being saved", alpha=0.3)+
  facet_grid(Sex~Year_lab)+
  theme_pander(base_size = 14)+
  theme(legend.position = "bottom",
        legend.title=element_text(size=12),
        plot.title = element_text(size = 12))+
  #strip.text.x = element_blank())+
  labs(y=expression(paste("Revivorship", ~ italic(l)[i]), x="Age"))+
  ggtitle(label="Revivorship, Sweden")


# Denmark

#survivorship
surv_den<-ggplot()+
  geom_line(data=res_long %>% filter(lambda!=0 & i=="New" & country=="Denmark" & Age>0),
            aes(Age,res, color=factor(i), group=factor(i)), color="black",linetype="dashed", size=1.6)+
  facet_grid(Sex~Year_lab)+
  geom_line(data=res_long%>% filter(lambda!=0 & i==c("0")& country=="Denmark"),
            aes(Age,res,
                color=factor(i), group=factor(i)),size=1.6)+
  scale_color_viridis_d(drop = TRUE, guide="none")+
  geom_area(data=res_long%>% filter(lambda!=0 & i==c("0")& country=="Denmark"),
            aes(Age,res,
                color=factor(i), group=factor(i), fill=factor(i)),position = "identity") +
  scale_fill_manual(values = alpha(rev(suf_palette("hanwell", n = 2,type = "discrete")), 0.5), guide="none")+
  facet_grid(Sex~Year_lab)+
  theme_pander(base_size = 14)+
  theme(legend.position = "none",
        plot.title = element_text(size = 12))+
  #strip.text.x = element_blank())+
  labs(y=expression(paste("Revivorship", ~ italic(l)[i]), x="Age"))+
  ggtitle(label="Survivorship, Denmark")

#revivorship
res_den<-ggplot()+
  geom_line(data=res_long%>% filter(lambda!=0 & i==c("1","2","3")& country=="Denmark"),
            aes(Age,res,
                color=factor(i), group=factor(i)),size=1.6)+
  scale_color_viridis_d(name="Number of times being saved")+
  geom_area(data=res_long%>% filter(lambda!=0 & i==c("1","2","3")& country=="Denmark"),
            aes(Age,res,
                color=factor(i), group=factor(i), fill=factor(i)),position = "identity") +
  scale_fill_viridis_d(name="Number of times being saved", alpha=0.3)+
  facet_grid(Sex~Year_lab)+
  theme_pander(base_size = 14)+
  theme(legend.position = "bottom",
        legend.title=element_text(size=12),
        plot.title = element_text(size = 12))+
  #strip.text.x = element_blank())+
  labs(y=expression(paste("Revivorship", ~ italic(l)[i]), x="Age"))+
  ggtitle(label="Revivorship, Denmark")

# France

#survivorship
surv_fra<-ggplot()+
  geom_line(data=res_long %>% filter(lambda!=0 & i=="New" & country=="France" & Age>0),
            aes(Age,res, color=factor(i), group=factor(i)), color="black",linetype="dashed", size=1.6)+
  facet_grid(Sex~Year_lab)+
  geom_line(data=res_long%>% filter(lambda!=0 & i==c("0")& country=="France"),
            aes(Age,res,
                color=factor(i), group=factor(i)),size=1.6)+
  scale_color_viridis_d(drop = TRUE, guide="none")+
  geom_area(data=res_long%>% filter(lambda!=0 & i==c("0")& country=="France"),
            aes(Age,res,
                color=factor(i), group=factor(i), fill=factor(i)),position = "identity") +
  scale_fill_manual(values = alpha(rev(suf_palette("hanwell", n = 2,type = "discrete")), 0.5), guide="none")+
  facet_grid(Sex~Year_lab)+
  theme_pander(base_size = 14)+
  theme(legend.position = "none",
        plot.title = element_text(size = 12))+
  #strip.text.x = element_blank())+
  labs(y=expression(paste("Revivorship", ~ italic(l)[i]), x="Age"))+
  ggtitle(label="Survivorship, France")

#revivorship
res_fra<-ggplot()+
  geom_line(data=res_long%>% filter(lambda!=0 & i==c("1","2","3")& country=="France"),
            aes(Age,res,
                color=factor(i), group=factor(i)),size=1.6)+
  scale_color_viridis_d(name="Number of times being saved")+
  geom_area(data=res_long%>% filter(lambda!=0 & i==c("1","2","3")& country=="France"),
            aes(Age,res,
                color=factor(i), group=factor(i), fill=factor(i)),position = "identity") +
  scale_fill_viridis_d(name="Number of times being saved", alpha=0.3)+
  facet_grid(Sex~Year_lab)+
  theme_pander(base_size = 14)+
  theme(legend.position = "bottom",
        legend.title=element_text(size=12),
        plot.title = element_text(size = 12))+
  #strip.text.x = element_blank())+
  labs(y=expression(paste("Revivorship", ~ italic(l)[i]), x="Age"))+
  ggtitle(label="Revivorship, France")


# Japan

#survivorship
surv_jap<-ggplot()+
  geom_line(data=res_long %>% filter(lambda!=0 & i=="New" & country=="Japan" & Age>0),
            aes(Age,res, color=factor(i), group=factor(i)), color="black",linetype="dashed", size=1.6)+
  facet_grid(Sex~Year_lab)+
  geom_line(data=res_long%>% filter(lambda!=0 & i==c("0")& country=="Japan"),
            aes(Age,res,
                color=factor(i), group=factor(i)),size=1.6)+
  scale_color_viridis_d(drop = TRUE, guide="none")+
  geom_area(data=res_long%>% filter(lambda!=0 & i==c("0")& country=="Japan"),
            aes(Age,res,
                color=factor(i), group=factor(i), fill=factor(i)),position = "identity") +
  scale_fill_manual(values = alpha(rev(suf_palette("hanwell", n = 2,type = "discrete")), 0.5), guide="none")+
  facet_grid(Sex~Year_lab)+
  theme_pander(base_size = 14)+
  theme(legend.position = "none",
        plot.title = element_text(size = 12))+
  #strip.text.x = element_blank())+
  labs(y=expression(paste("Revivorship", ~ italic(l)[i]), x="Age"))+
  ggtitle(label="Survivorship, Japan")

#revivorship
res_jap<-ggplot()+
  geom_line(data=res_long%>% filter(lambda!=0 & i==c("1","2","3")& country=="Japan"),
            aes(Age,res,
                color=factor(i), group=factor(i)),size=1.6)+
  scale_color_viridis_d(name="Number of times being saved")+
  geom_area(data=res_long%>% filter(lambda!=0 & i==c("1","2","3")& country=="Japan"),
            aes(Age,res,
                color=factor(i), group=factor(i), fill=factor(i)),position = "identity") +
  scale_fill_viridis_d(name="Number of times being saved", alpha=0.3)+
  facet_grid(Sex~Year_lab)+
  theme_pander(base_size = 14)+
  theme(legend.position = "bottom",
        legend.title=element_text(size=12),
        plot.title = element_text(size = 12))+
  #strip.text.x = element_blank())+
  labs(y=expression(paste("Revivorship", ~ italic(l)[i]), x="Age"))+
  ggtitle(label="Revivorship, Japan")

# USA

#survivorship
surv_usa<-ggplot()+
  geom_line(data=res_long %>% filter(lambda!=0 & i=="New" & country=="USA" & Age>0),
            aes(Age,res, color=factor(i), group=factor(i)), color="black",linetype="dashed", size=1.6)+
  facet_grid(Sex~Year_lab)+
  geom_line(data=res_long%>% filter(lambda!=0 & i==c("0")& country=="USA"),
            aes(Age,res,
                color=factor(i), group=factor(i)),size=1.6)+
  scale_color_viridis_d(drop = TRUE, guide="none")+
  geom_area(data=res_long%>% filter(lambda!=0 & i==c("0")& country=="USA"),
            aes(Age,res,
                color=factor(i), group=factor(i), fill=factor(i)),position = "identity") +
  scale_fill_manual(values = alpha(rev(suf_palette("hanwell", n = 2,type = "discrete")), 0.5), guide="none")+
  facet_grid(Sex~Year_lab)+
  theme_pander(base_size = 14)+
  theme(legend.position = "none",
        plot.title = element_text(size = 12))+
  #strip.text.x = element_blank())+
  labs(y=expression(paste("Revivorship", ~ italic(l)[i]), x="Age"))+
 ggtitle(label="Survivorship, USA")

#revivorship
res_usa<-ggplot()+
  geom_line(data=res_long%>% filter(lambda!=0 & i==c("1","2","3")& country=="USA"),
            aes(Age,res,
                color=factor(i), group=factor(i)),size=1.6)+
  scale_color_viridis_d(name="Number of times being saved")+
  geom_area(data=res_long%>% filter(lambda!=0 & i==c("1","2","3")& country=="USA"),
            aes(Age,res,
                color=factor(i), group=factor(i), fill=factor(i)),position = "identity") +
  scale_fill_viridis_d(name="Number of times being saved", alpha=0.3)+
  facet_grid(Sex~Year_lab)+
  theme_pander(base_size = 14)+
  theme(legend.position = "bottom",
        legend.title=element_text(size=12),
        plot.title = element_text(size = 12))+
  #strip.text.x = element_blank())+
  labs(y=expression(paste("Revivorship", ~ italic(l)[i]), x="Age"))+
  ggtitle(label="Revivorship, USA")


# ---------------------------------------------------------------------------------------------------#
# Combining plots for survival and revivorship: Figs 4-8
# ---------------------------------------------------------------------------------------------------#


# Sweden
fig_res_swe<- ggarrange(surv_swe,res_swe, labels = c("A", "B"), ncol = 1, font.label = list(size=12))

fig_res_swe

# Denmark
fig_res_den<- ggarrange(surv_den,res_den, labels = c("A", "B"), ncol = 1, font.label = list(size=12))

fig_res_den

# France
fig_res_fra<- ggarrange(surv_fra,res_fra, labels = c("A", "B"), ncol = 1, font.label = list(size=12))

fig_res_fra

# USA
fig_res_us<- ggarrange(surv_usa,res_usa, labels = c("A", "B"), ncol = 1, font.label = list(size=12))

fig_res_us

# Japan
fig_res_jap<- ggarrange(surv_jap,res_jap, labels = c("A", "B"), ncol = 1, font.label = list(size=12))

fig_res_jap


# -----------------------------------------------------------------------------------------------------#
# now estimating the tau´s, which are the number of life years expected to be lived in each revivorship
# state. We chose to approximate the integrals through Riemanns´sum, as they yielded small errors.
# -----------------------------------------------------------------------------------------------------#

tau_all<-function (lx) {
    m<-91
    a<-0
    b<-90
    w<-(b-a)/m
    x<-seq(a+(w/2),b-(w/2),w)
    h<-lx
    ex<-sum(h*w)
    return(ex)
}

# computing taus, difference in ex0 and proportion of revived persons

res_long2<-res_long %>%
  group_by(country,Year_lab,Sex,i) %>%
  mutate(tau=tau_all(res))

res_tau<-res_long2 %>%
 # filter(lx>0) %>%
  group_by(country,Year_lab,Sex,i) %>%
  summarise(tau) %>%
  slice(n()) %>%
  spread(i,tau) %>%
  ungroup()

write.table(res_tau, file="res_tau.csv", sep=",", row.names = F) # info for building Table 1.


##--------------------------------------------------------------------------------------------------------------------#
# Estimating e-dagger and entropy for revivorship states. I adapted this function from Alyson's formulae in
# https://github.com/alysonvanraalte/LifeIneq. See source codes for how I adapted the function to include
# the higher order terms with factorials. Some results were compared to Aburto, J.M., van Raalte, A.
# Lifespan Dispersion in Times of Life Expectancy Fluctuation: The Case of Central and Eastern Europe.
# Demography 55, 2071–2096 (2018). https://doi.org/10.1007/s13524-018-0729-9 and also to the historical data in
# Vaupel and Yashin (1987a) used for the US (from Faber, J. F. 1982. Life Tables for the United States: 1900-2050.).
#--------------------------------------------------------------------------------------------------------------------#

lt_H<-res_long2 %>%
  filter(i%in%c("0","1","2","3")) %>%
group_by(country,Year,Sex,i) %>%
  mutate(edagg_1=ineq_edag_i(Age,dx,lx,ex,ax,1)*-1,
         edagg_2=ineq_edag_i(Age,dx,lx,ex,ax,2),
         edagg_3=ineq_edag_i(Age,dx,lx,ex,ax,3)*-1,
         edagg_4=ineq_edag_i(Age,dx,lx,ex,ax,4),
         edagg_5=ineq_edag_i(Age,dx,lx,ex,ax,5)*-1,
         H_1=ineq_H_i(Age,dx,lx,ex,ax,1)*-1,
         H_2=ineq_H_i(Age,dx,lx,ex,ax,2),
         H_3=ineq_H_i(Age,dx,lx,ex,ax,3)*-1,
         H_4=ineq_H_i(Age,dx,lx,ex,ax,4),
         H_5=ineq_H_i(Age,dx,lx,ex,ax,5)*-1) %>%
  ungroup() %>%
  select(1:3,11,12,14,16,19:28) %>%
  filter(Age%in%0 & i==0)

lt_H_spans<-lt_H %>%
  mutate(H1_gain=ex*H_1,                      # this is actually equal to edagger. But I estimated to check
         H2_gain=ex*H_2,
         H3_gain=ex*H_3,
         H4_gain=ex*H_4,
         H5_gain=ex*H_5,
         H1_span=ex+H1_gain,
         H2_span=H1_span+H2_gain,
         H3_span=H2_span+H3_gain,
         H4_span=H3_span+H4_gain,
         H5_span=H4_span+H5_gain)

write.table(lt_H_spans, file="entropy2.csv", sep=",", row.names = F) # info for building

