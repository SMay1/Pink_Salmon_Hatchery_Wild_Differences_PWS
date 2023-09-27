#Analysis of Hatchery-wild data, tables, and plots for manuscript titled:
#"Phenotype-assortative dispersal and differences in fitness-linked traits between hatchery and wild pink salmon (Oncorhynchus gorbuscha) in Prince William Sound, Alaska"


#setwd to source file location
library(lubridate)## Manipulate dates 
library(tidyverse)  
library(DT)       ## To display pretty tables
library(ggpubr)
library(AICcmodavg) 
library(MuMIn)   #For dredge functionality
library(minqa)
library(nloptr)
library(sjPlot) #plot mixed models
library(merTools) #for predictInterval
library(RLRsim) #for likelihood ratio tests
library(glmmTMB)
library(lme4)
library(viridis)


###### Data sourced from ADF&G###
HW1320raw <- read.csv("../data/CSVStreamSpeci_2013-2020_9-21-22.csv")

#make dates work for plotting etc.
HW1320raw$SurveyDate<-as.Date(HW1320raw$SurveyDate, "%m/%d/%Y")  
HW1320raw$Year<-year(HW1320raw$SurveyDate)  
HW1320raw$DOY<- yday (HW1320raw$SurveyDate) 

#Group tributaries of Fitness Streams together in a new column "StreamGroup" for a whole stream
HW1320raw$StreamGroup <- NA
HW1320raw$StreamGroup <- as.character(HW1320raw$StreamName) #Pull raw streams into "StreamGroup for graphing later
HW1320raw[which(HW1320raw$StreamName == "Paddy C" | HW1320raw$StreamName == "Paddy Lower Right Trib" | HW1320raw$StreamName == "Paddy Left Trib"), 
]$StreamGroup <- "Paddy Creek"
HW1320raw[which(HW1320raw$StreamName == "Stockdale C" | HW1320raw$StreamName == "Stockdale Right Trib"), 
]$StreamGroup <- "Stockdale Creek"
HW1320raw[which(HW1320raw$StreamName == "Gilmour C" | HW1320raw$StreamName == "Gilmour Right Trib Below Lake"), 
]$StreamGroup <- "Gilmour Creek"
HW1320raw[which(HW1320raw$StreamName == "Hogan Bay"), 
]$StreamGroup <- "Hogan Bay Creek"
HW1320raw[which(HW1320raw$StreamName == "Erb C"), 
]$StreamGroup <- "Erb Creek"
HW1320raw[which(HW1320raw$ADFGStreamCode == "221-20-10200"),    #Two "Spring Creeks" so identify by AWC code. 
]$StreamGroup <- "Spring CreekF"

#make Fitness and Straying stream labels for plotting
HW1320raw$ADFGStreamCode <- as.character(HW1320raw$ADFGStreamCode)
HW1320raw$StreamType <- NA
HW1320raw[which(HW1320raw$ADFGStreamCode == "221-20-10200" | HW1320raw$ADFGStreamCode == "226-20-16040"  | 
                  HW1320raw$ADFGStreamCode == "227-20-17520"  | HW1320raw$ADFGStreamCode == "226-20-16010"  | 
                  HW1320raw$ADFGStreamCode == "227-20-17480-2" | HW1320raw$ADFGStreamCode == "226-30-16810" |  
                  HW1320raw$ADFGStreamCode == "226-20-16010-2"  | HW1320raw$ADFGStreamCode == "227-20-17480" |  
                  HW1320raw$ADFGStreamCode == "226-20-16010-3" | HW1320raw$ADFGStreamCode == "227-20-17520-2"), 
]$StreamType <- "Fitness"
HW1320raw[!c(HW1320raw$ADFGStreamCode == "221-20-10200" | HW1320raw$ADFGStreamCode == "226-20-16040"  | 
               HW1320raw$ADFGStreamCode == "227-20-17520"  | HW1320raw$ADFGStreamCode == "226-20-16010"  | 
               HW1320raw$ADFGStreamCode == "227-20-17480-2" | HW1320raw$ADFGStreamCode == "226-30-16810" |  
               HW1320raw$ADFGStreamCode == "226-20-16010-2"  | HW1320raw$ADFGStreamCode == "227-20-17480" |  
               HW1320raw$ADFGStreamCode == "226-20-16010-3" | HW1320raw$ADFGStreamCode == "227-20-17520-2"), 
]$StreamType <- "Straying"


#QAQC on all the little things
HW1320all <- filter(HW1320raw,  
                    Species == 440,                               #Pink salmon = 440 chum = 450
                    MarkStatusDescription == "OK",                #exclude "missing" and "no read" otoliths
                    MarkPresent == "Y" | MarkPresent == "N",      # Remove "null" and blank categories
                    Length > 200 & Length < 600,                  #Remove small outliers - same as ADFG protocol, if use 250 fixes plotting outliers
                    Sex == "M" | Sex == "F",                      #Remove unknown sex and blanks. 
                    NoHWAnalysis == "False") %>%                  #ture = deemed 'mixed up' in field/lab by QAQC         
  separate(ADFGStreamCode, c("pre", "mid", "end")) %>%            #split stream numbers by region, etc to...(next line)
  mutate(Region = case_when(pre == 221 | pre == 228 ~ "Eastern",  #designate eastern and western streams
                            pre != 221 | pre != 228 ~ "Western"),
         Lineage = as.factor(case_when(Year == 2013 | 
                                         Year == 2015 |           #delineate even and odd years
                                         Year == 2017 | 
                                         Year == 2019 ~ "Odd Years",
                                       Year == 2014 | 
                                         Year == 2016 | 
                                         Year == 2018 | 
                                         Year == 2020 ~ "Even Years")), 
         Sex = case_when(Sex == "F" ~ "Female",                   #rename for plotting ease
                         Sex == "M" ~ "Male"),
         MarkPresent = case_when(MarkPresent == "N" ~ "Wild",     #rename for plotting ease
                                 MarkPresent == "Y" ~ "Hatchery"))

#make things factors for plotting/modeling
HW1320all$StreamType <- as.factor (HW1320all$StreamType)
HW1320all$Sex <- as.factor(HW1320all$Sex) 
HW1320all$MarkPresent <- as.factor(HW1320all$MarkPresent)
HW1320all$StreamGroup <- as.factor(HW1320all$StreamGroup)


with(HW1320all, table(Lineage))
#Even Years  Odd Years of 2013 - 2018
#51194      60986 
#Even Years  Odd Years of 2013-2020
#81674     137159 


#Table of Length, SD, # oto's read of H and W origin from 2013-2020
TableLengths <- HW1320all %>%
  group_by(Year, StreamGroup, MarkPresent) %>%
  summarize(mean_length = mean(Length),
            sd_length = sd(Length),
            oto_count = n()) %>%
  mutate(combined_values = paste(mean_length, sd_length, oto_count, sep = ", ")) %>%
  dplyr::select(!c(mean_length:oto_count)) %>%
  pivot_wider(names_from = Year, values_from = combined_values) %>%
  arrange(StreamGroup, MarkPresent)

write.csv(TableLengths, "../figures/Figure1TableLengths.csv")



###### Figure 2: Boxplots of length and timing ####
######BOXPLOT LENGTH: facet by even and odd, combine with model plots?

LengthBox <- 
  ggplot(HW1320all, aes(x = MarkPresent, y = Length, fill = MarkPresent))+ #fill = Hatchery is interesting but busier to look at
  geom_boxplot(width = 0.9,outlier.shape = NA)+  #make line at mean and outliers instead of full whiskers
  facet_grid(Sex~Lineage, scales = "free_y", space = "free_y") +
  xlab("Origin") +
  ylab("Length (mm)") +
  theme_classic(base_size = 20)+ 
  theme (#axis.line = element_line(color="black"), 
    #legend.title=element_blank(),
    legend.position = "none",
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(margin = margin(r = 0)),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.spacing.x = unit(2, "mm"),
  )+
  lims(y=c(290,520))+
  scale_fill_brewer(palette="Paired")

#split by stream for supplemental table? Try it

# some possible supplemental plots ######
#histogram of fish length by year: possible supplemental plot?
ggplot(HW1320all, aes(x = Length, fill = Sex))+
  geom_histogram(binwidth = 10, position = "dodge")+
  facet_grid(Year~MarkPresent)+
  coord_cartesian(xlim =c(295, 550))+
  theme_bw()

#organize and arrange streams by longitude (east-west)
StreamLong <- HW1320all %>% group_by(StreamGroup) %>%  
  summarise(StreamLong = min(na.rm = T, Longitude)) %>% 
  arrange(StreamLong)

#SUPPLIMENTAL #1 LENGTH ALL streams, male and female length even and odd, New S1 as of 9/27/23
  S1Lengths<- HW1320all %>% mutate(StreamGroup = factor(StreamGroup, levels = StreamLong$StreamGroup)) %>%  
  #filter (Lineage == "Odd Years") %>% #arrange(StreamGroup) %>% 
  ggplot(aes(x = StreamGroup, y = Length, fill = MarkPresent))+  #fill = Hatchery or fill = MarkPresent both interesting!
  #geom_violin(draw_quantiles = c(.25, .50, .75))+ #violin plot is kinda cool
  geom_boxplot()+
  facet_grid(Sex~Lineage)+
  xlab("Origin") +
  ylab("Length (mm)") +
  theme_bw(base_size = 24)+ 
  theme (axis.line = element_line(color="black"), 
         legend.title=element_blank(),
         legend.position = "none",
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.x = element_text(angle=45, vjust=1, hjust=1))
ggsave("../figures/FigureS1Length.png", S1Lengths, units = "in", dpi = 600, width = 18, height = 10)

#Straying streams, male and female length 
HW1320all %>% mutate(StreamGroup = factor(StreamGroup, levels = StreamLong$StreamGroup)) %>%  
  filter (StreamType == "Straying") %>% #arrange(StreamGroup) %>% 
  ggplot(aes(x = StreamGroup, y = Length, fill = MarkPresent))+  #fill = Hatchery or fill = MarkPresent both interesting!
  #geom_violin(draw_quantiles = c(.25, .50, .75))+ #violin plot is kinda cool
  geom_boxplot()+
  facet_grid(Sex~Lineage)+
  xlab("Population") +
  ylab("Length (mm)") +
  theme_bw(base_size = 24)+ 
  theme (axis.line = element_line(color="black"), 
         legend.title=element_blank(),
         legend.position = "none",
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.x = element_text(angle=45, vjust=1, hjust=1))


##BOXPLOT TIMING #####
#Figure 2b? 
TimingBox <- 
  ggplot(HW1320all, aes(x = MarkPresent, y = DOY, fill = MarkPresent))+  #fill = Hatchery or fill = MarkPresent both interesting!
  geom_boxplot(width = 0.9,outlier.shape = NA)+
  facet_grid(Sex~Lineage, scales = "free_y", space = "free_y")+ 
  xlab("Origin") +
  ylab("Recovery Date (DOY)") +
  theme_classic(base_size = 20)+ 
  theme (legend.position = "none",
         strip.text.x = element_blank(),
         panel.spacing = unit(2, "mm"),
         axis.title.x = element_blank(),
         strip.background = element_blank()
  ) +
  lims(y=c(190,270))+
  scale_fill_manual(values = c("#B2DF8A", "#33A02C"))


# Calculate CV for each group
compute_CV <- function(data, group, value) {
  data %>%
    group_by_(.dots = group) %>%
    summarise(mean_value = mean(!!sym(value), na.rm = TRUE),
              sd_value = sd(!!sym(value), na.rm = TRUE)) %>%
    mutate(CV = (sd_value/mean_value)*100)
}

CV_Length <- compute_CV(HW1320all, c("MarkPresent", "Sex", "Lineage"), "Length")
CV_Timing <- compute_CV(HW1320all, c("MarkPresent", "Sex", "Lineage"), "DOY")


LengthBox <- 
  LengthBox +
  geom_text(data = CV_Length, aes(label = sprintf("%.2f%%", CV), y = Inf), vjust = 1.5, hjust = -0.3, size = 4)

TimingBox <- 
  TimingBox +
  geom_text(data = CV_Timing, aes(label = sprintf("%.2f%%", CV), y = Inf), vjust = 1.5, hjust = -0.3, size = 4)




FigureTwo <- ggpubr::ggarrange(LengthBox,TimingBox, nrow = 2, legend = "none", align = "v" )
ggsave("FigureTwo.png", FigureTwo, units = "in", dpi = 600, width = 8, height = 10)


#just for boxplots break out into stream supplemental model: ####
StreamLong <- HW1320all %>% group_by(StreamGroup) %>%  
  summarise(StreamLong = min(na.rm = T, Longitude)) %>% 
  arrange(StreamLong)
#SUPPLIMENTAL #2 TIMING All streams, male sand female, even and odd, DOY by stream, New S2 as of 9/27/23
S2DOY<-HW1320all %>% mutate(StreamGroup = factor(StreamGroup, levels = StreamLong$StreamGroup)) %>% 
  #filter (StreamType == "Straying") %>% 
  ggplot(aes(x = StreamGroup, y = DOY, fill = MarkPresent))+  #fill = Hatchery or fill = MarkPresent both interesting!
  #geom_violin(draw_quantiles = c(.25, .50, .75))+ #violin plot is kinda cool
  geom_boxplot()+
  facet_grid(Sex~Lineage)+
  xlab("Population") +
  ylab("DOY") +
  theme_bw(base_size = 24)+ 
  theme (axis.line = element_line(color="black"), 
         legend.title=element_blank(),
         legend.position = "none",
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.x = element_text(angle=45, vjust=1, hjust=1))
#check this with hatcheries split out - some H's have earlier fish, some later fish...make boxes proportionaly match the sample size?
ggsave("../figures/FigureS2DOY.png", S2DOY, dpi=600, width = 18, height = 10,units = "in")


###### Model analysis data prep ####

# center DOY, ALL data  
HW1320all$DOYcentered<-scale(HW1320all$DOY)
#Reinstate factor for random intercept testing

###Split data by even odd
HWeven <- HW1320all %>% filter( Lineage == "Even Years")  
HWodd <- HW1320all %>% filter( Lineage == "Odd Years") 
###### See if not centering and scaling runs models with convergence - make plotting easier
###Which are we using glmer or gmlmTMB



###### Length, timing analyses######

#########EVEN LENGTH analysis

even_length_global <- lmer(Length ~ MarkPresent + DOYcentered + Sex + 
                             MarkPresent:Sex + DOYcentered:Sex + 
                             (1|StreamGroup) + (1|Year), 
                           data = HWeven, na.action = "na.fail")
summary(even_length_global)
even_length_dredge <- dredge(even_length_global) 
even_length_best <- lmer(as.vector(Length) ~ MarkPresent + DOYcentered + Sex + 
                           (1|StreamGroup) + (1|Year), 
                         data = HWeven, na.action = "na.fail") #only best by 0.90 AIC, top 3 all under 2 AIC
summary(even_length_best)
plot_model(even_length_best, type = "pred") # with glmer: Could not compute variance-covariance matrix of predictions. no confidence intervals
qqnorm(main = "Normal Q-Q Even Length Best", residuals(even_length_best))
qqline(residuals(even_length_best))

AIC_table_even_length<-as.data.frame(even_length_dredge)[7:11]
AIC_table_even_length[,2:5]<-round(AIC_table_even_length[,2:5],2)
names(AIC_table_even_length[1])<-"K"
AIC_table_even_length$Model<-rownames(AIC_table_even_length)
for(i in 1:nrow(AIC_table_even_length)) {AIC_table_even_length$Model[i]<- as.character(formula(get.models(even_length_dredge,subset = T)[[i]]))[3]}
write.csv(AIC_table_even_length,"../figures/AIC_table_even_length.csv")


##ODD LENGTH analysis

odd_length_global <- lmer(Length ~ MarkPresent + DOYcentered + Sex + 
                            MarkPresent:Sex + DOYcentered:Sex + 
                            (1|StreamGroup) + (1|Year), 
                          data = HWodd, na.action = "na.fail", ) #does not converge with glmmTMB
summary(odd_length_global)
odd_length_dredge <- dredge(odd_length_global) 
odd_length_best <- lmer(Length ~ MarkPresent + DOYcentered + Sex + 
                          MarkPresent:Sex + DOYcentered:Sex + 
                          (1|StreamGroup) + (1|Year), 
                        data = HWodd, na.action = "na.fail", ) #best by 79 AIC
summary(odd_length_best)
plot_model(odd_length_best, type = "pred")
qqnorm(main = "Normal Q-Q Odd Length Best", residuals(odd_length_best))
qqline(residuals(odd_length_best))


AIC_table_odd_length<-as.data.frame(odd_length_dredge)[7:11]
AIC_table_odd_length[,2:5]<-round(AIC_table_odd_length[,2:5],2)
names(AIC_table_odd_length[1])<-"K"
AIC_table_odd_length$Model<-rownames(AIC_table_odd_length)
for(i in 1:nrow(AIC_table_odd_length)) {AIC_table_odd_length$Model[i]<- as.character(formula(get.models(odd_length_dredge,subset = T)[[i]]))[3]}
write.csv(AIC_table_odd_length,"AIC_table_odd_length.csv")

###TIMING analysis (swap length for timing...anything to change in how treating timing variable?)

###EVEN TIMING
even_timing_global <- lmer(DOYcentered ~ MarkPresent + Length + Sex + 
                             MarkPresent:Sex + Length:Sex + 
                             (1|StreamGroup) + (1|Year), 
                           data = HWeven, na.action = "na.fail")  
summary(even_timing_global)
even_timing_dredge <- dredge(even_timing_global)                  
even_timing_best <- lmer(DOYcentered ~ MarkPresent + Sex + 
                           MarkPresent:Sex + 
                           (1|StreamGroup) + (1|Year), 
                         data = HWeven, na.action = "na.fail")
summary(even_timing_best)
plot_model(even_timing_best, type = "pred")
qqnorm(main = "Normal Q-Q Even Timing Best", residuals(even_timing_best))
qqline(residuals(even_timing_best))

AIC_table_even_timing<-as.data.frame(even_timing_dredge)[7:11]
AIC_table_even_timing[,2:5]<-round(AIC_table_even_timing[,2:5],2)
names(AIC_table_even_timing[1])<-"K"
AIC_table_even_timing$Model<-rownames(AIC_table_even_timing)
for(i in 1:nrow(AIC_table_even_timing)) {AIC_table_even_timing$Model[i]<- as.character(formula(get.models(even_timing_dredge,subset = T)[[i]]))[3]}
write.csv(AIC_table_even_timing,"AIC_table_even_timing.csv")



####ODD TIMING
odd_timing_global <- lmer(DOYcentered ~ MarkPresent + Length + Sex +    
                            MarkPresent:Sex + Length:Sex + 
                            (1|StreamGroup) + (1|Year), 
                          data = HWodd, na.action = "na.fail")

summary(odd_timing_global)
odd_timing_dredge <- dredge(odd_timing_global) 
odd_timing_best <- lmer(DOYcentered ~ MarkPresent + Length + Sex + 
                          MarkPresent:Sex + Length:Sex + 
                          (1|StreamGroup) + (1|Year), 
                        data = HWodd, na.action = "na.fail")
plot_model(odd_timing_best, type = "pred")

summary(odd_timing_best)
qqnorm(main = "Normal Q-Q Odd Timing Best", residuals(odd_timing_best))
qqline(residuals(odd_timing_best))

AIC_table_odd_timing<-as.data.frame(odd_timing_dredge)[7:11]
AIC_table_odd_timing[,2:5]<-round(AIC_table_odd_timing[,2:5],2)
names(AIC_table_odd_timing[1])<-"K"
AIC_table_odd_timing$Model<-rownames(AIC_table_odd_timing)
for(i in 1:nrow(AIC_table_odd_timing)) {AIC_table_odd_timing$Model[i]<- as.character(formula(get.models(odd_timing_dredge,subset = T)[[i]]))[3]}
write.csv(AIC_table_odd_timing,"AIC_table_odd_timing.csv")


##Table 2 AIC Table #####
# Combine dataframes and select columns, reordering 'Model' column first
AIC_table_all <- bind_rows(
  AIC_table_even_length %>% mutate(Model = "Even Length"),
  AIC_table_odd_length %>% mutate(Model = "Odd Length"),
  AIC_table_even_timing %>% mutate(Model = "Even Timing"),
  AIC_table_odd_timing %>% mutate(Model = "Odd Timing")) 
#%>% select(Model, df, logLik, AICc, delta, weight)




#####  OLD Model diagnostics######### still looks fine for selected new models
## reisduals Vs fitted  -looks pretty much fine for whole model
plot(All_fixmods[[5]]) 
qqnorm(resid(All_fixmods[[5]]))
qqline(resid(All_fixmods[[5]]))

## standardized residuals versus fitted values by sex, lineage, markpresent, etc
plot(All_fixmods[[5]], resid(., scaled=TRUE) ~ fitted(.) | Lineage, abline = 0) #Looks good, although even years narrow range of sizes!
plot(All_fixmods[[5]], resid(., scaled=TRUE) ~ fitted(.) | DOYcentered, abline = 0) #Hard to tell, but looks fine? Some days VERY few points...
plot(All_fixmods[[5]], resid(., scaled=TRUE) ~ fitted(.) | MarkPresent, abline = 0) #Two VERY diff blobs in H fish - even odd years?

## Check out even vs odd for those groups
ODD <- filter(HW1320all, Lineage == "Odd Years")
EVEN <- filter(HW1320all, Lineage =="Even Years")
Odd.process <- lmer(Length ~  DOYcentered * MarkPresent * Sex + StreamGroup + (1|Year), data = ODD, REML = FALSE)
Even.process <- lmer(Length ~  DOYcentered * MarkPresent * Sex + StreamGroup + (1|Year), data = EVEN, REML = FALSE)
#yep, those hatchery blobs are VARY diffs in even/odd years
plot(Odd.process, resid(., scaled=TRUE) ~ fitted(.) | MarkPresent, abline = 0) #yep, those hatchery blobs are diffs in even/odd
plot(Even.process, resid(., scaled=TRUE) ~ fitted(.) | MarkPresent, abline = 0) #yep, those hatchery blobs are diffs in even/odd

## Standardized residuals vs fitted for random effect group (possible?)
## observed versus fitted values by Subject
plot(All_fixmods[[5]], Length ~ fitted(.) | StreamGroup, abline = c(0,1))  #I wonder what this looks like with random slopes?




###### Modeling and Model plots############################

##CREATE dataframde for LENGTH predictions
pred_length_even = data.frame(DOYcentered = rep(seq(min(HW1320all$DOYcentered),max(HW1320all$DOYcentered),length.out=7), 4), 
                              Lineage = 'Even Years',
                              MarkPresent = c(rep('Wild', length(min(HW1320all$DOYcentered):max(HW1320all$DOYcentered))*2), 
                                              rep('Hatchery', length(min(HW1320all$DOYcentered):max(HW1320all$DOYcentered))*2)),                       
                              Sex=c(rep('Male', length(min(HW1320all$DOYcentered):max(HW1320all$DOYcentered))), 
                                    rep('Female', length(min(HW1320all$DOYcentered):max(HW1320all$DOYcentered))), 
                                    rep('Male', length(min(HW1320all$DOYcentered):max(HW1320all$DOYcentered))), 
                                    rep('Female', length(min(HW1320all$DOYcentered):max(HW1320all$DOYcentered)))),
                              StreamGroup=NA, 
                              Year= NA)
pred_length_odd = data.frame(DOYcentered = rep(seq(min(HW1320all$DOYcentered),max(HW1320all$DOYcentered),length.out=7), 4), 
                             Lineage = 'Odd Years',
                             MarkPresent = c(rep('Wild', length(min(HW1320all$DOYcentered):max(HW1320all$DOYcentered))*2), 
                                             rep('Hatchery', length(min(HW1320all$DOYcentered):max(HW1320all$DOYcentered))*2)),                       
                             Sex=c(rep('Male', length(min(HW1320all$DOYcentered):max(HW1320all$DOYcentered))), 
                                   rep('Female', length(min(HW1320all$DOYcentered):max(HW1320all$DOYcentered))), 
                                   rep('Male', length(min(HW1320all$DOYcentered):max(HW1320all$DOYcentered))), 
                                   rep('Female', length(min(HW1320all$DOYcentered):max(HW1320all$DOYcentered)))),
                             StreamGroup=NA, 
                             Year= NA)                      

##EVEN 
pred_length_even$Length <- predict(even_length_best, newdata=pred_length_even, allow.new.levels = T, re.form = NA)
#pred_length_even$se.fit <- predict(even_length_best, newdata=pred_length_even, re.form = NA, se.fit = TRUE)$se.fit
##ODD##
pred_length_odd$Length <- predict(odd_length_best, newdata=pred_length_odd, allow.new.levels = T, re.form = NA)
#pred_length_odd$se.fit <- predict(odd_length_best, newdata=pred_length_odd, allow.new.levels = T, re.form = NA, se.fit = TRUE)$se.fit
#combine
pred_length <- rbind(pred_length_even, pred_length_odd)
#Create confidence intervals
pred_length$lwr <- pred_length$DOYcentered - 1.96*(pred_length$se.fit) 
pred_length$upr <- pred_length$DOYcentered + 1.96*(pred_length$se.fit)


pred_length$DOY = pred_length$DOYcentered*attr(HW1320all$DOYcentered, "scaled:scale") + attr(HW1320all$DOYcentered, "scaled:center")

####### Length Plot #####
##PLOT LENGTH
LengthPlot <- 
  ggplot(pred_length, aes(y = Length , x = DOY, se.fit, 
                          color = interaction(MarkPresent, Sex, sep=" "), 
                          linetype = interaction(MarkPresent, Sex, sep=" ")))+   
  geom_jitter(data = HW1320all, alpha = 0.15, size = .2)+
  geom_line(lwd = 1.5)+
  scale_color_manual(values= c("#A6CEE3", "#1F78B4","#A6CEE3", "#1F78B4"))+
  scale_linetype_manual(values = c("solid", "solid", "twodash", "twodash"))+
  coord_cartesian(xlim = c(200, 270), ylim = c(300, 525))+
  #geom_ribbon(aes(as.factor(DOY)), ymin = lwr, ymax = upr, alpha = .2)+
  facet_grid(~Lineage)+
  xlab("Recovery Date (DOY)") +
  ylab("Length (mm)") +
  #ggtitle("length") +
  theme_classic(base_size = 20)+ 
  theme (axis.line = element_line(color="black"), 
         legend.title=element_blank(),
         legend.position = "bottom",
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         strip.background = element_blank(),
  )
#ggsave("Lengthplottesting.png", LengthPlot, units = "in", dpi = 600, width = 10, height = 8)



###### Timing dataframe for modeling #####
#Create dataframe for predicting Timing data
pred_timing_even = data.frame(Length = rep(min(HW1320all$Length):max(HW1320all$Length), 4), 
                              Lineage = 'Even Years',
                              MarkPresent = c(rep('Wild', length(min(HW1320all$Length):max(HW1320all$Length))*2), 
                                              rep('Hatchery', length(min(HW1320all$Length):max(HW1320all$Length))*2)),                       
                              Sex=c(rep('Male', length(min(HW1320all$Length):max(HW1320all$Length))), 
                                    rep('Female', length(min(HW1320all$Length):max(HW1320all$Length))), 
                                    rep('Male', length(min(HW1320all$Length):max(HW1320all$Length))), 
                                    rep('Female', length(min(HW1320all$Length):max(HW1320all$Length)))),
                              StreamGroup=NA, 
                              Year= NA)
pred_timing_odd = data.frame(Length = rep(min(HW1320all$Length):max(HW1320all$Length), 4), 
                             Lineage = 'Odd Years',
                             MarkPresent = c(rep('Wild', length(min(HW1320all$Length):max(HW1320all$Length))*2), 
                                             rep('Hatchery', length(min(HW1320all$Length):max(HW1320all$Length))*2)),                       
                             Sex=c(rep('Male', length(min(HW1320all$Length):max(HW1320all$Length))), 
                                   rep('Female', length(min(HW1320all$Length):max(HW1320all$Length))), 
                                   rep('Male', length(min(HW1320all$Length):max(HW1320all$Length))), 
                                   rep('Female', length(min(HW1320all$Length):max(HW1320all$Length)))),
                             StreamGroup=NA, 
                             Year= NA)

#TIMING Predict data using model and fill empty data frame for TIMING
##EVEN 
pred_timing_even$DOY <- predict(even_timing_best, newdata=pred_timing_even, allow.new.levels = T, re.form = NA)
#pred_timing$se.fit <- predict(even_timing_best, newdata=pred_timing, re.form = NA, se.fit = TRUE)$se.fit
##ODD##
pred_timing_odd$DOY <- predict(odd_timing_best, newdata=pred_timing_odd, allow.new.levels = T, re.form = NA)
#pred_timing$se.fit <- predict(even_timing_best, newdata=pred_timing, re.form = NA, se.fit = TRUE)$se.fit
pred_timing <- rbind(pred_timing_even, pred_timing_odd)

#Create confidence intervals
pred_timing$lwr <- pred_timing$DOY - 1.96*(pred_timing$se.fit) 
pred_timing$upr <- pred_timing$DOY + 1.96*(pred_timing$se.fit)

pred_timing$DOY = pred_timing$DOY*attr(HW1320all$DOYcentered, "scaled:scale") + attr(HW1320all$DOYcentered, "scaled:center")
####### PLOT TIMING #####


TimingPlot <- 
  ggplot(pred_timing, aes(y = DOY, x = Length, se.fit, color = interaction(MarkPresent, Sex, sep=" "), linetype = interaction(MarkPresent, Sex, sep=" ")))+   
  geom_jitter(data = HW1320all, alpha = 0.15, size = .2)+
  geom_line(lwd = 1.5)+
  scale_color_manual(values= c("#B2DF8A", "#33A02C","#B2DF8A", "#33A02C"))+
  scale_linetype_manual(values = c("solid", "solid", "twodash", "twodash"))+
  coord_cartesian(ylim = c(200, 270), xlim = c(300, 525))+
  #geom_ribbon(aes(as.factor(DOY)), ymin = lwr, ymax = upr, alpha = .2)+
  facet_grid(~Lineage)+
  xlab("Length (mm)") +
  ylab("Recovery Date (DOY)") +
  theme_classic(base_size = 20)+ 
  theme (axis.line = element_line(color="black"), 
         legend.title=element_blank(),
         legend.position = "bottom",
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         strip.background = element_blank(),
         strip.text = element_blank())
#ggsave("Timingplottesting.png", TimingPlot, units = "in", dpi = 600, width = 10, height = 8)


####### Figure Three ####
FigureThree <- ggpubr::ggarrange(LengthPlot, TimingPlot, nrow = 2, align = "v")
ggsave("FigureThreePoints.png", FigureThree, units = "in", dpi = 600, width = 9, height = 10)






# EXTRA STUFF #####
#Make centered DOY from predicted
#pred_All$DOYcentered<-scale(pred_All$DOY, center = TRUE, scale = TRUE)

#Constrain dates to actual dates by subsetting
Subset<- HW1320all %>%   #find min and max dates from original data
  filter(Hatchery != "SOLOMON GULCH") %>% #Why take solomon gulch out? (no difference in dates...)
  group_by(Lineage, MarkPresent)%>%  #Year, MarkPresent, Sex for data below
  count(min(DOY), max(DOY))


a <- pred_All %>%
  filter(Lineage == "Even Years" & MarkPresent == "Hatchery") %>% filter(DOY >= 205 & DOY <= 261) 
b <- pred_All %>%
  filter(Lineage == "Even Years" & MarkPresent == "Wild") %>% filter(DOY >= 200 & DOY <= 262)
c <- pred_All %>%
  filter(Lineage == "Odd Years" & MarkPresent == "Hatchery") %>% filter(DOY >= 199 & DOY <= 267) 
d <- pred_All %>%
  filter(Lineage == "Odd Years" & MarkPresent == "Wild") %>% filter(DOY >= 191 & DOY <= 267)

pred_All <- rbind(a,b,c,d)


#Fited model Plots

###Length plot!
#HW1320all4way<-
ggplot(pred_length, aes(as.factor(DOY), se.fit, col = MarkPresent, linetype = Sex))+   
  geom_line(size = 1.5)+
  #geom_ribbon(aes(as.factor(DOY)), ymin = lwr, ymax = upr, alpha = .2)+
  facet_grid(~Lineage)+
  xlab("Recovery Date") +
  ylab("Length (mm)") +
  #ggtitle("All Samples, 4-way, Predicted") +
  theme_bw(base_size = 24)+ 
  theme (axis.line = element_line(color="black"), 
         legend.title=element_blank(),
         #legend.position = "none",
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())+
  scale_linetype_manual(values = c(6,6,1,1), name ="", labels = c("Hatchery Female", "Hatchery Male", "Wild Female", "Wild Male")) +
  scale_color_manual(values = c("grey46", "black","grey46", "black"), name = "", labels = c("Hatchery Female", "Hatchery Male", "Wild Female", "Wild Male")) +
  scale_x_discrete(breaks = c(206, 232, 258), labels = c ("Jul 25", "Aug 20","Sept 15")) 
#ggsave("HW1320allPredict4-way.jpeg", HW1320all4way, height = 7, width = 9, units = "in", dpi = 300)

##### Try this with HW data points overlayed.  In manuscript as of 4/1/20
color<-ggplot()+ 
  geom_jitter(data = HW1320all, aes(as.factor(DOY), Length, group = MarkPresent:Sex, col = MarkPresent:Sex), alpha = 1/5)+
  geom_ribbon(data = pred_All, aes(as.factor(DOY), ymin = lwr, ymax = upr, fill = MarkPresent:Sex,group = MarkPresent:Sex), alpha = .4)+
  geom_line(data = pred_All, aes(as.factor(DOY), fit, group = MarkPresent:Sex, col = MarkPresent:Sex, linetype = MarkPresent:Sex), size = 1.5)+
  facet_wrap(~Lineage)+
  xlab("Recovery Date") +
  ylab("Length (mm)") +
  #ggtitle("All Samples, 4-way, Predicted") +
  theme_bw(base_size = 24)+ 
  theme (axis.line = element_line(color="black"), 
         legend.title=element_blank(),
         legend.position = "bottom",
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())+
  scale_linetype_manual(values = c(1,6,1,6), name ="", labels = c("Hatchery Female", "Hatchery Male", "Wild Female", "Wild Male")) +
  scale_color_viridis(discrete = TRUE, name = "", labels = c("Hatchery Female", "Hatchery Male", "Wild Female", "Wild Male")) +
  scale_fill_viridis(discrete = TRUE, name = "", labels = c("Hatchery Female", "Hatchery Male", "Wild Female", "Wild Male"))+
  scale_x_discrete(breaks = c(206, 232, 258), labels = c ("Jul 25", "Aug 20","Sept 15")) +
  coord_cartesian(ylim = c(325, 525))
ggsave("HW1320allPredict4-waycolor.jpeg", color, height = 7, width = 10, units = "in", dpi = 600)        


#######  Figure 4:  Size and timing selective sorting by hatchery strays######

Stream_averagesize <- HW1320all %>%  group_by (Year, StreamGroup, MarkPresent) %>% 
  summarise(stream_average_size = mean(Length),n=n())
All_averagesize <- HW1320all %>%  group_by (Year, MarkPresent) %>% 
  summarise(all_average_size = mean(Length))

AveSizeDiff <- Stream_averagesize %>% left_join(All_averagesize, by = c("Year", "MarkPresent")) %>% 
  mutate(size_diff = stream_average_size - all_average_size) %>% 
  dplyr::select(Year, StreamGroup, MarkPresent, size_diff) %>% 
  group_by(Year, StreamGroup) %>% 
  pivot_wider(names_from = "MarkPresent", values_from = "size_diff") %>% 
  drop_na() %>% 
  mutate(Lineage = ifelse(Year%%2, "Odd years", "Even years"),
         Lineage = as.factor(Lineage), StreamGroup=as.factor(StreamGroup),
         Year = as.factor(Year))


GlobalAveDiff <- lmer(Hatchery~Wild * Lineage + (1|Year), data = AveSizeDiff, na.action = "na.fail")
summary((GlobalAveDiff))
dredge(GlobalAveDiff)
plot_model(GlobalAveDiff, type = "pred")
hist(AveSizeDiff$Hatchery, breaks = 20)
AveDiffDredge<-dredge(GlobalAveDiff)

ggplot(AveSizeDiff, aes(Wild, Hatchery)) +
  geom_point()+
  facet_wrap(~Year)+
  geom_smooth(method="lm")+
  geom_abline(slope=1)

Pred_Size_Diff <- data.frame(Wild = rep(seq(min(AveSizeDiff$Wild),max(AveSizeDiff$Wild),length.out=1000),2), 
                             Lineage = rep(c('Even years',"Odd years"),each=1000), 
                             Year= NA)
Pred_Size_Diff$Hatchery<-predict(GlobalAveDiff,newdata = Pred_Size_Diff,allow.new.levels=T)
Pred_Size_Diff$se<-predict(GlobalAveDiff,newdata = Pred_Size_Diff,se.fit = T)$se.fit

AveSizeDiff %>% ggplot(aes(x=Hatchery, y=Wild, color = Lineage))+
  geom_point(alpha=0.2)+
  geom_ribbon(data=Pred_Size_Diff,aes(ymin=Wild-se,ymax=Wild+se,fill=Lineage),alpha=0.2)+
  geom_line(data=Pred_Size_Diff,inherit.aes = T)+
  ylim(c(-40,40))+
  xlim(c(-40,40))+
  theme_classic()

library(merTools)
##Even Size:
AveSizeDiff_even<- AveSizeDiff %>% filter(Lineage=="Even years")
GlobalAveDiff_even <- lmer(Hatchery~Wild + (1|Year), data = AveSizeDiff_even, na.action = "na.fail")
AveDiff_even_dredge<-dredge(GlobalAveDiff_even)
GlobalAveDiff_even <- lmer(Hatchery~Wild + (1|Year), data = AveSizeDiff_even, na.action = "na.fail")
BestAveDiff_even <- lmer(Hatchery~1 + (1|Year), data = AveSizeDiff_even, na.action = "na.fail")
r.squaredGLMM(BestAveDiff_even)
Pred_Size_Diff_even <- data.frame(Wild = seq(min(AveSizeDiff_even$Wild),max(AveSizeDiff_even$Wild),length.out=1000), 
                                  Lineage = "Even years", 
                                  Year= factor(NA, levels = levels(AveSizeDiff_even$Year)))
Pred_Size_Diff_even$Hatchery <- predictSE(GlobalAveDiff_even, newdata = Pred_Size_Diff_even, 
                                          re.form = ~0, estimate.se = F)$fit
Pred_Size_Diff_even$se <- predictSE(GlobalAveDiff_even, newdata = Pred_Size_Diff_even, 
                                    re.form = ~0, estimate.se = TRUE)$se.fit


##Odd Size:
AveSizeDiff_odd<- AveSizeDiff %>% filter(Lineage=="Odd years")
GlobalAveDiff_odd <- lmer(Hatchery~Wild + (1|Year), data = AveSizeDiff_odd, na.action = "na.fail")
dredge(GlobalAveDiff_odd)
Pred_Size_Diff_odd <- data.frame(Wild = seq(min(AveSizeDiff_odd$Wild),max(AveSizeDiff_odd$Wild),length.out=1000), 
                                 Lineage = "Odd years", 
                                 Year= factor(NA, levels = levels(AveSizeDiff_odd$Year)))
Pred_Size_Diff_odd$Hatchery <- predictSE(GlobalAveDiff_odd, newdata = Pred_Size_Diff_odd, 
                                         re.form = ~0, estimate.se = TRUE)$fit
Pred_Size_Diff_odd$se <- predictSE(GlobalAveDiff_odd, newdata = Pred_Size_Diff_odd, 
                                   re.form = ~0, estimate.se = TRUE)$se.fit


Pred_Size_Diff <- rbind(Pred_Size_Diff_odd,Pred_Size_Diff_even)
#Pred_Size_Diff<-rbind(Pred_Size_Diff_even,Pred_Size_Diff_odd)
#Size assortative plot ####
AveSizeDiffPlot <- ggplot(AveSizeDiff, aes(x=Wild, y=Hatchery, color = Lineage))+
  geom_point(alpha = .6, size = 2)+
  geom_ribbon(data=Pred_Size_Diff,aes(ymin=Hatchery-se*1.96,ymax=Hatchery+se*1.96,fill = Lineage),alpha=0.15, linetype = 2)+
  geom_line(data=Pred_Size_Diff)+
  scale_color_manual(values = c("deepskyblue4","sienna1"))+
  #ylim(c(-40,40))+
  #xlim(c(-40,40))+
  theme_classic(base_size = 14)+
  theme(legend.position = "none")+
  labs(y=expression(atop(mu[PWS] - mu[stream], "body length of hatchery fish (mm)")),
       x=expression(atop(mu[PWS] - mu[stream], "body length of wild fish (mm)")))

##timing####
Stream_averagetime <- HW1320all %>%  group_by (Year, StreamGroup, MarkPresent) %>% 
  summarise(stream_average_time = mean(DOY))
All_averagetime <- HW1320all %>%  group_by (Year, MarkPresent) %>% 
  summarise(all_average_time = mean(DOY))

AveTimeDiff <- Stream_averagetime %>% left_join(All_averagetime, by = c("Year", "MarkPresent")) %>% 
  mutate(time_diff = stream_average_time - all_average_time) %>% 
  dplyr::select(Year, StreamGroup, MarkPresent, time_diff) %>% 
  group_by(Year, StreamGroup) %>% 
  pivot_wider(names_from = "MarkPresent", values_from = "time_diff") %>% 
  drop_na()%>% 
  mutate(Lineage = ifelse(Year%%2, "Odd years", "Even years"))

AveTimeDiff %>% ggplot(aes(x=Hatchery, y=Wild))+
  geom_point()

GlobalAveTimeDiff <- lmer(Hatchery~Wild + Lineage + (1|Year), data = AveTimeDiff, na.action = "na.fail")
summary((GlobalAveTimeDiff))

plot_model(GlobalAveTimeDiff, type = "pred")
hist(AveTimeDiff$Hatchery, breaks = 20)
dredge(GlobalAveTimeDiff)

BestAveTimeDiff <- glmmTMB(Hatchery~Wild + (1|Year), data = AveTimeDiff, na.action = "na.fail")

Pred_Time_Diff <- data.frame(Wild = rep(seq(min(AveTimeDiff$Wild),max(AveTimeDiff$Wild),length.out=1000),2), 
                             Lineage = rep(c('Even years',"Odd years"),each=1000), 
                             Year= NA)
Pred_Time_Diff$Hatchery<-predict(BestAveTimeDiff,newdata = Pred_Time_Diff)
Pred_Time_Diff$se<-predict(BestAveTimeDiff,newdata = Pred_Time_Diff,se.fit = T)$se.fit


#Timing split for plotting
##Even Timing:
AveTimeDiff_even<- AveTimeDiff %>% filter(Lineage=="Even years")
GlobalAveTimeDiff_even <- lmer(Hatchery~Wild + (1|Year), data = AveTimeDiff_even, na.action = "na.fail")
dredge(GlobalAveTimeDiff_even)
summary(GlobalAveTimeDiff_even)
Pred_Time_Diff_even <- data.frame(Wild = seq(min(AveTimeDiff_even$Wild),max(AveTimeDiff_even$Wild),length.out=1000), 
                                  Lineage = "Even years", 
                                  Year= factor(NA, levels = levels(AveTimeDiff_even$Year)))
Pred_Time_Diff_even$Hatchery <- predictSE(GlobalAveTimeDiff_even, newdata = Pred_Time_Diff_even, 
                                          re.form = ~0, estimate.se = TRUE)$fit
Pred_Time_Diff_even$se <- predictSE(GlobalAveTimeDiff_even, newdata = Pred_Time_Diff_even, 
                                    re.form = ~0, estimate.se = TRUE)$se.fit

##Odd Timing:
AveTimeDiff_odd<- AveTimeDiff %>% filter(Lineage=="Odd years")
GlobalAveTimeDiff_odd <- lmer(Hatchery~Wild + (1|Year), data = AveTimeDiff_odd, na.action = "na.fail")
dredge(GlobalAveTimeDiff_odd)
summary(GlobalAveTimeDiff_odd)
Pred_Time_Diff_odd <- data.frame(Wild = seq(min(AveTimeDiff_odd$Wild),max(AveTimeDiff_odd$Wild),length.out=1000), 
                                 Lineage = "Odd years", 
                                 Year= factor(NA, levels = levels(AveTimeDiff_odd$Year)))
Pred_Time_Diff_odd$Hatchery <- predictSE(GlobalAveTimeDiff_odd, newdata = Pred_Time_Diff_odd, 
                                         re.form = ~0, estimate.se = TRUE)$fit
Pred_Time_Diff_odd$se <- predictSE(GlobalAveTimeDiff_odd, newdata = Pred_Time_Diff_odd, 
                                   re.form = ~0, estimate.se = TRUE)$se.fit

Pred_Time_Diff<-rbind(Pred_Time_Diff_even,Pred_Time_Diff_odd)
#Timing assortative plot ####
AveTimeDiffPlot <- ggplot(AveTimeDiff, aes(x=Wild, y=Hatchery, color = Lineage))+
  geom_point(alpha = .6, size = 2)+
  geom_ribbon(data=Pred_Time_Diff,aes(ymin=Hatchery-se*1.96,ymax=Hatchery+se*1.96,fill=Lineage),alpha=0.15, linetype = 2)+
  geom_line(data=Pred_Time_Diff)+
  scale_color_manual(values = c("sienna1", "deepskyblue4"))+
  #ylim(c(-40,40))+
  #xlim(c(-40,40))+
  theme_classic(base_size = 14)+
  theme(legend.title = element_blank(), 
        legend.position = c(0.8, 0.2))+
  labs(y=expression(atop(mu[PWS] - mu[stream], "recovery date of hatchery fish (DOY)")),
       x=expression(atop(mu[PWS] - mu[stream], "recovery date of wild fish (DOY)")))

FigureFour <- ggpubr::ggarrange(AveSizeDiffPlot, AveTimeDiffPlot, nrow = 2, align = "v")
ggsave("../figures/FigureFour.png", FigureFour, units = "in", dpi = 600, width = 6, height = 8)




#####results writing part####
#predicted lengths for the results writing part: 
#Compare Hatchery to wild all years sexes combined:
pred_All %>%
  group_by(MarkPresent) %>%         #OR group by MarkPresent or Sex to specifically compare those. 
  summarise(mean_length = mean(fit), 
            sd_length = sd(fit),
            CV = (sd(fit)/mean(fit))*100,  #Higher CV's indicate more variation around the mean; ie the "level of variation relative to the mean"
            count = length(fit))
t.test(fit~MarkPresent, data = pred_All)  #t = 4.985, df = 1066.4, p-value = 7.229e-07
#Over all years wild fish were x bigger than hatchery:
abs((mean(subset(pred_All, MarkPresent == "Wild")$fit) - mean(subset(pred_All, MarkPresent == "Hatchery")$fit))/((mean(subset(pred_All, MarkPresent == "Wild")$fit) + mean(subset(pred_All, MarkPresent == "Hatchery")$fit))/2))*100 
mean(subset(pred_All, MarkPresent == "Wild")$fit) - mean(subset(pred_All, MarkPresent == "Hatchery")$fit)
# = 1.08 % bigger wild fish than hatchery fish, or 4.4 mm larger

#Compare even to odd year fish:
pred_All %>% 
  group_by(Lineage) %>% 
  summarise(mean_length = mean(fit), 
            sd_length = sd(fit), 
            CV = (sd(fit)/mean(fit))*100,
            count = length(fit))
#Even year fish were x% bigger than odd: 
abs((mean(subset(pred_All, Lineage == "Even Years")$fit) - mean(subset(pred_All, Lineage == "Odd Years")$fit))/((mean(subset(pred_All, Lineage == "Even Years")$fit) + mean(subset(pred_All, Lineage == "Odd Years")$fit))/2))*100 
mean(subset(pred_All, Lineage == "Even Years")$fit) - mean(subset(pred_All, Lineage == "Odd Years")$fit)
#= 4.6% bigger even year fish, or 18.9 mm larger even year fish

#compare groups of hatchery and wild between even and odd years
#For HW1320 data
HW1320all %>%
  group_by(Lineage, MarkPresent) %>%      
  summarise(mean_length = mean(Length), 
            sd_length = sd(Length), 
            CV = (sd(Length)/mean(Length))*100,
            count = length(Length))

pred_All %>%
  group_by(Lineage, MarkPresent) %>%      
  summarise(mean_length = mean(fit), 
            sd_length = sd(fit), 
            CV = (sd(fit)/mean(fit))*100,
            count = length(fit))
#in Even years, hatchery fish 1.2% , or 5.1 mm longer than wild fish:
abs((mean(subset(pred_All, Lineage == "Even Years" & MarkPresent == "Wild")$fit) - mean(subset(pred_All, Lineage == "Even Years" & MarkPresent == "Hatchery")$fit))/((mean(subset(pred_All, Lineage == "Even Years"& MarkPresent == "Wild")$fit) + mean(subset(pred_All, Lineage == "Even Years"& MarkPresent == "Hatchery")$fit))/2))*100  
mean(subset(pred_All, Lineage == "Even Years" & MarkPresent == "Hatchery")$fit) - mean(subset(pred_All, Lineage == "Even Years" & MarkPresent == "Wild")$fit) 
#3.1% longer Wild fish in odd years, or 12.7 mm longer wild fish in odd years
abs((mean(subset(pred_All, Lineage == "Odd Years" & MarkPresent == "Wild")$fit) - mean(subset(pred_All, Lineage == "Odd Years" & MarkPresent == "Hatchery")$fit))/((mean(subset(pred_All, Lineage == "Odd Years"& MarkPresent == "Wild")$fit) + mean(subset(pred_All, Lineage == "Odd Years"& MarkPresent == "Hatchery")$fit))/2))*100  
mean(subset(pred_All, Lineage == "Odd Years" & MarkPresent == "Wild")$fit) - mean(subset(pred_All, Lineage == "Odd Years" & MarkPresent == "Hatchery")$fit)

#Males vs females size in even years
#in Even years, female fish 0.5% , or 2.3 mm longer than male fish
abs((mean(subset(pred_All, Lineage == "Even Years" & Sex == "Male")$fit) - mean(subset(pred_All, Lineage == "Even Years" & Sex == "Female")$fit))/((mean(subset(pred_All, Lineage == "Even Years"& Sex == "Male")$fit) + mean(subset(pred_All, Lineage == "Even Years"& Sex == "Female")$fit))/2))*100
mean(subset(pred_All, Lineage == "Even Years" & Sex == "Male")$fit) - mean(subset(pred_All, Lineage == "Even Years" & Sex == "Female")$fit) 
#in Odd years, male fish 1.6% , or 6.5 mm longer than female fish
abs((mean(subset(pred_All, Lineage == "Odd Years" & Sex == "Male")$fit) - mean(subset(pred_All, Lineage == "Odd Years" & Sex == "Female")$fit))/((mean(subset(pred_All, Lineage == "Odd Years"& Sex == "Male")$fit) + mean(subset(pred_All, Lineage == "Odd Years"& Sex == "Female")$fit))/2))*100
mean(subset(pred_All, Lineage == "Odd Years" & Sex == "Male")$fit) - mean(subset(pred_All, Lineage == "Odd Years" & Sex == "Female")$fit) 

f<- pred_All %>% filter(Lineage == "Even Years", Sex == "Female", MarkPresent == "Hatchery")
g<- pred_All %>% filter(Lineage == "Even Years", Sex != "Female" | MarkPresent != "Hatchery")
#in Even years, female fish 1.8% , or 7.7 mm longer than all other fish
abs((mean(g$fit) - mean(f$fit))/((mean(g$fit) + mean(f$fit))/2))*100
mean(g$fit) - mean(f$fit) 
g %>%   
  group_by(Lineage, MarkPresent, Sex)%>%
  count()















## EXTRA STUFF #####
#for the Freqpoly plot (no longer in manuscript) looking at salmon return dates by hatchery
Test <- filter(HW1320raw,
               Species == 440,
               MarkStatusDescription == "OK",                #exclude "missing" and "no read" oto's
               MarkPresent == "Y" | MarkPresent == "N",      # Remove "null" and blank categories in column
               
               Length > 200,                                #Remove small outliers - same as ADFG protocol
               Sex == "M" | Sex == "F",                        #Remove unknown sex and blanks
               NoHWAnalysis == "False") %>%                       
  separate(ADFGStreamCode, c("pre", "mid", "end")) %>%     #split stream numbers by region, etc
  mutate(Region = as.factor(case_when(pre == 221 | pre == 228 ~ "Eastern",  #designate eastern and western streams
                                      pre != 221 | pre != 228 ~ "Western")),
         Lineage = as.factor(case_when(Year == 2013 | Year == 2015 | Year == 2017 | Year == 2019 ~ "Odd Years",
                                       Year == 2014 | Year == 2016 | Year == 2018 | Year == 2020 ~ "Even Years")), 
         Sex = as.factor(case_when(Sex == "F" ~ "Female",
                                   Sex == "M" ~ "Male")),
         MarkPresent = as.factor(case_when(MarkPresent == "N" ~ "Wild",
                                           MarkPresent == "Y" ~ "Hatchery")),
         Origin = as.factor(case_when(Hatchery == "SOLOMON GULCH" ~ "SGH",
                                      MarkPresent == "Wild" ~ "Wild",
                                      Hatchery != "SOLOMON GULCH" | Hatchery != "" ~ "AKF, CCH, WNB")))  

