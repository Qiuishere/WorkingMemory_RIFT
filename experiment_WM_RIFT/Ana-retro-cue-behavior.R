
## 0. load all the functions we need --------------------

rm(list = ls())

library( pacman)
p_load(bruceR, afex, tidyverse, rstudioapi, ggpubr, psycho, install = TRUE)

## plots
YB = c("#00A19D", "#FFB344")
GP  = c("#4A996B", "#677FA3", "#636395", "#DB6B97")
POP = c("#A8DACD", "#EB68A0", "#7A4D9F", "#22235F")
Grandient.P = c('#484c7f', "#7A4D9F",'#ac8daf','#ddb6c6', '#f1d4d4')
Grad.orange = c('#ffd7b5','#ffb38a', '#ff9248', '#ff6700')
Grad.blue = c('#bbe4e9','#79c2d0', '#53a8b6', '#54878f')
BW = c('black','white')

Sizes = data.frame(matrix(ncol = 1, nrow = 1))
Sizes$point = 2
Sizes$line = 0.6
Sizes$errorbar = 2
Sizes$titletext = 12
Sizes$axistext = 10
Sizes$texttext = 9
Sizes$striptext = 16

mycurve = list(
  stat_summary(geom = "point", fun = "mean", size = 2.5, alpha = 0.6),
  stat_summary(geom = "errorbar", position = position_dodge(0.7),  fun.data = "mean_cl_boot", width = 0.1, size = 0.5),
  stat_summary(geom = 'line', fun = 'mean', size = 1),
  scale_color_manual(values =Grandient.P),  
  theme_gray(),
  theme(strip.text.x = element_text(size = 14, color = "black", face = "bold"), 
        strip.text.y = element_text(size = 14, color = "black", face = "bold"),
        strip.background = element_blank(),
        panel.background = element_rect(fill = 'grey95',color = 'white'),
        panel.grid.major.x = element_line(linetype = 'solid', color='grey90'),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(linetype = 'solid',color='grey90'),
        panel.grid.minor.y = element_line(linetype = 'solid',color='grey90'),
        plot.title = element_text(hjust = 0.5, size = 18, color = "black",face = "bold"), 
        axis.title.x = element_text(size = 12, color = "black", face = "plain"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"), 
        axis.text.x = element_text(size = 12, color = "black", face = "bold"), 
        axis.text.y = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 12, face = "bold"),legend.title = element_blank(),legend.position = "right", strip.placement = "outside",
  )
)

barWidth = 0.9
mybar.2 = list(
  stat_summary(geom = "bar", fun = "mean", position = position_dodge(), color='black'),
  geom_line(aes(group = subject),position = position_dodge(0), color = 'grey90', linetype = 'solid'),
  stat_summary(geom = "errorbar", position = position_dodge(0.9),  fun.data = "mean_cl_boot", color = 'grey20', width = 0.2, size = 1),
  geom_point(aes(group = stimulus), position = position_dodge(0.9), color = 'black', fill = 'grey80'),
  theme_gray(), scale_fill_manual(values = POP),  
  geom_hline(yintercept = 0, size = 1.2), 
  theme(strip.text.x = element_text(size = 14, color = "black", face = "bold"), 
        strip.text.y = element_text(size = 14, color = "black", face = "bold"),
        strip.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "white",linetype = 'solid'),
        panel.grid.minor.y = element_line(color = "white",linetype = 'dashed'),
        plot.title = element_text(hjust = 0.5, size = 18, color = "black",face = "bold"), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 12, color = "black", face = "bold"), 
        axis.title.y = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "none",
  )
)

mybar.3 = list(
  stat_summary(geom = "bar", fun = "mean", position = position_dodge()),
 # geom_line(aes(group = subject),position = position_dodge(0), color = 'grey90', linetype = 'solid'),
  stat_summary(geom = "errorbar", position = position_dodge(0.9),  fun.data = "mean_cl_boot", color = 'grey40', width = 0.3, size = 1.5),
  theme_gray(), scale_color_manual(values = POP),  scale_fill_manual(values = POP),  
  geom_hline(yintercept = 0, size = 2), 
  theme(strip.text.x = element_text(size = 14, color = "black", face = "bold"), 
        strip.text.y = element_text(size = 14, color = "black", face = "bold"),
        strip.background = element_blank(),
        panel.background = element_rect(color = 'white', fill = 'white'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "white",linetype = 'solid'),
        panel.grid.minor.y = element_line(color = "white",linetype = 'dashed'),
        plot.title = element_text(hjust = 0.5, size = 18, color = "black",face = "bold"), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 11, color = "black", face = "bold"), 
        axis.title.y = element_text(size = 14, color = "black", face = "bold"), 
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "none",
  )
)


## 1. load data ------------------------
currentFolder = dirname(getSourceEditorContext()$path)
dataFolder = paste(currentFolder, "/", sep = "")
fileName = "Data_Demo/*/Test/ResultTable*.csv"
fileList = Sys.glob(file.path(dataFolder, fileName))

D.raw = map_df(fileList, read_csv, col_types = cols(subjectid = col_factor(), stimulus = col_factor(), targetId = col_factor(), view1 = col_factor(), view2 = col_factor()))

cutoff.absError = 30
cutoff.initiationTime = 3000
cutoff.rt = 7

D = D.raw %>% 
  rename(subject = subjectid) %>% 
  filter(!is.na(Error)) %>% 
  filter(rt < cutoff.rt)%>%
  mutate(absError = abs(Error),
         angleUncued = if_else(targetId==1, angle2, angle1),
         ErrorUncued = resp - angleUncued,
         stimulus = factor(stimulus, levels = c('Body','Bar')),
         )
 # filter(absError < cutoff.absError) %>% 
 # mutate(Error = -Error) # now negative error means tend to adjust to lower, for clearer presentation
# note that error is referred to the screen, thus 

Nsub.raw = length(unique(D$subject))

time = D.raw %>% 
  group_by(subject) %>% 
  dplyr::summarise(totaltime = max(adjustOnset)/60)
aveTime = mean(time$totaltime)

## 2.  screen the bad data  ----------------------------------------------------



D %>%#filter(subject=='l7tlsvp') %>% 
  ggplot(aes(x = stimulus, y = Error, fill = stimulus, color = 'black')) +
  mybar.2 +  
  labs(y = 'Lower <--- Error---> Higher') + scale_y_reverse()+
  facet_grid(~subject, ncol = 2)


## 3. some inspections ---------------------------------------------------------
# see if people tend to underestimate the adjustment distance: bised towrds the startAngle
D %>% ggplot(aes(x = angleProbe- angleTarget, y = Error))+
  geom_point(alpha = 0.05)+
  geom_smooth(method=lm) + #  (by default includes 95% confidence region)
  coord_fixed(ratio=1)+
  facet_wrap(~subject, ncol = 2)

D %>% ggplot(aes(x= angleProbe, y = resp)) +
  geom_point()+
  geom_smooth(method=lm) +
  facet_wrap(~subject,ncol = 2)

# initiation time distribution
D %>% group_by(subject) %>% 
  ggplot( aes(x= initiationTime)) +
  geom_histogram(alpha = 0.6, binwidth = 400, position = 'identity') +
  facet_grid(~subject)

D = D %>% group_by(subject) %>%
  mutate(medInitTime = median(initiationTime)) %>% ungroup %>% 
  mutate(longTime = if_else(initiationTime > medInitTime, 1, 0))

f.tocued = D %>% ggplot(aes(x = angleTarget, y = resp))+
  geom_point()+
  geom_smooth(method=lm)+
  facet_grid(~subject)

f.touncued = D %>% ggplot(aes(x = angleUncued, y = resp))+
  geom_point()+
  geom_smooth(method=lm)+
  facet_grid(~subject)

hi = ggarrange(f.tocued, f.touncued,  ncol = 2)
plot(hi)


# attraction to the uncued target, separated by order to see order effect
f.angle1 = D %>% filter(targetId==2) %>% 
  ggplot(aes(x= angle1, y = resp)) +
  geom_point()+
  geom_smooth(method=lm) 

f.angle2 = D %>% filter(targetId==1) %>% 
  ggplot(aes(x= angle2, y = resp)) +
  geom_point()+
  geom_smooth(method=lm)
recencyEffect = ggarrange(f.angle1, f.angle2,  ncol = 2)
plot(recencyEffect)


## 4. plot results -------------------------------------------------------------


 ## group plot, by condition ---------------------------------------------------

ErrorbyCondi = D %>% 
  group_by(subject,stimulus) %>% 
  dplyr::summarise(SD = sd(Error,na.rm = TRUE), Error = mean(Error), ErrorUncued = mean(ErrorUncued), RT = mean(rt,na.rm=TRUE), Count = n(), nNA = sum(is.na(Error))) %>% 
  ungroup() 

ErrorbyCondi %>% ggplot(aes(x = stimulus, y = Error, fill = stimulus)) +
  mybar.2 + scale_y_reverse() +
  labs(y = 'Lower <--- Error---> Higher') 

ErrorbyCondi %>% ggplot(aes(x = stimulus, y = ErrorUncued, fill = stimulus)) +
  mybar.2 + scale_y_reverse() +
  labs(y = 'Lower <--- Error---> Higher') 

ErrorbyCondi %>% ggplot(aes(x = stimulus, y = SD, fill = stimulus)) +
  mybar.2 +
  labs(y = 'STD') 

model = bruceR::MANOVA(ErrorbyCondi,
                       subID = "subject",
                       dv = 'Error',
                       within = c( "stimulus"),
                       sph.correction = "GG")


stimulus = ErrorbyCondi %>%  dplyr::ungroup() %>%  
  dplyr::select(stimulus, Error, subject) %>%  
  tidyr::pivot_wider(names_from = c("stimulus"), values_from = "Error", names_sep = '.') %>% 
  dplyr::select(-'subject') 
TTEST(stimulus, c('Upright', 'Inverted'), paired = TRUE) 


D %>%
  ggplot(aes(x = realDelay, y = Error, color = isi)) +
  geom_point(aes(group = isi), position = position_jitter(width = 0.8))+
  labs(y = 'Lower <--- Error---> Higher') + scale_y_reverse() +
  facet_wrap(~stimulus)

## index of Gravity bias and biomechanical constraints--------------------------


biasIndex = ErrorbyCondi %>%
  dplyr::select(c(subject,blurLevel,armHeight,direction,Error)) %>%
  tidyr::pivot_wider(names_from = c(armHeight,direction), values_from = Error, names_sep = '.') %>% 
  mutate(Gravity = (Upper.Front + Upper.Back + Lower.Front + Lower.Back)/4 ,
         Biomechanical = (Upper.Front - Upper.Back - Lower.Front + Lower.Back)/2)



f.bias.g = biasIndex %>% 
  ggplot(aes(x = blurLevel, y = Gravity, fill = blurLevel)) +
  mybar.3+scale_color_manual(values = Grad.orange)+  scale_fill_manual(values = Grad.orange)+
  geom_point(aes(group = blurLevel), position = position_dodge(0.9), color = 'black', fill = 'grey80')+
  labs( y = 'Bias Index') #title = 'Gravity',

f.bias.b = biasIndex %>% 
  ggplot(aes(x = blurLevel, y = Biomechanical, fill = blurLevel)) +
  mybar.3+scale_color_manual(values = Grad.blue)+  scale_fill_manual(values = Grad.blue)+
  geom_point(aes(group = blurLevel), position = position_dodge(0.9), color = 'black', fill = 'grey80')+
  labs(y = 'Bias Index') #title = 'Biomechanical\n Constraints',

hi = ggarrange(f.bias.g, f.bias.b,  nrow = 1)
plot(hi) #(555 * 555)

fileName = paste0(dataFolder, "/blurLevel.eps")
w = 14/2.54
aspect_ratio <- 0.7
ggsave(fileName, hi, height = w* aspect_ratio, width = w, dpi = "retina")

## statistics ------------------------------------------------------------------
model = bruceR::MANOVA(biasIndex,
                       subID = "subject",
                       dv = 'Biomechanical',
                       within = c( "blurLevel"),
                       sph.correction = "HF")

EMMEANS(
  model,
  effect = "blurLevel",
  by = NULL,
  contrast = "pairwise",
  p.adjust = "bonferroni",
  model.type = "multivariate"
)



# BF analysis
result.G = anovaBF(Gravity ~ blurLevel + subject, data = biasIndex, whichRandom = 'subject', whichModels = 'bottom')
summary(result.G)
result.B = anovaBF(Biomechanical ~ blurLevel + subject, data = biasIndex, whichRandom = 'subject', whichModels = 'bottom')
summary(result.B)


# bias for each condition against 0
biasIndex.wide = biasIndex %>% 
  dplyr::select(c(subject,blurLevel,Biomechanical)) %>% 
  pivot_wider(names_from = 'blurLevel', values_from = 'Biomechanical')

variables = colnames(biasIndex.wide)
stat.error = map(variables[2:length(variables)], function(x) TTEST(biasIndex.wide, x,  test.value = 0) )
stat.error = do.call(rbind, stat.error)


stat.error = stat.error %>% 
  mutate(annotation = ifelse(pval< 0.001,'***', ifelse(pval < 0.01, '**', ifelse(pval < 0.05, '*','n.s.'))))


## fit the indices into linear curve, analyze the slopes -----------------------

temp = biasIndex %>%
  mutate(blurRadius = case_when(blurLevel== 'No blur'~ 0, blurLevel== 'Weak'~ 1, blurLevel== 'Medium'~ 2,blurLevel== 'Strong'~ 3,)) %>% 
  group_by(subject) %>% 
  nest() 

slopes =  temp %>% mutate(model1 = lapply(data, function(df) lm(Gravity ~ blurRadius, data = df) %>% tidy())) %>% 
  unnest(model1) %>%
  mutate(Gravity = estimate) %>% 
  filter(term != "(Intercept)") %>% #remove the intercepts
  dplyr::select(-c(data:p.value))

slopes2 = temp %>% mutate(model2 = lapply(data, function(df) lm(Biomechanical ~ blurRadius, data = df) %>% tidy())) %>% ## do LM to each nested DF
  unnest(model2)  %>% 
  filter(term != "(Intercept)")

slopes$Biomechanical = slopes2$estimate
variables = colnames(slopes)
stat.error = map(variables[2:length(variables)], function(x) TTEST(slopes, x,  test.value = 0) )

slopes %>% pivot_longer( cols = c(Gravity, Biomechanical), names_to = "bias", values_to = 'slope') %>% 
  ggplot(aes(x = bias, y = slope, fill = bias)) +
  mybar.3 + 
  geom_point(position = position_jitter(width = 0.1))+
  labs(y = 'Slope of fitting')

## demography ------------------------------------------------------------------
#excSub = unique(excId.badPerf$subject)
prolific = D.raw %>%
  filter(!is.na(prolificID)) %>%
  #  filter(subject %in% excSub) %>% 
  select(c(subject,prolificID))

age = D.raw %>% 
  filter(!is.na(age)) %>% 
  mutate(age = as.numeric(parse_number(age))) %>% 
  dplyr::summarise(meanAge = mean(age, na.rm = TRUE), rangeAge = range(age, na.rm = TRUE) )

sex = D.raw %>% 
  filter(!is.na(sex)) %>% 
  mutate(sex = str_remove(sex,'<b>'),
         sex = str_remove(sex,'</b>'),
         sex = as.factor(sex)) %>% 
  group_by(sex) %>%
  tally()

hand = D.raw %>% 
  filter(!is.na(hand)) %>% 
  group_by(hand) %>%
  tally()


