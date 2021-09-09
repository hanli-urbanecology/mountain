install.packages('corrplot')
install.packages('wesanderson')
install.packages('deldir')
install.packages('concaveman')

library(ggplot2)
library(mgcv)
library(tidyr)
library(corrplot)
library(vegan)
library(ggforce)
library(deldir)
library(concaveman)
library(RColorBrewer)
library(wesanderson)


### site feature ###
site <- read.csv("net_near.csv")
site_cor <- cor(site[, 5:13])
corrplot(site_cor, method = "pie")

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

p.mat <- cor.mtest(site[, 5:13])

corrplot(site_cor, type="upper", order="hclust", 
         p.mat = p.mat, sig.level = 0.01)

###############################################################################

capture <- read.csv("capture3.csv")

#######nmds####################
bat <- as.matrix(capture[, c(28:32, 34:38)])
bat.mds <- metaMDS(bat, k=2, distance  = "euclidean", trymax = 10000) 


bat.mds_capture <- bat.mds
bat.mds_capture
plot(bat.mds_capture, type = "t")

data.scores <- as.data.frame(scores(bat.mds_capture))  
data.scores$X <- as.numeric(rownames(data.scores))

species.scores <- as.data.frame(scores(bat.mds_capture, "species"))  
species.scores$species <- rownames(species.scores)
write.csv(species.scores, "species_scores_capture.csv")


capture_graph <- merge(data.scores, capture)
capture_graph$year1 <- as.factor(capture_graph$year)
write.csv(capture_graph, "capture_graph.csv")

species.scores <- read.csv("species_scores_capture1.csv")
capture_graph <- read.csv("capture_graph1.csv")

capture_nmds_graph <-ggplot() +
  #geom_delaunay_tile(data=capture_graph, aes(x=NMDS1, y=NMDS2), alpha = 0.3) + 
  #geom_delaunay_segment2(data=capture_graph, aes(x=NMDS1, y=NMDS2, colour = wns, group = -1), size = 2,
  #                       lineend = 'round') +
  
  #data= subset(capture_graph, #NMDS1<0.2&NMDS1>-0.15&NMDS2<0.2&NMDS2>-0.2
  #                  ), aes(x=NMDS1, y=NMDS2, fill=wns)) + 
  geom_mark_ellipse(data=capture_graph, aes(x=NMDS1, y=NMDS2,color = wns, fill = wns), alpha =0.2, 
                    expand = unit(1, "mm")) + 
  geom_text(data= species.scores, aes(x=NMDS1,y=NMDS2, label=species, family="serif", fontface = "bold"), size=4, alpha=1)+
  #geom_polygon(group = tran_graph$year1, fill = tran_graph$year1, alpha = 0.5) + 
  geom_point(data=capture_graph, aes(x=NMDS1, y=NMDS2, color=wns, fill=wns, shape = wns), size = 3, alpha =0.3) +
  xlim(-0.9, 0.65) +
  ylim(-0.75, 0.75) +
  scale_fill_manual(values=wes_palette(n=2, name = "Royal1"), name="Site", 
                    breaks = c("pre", "post"), labels=c("Pre-WNS", "Post-WNS")) +
  scale_color_manual(values=wes_palette(n=2, name = "Royal1"), name="Site",  
                     breaks = c("pre", "post"), labels=c("Pre-WNS", "Post-WNS")) +
  scale_shape_discrete(name="Site", 
                     breaks = c("pre", "post"), labels=c("Pre-WNS", "Post-WNS")) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(family="serif", size =12, face = "bold"),
    legend.key.size = unit(1, "cm"),
    legend.position = "right",
    legend.title = element_text(family="serif",size=12, face = "bold"),
    legend.text = element_text(family="serif",size=10),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(family="serif", size=12), 
    axis.text = element_text(family="serif", size =12), 
  )

capture_nmds_graph

ggsave("capture_nmds_graph.png",plot=capture_nmds_graph, dpi=300, dev='png', height=14, width=22, units="cm")  



#### MRPP #####
capture$year1 <- as.factor(capture$year)
capture$year2 <- capture$year1
capture[capture$year1=="2000",]$year2 <- "2003"
capture[capture$year1=="2001",]$year2 <- "2003"
capture[capture$year1=="2002",]$year2 <- "2003"

capture_mrpp <- mrpp(dat= capture[, c(28:38)], 
                  grouping = capture$wns)

capture_mrpp

prewns <- subset(capture, wns=="pre")
capture_mrpp_pre <- mrpp(dat= prewns[, c(28:38)], 
                     grouping = prewns$year2)

capture_mrpp_pre

postwns <- subset(capture, wns=="post")
capture_mrpp_post <- mrpp(dat= postwns[, c(28:38)], 
                     grouping = postwns$year1)

capture_mrpp_post


#######graphs#####################
bat1 <- subset(capture, known >0)
wide <- bat1[, c(5:7, 11:24)]
long <- gather(wide, species, capture, CORA:PESU)

bat11 <- rbind(subset(long, species=="CORA" & capture<0.4), subset(long, species=="EPFU" & capture<0.92))
bat11 <- rbind(bat11, subset(long, species=="LANO" & capture <0.82))
bat11 <- rbind(bat11, subset(long, species=="LABO" & capture <0.80))
bat11 <- rbind(bat11, subset(long, species=="LACI" & capture <0.15))
bat11 <- rbind(bat11, subset(long, species=="MYGR" & capture <0.99))
bat11 <- rbind(bat11, subset(long, species=="MYLE" & capture <0.26))
bat11 <- rbind(bat11, subset(long, species=="MYLU" & capture <2.16))
bat11 <- rbind(bat11, subset(long, species=="MYSE" & capture <0.78))
bat11 <- rbind(bat11, subset(long, species=="MYSO" & capture <0.26))
bat11 <- rbind(bat11, subset(long, species=="PESU" & capture <0.59))

 
capture <- ggplot(data=bat11, aes(x=year, y=capture))+
  geom_point() +
  geom_smooth(method = "gam") +
  geom_vline(aes(xintercept=2010.5), linetype="dashed", color="gold", size=1) +
  xlab("Year") +
  ylab("Capture (number per net hour)") +
  facet_wrap(~species, nrow=3, scales = "free_y") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(family="serif", size =20, face = "bold"),
    legend.key.size = unit(2, "cm"),
    legend.position = c(0.85, 0.2),
    legend.title = element_text(family="serif",size=20),
    legend.text = element_text(family="serif",size=18),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(family="serif", size=20), 
    axis.text = element_text(family="serif", size =20), 
    axis.text.x = element_text(family="serif", size =18, )
  )
 
ggsave("capture.png",plot=capture, dpi=500, dev='png', height=30, width=48, units="cm")

###################gam###########################

#######loops########
sitelist <- names(capture)[16:26]

batlist1 <- names(capture)[c(27, 35)]
batlist2 <- names(capture)[c(28, 32, 34, 37)]
batlist3 <- names(capture)[c(29, 30, 31, 33, 36, 38)]


######year*wns#######

models <- lapply(batlist1, function(x) {
  gam(substitute(i ~ s(julian, k=7) + year*wns, list(i = as.name(x))), family="nb", method = "REML", data = capture)
})

models <- lapply(batlist2, function(x) {
  gam(substitute(i ~ s(julian, k=7) + year*wns, list(i = as.name(x))), family="quasibinomial", method = "REML", data = capture)
})

models <- lapply(batlist3, function(x) {
  gam(substitute(i ~ s(julian, k=7) + year*wns, list(i = as.name(x))), family="quasipoisson", method = "REML", data = capture)
})

model_result <- lapply(models, summary)
sink("capture_result/all/test.txt")
print(model_result)
sink()

### species without wns difference ###

### pesu ###
#### year*wns significant####
fit <- gam(data=subset(capture, wns=="pre"), PESU ~ s(julian, k=7) + year, family = "quasipoisson",  method = "REML")
summary(fit)

fit <- gam(data=subset(capture, wns=="post"), PESU ~ s(julian, k=7) + year, family = "quasipoisson",  method = "REML")
summary(fit)

models <- lapply(sitelist, function(x) {
  gam(substitute(PESU ~ s(julian, k=7) + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="pre"))
})

models <- lapply(sitelist, function(x) {
  gam(substitute(PESU ~ s(julian, k=7) + year * i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="post"))
})

models <- lapply(sitelist, function(x) {
  gam(substitute(PESU ~ s(julian, k=7) + year + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="post"))
})

models <- lapply(sitelist, function(x) {
  gam(substitute(PESU ~ s(julian, k=7) + s(year, k=7) + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="post"))
})

model_result <- lapply(models, summary)
sink("result/pesu/test.txt")
print(model_result)
sink()


### myso ###
#### year*wns significant####
fit <- gam(data=subset(capture, wns=="pre"), MYSO ~ s(julian, k=7) + year, family = "quasibinomial",  method = "REML")
summary(fit)

fit <- gam(data=subset(capture, wns=="post"), MYSO ~ s(julian, k=7) + year, family = "quasibinomial",  method = "REML")
summary(fit)

models <- lapply(sitelist, function(x) {
  gam(substitute(MYSO ~ s(julian, k=7) + year * i, list(i = as.name(x))), family="quasibinomial", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/MYSO/MYSO-1-pre-yearXsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYSO ~ s(julian, k=7) + year + i, list(i = as.name(x))), family="quasibinomial", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/MYSO/MYSO-2-pre-year+site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYSO ~ s(julian, k=7) + s(year, k=7) + i, list(i = as.name(x))), family="quasibinomial", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/MYSO/MYSO-3-pre-gamsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYSO ~ s(julian, k=7) + year * i, list(i = as.name(x))), family="quasibinomial", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/MYSO/MYSO-4-post-yearXsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYSO ~ s(julian, k=7) + year + i, list(i = as.name(x))), family="quasibinomial", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/MYSO/MYSO-5-post-year+site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYSO ~ s(julian, k=7) + s(year, k=7) + i, list(i = as.name(x))), family="quasibinomial", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/MYSO/MYSO-6-post-gamsite.txt")
print(model_result)
sink()

### MYSE ###
#### year*wns significant####
fit <- gam(data=subset(capture, wns=="pre"), MYSE ~ s(julian, k=7) + year, family = "quasipoisson",  method = "REML")
summary(fit)

fit <- gam(data=subset(capture, wns=="post"), MYSE ~ s(julian, k=7) + year, family = "quasipoisson",  method = "REML")
summary(fit)

models <- lapply(sitelist, function(x) {
  gam(substitute(MYSE ~ s(julian, k=7) +  i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/MYSE/MYSE-0-pre-site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYSE ~ s(julian, k=7) + year * i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/MYSE/MYSE-1-pre-yearXsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYSE ~ s(julian, k=7) + year + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/MYSE/MYSE-2-pre-year+site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYSE ~ s(julian, k=7) + s(year, k=7) + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/MYSE/MYSE-3-pre-gamsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYSE ~ s(julian, k=7) + year * i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/MYSE/MYSE-4-post-yearXsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYSE ~ s(julian, k=7) + year + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/MYSE/MYSE-5-post-year+site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYSE ~ s(julian, k=7) + s(year, k=7) + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/MYSE/MYSE-6-post-gamsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYSE ~ s(julian, k=7) + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/MYSE/MYSE-7-post-site.txt")
print(model_result)
sink()

### MYLU ###
#### year*wns significant####
fit <- gamm(data=subset(capture, wns=="pre"), MYLU ~ s(julian, k=7) + year, random=list(site=~1), family = "nb",  method = "REML")
lapply(fit, summary)

fit <- gamm(data=subset(capture, wns=="post"), MYLU ~ s(julian, k=7) + year, random=list(site=~1), family = "nb",  method = "REML")
lapply(fit, summary)

fit <- gam(data=subset(capture, wns=="pre"), MYLU ~ s(julian, k=7) + year, family = "nb",  method = "REML")
summary(fit)

fit <- gam(data=subset(capture, wns=="post"), MYLU ~ s(julian, k=7) + year, family = "nb",  method = "REML")
summary(fit)

models <- lapply(sitelist, function(x) {
  gam(substitute(MYLU ~ s(julian, k=7) +  i, list(i = as.name(x))), family="nb", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/MYLU/MYLU-0-pre-site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYLU ~ s(julian, k=7) + year * i, list(i = as.name(x))), family="nb", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/MYLU/MYLU-1-pre-yearXsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYLU ~ s(julian, k=7) + year + i, list(i = as.name(x))), family="nb", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/MYLU/MYLU-2-pre-year+site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYLU ~ s(julian, k=7) + s(year, k=7) + i, list(i = as.name(x))), family="nb", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/MYLU/MYLU-3-pre-gamsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYLU ~ s(julian, k=7) + year * i, list(i = as.name(x))), family="nb", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/MYLU/MYLU-4-post-yearXsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYLU ~ s(julian, k=7) + year + i, list(i = as.name(x))), family="nb", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/MYLU/MYLU-5-post-year+site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYLU ~ s(julian, k=7) + s(year, k=7) + i, list(i = as.name(x))), family="nb", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/MYLU/MYLU-6-post-gamsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYLU ~ s(julian, k=7) + i, list(i = as.name(x))), family="nb", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/MYLU/MYLU-7-post-site.txt")
print(model_result)
sink()

### MYLE ###
#### year*wns significant####
fit <- gamm(data=subset(capture, wns=="pre"), MYLE ~ s(julian, k=7) + year, random=list(site=~1), family = "quasibinomial",  method = "REML")
lapply(fit, summary)

fit <- gamm(data=subset(capture, wns=="post"), MYLE ~ s(julian, k=7) + year, random=list(site=~1), family = "quasibinomial",  method = "REML")
lapply(fit, summary)

fit <- gam(data=subset(capture, wns=="pre"), MYLE ~ s(julian, k=7) + year, family = "quasibinomial",  method = "REML")
summary(fit)

fit <- gam(data=subset(capture, wns=="post"), MYLE ~ s(julian, k=7) + year, family = "quasibinomial",  method = "REML")
summary(fit)

models <- lapply(sitelist, function(x) {
  gam(substitute(MYLE ~ s(julian, k=7) +  i, list(i = as.name(x))), family="quasibinomial", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/MYLE/MYLE-0-pre-site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYLE ~ s(julian, k=7) + year * i, list(i = as.name(x))), family="quasibinomial", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/MYLE/MYLE-1-pre-yearXsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYLE ~ s(julian, k=7) + year + i, list(i = as.name(x))), family="quasibinomial", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/MYLE/MYLE-2-pre-year+site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYLE ~ s(julian, k=7) + s(year, k=7) + i, list(i = as.name(x))), family="quasibinomial", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/MYLE/MYLE-3-pre-gamsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYLE ~ s(julian, k=7) + year * i, list(i = as.name(x))), family="quasibinomial", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/MYLE/MYLE-4-post-yearXsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYLE ~ s(julian, k=7) + year + i, list(i = as.name(x))), family="quasibinomial", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/MYLE/MYLE-5-post-year+site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYLE ~ s(julian, k=7) + s(year, k=7) + i, list(i = as.name(x))), family="quasibinomial", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/MYLE/MYLE-6-post-gamsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYLE ~ s(julian, k=7) + i, list(i = as.name(x))), family="quasibinomial", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/MYLE/MYLE-7-post-site.txt")
print(model_result)
sink()

### MYGR ###
#### year*wns significant####
fit <- gamm(data=subset(capture, wns=="pre"), MYGR ~ s(julian, k=7) + year, random=list(site=~1), family = "quasipoisson",  method = "REML")
lapply(fit, summary)

fit <- gamm(data=subset(capture, wns=="post"), MYGR ~ s(julian, k=7) + year, random=list(site=~1), family = "quasipoisson",  method = "REML")
lapply(fit, summary)

fit <- gam(data=subset(capture, wns=="pre"), MYGR ~ s(julian, k=7) + year, family = "quasipoisson",  method = "REML")
summary(fit)

fit <- gam(data=subset(capture, wns=="post"), MYGR ~ s(julian, k=7) + year, family = "quasipoisson",  method = "REML")
summary(fit)

models <- lapply(sitelist, function(x) {
  gam(substitute(MYGR ~ s(julian, k=7) +  i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/MYGR/MYGR-0-pre-site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYGR ~ s(julian, k=7) + year * i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/MYGR/MYGR-1-pre-yearXsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYGR ~ s(julian, k=7) + year + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/MYGR/MYGR-2-pre-year+site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYGR ~ s(julian, k=7) + s(year, k=7) + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/MYGR/MYGR-3-pre-gamsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYGR ~ s(julian, k=7) + year * i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/MYGR/MYGR-4-post-yearXsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYGR ~ s(julian, k=7) + year + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/MYGR/MYGR-5-post-year+site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYGR ~ s(julian, k=7) + s(year, k=7) + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/MYGR/MYGR-6-post-gamsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(MYGR ~ s(julian, k=7) + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/MYGR/MYGR-7-post-site.txt")
print(model_result)
sink()


### LACI ###
#### year*wns NOT significant####
fit <- gamm(data=capture, LACI ~ s(julian, k=7) + year, random=list(site=~1), family = "quasibinomial",  method = "REML")
lapply(fit, summary)

fit <- gamm(data=subset(capture, wns=="pre"), LACI ~ s(julian, k=7) + year, random=list(site=~1), family = "quasibinomial",  method = "REML")
lapply(fit, summary)

fit <- gamm(data=subset(capture, wns=="post"), LACI ~ s(julian, k=7) + year, random=list(site=~1), family = "quasibinomial",  method = "REML")
lapply(fit, summary)

fit <- gam(data=subset(capture, wns=="pre"), LACI ~ s(julian, k=7) + year, family = "quasibinomial",  method = "REML")
summary(fit)

fit <- gam(data=subset(capture, wns=="post"), LACI ~ s(julian, k=7) + year, family = "quasibinomial",  method = "REML")
summary(fit)

models <- lapply(sitelist, function(x) {
  gam(substitute(LACI ~ s(julian, k=7) +  i, list(i = as.name(x))), family="quasibinomial", method = "REML", data=capture)
})
model_result <- lapply(models, summary)
sink("result/LACI/LACI-0-ALL-site.txt")
print(model_result)
sink()

### LABO ###
#### year*wns significant####
fit <- gamm(data=subset(capture, wns=="pre"), LABO ~ s(julian, k=7) + year, random=list(site=~1), family = "quasipoisson",  method = "REML")
lapply(fit, summary)

fit <- gamm(data=subset(capture, wns=="post"), LABO ~ s(julian, k=7) + year, random=list(site=~1), family = "quasipoisson",  method = "REML")
lapply(fit, summary)

fit <- gam(data=subset(capture, wns=="pre"), LABO ~ s(julian, k=7) + year, family = "quasipoisson",  method = "REML")
summary(fit)

fit <- gam(data=subset(capture, wns=="post"), LABO ~ s(julian, k=7) + year, family = "quasipoisson",  method = "REML")
summary(fit)

models <- lapply(sitelist, function(x) {
  gam(substitute(LABO ~ s(julian, k=7) +  i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/LABO/LABO-0-pre-site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(LABO ~ s(julian, k=7) + year * i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/LABO/LABO-1-pre-yearXsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(LABO ~ s(julian, k=7) + year + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/LABO/LABO-2-pre-year+site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(LABO ~ s(julian, k=7) + s(year, k=7) + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/LABO/LABO-3-pre-gamsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(LABO ~ s(julian, k=7) + year * i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/LABO/LABO-4-post-yearXsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(LABO ~ s(julian, k=7) + year + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/LABO/LABO-5-post-year+site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(LABO ~ s(julian, k=7) + s(year, k=7) + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/LABO/LABO-6-post-gamsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(LABO ~ s(julian, k=7) + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/LABO/LABO-7-post-site.txt")
print(model_result)
sink()

### LANO ###
#### year*wns NOT significant####
fit <- gamm(data=capture, LANO ~ s(julian, k=7) + year, random=list(site=~1), family = "quasipoisson",  method = "REML")
lapply(fit, summary)

fit <- gamm(data=subset(capture, wns=="pre"), LANO ~ s(julian, k=7) + year, random=list(site=~1), family = "quasipoisson",  method = "REML")
lapply(fit, summary)

fit <- gamm(data=subset(capture, wns=="post"), LANO ~ s(julian, k=7) + year, random=list(site=~1), family = "quasipoisson",  method = "REML")
lapply(fit, summary)

fit <- gam(data=subset(capture, wns=="pre"), LANO ~ s(julian, k=7) + year, family = "quasipoisson",  method = "REML")
summary(fit)

fit <- gam(data=subset(capture, wns=="post"), LANO ~ s(julian, k=7) + year, family = "quasipoisson",  method = "REML")
summary(fit)

models <- lapply(sitelist, function(x) {
  gam(substitute(LANO ~ s(julian, k=7) +  i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/LANO/LANO-0-pre-site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(LANO ~ s(julian, k=7) + year * i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/LANO/LANO-1-pre-yearXsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(LANO ~ s(julian, k=7) + year + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/LANO/LANO-2-pre-year+site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(LANO ~ s(julian, k=7) + s(year, k=7) + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/LANO/LANO-3-pre-gamsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(LANO ~ s(julian, k=7) + year * i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/LANO/LANO-4-post-yearXsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(LANO ~ s(julian, k=7) + year + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/LANO/LANO-5-post-year+site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(LANO ~ s(julian, k=7) + s(year, k=7) + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/LANO/LANO-6-post-gamsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(LANO ~ s(julian, k=7) + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/LANO/LANO-7-post-site.txt")
print(model_result)
sink()

### EPFU ###
#### year*wns NOT significant####
fit <- gamm(data=capture, EPFU ~ s(julian, k=7) + year, random=list(site=~1), family = "quasipoisson",  method = "REML")
lapply(fit, summary)

fit <- gamm(data=subset(capture, wns=="pre"), EPFU ~ s(julian, k=7) + year, random=list(site=~1), family = "quasipoisson",  method = "REML")
lapply(fit, summary)

fit <- gamm(data=subset(capture, wns=="post"), EPFU ~ s(julian, k=7) + year, random=list(site=~1), family = "quasipoisson",  method = "REML")
lapply(fit, summary)

fit <- gam(data=subset(capture, wns=="pre"), EPFU ~ s(julian, k=7) + year, family = "quasipoisson",  method = "REML")
summary(fit)

fit <- gam(data=subset(capture, wns=="post"), EPFU ~ s(julian, k=7) + year, family = "quasipoisson",  method = "REML")
summary(fit)

models <- lapply(sitelist, function(x) {
  gam(substitute(EPFU ~ s(julian, k=7) +  i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/EPFU/EPFU-0-pre-site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(EPFU ~ s(julian, k=7) + year * i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/EPFU/EPFU-1-pre-yearXsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(EPFU ~ s(julian, k=7) + year + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/EPFU/EPFU-2-pre-year+site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(EPFU ~ s(julian, k=7) + s(year, k=7) + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/EPFU/EPFU-3-pre-gamsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(EPFU ~ s(julian, k=7) + year * i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/EPFU/EPFU-4-post-yearXsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(EPFU ~ s(julian, k=7) + year + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/EPFU/EPFU-5-post-year+site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(EPFU ~ s(julian, k=7) + s(year, k=7) + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/EPFU/EPFU-6-post-gamsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(EPFU ~ s(julian, k=7) + i, list(i = as.name(x))), family="quasipoisson", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/EPFU/EPFU-7-post-site.txt")
print(model_result)
sink()

### CORA ###
#### year*wns NOT significant####
fit <- gamm(data=capture, CORA ~ s(julian, k=7) + year, random=list(site=~1), family = "quasibinomial",  method = "REML")
lapply(fit, summary)

fit <- gamm(data=subset(capture, wns=="pre"), CORA ~ s(julian, k=7) + year, random=list(site=~1), family = "quasibinomial",  method = "REML")
lapply(fit, summary)

fit <- gamm(data=subset(capture, wns=="post"), CORA ~ s(julian, k=7) + year, random=list(site=~1), family = "quasibinomial",  method = "REML")
lapply(fit, summary)

fit <- gam(data=capture, CORA ~ s(julian, k=7) + year, family = "quasibinomial",  method = "REML")
summary(fit)

fit <- gam(data=subset(capture, wns=="pre"), CORA ~ s(julian, k=7) + year, family = "quasibinomial",  method = "REML")
summary(fit)

fit <- gam(data=subset(capture, wns=="post"), CORA ~ s(julian, k=7) + year, family = "quasibinomial",  method = "REML")
summary(fit)

models <- lapply(sitelist, function(x) {
  gam(substitute(CORA ~ s(julian, k=7) +  i, list(i = as.name(x))), family="quasibinomial", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/CORA/CORA-0-pre-site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(CORA ~ s(julian, k=7) + year * i, list(i = as.name(x))), family="quasibinomial", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/CORA/CORA-1-pre-yearXsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(CORA ~ s(julian, k=7) + year + i, list(i = as.name(x))), family="quasibinomial", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/CORA/CORA-2-pre-year+site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(CORA ~ s(julian, k=7) + s(year, k=7) + i, list(i = as.name(x))), family="quasibinomial", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/CORA/CORA-3-pre-gamsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(CORA ~ s(julian, k=7) + year * i, list(i = as.name(x))), family="quasibinomial", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/CORA/CORA-4-post-yearXsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(CORA ~ s(julian, k=7) + year + i, list(i = as.name(x))), family="quasibinomial", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/CORA/CORA-5-post-year+site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(CORA ~ s(julian, k=7) + s(year, k=7) + i, list(i = as.name(x))), family="quasibinomial", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/CORA/CORA-6-post-gamsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(CORA ~ s(julian, k=7) + i, list(i = as.name(x))), family="quasibinomial", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/CORA/CORA-7-post-site.txt")
print(model_result)
sink()

### all ###
#### year*wns NOT significant####
fit <- gamm(data=capture, all ~ s(julian, k=7) + year, random=list(site=~1), family = "nb",  method = "REML")
lapply(fit, summary)

fit <- gamm(data=subset(capture, wns=="pre"), all ~ s(julian, k=7) + year, random=list(site=~1), family = "nb",  method = "REML")
lapply(fit, summary)

fit <- gamm(data=subset(capture, wns=="post"), all ~ s(julian, k=7) + year, random=list(site=~1), family = "nb",  method = "REML")
lapply(fit, summary)

fit <- gam(data=capture, all ~ s(julian, k=7) + year, family = "nb",  method = "REML")
summary(fit)

fit <- gam(data=subset(capture, wns=="pre"), all ~ s(julian, k=7) + year, family = "nb",  method = "REML")
summary(fit)

fit <- gam(data=subset(capture, wns=="post"), all ~ s(julian, k=7) + year, family = "nb",  method = "REML")
summary(fit)

models <- lapply(sitelist, function(x) {
  gam(substitute(all ~ s(julian, k=7) +  i, list(i = as.name(x))), family="nb", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/all/all-0-pre-site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(all ~ s(julian, k=7) + year * i, list(i = as.name(x))), family="nb", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/all/all-1-pre-yearXsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(all ~ s(julian, k=7) + year + i, list(i = as.name(x))), family="nb", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/all/all-2-pre-year+site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(all ~ s(julian, k=7) + s(year, k=7) + i, list(i = as.name(x))), family="nb", method = "REML", data=subset(capture, wns=="pre"))
})
model_result <- lapply(models, summary)
sink("result/all/all-3-pre-gamsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(all ~ s(julian, k=7) + year * i, list(i = as.name(x))), family="nb", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/all/all-4-post-yearXsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(all ~ s(julian, k=7) + year + i, list(i = as.name(x))), family="nb", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/all/all-5-post-year+site.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(all ~ s(julian, k=7) + s(year, k=7) + i, list(i = as.name(x))), family="nb", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/all/all-6-post-gamsite.txt")
print(model_result)
sink()

models <- lapply(sitelist, function(x) {
  gam(substitute(all ~ s(julian, k=7) + i, list(i = as.name(x))), family="nb", method = "REML", data=subset(capture, wns=="post"))
})
model_result <- lapply(models, summary)
sink("result/all/all-7-post-site.txt")
print(model_result)
sink()
