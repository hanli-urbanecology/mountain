library(ggplot2)
library(mgcv)
library(tidyr)
library(corrplot)
library(vegan)
library(ggforce)
library("RColorBrewer")

require(RColorBrewer)
install.packages('ggforce')

tran <- read.csv("tran.csv")
tran_site <- read.csv("tran_near.csv")

tran <- merge(tran, tran_site, by.x = "site", by.y = "route")
write.csv(tran, "tran2.csv")

tran <-read.csv("tran2.csv")

#### nmds test ####
bat <- as.matrix(tran[, c(10, 11, 12, 13, 14, 16, 19, 20, 21)])
bat.mds <- metaMDS(bat, k=2, distance  = "euclidean", trymax = 10000) 
bat.mds_tran <- bat.mds
bat.mds_tran
plot(bat.mds_tran, type = "t")

data.scores <- as.data.frame(scores(bat.mds_tran))  
data.scores$X <- as.numeric(rownames(data.scores))

tran_graph <- merge(data.scores, tran)
write.csv(tran_graph, "tran_graph.csv")


ggplot(data= subset(tran_graph, NMDS1<0.2&NMDS1>-0.15&NMDS2<0.2&NMDS2>-0.2), aes(x=NMDS1, y=NMDS2)) + 
  geom_mark_ellipse(aes(color = year1)) + 
  #geom_polygon(group = tran_graph$year1, fill = tran_graph$year1, alpha = 0.5) + 
  geom_point(aes(color=year1)) +
  xlim(-0.18, 0.21)+
  ylim(-0.22, 0.22)+
  scale_fill_brewer(palette="RdGy")+
  scale_color_brewer(palette="RdGy")


#### MRPP #####
tran$year1 <- as.factor(tran$year)
tran_mrpp <- mrpp(dat= tran[, c(10, 11, 12, 13, 14, 16, 19, 20, 21)], 
                  grouping = tran$year1)

#### transect characteristics #####

elev_long <- read.csv("elev_long.csv")
elev_mean <- aggregate(elev_long[3], by=list(elev_long$route), FUN = mean)
elev_max <- aggregate(elev~route, elev_long, max)
elev_min <- aggregate(elev_long$elev, by=list(elev_long$route), FUN = min)

names(elev_mean)[names(elev_mean) == "Group.1"] <- "route"
names(elev_mean)[names(elev_mean) == "elev"] <- "elevmean"
names(elev_max)[names(elev_max)=="elev"] <- "elevmax"
names(elev_min)[names(elev_min)=="x"] <-"elevmin"
names(elev_min)[names(elev_min)=="Group.1"] <- "route"

elev <- merge(elev_max, elev_mean, "route")
elev <- merge(elev_min, elev, "route")
elev$elevchange <- elev$elevmax-elev$elevmin

write.csv(elev, "elev.csv")

### transect feature ###
site <- read.csv("tran_near.csv")
site_cor <- cor(site[, 2:8])
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

p.mat <- cor.mtest(site[, 2:8])

corrplot(site_cor, type="upper", order="hclust", 
         p.mat = p.mat, sig.level = 0.01)

###############################################################################



##total plot##
ggplot(data=tran, aes(x=year, y=total)) +
  geom_point(size=1, ) + 
  geom_smooth(method = "gam", formula = y ~ s(x, k=7)) #+scale_y_continuous(limits = c(0, 150))#

## species plot##
wide <- tran[, c(1:22, 24)]
long <- gather(wide, species, transect, CORA:total)
bat11 <- rbind(subset(long, species=="EPFU"), subset(long, species=="LACI"))
bat11 <- rbind(bat11, subset(long, species=="LACI"))
bat11 <- rbind(bat11, subset(long, species=="LANO"))
bat11 <- rbind(bat11, subset(long, species=="NYHU"))
bat11 <- rbind(bat11, subset(long, species=="PESU"))
bat11 <- rbind(bat11, subset(long, species=="TABR"))
bat11 <- rbind(bat11, subset(long, species=="MYLU"))
bat11 <- rbind(bat11, subset(long, species=="MYGR"))
bat11 <- rbind(bat11, subset(long, species=="MYspp"))

transect <- ggplot(data=bat11, aes(x=year, y=transect))+
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k=4)) +
  xlab("Year") +
  ylab("Bat pass (number per transect)") +
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
    axis.text = element_text(family="serif", size =18), 
    axis.text.x = element_text(family="serif", size =18)
  )

ggsave("transect.png",plot=transect, dpi=500, dev='png', height=30, width=48, units="cm")


################################################################################################
#####GAM#####

tran <-read.csv("tran2.csv")

#######loops########
sitelist <- names(tran)[26:34]
batlist <- names(tran)[c(10:14, 16, 19:21, 23, 25)]

#####gamm for year#####
models <- lapply(batlist, function(x) {
  lapply(gamm(data = tran, substitute(i ~ s(julian, k=4) + year, list(i = as.name(x))), random=list(site=~1), family="nb", method = "REML"), summary)
})

sink("transect_result/gamm.txt")
print(models)
sink()

####gam for site#####

###EPFU, year not significant###
models <- lapply(sitelist, function(x) {
  summary(gam(substitute(EPFU ~ s(julian, k=4) + i, list(i = as.name(x))), family="nb", method = "REML", data=tran))
})
sink("transect_result/EPFU.txt")
print(models)
sink()

###LABO, year not significant###
models <- lapply(sitelist, function(x) {
  summary(gam(substitute(LABO ~ s(julian, k=4) + i, list(i = as.name(x))), family="nb", method = "REML", data=tran))
})
sink("transect_result/LABO.txt")
print(models)
sink()

###LACI, year  significant###
models <- lapply(sitelist, function(x) {
  summary(gam(substitute(LACI ~ s(julian, k=4) + year*i, list(i = as.name(x))), family="nb", method = "REML", data=tran))
})
sink("transect_result/LACI1.txt")
print(models)
sink()

models <- lapply(sitelist, function(x) {
  summary(gam(substitute(LACI ~ s(julian, k=4) + s(year, k=4) + i, list(i = as.name(x))), family="nb", method = "REML", data=tran))
})
sink("transect_result/LACI2.txt")
print(models)
sink()

###LANO, year  significant###
models <- lapply(sitelist, function(x) {
  summary(gam(substitute(LANO ~ s(julian, k=4) + year*i, list(i = as.name(x))), family="nb", method = "REML", data=tran))
})
sink("transect_result/LANO1.txt")
print(models)
sink()

models <- lapply(sitelist, function(x) {
  summary(gam(substitute(LANO ~ s(julian, k=4) + s(year, k=4) + i, list(i = as.name(x))), family="nb", method = "REML", data=tran))
})
sink("transect_result/LANO2.txt")
print(models)
sink()

###MYGR, year not significant###
models <- lapply(sitelist, function(x) {
  summary(gam(substitute(MYGR ~ s(julian, k=4) + i, list(i = as.name(x))), family="nb", method = "REML", data=tran))
})
sink("transect_result/MYGR.txt")
print(models)
sink()

###MYLU, year  significant###
models <- lapply(sitelist, function(x) {
  summary(gam(substitute(MYLU ~ s(julian, k=4) + year*i, list(i = as.name(x))), family="nb", method = "REML", data=tran))
})
sink("transect_result/MYLU1.txt")
print(models)
sink()

models <- lapply(sitelist, function(x) {
  summary(gam(substitute(MYLU ~ s(julian, k=4) + s(year, k=4) + i, list(i = as.name(x))), family="nb", method = "REML", data=tran))
})
sink("transect_result/MYLU2.txt")
print(models)
sink()

###NYHU, year not significant###
models <- lapply(sitelist, function(x) {
  summary(gam(substitute(NYHU ~ s(julian, k=4) + i, list(i = as.name(x))), family="nb", method = "REML", data=tran))
})
sink("transect_result/NYHU.txt")
print(models)
sink()

###PESU, year  significant###
models <- lapply(sitelist, function(x) {
  summary(gam(substitute(PESU ~ s(julian, k=4) + year*i, list(i = as.name(x))), family="nb", method = "REML", data=tran))
})
sink("transect_result/PESU1.txt")
print(models)
sink()

models <- lapply(sitelist, function(x) {
  summary(gam(substitute(PESU ~ s(julian, k=4) + s(year, k=4) + i, list(i = as.name(x))), family="nb", method = "REML", data=tran))
})
sink("transect_result/PESU2.txt")
print(models)
sink()

###TABR, year not significant###
models <- lapply(sitelist, function(x) {
  summary(gam(substitute(TABR ~ s(julian, k=4) + i, list(i = as.name(x))), family="nb", method = "REML", data=tran))
})
sink("transect_result/TABR.txt")
print(models)
sink()


###total, year  significant###
models <- lapply(sitelist, function(x) {
  summary(gam(substitute(total ~ s(julian, k=4) + year*i, list(i = as.name(x))), family="nb", method = "REML", data=tran))
})
sink("transect_result/total1.txt")
print(models)
sink()

models <- lapply(sitelist, function(x) {
  summary(gam(substitute(total ~ s(julian, k=4) + s(year, k=4) + i, list(i = as.name(x))), family="nb", method = "REML", data=tran))
})
sink("transect_result/total2.txt")
print(models)
sink()