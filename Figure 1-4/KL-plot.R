
# 22.05.2022

library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)

Model1 = "Model1"
Data1 = read.csv(paste(Model1,".csv",sep=""))
#Data1 = Data1[-(1:300),]

Model2 = "Model2"
Data2 = read.csv(paste(Model2,".csv",sep=""))
#Data2 = Data2[-(1:300),]

Model3 = "Model3"
Data3 = read.csv(paste(Model3,".csv",sep=""))
#Data3 = Data3[-(1:300),]

Model4 = "Model4"
Data4 = read.csv(paste(Model4,".csv",sep=""))
#Data4 = Data4[-(1:300),]

Model5 = "Model5"
Data5 = read.csv(paste(Model5,".csv",sep=""))
#Data5 = Data5[-(1:300),]

Model6 = "Model6"
Data6 = read.csv(paste(Model6,".csv",sep=""))
#Data6 = Data6[-(1:300),]


KLData = rbind(Data1[,c(2,3,6)],Data2[,c(2,3,6)],Data3[,c(2,3,6)],Data4[,c(2,3,6)],Data5[,c(2,3,6)],Data6[,c(2,3,6)])

colnames(KLData) = c("Method","Dimention","LossVal")
KLData["ModelName"] = c(rep("Model1",3600),rep("Model2",3600),rep("Model3",3600),rep("Model4",3600),
                        rep("Model5",3600),rep("Model6",3600))


detach(KLData)

###########################################

names = c(Model1,Model2,Model3,Model4,Model5,Model6)

plist = list()
plist[]

limits = list()

hila = matrix(1:21600,ncol=6)

se = function(x){ sd(x)/sqrt(length(x))}

d1 = KLData[hila[,1],]
d2 = KLData[hila[,2],]
d3 = KLData[hila[,3],]
d4 = KLData[hila[,4],]
d5 = KLData[hila[,5],]
d6 = KLData[hila[,6],]

###################################################

d1_mean = as.data.frame(matrix(0,12,4))
colnames(d1_mean) = c("Method","20","50","100")

d1_se = as.data.frame(matrix(0,12,4))
colnames(d1_se) = c("Method","20","50","100")

d1_mean$Method = factor(c("GLASSO","GLASSO-I","GLASSO-vI"
       ,"GELNET","GELNET-I","GELNET-vI","ROPE","ROPE-I"," ROPE-vI","GSOS","GSOS-I","GSOS-vI"),
       levels=c("GLASSO","GLASSO-I","GLASSO-vI"
       ,"GELNET","GELNET-I","GELNET-vI","ROPE","ROPE-I"," ROPE-vI","GSOS","GSOS-I","GSOS-vI"))

d1_mean[1,2:4] = c(mean(d1[1:100,3]),mean(d1[101:200,3]),mean(d1[201:300,3]))
d1_mean[2,2:4] = c(mean(d1[301:400,3]),mean(d1[401:500,3]),mean(d1[501:600,3]))
d1_mean[3,2:4] = c(mean(d1[601:700,3]),mean(d1[701:800,3]),mean(d1[801:900,3]))
d1_mean[4,2:4] = c(mean(d1[901:1000,3]),mean(d1[1001:1100,3]),mean(d1[1101:1200,3]))
d1_mean[5,2:4] = c(mean(d1[1201:1300,3]),mean(d1[1301:1400,3]),mean(d1[1401:1500,3]))
d1_mean[6,2:4] = c(mean(d1[1501:1600,3]),mean(d1[1601:1700,3]),mean(d1[1701:1800,3]))
d1_mean[7,2:4] = c(mean(d1[1801:1900,3]),mean(d1[1901:2000,3]),mean(d1[2001:2100,3]))
d1_mean[8,2:4] = c(mean(d1[2101:2200,3]),mean(d1[2201:2300,3]),mean(d1[2301:2400,3]))
d1_mean[9,2:4] = c(mean(d1[2401:2500,3]),mean(d1[2501:2600,3]),mean(d1[2601:2700,3]))
d1_mean[10,2:4] = c(mean(d1[2701:2800,3]),mean(d1[2801:2900,3]),mean(d1[2901:3000,3]))
d1_mean[11,2:4] = c(mean(d1[3001:3100,3]),mean(d1[3101:3200,3]),mean(d1[3201:3300,3]))
d1_mean[12,2:4] = c(mean(d1[3301:3400,3]),mean(d1[3401:3500,3]),mean(d1[3501:3600,3]))

d1_se$Method = factor(c("GLASSO","GLASSO-I","GLASSO-vI"
       ,"GELNET","GELNET-I","GELNET-vI","ROPE","ROPE-I"," ROPE-vI","GSOS","GSOS-I","GSOS-vI"),
       levels=c("GLASSO","GLASSO-I","GLASSO-vI"
       ,"GELNET","GELNET-I","GELNET-vI","ROPE","ROPE-I"," ROPE-vI","GSOS","GSOS-I","GSOS-vI"))

d1_se[1,2:4] = c(se(d1[1:100,3]),se(d1[101:200,3]),se(d1[201:300,3]))
d1_se[2,2:4] = c(se(d1[301:400,3]),se(d1[401:500,3]),se(d1[501:600,3]))
d1_se[3,2:4] = c(se(d1[601:700,3]),se(d1[701:800,3]),se(d1[801:900,3]))
d1_se[4,2:4] = c(se(d1[901:1000,3]),se(d1[1001:1100,3]),se(d1[1101:1200,3]))
d1_se[5,2:4] = c(se(d1[1201:1300,3]),se(d1[1301:1400,3]),se(d1[1401:1500,3]))
d1_se[6,2:4] = c(se(d1[1501:1600,3]),se(d1[1601:1700,3]),se(d1[1701:1800,3]))
d1_se[7,2:4] = c(se(d1[1801:1900,3]),se(d1[1901:2000,3]),se(d1[2001:2100,3]))
d1_se[8,2:4] = c(se(d1[2101:2200,3]),se(d1[2201:2300,3]),se(d1[2301:2400,3]))
d1_se[9,2:4] = c(se(d1[2401:2500,3]),se(d1[2501:2600,3]),se(d1[2601:2700,3]))
d1_se[10,2:4] = c(se(d1[2701:2800,3]),se(d1[2801:2900,3]),se(d1[2901:3000,3]))
d1_se[11,2:4] = c(se(d1[3001:3100,3]),se(d1[3101:3200,3]),se(d1[3201:3300,3]))
d1_se[12,2:4] = c(se(d1[3301:3400,3]),se(d1[3401:3500,3]),se(d1[3501:3600,3]))

d1f_mean = melt(d1_mean, id.vars=c("Method"), variable.name = "Dimension", value.name="LossMean")

d1f_se = melt(d1_se, id.vars=c("Method"), variable.name = "Dimension", value.name="LossMean")

limits1 = aes(ymax = d1f_mean[,"LossMean"] + d1f_se[,"LossMean"], ymin=d1f_mean[,"LossMean"] - d1f_se[,"LossMean"])

p1 =  ggplot(d1f_mean, aes(Dimension, LossMean, fill = Method)) + 
  geom_bar(position="dodge", stat="identity") + 
  geom_pointrange(limits1, position=position_dodge(.9),fatten=0.1) +
  #scale_fill_brewer(palette = "Set3") + 
  scale_fill_manual(values = c("#9DAEE1FF", "#879DE1FF","#718DE1FF",
     "#6584E1FF","#4F73E1FF", "#436BE1FF", "#2D5AE1FF","#1649E1FF","#0038E1FF","#FF7F50", "#A52A2A", "#800000")) +
  theme_gray() +
  labs(y = "Loss mean", x = "Dimension", fill = NULL) +
  theme(legend.position='none') +
  ggtitle("Model 1") +
  coord_cartesian(ylim=c(0,30))

plist[[names[1]]] = p1

###########################

d2_mean = as.data.frame(matrix(0,12,4))
colnames(d2_mean) = c("Method","20","50","100")

d2_se = as.data.frame(matrix(0,12,4))
colnames(d2_se) = c("Method","20","50","100")

d2_mean$Method = factor(c("GLASSO","GLASSO-I","GLASSO-vI"
       ,"GELNET","GELNET-I","GELNET-vI","ROPE","ROPE-I"," ROPE-vI","GSOS","GSOS-I","GSOS-vI"),
       levels=c("GLASSO","GLASSO-I","GLASSO-vI"
       ,"GELNET","GELNET-I","GELNET-vI","ROPE","ROPE-I"," ROPE-vI","GSOS","GSOS-I","GSOS-vI"))

d2_mean[1,2:4] = c(mean(d2[1:100,3]),mean(d2[101:200,3]),mean(d2[201:300,3]))
d2_mean[2,2:4] = c(mean(d2[301:400,3]),mean(d2[401:500,3]),mean(d2[501:600,3]))
d2_mean[3,2:4] = c(mean(d2[601:700,3]),mean(d2[701:800,3]),mean(d2[801:900,3]))
d2_mean[4,2:4] = c(mean(d2[901:1000,3]),mean(d2[1001:1100,3]),mean(d2[1101:1200,3]))
d2_mean[5,2:4] = c(mean(d2[1201:1300,3]),mean(d2[1301:1400,3]),mean(d2[1401:1500,3]))
d2_mean[6,2:4] = c(mean(d2[1501:1600,3]),mean(d2[1601:1700,3]),mean(d2[1701:1800,3]))
d2_mean[7,2:4] = c(mean(d2[1801:1900,3]),mean(d2[1901:2000,3]),mean(d2[2001:2100,3]))
d2_mean[8,2:4] = c(mean(d2[2101:2200,3]),mean(d2[2201:2300,3]),mean(d2[2301:2400,3]))
d2_mean[9,2:4] = c(mean(d2[2401:2500,3]),mean(d2[2501:2600,3]),mean(d2[2601:2700,3]))
d2_mean[10,2:4] = c(mean(d2[2701:2800,3]),mean(d2[2801:2900,3]),mean(d2[2901:3000,3]))
d2_mean[11,2:4] = c(mean(d2[3001:3100,3]),mean(d2[3101:3200,3]),mean(d2[3201:3300,3]))
d2_mean[12,2:4] = c(mean(d2[3301:3400,3]),mean(d2[3401:3500,3]),mean(d2[3501:3600,3]))


d2_se$Method = factor(c("GLASSO","GLASSO-I","GLASSO-vI"
       ,"GELNET","GELNET-I","GELNET-vI","ROPE","ROPE-I"," ROPE-vI","GSOS","GSOS-I","GSOS-vI"),
       levels=c("GLASSO","GLASSO-I","GLASSO-vI"
       ,"GELNET","GELNET-I","GELNET-vI","ROPE","ROPE-I"," ROPE-vI","GSOS","GSOS-I","GSOS-vI"))

d2_se[1,2:4] = c(se(d2[1:100,3]),se(d2[101:200,3]),se(d2[201:300,3]))
d2_se[2,2:4] = c(se(d2[301:400,3]),se(d2[401:500,3]),se(d2[501:600,3]))
d2_se[3,2:4] = c(se(d2[601:700,3]),se(d2[701:800,3]),se(d2[801:900,3]))
d2_se[4,2:4] = c(se(d2[901:1000,3]),se(d2[1001:1100,3]),se(d2[1101:1200,3]))
d2_se[5,2:4] = c(se(d2[1201:1300,3]),se(d2[1301:1400,3]),se(d2[1401:1500,3]))
d2_se[6,2:4] = c(se(d2[1501:1600,3]),se(d2[1601:1700,3]),se(d2[1701:1800,3]))
d2_se[7,2:4] = c(se(d2[1801:1900,3]),se(d2[1901:2000,3]),se(d2[2001:2100,3]))
d2_se[8,2:4] = c(se(d2[2101:2200,3]),se(d2[2201:2300,3]),se(d2[2301:2400,3]))
d2_se[9,2:4] = c(se(d2[2401:2500,3]),se(d2[2501:2600,3]),se(d2[2601:2700,3]))
d2_se[10,2:4] = c(se(d2[2701:2800,3]),se(d2[2801:2900,3]),se(d2[2901:3000,3]))
d2_se[11,2:4] = c(se(d2[3001:3100,3]),se(d2[3101:3200,3]),se(d2[3201:3300,3]))
d2_se[12,2:4] = c(se(d2[3301:3400,3]),se(d2[3401:3500,3]),se(d2[3501:3600,3]))



d2f_mean = melt(d2_mean, id.vars=c("Method"), variable.name = "Dimension", value.name="LossMean")

d2f_se = melt(d2_se, id.vars=c("Method"), variable.name = "Dimension", value.name="LossMean")

limits2 = aes(ymax = d2f_mean[,"LossMean"] + d2f_se[,"LossMean"], ymin=d2f_mean[,"LossMean"] - d2f_se[,"LossMean"])

p2 =  ggplot(d2f_mean, aes(Dimension, LossMean, fill = Method)) + 
  geom_bar(position="dodge", stat="identity") + 
  geom_pointrange(limits2, position=position_dodge(.9),fatten=0.1) +
  #scale_fill_brewer(palette = "Set3") + 
  scale_fill_manual(values = c("#9DAEE1FF", "#879DE1FF","#718DE1FF",
     "#6584E1FF","#4F73E1FF", "#436BE1FF", "#2D5AE1FF","#1649E1FF","#0038E1FF","#FF7F50", "#A52A2A", "#800000")) +
  
  theme_gray() +
  labs(y = " ", x = "Dimension", fill = NULL) + 
  theme(legend.position='none') +
  ggtitle("Model 2") +
  coord_cartesian(ylim=c(0,30))


plist[[names[2]]] = p2

###########################

d3_mean = as.data.frame(matrix(0,12,4))
colnames(d3_mean) = c("Method","20","50","100")

d3_se = as.data.frame(matrix(0,12,4))
colnames(d3_se) = c("Method","20","50","100")

d3_mean$Method = factor(c("GLASSO","GLASSO-I","GLASSO-vI"
       ,"GELNET","GELNET-I","GELNET-vI","ROPE","ROPE-I"," ROPE-vI","GSOS","GSOS-I","GSOS-vI"),
       levels=c("GLASSO","GLASSO-I","GLASSO-vI"
       ,"GELNET","GELNET-I","GELNET-vI","ROPE","ROPE-I"," ROPE-vI","GSOS","GSOS-I","GSOS-vI"))

d3_mean[1,2:4] = c(mean(d3[1:100,3]),mean(d3[101:200,3]),mean(d3[201:300,3]))
d3_mean[2,2:4] = c(mean(d3[301:400,3]),mean(d3[401:500,3]),mean(d3[501:600,3]))
d3_mean[3,2:4] = c(mean(d3[601:700,3]),mean(d3[701:800,3]),mean(d3[801:900,3]))
d3_mean[4,2:4] = c(mean(d3[901:1000,3]),mean(d3[1001:1100,3]),mean(d3[1101:1200,3]))
d3_mean[5,2:4] = c(mean(d3[1201:1300,3]),mean(d3[1301:1400,3]),mean(d3[1401:1500,3]))
d3_mean[6,2:4] = c(mean(d3[1501:1600,3]),mean(d3[1601:1700,3]),mean(d3[1701:1800,3]))
d3_mean[7,2:4] = c(mean(d3[1801:1900,3]),mean(d3[1901:2000,3]),mean(d3[2001:2100,3]))
d3_mean[8,2:4] = c(mean(d3[2101:2200,3]),mean(d3[2201:2300,3]),mean(d3[2301:2400,3]))
d3_mean[9,2:4] = c(mean(d3[2401:2500,3]),mean(d3[2501:2600,3]),mean(d3[2601:2700,3]))
d3_mean[10,2:4] = c(mean(d3[2701:2800,3]),mean(d3[2801:2900,3]),mean(d3[2901:3000,3]))
d3_mean[11,2:4] = c(mean(d3[3001:3100,3]),mean(d3[3101:3200,3]),mean(d3[3201:3300,3]))
d3_mean[12,2:4] = c(mean(d3[3301:3400,3]),mean(d3[3401:3500,3]),mean(d3[3501:3600,3]))


d3_se$Method = factor(c("GLASSO","GLASSO-I","GLASSO-vI"
       ,"GELNET","GELNET-I","GELNET-vI","ROPE","ROPE-I"," ROPE-vI","GSOS","GSOS-I","GSOS-vI"),
       levels=c("GLASSO","GLASSO-I","GLASSO-vI"
       ,"GELNET","GELNET-I","GELNET-vI","ROPE","ROPE-I"," ROPE-vI","GSOS","GSOS-I","GSOS-vI"))

d3_se[1,2:4] = c(se(d3[1:100,3]),se(d3[101:200,3]),se(d3[201:300,3]))
d3_se[2,2:4] = c(se(d3[301:400,3]),se(d3[401:500,3]),se(d3[501:600,3]))
d3_se[3,2:4] = c(se(d3[601:700,3]),se(d3[701:800,3]),se(d3[801:900,3]))
d3_se[4,2:4] = c(se(d3[901:1000,3]),se(d3[1001:1100,3]),se(d3[1101:1200,3]))
d3_se[5,2:4] = c(se(d3[1201:1300,3]),se(d3[1301:1400,3]),se(d3[1401:1500,3]))
d3_se[6,2:4] = c(se(d3[1501:1600,3]),se(d3[1601:1700,3]),se(d3[1701:1800,3]))
d3_se[7,2:4] = c(se(d3[1801:1900,3]),se(d3[1901:2000,3]),se(d3[2001:2100,3]))
d3_se[8,2:4] = c(se(d3[2101:2200,3]),se(d3[2201:2300,3]),se(d3[2301:2400,3]))
d3_se[9,2:4] = c(se(d3[2401:2500,3]),se(d3[2501:2600,3]),se(d3[2601:2700,3]))
d3_se[10,2:4] = c(se(d3[2701:2800,3]),se(d3[2801:2900,3]),se(d3[2901:3000,3]))
d3_se[11,2:4] = c(se(d3[3001:3100,3]),se(d3[3101:3200,3]),se(d3[3201:3300,3]))
d3_se[12,2:4] = c(se(d3[3301:3400,3]),se(d3[3401:3500,3]),se(d3[3501:3600,3]))


d3f_mean = melt(d3_mean, id.vars=c("Method"), variable.name = "Dimension", value.name="LossMean")

d3f_se = melt(d3_se, id.vars=c("Method"), variable.name = "Dimension", value.name="LossMean")

limits3 = aes(ymax = d3f_mean[,"LossMean"] + d3f_se[,"LossMean"], ymin=d3f_mean[,"LossMean"] - d3f_se[,"LossMean"])

p3 =  ggplot(d3f_mean, aes(Dimension, LossMean, fill = Method)) + 
  geom_bar(position="dodge", stat="identity") + 
  geom_pointrange(limits3, position=position_dodge(.9),fatten=0.1) +
  #scale_fill_brewer(palette = "Set3") + 
  scale_fill_manual(values = c("#9DAEE1FF", "#879DE1FF","#718DE1FF",
     "#6584E1FF","#4F73E1FF", "#436BE1FF", "#2D5AE1FF","#1649E1FF","#0038E1FF","#FF7F50", "#A52A2A", "#800000")) +
  theme_gray() +
  labs(y = " ", x = "Dimension", fill = NULL) + 
  theme(legend.position='none') +
  ggtitle("Model 3") +
  coord_cartesian(ylim=c(0,15))


plist[[names[3]]] = p3

###########################

d4_mean = as.data.frame(matrix(0,12,4))
colnames(d4_mean) = c("Method","20","50","100")

d4_se = as.data.frame(matrix(0,12,4))
colnames(d4_se) = c("Method","20","50","100")

d4_mean$Method = factor(c("GLASSO","GLASSO-I","GLASSO-vI"
       ,"GELNET","GELNET-I","GELNET-vI","ROPE","ROPE-I"," ROPE-vI","GSOS","GSOS-I","GSOS-vI"),
       levels=c("GLASSO","GLASSO-I","GLASSO-vI"
       ,"GELNET","GELNET-I","GELNET-vI","ROPE","ROPE-I"," ROPE-vI","GSOS","GSOS-I","GSOS-vI"))

d4_mean[1,2:4] = c(mean(d4[1:100,3]),mean(d4[101:200,3]),mean(d4[201:300,3]))
d4_mean[2,2:4] = c(mean(d4[301:400,3]),mean(d4[401:500,3]),mean(d4[501:600,3]))
d4_mean[3,2:4] = c(mean(d4[601:700,3]),mean(d4[701:800,3]),mean(d4[801:900,3]))
d4_mean[4,2:4] = c(mean(d4[901:1000,3]),mean(d4[1001:1100,3]),mean(d4[1101:1200,3]))
d4_mean[5,2:4] = c(mean(d4[1201:1300,3]),mean(d4[1301:1400,3]),mean(d4[1401:1500,3]))
d4_mean[6,2:4] = c(mean(d4[1501:1600,3]),mean(d4[1601:1700,3]),mean(d4[1701:1800,3]))
d4_mean[7,2:4] = c(mean(d4[1801:1900,3]),mean(d4[1901:2000,3]),mean(d4[2001:2100,3]))
d4_mean[8,2:4] = c(mean(d4[2101:2200,3]),mean(d4[2201:2300,3]),mean(d4[2301:2400,3]))
d4_mean[9,2:4] = c(mean(d4[2401:2500,3]),mean(d4[2501:2600,3]),mean(d4[2601:2700,3]))
d4_mean[10,2:4] = c(mean(d4[2701:2800,3]),mean(d4[2801:2900,3]),mean(d4[2901:3000,3]))
d4_mean[11,2:4] = c(mean(d4[3001:3100,3]),mean(d4[3101:3200,3]),mean(d4[3201:3300,3]))
d4_mean[12,2:4] = c(mean(d4[3301:3400,3]),mean(d4[3401:3500,3]),mean(d4[3501:3600,3]))


d4_se$Method = factor(c("GLASSO","GLASSO-I","GLASSO-vI"
       ,"GELNET","GELNET-I","GELNET-vI","ROPE","ROPE-I"," ROPE-vI","GSOS","GSOS-I","GSOS-vI"),
       levels=c("GLASSO","GLASSO-I","GLASSO-vI"
       ,"GELNET","GELNET-I","GELNET-vI","ROPE","ROPE-I"," ROPE-vI","GSOS","GSOS-I","GSOS-vI"))

d4_se[1,2:4] = c(se(d4[1:100,3]),se(d4[101:200,3]),se(d4[201:300,3]))
d4_se[2,2:4] = c(se(d4[301:400,3]),se(d4[401:500,3]),se(d4[501:600,3]))
d4_se[3,2:4] = c(se(d4[601:700,3]),se(d4[701:800,3]),se(d4[801:900,3]))
d4_se[4,2:4] = c(se(d4[901:1000,3]),se(d4[1001:1100,3]),se(d4[1101:1200,3]))
d4_se[5,2:4] = c(se(d4[1201:1300,3]),se(d4[1301:1400,3]),se(d4[1401:1500,3]))
d4_se[6,2:4] = c(se(d4[1501:1600,3]),se(d4[1601:1700,3]),se(d4[1701:1800,3]))
d4_se[7,2:4] = c(se(d4[1801:1900,3]),se(d4[1901:2000,3]),se(d4[2001:2100,3]))
d4_se[8,2:4] = c(se(d4[2101:2200,3]),se(d4[2201:2300,3]),se(d4[2301:2400,3]))
d4_se[9,2:4] = c(se(d4[2401:2500,3]),se(d4[2501:2600,3]),se(d4[2601:2700,3]))
d4_se[10,2:4] = c(se(d4[2701:2800,3]),se(d4[2801:2900,3]),se(d4[2901:3000,3]))
d4_se[11,2:4] = c(se(d4[3001:3100,3]),se(d4[3101:3200,3]),se(d4[3201:3300,3]))
d4_se[12,2:4] = c(se(d4[3301:3400,3]),se(d4[3401:3500,3]),se(d4[3501:3600,3]))


d4f_mean = melt(d4_mean, id.vars=c("Method"), variable.name = "Dimension", value.name="LossMean")

d4f_se = melt(d4_se, id.vars=c("Method"), variable.name = "Dimension", value.name="LossMean")

limits4 = aes(ymax = d4f_mean[,"LossMean"] + d4f_se[,"LossMean"], ymin=d4f_mean[,"LossMean"] - d4f_se[,"LossMean"])

p4 =  ggplot(d4f_mean, aes(Dimension, LossMean, fill = Method)) + 
  geom_bar(position="dodge", stat="identity") + 
  geom_pointrange(limits4, position=position_dodge(.9),fatten=0.1) +
  #scale_fill_brewer(palette = "Set3") + 
  scale_fill_manual(values = c("#9DAEE1FF", "#879DE1FF","#718DE1FF",
     "#6584E1FF","#4F73E1FF", "#436BE1FF", "#2D5AE1FF","#1649E1FF","#0038E1FF","#FF7F50", "#A52A2A", "#800000")) +
  theme_gray() +
  labs(y = "Loss mean", x = "Dimension", fill = NULL) + 
  theme(legend.position='none') +
  ggtitle("Model 4") +
  coord_cartesian(ylim=c(0,15))

plist[[names[4]]] = p4

###########################

d5_mean = as.data.frame(matrix(0,12,4))
colnames(d5_mean) = c("Method","20","50","100")

d5_se = as.data.frame(matrix(0,12,4))
colnames(d5_se) = c("Method","20","50","100")

d5_mean$Method = factor(c("GLASSO","GLASSO-I","GLASSO-vI"
       ,"GELNET","GELNET-I","GELNET-vI","ROPE","ROPE-I"," ROPE-vI","NEW","NEW-I","NEW-vI"),
       levels=c("GLASSO","GLASSO-I","GLASSO-vI"
       ,"GELNET","GELNET-I","GELNET-vI","ROPE","ROPE-I"," ROPE-vI","NEW","NEW-I","NEW-vI"))

d5_mean[1,2:4] = c(mean(d5[1:100,3]),mean(d5[101:200,3]),mean(d5[201:300,3]))
d5_mean[2,2:4] = c(mean(d5[301:400,3]),mean(d5[401:500,3]),mean(d5[501:600,3]))
d5_mean[3,2:4] = c(mean(d5[601:700,3]),mean(d5[701:800,3]),mean(d5[801:900,3]))
d5_mean[4,2:4] = c(mean(d5[901:1000,3]),mean(d5[1001:1100,3]),mean(d5[1101:1200,3]))
d5_mean[5,2:4] = c(mean(d5[1201:1300,3]),mean(d5[1301:1400,3]),mean(d5[1401:1500,3]))
d5_mean[6,2:4] = c(mean(d5[1501:1600,3]),mean(d5[1601:1700,3]),mean(d5[1701:1800,3]))
d5_mean[7,2:4] = c(mean(d5[1801:1900,3]),mean(d5[1901:2000,3]),mean(d5[2001:2100,3]))
d5_mean[8,2:4] = c(mean(d5[2101:2200,3]),mean(d5[2201:2300,3]),mean(d5[2301:2400,3]))
d5_mean[9,2:4] = c(mean(d5[2401:2500,3]),mean(d5[2501:2600,3]),mean(d5[2601:2700,3]))
d5_mean[10,2:4] = c(mean(d5[2701:2800,3]),mean(d5[2801:2900,3]),mean(d5[2901:3000,3]))
d5_mean[11,2:4] = c(mean(d5[3001:3100,3]),mean(d5[3101:3200,3]),mean(d5[3201:3300,3]))
d5_mean[12,2:4] = c(mean(d5[3301:3400,3]),mean(d5[3401:3500,3]),mean(d5[3501:3600,3]))


d5_se$Method = factor(c("GLASSO","GLASSO-I","GLASSO-vI"
       ,"GELNET","GELNET-I","GELNET-vI","ROPE","ROPE-I"," ROPE-vI","GSOS","GSOS-I","GSOS-vI"),
       levels=c("GLASSO","GLASSO-I","GLASSO-vI"
       ,"GELNET","GELNET-I","GELNET-vI","ROPE","ROPE-I"," ROPE-vI","GSOS","GSOS-I","GSOS-vI"))

d5_se[1,2:4] = c(se(d5[1:100,3]),se(d5[101:200,3]),se(d5[201:300,3]))
d5_se[2,2:4] = c(se(d5[301:400,3]),se(d5[401:500,3]),se(d5[501:600,3]))
d5_se[3,2:4] = c(se(d5[601:700,3]),se(d5[701:800,3]),se(d5[801:900,3]))
d5_se[4,2:4] = c(se(d5[901:1000,3]),se(d5[1001:1100,3]),se(d5[1101:1200,3]))
d5_se[5,2:4] = c(se(d5[1201:1300,3]),se(d5[1301:1400,3]),se(d5[1401:1500,3]))
d5_se[6,2:4] = c(se(d5[1501:1600,3]),se(d5[1601:1700,3]),se(d5[1701:1800,3]))
d5_se[7,2:4] = c(se(d5[1801:1900,3]),se(d5[1901:2000,3]),se(d5[2001:2100,3]))
d5_se[8,2:4] = c(se(d5[2101:2200,3]),se(d5[2201:2300,3]),se(d5[2301:2400,3]))
d5_se[9,2:4] = c(se(d5[2401:2500,3]),se(d5[2501:2600,3]),se(d5[2601:2700,3]))
d5_se[10,2:4] = c(se(d5[2701:2800,3]),se(d5[2801:2900,3]),se(d5[2901:3000,3]))
d5_se[11,2:4] = c(se(d5[3001:3100,3]),se(d5[3101:3200,3]),se(d5[3201:3300,3]))
d5_se[12,2:4] = c(se(d5[3301:3400,3]),se(d5[3401:3500,3]),se(d5[3501:3600,3]))


d5f_mean = melt(d5_mean, id.vars=c("Method"), variable.name = "Dimension", value.name="LossMean")

d5f_se = melt(d5_se, id.vars=c("Method"), variable.name = "Dimension", value.name="LossMean")

limits5 = aes(ymax = d5f_mean[,"LossMean"] + d5f_se[,"LossMean"], ymin=d5f_mean[,"LossMean"] - d5f_se[,"LossMean"])

p5 =  ggplot(d5f_mean, aes(Dimension, LossMean, fill = Method)) + 
  geom_bar(position="dodge", stat="identity") + 
  geom_pointrange(limits5, position=position_dodge(.9),fatten=0.1) +
  #scale_fill_brewer(palette = "Set3") + 
  scale_fill_manual(values = c("#9DAEE1FF", "#879DE1FF","#718DE1FF",
     "#6584E1FF","#4F73E1FF", "#436BE1FF", "#2D5AE1FF","#1649E1FF","#0038E1FF","#FF7F50", "#A52A2A", "#800000")) +
  theme_gray() +
  labs(y = " ", x = "Dimension", fill = NULL) + 
  theme(legend.position='none') +
  ggtitle("Model 5") +
  coord_cartesian(ylim=c(0,15))

plist[[names[5]]] = p5

###########################

d6_mean = as.data.frame(matrix(0,12,4))
colnames(d6_mean) = c("Method","20","50","100")

d6_se = as.data.frame(matrix(0,12,4))
colnames(d6_se) = c("Method","20","50","100")

d6_mean$Method = factor(c("GLASSO","GLASSO-I","GLASSO-vI"
       ,"GELNET","GELNET-I","GELNET-vI","ROPE","ROPE-I"," ROPE-vI","GSOS","GSOS-I","GSOS-vI"),
       levels=c("GLASSO","GLASSO-I","GLASSO-vI"
       ,"GELNET","GELNET-I","GELNET-vI","ROPE","ROPE-I"," ROPE-vI","GSOS","GSOS-I","GSOS-vI"))

d6_mean[1,2:4] = c(mean(d6[1:100,3]),mean(d6[101:200,3]),mean(d6[201:300,3]))
d6_mean[2,2:4] = c(mean(d6[301:400,3]),mean(d6[401:500,3]),mean(d6[501:600,3]))
d6_mean[3,2:4] = c(mean(d6[601:700,3]),mean(d6[701:800,3]),mean(d6[801:900,3]))
d6_mean[4,2:4] = c(mean(d6[901:1000,3]),mean(d6[1001:1100,3]),mean(d6[1101:1200,3]))
d6_mean[5,2:4] = c(mean(d6[1201:1300,3]),mean(d6[1301:1400,3]),mean(d6[1401:1500,3]))
d6_mean[6,2:4] = c(mean(d6[1501:1600,3]),mean(d6[1601:1700,3]),mean(d6[1701:1800,3]))
d6_mean[7,2:4] = c(mean(d6[1801:1900,3]),mean(d6[1901:2000,3]),mean(d6[2001:2100,3]))
d6_mean[8,2:4] = c(mean(d6[2101:2200,3]),mean(d6[2201:2300,3]),mean(d6[2301:2400,3]))
d6_mean[9,2:4] = c(mean(d6[2401:2500,3]),mean(d6[2501:2600,3]),mean(d6[2601:2700,3]))
d6_mean[10,2:4] = c(mean(d6[2701:2800,3]),mean(d6[2801:2900,3]),mean(d6[2901:3000,3]))
d6_mean[11,2:4] = c(mean(d6[3001:3100,3]),mean(d6[3101:3200,3]),mean(d6[3201:3300,3]))
d6_mean[12,2:4] = c(mean(d6[3301:3400,3]),mean(d6[3401:3500,3]),mean(d6[3501:3600,3]))


d6_se$Method = factor(c("GLASSO","GLASSO-I","GLASSO-vI"
       ,"GELNET","GELNET-I","GELNET-vI","ROPE","ROPE-I"," ROPE-vI","GSOS","GSOS-I","GSOS-vI"),
       levels=c("GLASSO","GLASSO-I","GLASSO-vI"
       ,"GELNET","GELNET-I","GELNET-vI","ROPE","ROPE-I"," ROPE-vI","GSOS","GSOS-I","GSOS-vI"))

d6_se[1,2:4] = c(se(d6[1:100,3]),se(d6[101:200,3]),se(d6[201:300,3]))
d6_se[2,2:4] = c(se(d6[301:400,3]),se(d6[401:500,3]),se(d6[501:600,3]))
d6_se[3,2:4] = c(se(d6[601:700,3]),se(d6[701:800,3]),se(d6[801:900,3]))
d6_se[4,2:4] = c(se(d6[901:1000,3]),se(d6[1001:1100,3]),se(d6[1101:1200,3]))
d6_se[5,2:4] = c(se(d6[1201:1300,3]),se(d6[1301:1400,3]),se(d6[1401:1500,3]))
d6_se[6,2:4] = c(se(d6[1501:1600,3]),se(d6[1601:1700,3]),se(d6[1701:1800,3]))
d6_se[7,2:4] = c(se(d6[1801:1900,3]),se(d6[1901:2000,3]),se(d6[2001:2100,3]))
d6_se[8,2:4] = c(se(d6[2101:2200,3]),se(d6[2201:2300,3]),se(d6[2301:2400,3]))
d6_se[9,2:4] = c(se(d6[2401:2500,3]),se(d6[2501:2600,3]),se(d6[2601:2700,3]))
d6_se[10,2:4] = c(se(d6[2701:2800,3]),se(d6[2801:2900,3]),se(d6[2901:3000,3]))
d6_se[11,2:4] = c(se(d6[3001:3100,3]),se(d6[3101:3200,3]),se(d6[3201:3300,3]))
d6_se[12,2:4] = c(se(d6[3301:3400,3]),se(d6[3401:3500,3]),se(d6[3501:3600,3]))


d6f_mean = melt(d6_mean, id.vars=c("Method"), variable.name = "Dimension", value.name="LossMean")

d6f_se = melt(d6_se, id.vars=c("Method"), variable.name = "Dimension", value.name="LossMean")

limits6 = aes(ymax = d6f_mean[,"LossMean"] + d6f_se[,"LossMean"], ymin=d6f_mean[,"LossMean"] - d6f_se[,"LossMean"])

p6 =  ggplot(d6f_mean, aes(Dimension, LossMean, fill = Method)) + 
  geom_bar(position="dodge", stat="identity") + 
  geom_pointrange(limits6, position=position_dodge(.9),fatten=0.1) +
  #scale_fill_brewer(palette = "Set3") + 
  scale_fill_manual(values = c("#9DAEE1FF", "#879DE1FF","#718DE1FF",
     "#6584E1FF","#4F73E1FF", "#436BE1FF", "#2D5AE1FF","#1649E1FF","#0038E1FF","#FF7F50", "#A52A2A", "#800000")) +
  theme_gray() +
  labs(y = " ", x = "Dimension", fill = NULL) + 
  ggtitle("Model 6") +
  #theme(legend.position='none') +
  theme(legend.position = "bottom", legend.text = element_text(size=7),legend.box = "horizontal") + 
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) + 
  coord_cartesian(ylim=c(0,15))

plist[[names[6]]] = p6

###########################

#do.call("grid.arrange", c(plist, ncol=2))

main=textGrob("KL",gp=gpar(fontsize=14,font=2))

g_legend<-function(a.gplot){
  
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  
  return(legend)
  
}

mylegend<-g_legend(p6)

grid.arrange(arrangeGrob(p1,p2,p3,p4,p5,p6 + theme(legend.position='none'),nrow=2,ncol=3),
             mylegend, nrow=2,heights=c(10,1),top=main)


#########################################

height = 6.42*0.8
width = 9.29*0.8

#height = 5.42
#width = 6.5

pdf("KL.pdf", height=height, width=width)
grid.arrange(arrangeGrob(p1,p2,p3,p4,p5,p6 + theme(legend.position='none'),nrow=2,ncol=3),
             mylegend, nrow=2,heights=c(10,1),top=main)
dev.off()