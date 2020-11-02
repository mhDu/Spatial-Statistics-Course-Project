##################################################
## Project: 空间统计基础除函数拟合代码
## Script purpose: 《空间统计基础》疫情不均衡研究，实现数据分析与绘图
## Date: 2020-07
## Author: Minghao Du
## email: 0108170318@csu.edu.cn / 695948191@qq.com
## Code：R 4.0
##################################################

rm(list = ls())

path = '/Users/Minghao/课程/2020_空间统计基础/数据'
setwd(path)

#####Packages#####
library(bibliometrix)
library(ggplot2)
library(readxl)
library(cowplot)
library(reshape2)
library(matrixStats)
library(scales)
library(car)


#####Load the data#####
# 从PubMed上下载的数据
# 问题一
file_cov <- c('1949-2007.txt', '2008-2019.txt', '2020.txt')
file_sars <- 'SARS.txt'
file_mers <- 'MERS.txt'
file_cov19 <- 'covid-19.txt'

## Excel数据
# 问题二
Sars_in <- read_excel('data.xlsx')[2:6, 1:2]
Mers_in <- read_excel('data.xlsx')[2:33, 4:5]
H1N1_in <- read_excel('data.xlsx')[2:24, 7:8]
Ebola_in <- read_excel('data.xlsx')[2:15, 10:11]
Sars_de <- read_excel('data.xlsx')[2:17, 13:14]
Mers_de <- read_excel('data.xlsx')[2:5, 16:17]
H1N1_de <- read_excel('data.xlsx')[2:10, 19:20]
Ebola_de <- read_excel('data.xlsx')[2:6, 22:23]
Cov19_in <- read_excel('data.xlsx')[2:6, 25:26]

# 问题三
cov_type1 <- read_excel('data.xlsx', sheet=2)[1:67, 1:3]
cov_type2 <- read_excel('data.xlsx', sheet=2)[1:67, 4:6]
cov_type3 <- read_excel('data.xlsx', sheet=2)[1:67, 7:9]
cov_type4 <- read_excel('data.xlsx', sheet=2)[1:67, 10:12]
cov_type5 <- read_excel('data.xlsx', sheet=2)[1:67, 13:15]
cov_type6 <- read_excel('data.xlsx', sheet=2)[1:67, 16:18]
cov_type7 <- read_excel('data.xlsx', sheet=2)[1:67, 19:21]
cov_taxa <- read_excel('data.xlsx', sheet=3)

# 问题四
cov19_data <- read_excel('data.xlsx', sheet=4, col_names=TRUE)[2:125, 1:4]
sars_data <- read_excel('data.xlsx', sheet=4, col_names=TRUE)[2:31, 6:9]
mers_data <- read_excel('data.xlsx', sheet=4, col_names=TRUE)[2:28, 11:14]

############################问题一############################
#转换格式，方便使用bibliometrix包进行分析
M_cov <- convert2df(file_cov, dbsource = "pubmed", format = "pubmed")
M_sars <- convert2df(file_sars, dbsource = "pubmed", format = "pubmed")
M_mers <- convert2df(file_mers, dbsource = "pubmed", format = "pubmed")
#M_cov19 <- convert2df(file_cov19, dbsource = "pubmed", format = "pubmed")

# 简单的分析
results_cov <- biblioAnalysis(M_cov, sep = ";")
results_sars <- biblioAnalysis(M_sars, sep = ";")
results_mers <- biblioAnalysis(M_mers, sep = ";")
#results_cov19 <- biblioAnalysis(M_cov19, sep = ";")

# 把year从上面的分析结果中提取出来
years_cov<- as.data.frame(table(results_cov$Years))
years_sars <- as.data.frame(table(results_sars$Years))
years_mers <- as.data.frame(table(results_mers$Years))

# 把几个数据汇总到一起，方便下面计算剩余数据
years <- merge(years_cov, years_mers, by='Var1', all=T)
years <- merge(years, years_sars, by='Var1', all=T)
years[is.na(years)] <- 0

# 给类别增加上标签，方便叠置直方图
years_sars$veg <- 'sars'
years_mers$veg <- 'mers'
years_cov$veg <- 'all'

# 计算剩余数据
years_re <- years$Freq.x-years$Freq.y-years$Freq
# 把全部数据替换为剩余数据
years_cov$Freq <- years_re
# 把三个数据放在一起，方便下面绘图
years_row <- rbind(years_cov, years_mers, years_sars)
# 绘制直方图
ggplot(years_row, aes(x=Var1, y=Freq, fill=veg)) +
  geom_bar(stat='identity') +
  guides(fill=guide_legend(title=NULL)) +
  scale_fill_discrete(labels=c('residual', 'mers', 'sars')) +
  labs(x='Year', y='Number') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


############################问题二############################
# 把从Excel导入的数据转化为数值型
Sars_in <- as.data.frame(lapply(Sars_in,as.numeric))
Mers_in <- as.data.frame(lapply(Mers_in,as.numeric))
H1N1_in <- as.data.frame(lapply(H1N1_in,as.numeric))
Ebola_in <- as.data.frame(lapply(Ebola_in,as.numeric))
Cov19_in <- as.data.frame(lapply(Cov19_in,as.numeric))
Sars_de <- as.data.frame(lapply(Sars_de,as.numeric))
Mers_de <- as.data.frame(lapply(Mers_de,as.numeric))
H1N1_de <- as.data.frame(lapply(H1N1_de,as.numeric))
Ebola_de <- as.data.frame(lapply(Ebola_de,as.numeric))

#####增长模型#####

##SARS
# 把Sars第一行的数据去除
Sars_in_fit <- Sars_in[-1, ]
names(Sars_in_fit) <- c('t', 'n')
Sars_in_fit$t <- c(1:4)
# 估算Logistic的参数
coef(lm(logit(Sars_in_fit$n/460)~Sars_in_fit$t, data=Sars_in_fit))
# Fit the model
Sars_in_model<-nls(Sars_in_fit$n~phi1/(1+exp(-(phi2+phi3*Sars_in_fit$t))),
                   start=list(phi1=460,phi2=-7.63,phi3=1.97),data=Sars_in_fit,trace=TRUE)
# Model fitting summary
summary(Sars_in_model)
#set parameters
Sars_in_phi1 <- coef(Sars_in_model)[1]
Sars_in_phi2 <- coef(Sars_in_model)[2]
Sars_in_phi3 <- coef(Sars_in_model)[3]
Sars_in_t <- seq(0,4,0.01) #construct a range of x values bounded by the data
Sars_in_n <- Sars_in_phi1/(1+exp(-(Sars_in_phi2+Sars_in_phi3*Sars_in_t))) #predicted mass
Sars_in_predict<-data.frame(Sars_in_t,Sars_in_n) #create the prediction data frame#And add a nice plot 
Sars_in_plot <- ggplot(data=Sars_in,aes(x=1:5,y=...2)) +
  geom_point(color='blue',size=5) +
  geom_line(size=1) +
  theme_bw() +
  labs(x='Month',y='Number')+
  scale_x_continuous(labels=c('2003/01', '2003/02', '2003/03', '2003/04', '2003/05')) +
  theme(axis.text=element_text(size=18),axis.title=element_text(size=24,),
        axis.text.x = element_text(angle = 45, hjust=1)) +
  geom_line(data=Sars_in_predict,aes(x=Sars_in_t+1,y=Sars_in_n), size=1, color='red')

##MERS
# 把Mers第一行的数据去除
Mers_in_fit <- Mers_in[-1, ]
Mers_in_fit <- Mers_in_fit[-31, ]
names(Mers_in_fit) <- c('t', 'n')
Mers_in_fit$t <- c(1:30)
# 估算Logistic的参数
coef(lm(logit(Mers_in_fit$n/34)~Mers_in_fit$t, data=Mers_in_fit))
# Fit the model
Mers_in_model<-nls(Mers_in_fit$n~phi1/(1+exp(-(phi2+phi3*Mers_in_fit$t))),
                   start=list(phi1=34.95,phi2=-2.29,phi3=0.33),data=Mers_in_fit,trace=TRUE)
# Model fitting summary
summary(Mers_in_model)
#set parameters
Mers_in_phi1 <- coef(Mers_in_model)[1]
Mers_in_phi2 <- coef(Mers_in_model)[2]
Mers_in_phi3 <- coef(Mers_in_model)[3]
Mers_in_t <- seq(0,31,0.01) #construct a range of x values bounded by the data
Mers_in_n <- Mers_in_phi1/(1+exp(-(Mers_in_phi2+Mers_in_phi3*Mers_in_t))) #predicted mass
Mers_in_predict<-data.frame(Mers_in_t,Mers_in_n) #create the prediction data frame#And add a nice plot
Mers_in_plot <- ggplot(data=Mers_in,aes(x=1:32,y=...5)) +
  geom_point(color='blue',size=5) +
  geom_line(size=1) +
  theme_bw() +
  labs(x='Month',y='Number')+
  scale_x_continuous(breaks=seq(1,32,1),
    labels=c('2013/01', '2013/02', '2013/03', '2013/04', '2013/05', '2013/06', '2013/07', 
             '2013/08', '2013/09', '2013/10', '2013/11', '2013/12', '2014/01', '2014/02', 
             '2014/03', '2014/04', '2014/05', '2014/06', '2014/07', '2014/08', '2014/09', 
             '2014/10', '2014/11', '2014/12', '2015/01', '2015/02', '2015/03', '2015/04', 
             '2015/05', '2015/06', '2015/07', '2015/08')) +
  theme(axis.text=element_text(size=18),axis.title=element_text(size=24), 
        axis.text.x = element_text(angle = 45, hjust=1)) +
  geom_line(data=Mers_in_predict,aes(x=Mers_in_t+1,y=Mers_in_n), size=1, color='red')

##H1N1
# 把H1N1最后一行的数据去除
H1N1_in_fit <- H1N1_in[-23,]
names(H1N1_in_fit) <- c('t', 'n')
H1N1_in_fit$t <- c(1:22)
# 估算Logistic的参数
coef(lm(logit(H1N1_in_fit$n/300)~H1N1_in_fit$t, data=H1N1_in_fit))
# Fit the model
H1N1_in_model<-nls(H1N1_in_fit$n~phi1/(1+exp(-(phi2+phi3*H1N1_in_fit$t))),
                   start=list(phi1=282,phi2=-2.2,phi3=0.65),data=H1N1_in_fit,trace=TRUE)
# Model fitting summary
summary(H1N1_in_model)
#set parameters
H1N1_in_phi1 <- coef(H1N1_in_model)[1]
H1N1_in_phi2 <- coef(H1N1_in_model)[2]
H1N1_in_phi3 <- coef(H1N1_in_model)[3]
H1N1_in_t <- seq(1,23,0.01) #construct a range of x values bounded by the data
H1N1_in_n <- H1N1_in_phi1/(1+exp(-(H1N1_in_phi2+H1N1_in_phi3*H1N1_in_t))) #predicted mass
H1N1_in_predict<-data.frame(H1N1_in_t,H1N1_in_n) #create the prediction data frame#And add a nice plot
H1N1_in_plot <- ggplot(data=H1N1_in,aes(x=1:23,y=...8)) +
  geom_point(color='blue',size=5) +
  geom_line(size=1) +
  theme_bw() +
  labs(x='Month',y='Number')+
  scale_x_continuous(breaks=seq(1,23,1),
                     labels=c('2009/03', '2009/04', '2009/05', '2009/06', '2009/07', '2009/08', '2009/09', 
                              '2009/10', '2009/11', '2009/12', '2010/01', '2010/02', '2010/03', '2010/04', 
                              '2010/05', '2010/06', '2010/07', '2010/08', '2010/09', '2010/10', '2010/11', 
                              '2010/12', '2011/01')) +
  theme(axis.text=element_text(size=18),axis.title=element_text(size=24),
        axis.text.x = element_text(angle = 45, hjust=1)) +
  geom_line(data=H1N1_in_predict,aes(x=H1N1_in_t,y=H1N1_in_n), size=1, color='red')   

##Ebola
# 选取Ebola的第4到13行
Ebola_in_fit <- Ebola_in[4:13,]
names(Ebola_in_fit) <- c('t', 'n')
Ebola_in_fit$t <- c(1:10)
# 估算Logistic的参数
coef(lm(logit(Ebola_in_fit$n/250)~Ebola_in_fit$t, data=Ebola_in_fit))
# Fit the model
Ebola_in_model<-nls(Ebola_in_fit$n~phi1/(1+exp(-(phi2+phi3*Ebola_in_fit$t))),
                    start=list(phi1=271.97,phi2=-10.46,phi3=1.54),data=Ebola_in_fit,trace=TRUE)
# Model fitting summary
summary(Ebola_in_model)
#set parameters
Ebola_in_phi1 <- coef(Ebola_in_model)[1]
Ebola_in_phi2 <- coef(Ebola_in_model)[2]
Ebola_in_phi3 <- coef(Ebola_in_model)[3]
Ebola_in_t <- seq(-2,11,0.01) #construct a range of x values bounded by the data
Ebola_in_n <- Ebola_in_phi1/(1+exp(-(Ebola_in_phi2+Ebola_in_phi3*Ebola_in_t))) #predicted mass
Ebola_in_predict<-data.frame(Ebola_in_t,Ebola_in_n) #create the prediction data frame#And add a nice plot
Ebola_in_plot <- ggplot(data=Ebola_in,aes(x=1:14,y=...11)) +
  geom_point(color='blue',size=5) +
  geom_line(size=1) +
  theme_bw() +
  labs(x='Month',y='Number')+
  scale_x_continuous(breaks=seq(1,14,1),
                     labels=c('2013/12', '2014/01', '2014/02', '2014/03', '2014/04', '2014/05', '2014/06', 
                              '2014/07', '2014/08', '2014/09', '2014/10', '2014/11', '2014/12', '2015/01')) +
  theme(axis.text=element_text(size=18),axis.title=element_text(size=24),
        axis.text.x = element_text(angle = 45, hjust=1)) +
  geom_line(data=Ebola_in_predict,aes(x=Ebola_in_t+3,y=Ebola_in_n), size=1, color='red')

##COVID-19
Cov19_in_fit <- Cov19_in
names(Cov19_in_fit) <- c('t', 'n')
Cov19_in_fit$t <- c(1:5)
# 估算Logistic的参数
coef(lm(logit(Cov19_in_fit$n/15000)~Cov19_in_fit$t, data=Cov19_in_fit))
# Fit the model
Cov19_in_model<-nls(Cov19_in_fit$n~phi1/(1+exp(-(phi2+phi3*Cov19_in_fit$t))),
                    start=list(phi1=10970,phi2=-7.1,phi3=1.9),data=Cov19_in_fit,trace=TRUE)
# Model fitting summary
summary(Cov19_in_model)
#set parameters
Cov19_in_phi1 <- coef(Cov19_in_model)[1]
Cov19_in_phi2 <- coef(Cov19_in_model)[2]
Cov19_in_phi3 <- coef(Cov19_in_model)[3]
Cov19_in_t <- seq(1,5,0.01) #construct a range of x values bounded by the data
Cov19_in_n <- Cov19_in_phi1/(1+exp(-(Cov19_in_phi2+Cov19_in_phi3*Cov19_in_t))) #predicted mass
Cov19_in_predict<-data.frame(Cov19_in_t,Cov19_in_n) #create the prediction data frame#And add a nice plot
Cov19_in_plot <- ggplot(data=Cov19_in,aes(x=1:5,y=...26)) +
  geom_point(color='blue',size=5) +
  geom_line(size=1) +
  theme_bw() +
  labs(x='Month',y='Number')+
  scale_x_continuous(labels=c('2020/01', '2020/02', '2020/03', '2020/04', '2020/05')) +
  theme(axis.text=element_text(size=18),axis.title=element_text(size=24),
        axis.text.x = element_text(angle = 45, hjust=1)) +
  geom_line(data=Cov19_in_predict,aes(x=Cov19_in_t,y=Cov19_in_n), size=1, color='red')

#####衰退模型#####

##SARS
# Sars选取前9年数据
names(Sars_de) <- c('t', 'n')
Sars_de_fit <- Sars_de[1:9,]
Sars_de_fit$t <- c(1:9)
# Fit the model
Sars_de_model<-nls(Sars_de_fit$n~phi1*exp(-phi2*Sars_de_fit$t),
                   start=list(phi1=2139, phi2=0.34),data=Sars_de_fit,trace=TRUE)
# Model fitting summary
summary(Sars_de_model)
#set parameters
Sars_de_phi1 <- coef(Sars_de_model)[1]
Sars_de_phi2 <- coef(Sars_de_model)[2]
Sars_de_t <- seq(1,16,0.01) #construct a range of x values bounded by the data
Sars_de_n <- Sars_de_phi1*exp(-Sars_de_phi2*Sars_de_t) #predicted mass
Sars_de_predict <- data.frame(Sars_de_t,Sars_de_n) #create the prediction data frame#And add a nice plot 
Sars_de_plot <- ggplot(data=Sars_de,aes(x=1:16,y=n)) +
  geom_point(color='blue',size=5) +
  geom_line(size=1) +
  theme_bw() +
  labs(x='Year',y='Number')+
  scale_x_continuous(breaks=seq(1,16,1),
                     labels=c('2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010', 
                              '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018')) +
  theme(axis.text=element_text(size=18),axis.title=element_text(size=24),
        axis.text.x = element_text(angle = 45, hjust=1)) +
  geom_line(data=Sars_de_predict,aes(x=Sars_de_t,y=Sars_de_n), size=1, color='red')

##MERS
names(Mers_de) <- c('t', 'n')
Mers_de_fit <- Mers_de
Mers_de_fit$t <- c(1:4)
# Fit the model
Mers_de_model<-nls(Mers_de_fit$n~phi1*exp(-phi2*Mers_de_fit$t),
                   start=list(phi1=518, phi2=0.18),data=Mers_de_fit,trace=TRUE)
# Model fitting summary
summary(Mers_de_model)
#set parameters
Mers_de_phi1 <- coef(Mers_de_model)[1]
Mers_de_phi2 <- coef(Mers_de_model)[2]
Mers_de_t <- seq(1,4,0.01) #construct a range of x values bounded by the data
Mers_de_n <- Mers_de_phi1*exp(-Mers_de_phi2*Mers_de_t) #predicted mass
Mers_de_predict <- data.frame(Mers_de_t,Mers_de_n) #create the prediction data frame#And add a nice plot 
Mers_de_plot <- ggplot(data=Mers_de,aes(x=1:4,y=n)) +
  geom_point(color='blue',size=5) +
  geom_line(size=1) +
  theme_bw() +
  labs(x='Year',y='Number')+
  ylim(0,450) +
  scale_x_continuous(labels=c('2015', '2016', '2017', '2018')) +
  theme(axis.text=element_text(size=18),axis.title=element_text(size=24),
        axis.text.x = element_text(angle = 45, hjust=1)) +
  geom_line(data=Mers_de_predict,aes(x=Mers_de_t,y=Mers_de_n), size=1, color='red')

##H1N1
names(H1N1_de) <- c('t', 'n')
H1N1_de_fit <- H1N1_de
H1N1_de_fit$t <- c(1:9)
# Fit the model
H1N1_de_model<-nls(H1N1_de_fit$n~phi1*exp(-phi2*H1N1_de_fit$t),
                   start=list(phi1=3361, phi2=0.23),data=H1N1_de_fit,trace=TRUE)
# Model fitting summary
summary(H1N1_de_model)
#set parameters
H1N1_de_phi1 <- coef(H1N1_de_model)[1]
H1N1_de_phi2 <- coef(H1N1_de_model)[2]
H1N1_de_t <- seq(1,9,0.01) #construct a range of x values bounded by the data
H1N1_de_n <- H1N1_de_phi1*exp(-H1N1_de_phi2*H1N1_de_t) #predicted mass
H1N1_de_predict <- data.frame(H1N1_de_t,H1N1_de_n) #create the prediction data frame#And add a nice plot 
H1N1_de_plot <- ggplot(data=H1N1_de,aes(x=1:9,y=n)) +
  geom_point(color='blue',size=5) +
  geom_line(size=1) +
  theme_bw() +
  labs(x='Year',y='Number')+
  ylim(0, 3000) +
  scale_x_continuous(breaks=seq(1,9,1),
                     labels=c('2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018', '2019')) +
  theme(axis.text=element_text(size=18),axis.title=element_text(size=24),
        axis.text.x = element_text(angle = 45, hjust=1)) +
  geom_line(data=H1N1_de_predict,aes(x=H1N1_de_t,y=H1N1_de_n), size=1, color='red')

##EBOLA
names(Ebola_de) <- c('t', 'n')
Ebola_de_fit <- Ebola_de
Ebola_de_fit$t <- c(1:5)
# Fit the model
Ebola_de_model<-nls(Ebola_de_fit$n~phi1*exp(-phi2*Ebola_de_fit$t),
                    start=list(phi1=2328, phi2=0.31),data=Ebola_de_fit,trace=TRUE)
# Model fitting summary
summary(Ebola_de_model)
#set parameters
Ebola_de_phi1 <- coef(Ebola_de_model)[1]
Ebola_de_phi2 <- coef(Ebola_de_model)[2]
Ebola_de_t <- seq(1,5,0.01) #construct a range of x values bounded by the data
Ebola_de_n <- Ebola_de_phi1*exp(-Ebola_de_phi2*Ebola_de_t) #predicted mass
Ebola_de_predict <- data.frame(Ebola_de_t,Ebola_de_n) #create the prediction data frame#And add a nice plot 
Ebola_de_plot <- ggplot(data=Ebola_de,aes(x=1:5,y=n)) +
  geom_point(color='blue',size=5) +
  geom_line(size=1) +
  theme_bw() +
  labs(x='Year',y='Number')+
  ylim(0, 1800) +
  scale_x_continuous(labels=c('2015', '2016', '2017', '2018', '2019')) +
  theme(axis.text=element_text(size=18),axis.title=element_text(size=24),
        axis.text.x = element_text(angle = 45, hjust=1)) +
  geom_line(data=Ebola_de_predict,aes(x=Ebola_de_t,y=Ebola_de_n), size=1, color='red')

#####把所有的图集中#####
plot_grid(Sars_in_plot, Mers_in_plot, H1N1_in_plot, Ebola_in_plot, Cov19_in_plot,
          Sars_de_plot, Mers_de_plot, H1N1_de_plot, Ebola_de_plot,
          nrow=3, ncol=3, labels='auto', label_size=50)


############################问题三############################
# 修改从Excel导入的数据列名
names(cov_type1) <- c('year', 'number', 'type')
names(cov_type2) <- c('year', 'number', 'type')
names(cov_type3) <- c('year', 'number', 'type')
names(cov_type4) <- c('year', 'number', 'type')
names(cov_type5) <- c('year', 'number', 'type')
names(cov_type6) <- c('year', 'number', 'type')
names(cov_type7) <- c('year', 'number', 'type')
# 把7个数据放在一起，方便下面绘图
cov_type <- rbind(cov_type1, cov_type2, cov_type3, cov_type4, cov_type5, cov_type6, cov_type7)
# 绘制直方图
ggplot(cov_type, aes(x=year, y=number, fill=type)) +
  geom_bar(stat='identity') +
  guides(fill=guide_legend(title=NULL)) +
  scale_x_continuous(breaks=seq(1953,2020,1)) +
  labs(x='Year', y='Number') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 计算变异系数
cov_taxa[cov_taxa == 0] <- NA
cov_taxa$mean <- rowMeans(cov_taxa[-1], na.rm=TRUE)
cov_taxa$variance <- apply(cov_taxa[2:8], 1, sd, na.rm=TRUE)
cov_taxa[is.na(cov_taxa)] <- 0
cov_taxa[2, 11] <- 0
cov_taxa$cv <- cov_taxa$variance/cov_taxa$mean
# 变异系数曲线
ggplot(cov_taxa, aes(x=year, y=cv)) +
  geom_line(color='red', size=1) +
  guides(fill=guide_legend(title=NULL)) +
  scale_x_continuous(breaks=seq(1953,2020,1)) +
  labs(x='Year', y='Number') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))   


############################问题四############################
#####所有国家冠状病毒论文概况#####
# 把从PubMed上面下载的数据导入country
Mccn_cov <- metaTagExtraction(M_cov, Field = "AU_CO", sep = ";")
# 把country提取出来
country_cov <- Mccn_cov$AU_CO
# 把country字段分割
country_cov <- strsplit(country_cov, ';')
# 把分割后的list拆分
country_cov <- unlist(country_cov, use.name=FALSE)
# 变为数据框
country_cov <- as.data.frame(table(country_cov))
# 把数据框排序
country_cov <- country_cov[order(-country_cov$Freq),]
# 修改行名
rownames(country_cov) <- c(1:137)
# 选择频数大于1000的行
country_cov_sub <- subset(country_cov, Freq>1000)
# 绘制频数直方图
ggplot(country_cov_sub, aes(x=reorder(country_cov, -Freq), y=Freq)) + # reoder是为了排序
  geom_bar(stat='identity') +
  labs(x='Country', y='Number') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# Calculate the quantile
country_cov_quantile <- quantile(country_cov$Freq, seq(0, 1, 0.05))
# Calculate the mean
country_cov_mean <- mean(country_cov$Freq)
# Calculate the standard deviation
country_cov_sd <- sd(country_cov$Freq)
# Calculate the Coefficient of Variation
country_cov_cv <-  country_cov_sd/country_cov_mean


#####计算e#####
# 修改列名
names(cov19_data) <- c('country', 'paper', 'diagnose', 'gdp')
names(sars_data) <- c('country', 'paper', 'diagnose', 'gdp')
names(mers_data) <- c('country', 'paper', 'diagnose', 'gdp')
# 把从Excel导入的数据转化为数值型
cov19_data$paper <- as.numeric(cov19_data$paper)
cov19_data$diagnose <- as.numeric(cov19_data$diagnose)
cov19_data$gdp <- as.numeric(cov19_data$gdp)
mers_data$paper <- as.numeric(mers_data$paper)
mers_data$diagnose <- as.numeric(mers_data$diagnose)
mers_data$gdp <- as.numeric(mers_data$gdp)
sars_data$paper <- as.numeric(sars_data$paper)
sars_data$diagnose <- as.numeric(sars_data$diagnose)
sars_data$gdp <- as.numeric(sars_data$gdp)

## Mers
# Covert quantitative data to ranking data
mers_paper_rank <- rank(mers_data$paper, ties.method='average')
mers_dia_rank <- rank(mers_data$diagnose, ties.method='average')
mers_gdp_rank <- rank(mers_data$gdp, ties.method='average')
# Add the ranking data to data frame
mers_data$dia_rank <- mers_dia_rank
mers_data$paper_rank <- mers_paper_rank
mers_data$gdp_rank <- mers_gdp_rank
## 计算spearman系数
# mers 诊断数与论文数
cor.test(mers_data$diagnose, mers_data$paper, method='spearman')
cor.test(mers_data$diagnose, mers_data$paper, method='kendall')
# mers gpd与论文数
cor.test(mers_data$gdp, mers_data$paper, method='spearman')
cor.test(mers_data$gdp, mers_data$paper, method='kendall')
# Calculate the e
e_mers <- rep(1, 27)
for (i in 1:27) {
  if(mers_dia_rank[i] > mers_paper_rank[i]){
    e_mers[i] <- (mers_paper_rank[i]-mers_dia_rank[i])/sqrt(mers_dia_rank[i])
  }else if(mers_dia_rank[i] == mers_paper_rank[i]){
    e_mers[i] <- 0
  }else{
    e_mers[i] <- (mers_paper_rank[i]-mers_dia_rank[i])/sqrt(mers_paper_rank[i])
  }
}
mers_data$e <- e_mers
sum(e_mers)
# gdp 的e
e_mers_gdp <- rep(1, 27)
for (i in 1:27) {
  if(mers_gdp_rank[i] > mers_paper_rank[i]){
    e_mers_gdp[i] <- (mers_paper_rank[i]-mers_gdp_rank[i])/sqrt(mers_gdp_rank[i])
  }else if(mers_gdp_rank[i] == mers_paper_rank[i]){
    e_mers_gdp[i] <- 0
  }else{
    e_mers_gdp[i] <- (mers_paper_rank[i]-mers_gdp_rank[i])/sqrt(mers_paper_rank[i])
  }
}
mers_data$e_gdp <- e_mers_gdp
sum(e_mers_gdp)
# add e=3 curve
mers_d_3 <- seq(1,30,0.1)
mers_p_3 <- (mers_d_3 - 3*sqrt(mers_d_3))
mers_3 <- data.frame(mers_d_3, mers_p_3)
mers_p_3_a <- seq(1,30,0.1)
mers_d_3_a <- (mers_p_3_a - 3*sqrt(mers_p_3_a))
mers_3_a <- data.frame(mers_d_3_a, mers_p_3_a)
# plot scatter plot
sp_mers <- ggplot(data=mers_data, aes(x=dia_rank, y=paper_rank)) +
  geom_point(col = ifelse(mers_data$e>0, 'red', ifelse(mers_data$e==0, 'black', 'blue')), size=3) +
  geom_abline(size=2) +
  theme_classic() +
  scale_x_continuous(limit=c(0,30)) +
  scale_y_continuous(limit=c(0,30)) +
  geom_line(data=mers_3, aes(x=mers_d_3, y=mers_p_3), size=2) +
  geom_line(data=mers_3_a, aes(x=mers_d_3_a, y=mers_p_3_a), size=2) +
  coord_fixed(ratio = 1) +
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=34,face="bold"))
# gdp
sp_mers_gdp <- ggplot(data=mers_data, aes(x=gdp_rank, y=paper_rank)) +
  geom_point(col = ifelse(mers_data$e_gdp>0, 'red', ifelse(mers_data$e_gdp==0, 'black', 'blue')), size=3) +
  geom_abline(size=2) +
  theme_classic() +
  scale_x_continuous(limit=c(0,30)) +
  scale_y_continuous(limit=c(0,30)) +
  geom_line(data=mers_3, aes(x=mers_d_3, y=mers_p_3), size=2) +
  geom_line(data=mers_3_a, aes(x=mers_d_3_a, y=mers_p_3_a), size=2) +
  coord_fixed(ratio = 1) +
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=34,face="bold"))
# Plot PDF of e
e_mers_pdf <- ggplot(mers_data, aes(x=e_gdp)) +
  geom_density(color='red', size=2) +
  geom_density(aes(x=e), size=2) +
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=34,face="bold"))
# Calculate the quantile
country_mers_quantile <- quantile(mers_data$paper, seq(0, 1, 0.05))
# Calculate the mean
country_mers_mean <- mean(mers_data$paper)
# Calculate the standard deviation
country_mers_sd <- sd(mers_data$paper)
# Calculate the Coefficient of Variation
country_mers_cv <-  country_mers_sd/country_mers_mean


# 患病数大于三例的
mers_data_3 <- subset(mers_data, diagnose>3)
mers_paper_rank_3 <- rank(mers_data_3$paper, ties.method='average')
mers_dia_rank_3 <- rank(mers_data_3$diagnose, ties.method='average')
mers_gdp_rank_3 <- rank(mers_data_3$gdp, ties.method='average')
# Add the ranking data to data frame
mers_data_3$dia_rank <- mers_dia_rank_3
mers_data_3$paper_rank <- mers_paper_rank_3
mers_data_3$gdp_rank <- mers_gdp_rank_3
## 计算spearman系数
# mers 诊断数与论文数
cor.test(mers_data_3$diagnose, mers_data_3$paper, method='spearman')
cor.test(mers_data_3$diagnose, mers_data_3$paper, method='kendall')
# mers gdp与论文数
cor.test(mers_data_3$gdp, mers_data_3$paper, method='spearman')
cor.test(mers_data_3$gdp, mers_data_3$paper, method='kendall')
# Calculate the e
e_mers_3 <- rep(1, 9)
for (i in 1:9) {
  if(mers_dia_rank_3[i] > mers_paper_rank_3[i]){
    e_mers_3[i] <- (mers_paper_rank_3[i]-mers_dia_rank_3[i])/sqrt(mers_dia_rank_3[i])
  }else if(mers_dia_rank_3[i] == mers_paper_rank_3[i]){
    e_mers_3[i] <- 0
  }else{
    e_mers_3[i] <- (mers_paper_rank_3[i]-mers_dia_rank_3[i])/sqrt(mers_paper_rank_3[i])
  }
}
mers_data_3$e <- e_mers_3
sum(e_mers_3)
# e_gdp
e_mers_gdp_3 <- rep(1, 9)
for (i in 1:9) {
  if(mers_gdp_rank_3[i] > mers_paper_rank_3[i]){
    e_mers_gdp_3[i] <- (mers_paper_rank_3[i]-mers_gdp_rank_3[i])/sqrt(mers_gdp_rank_3[i])
  }else if(mers_dia_rank_3[i] == mers_paper_rank_3[i]){
    e_mers_gdp_3[i] <- 0
  }else{
    e_mers_gdp_3[i] <- (mers_paper_rank_3[i]-mers_gdp_rank_3[i])/sqrt(mers_paper_rank_3[i])
  }
}
mers_data_3$e_gdp <- e_mers_gdp_3
sum(e_mers_gdp_3)
# add e=3 curve
mers_d_1 <- seq(1,9,0.1)
mers_p_1 <- (mers_d_1 - sqrt(mers_d_1))
mers_1 <- data.frame(mers_d_1, mers_p_1)
mers_p_1_a <- seq(1,9,0.1)
mers_d_1_a <- (mers_p_1_a - sqrt(mers_p_1_a))
mers_1_a <- data.frame(mers_d_1_a, mers_p_1_a)
# plot scatter plot
sp_mers_3 <- ggplot(data=mers_data_3, aes(x=dia_rank, y=paper_rank)) +
  geom_abline(size=2) +
  theme_classic() +
  scale_x_continuous(limit=c(0,9)) +
  scale_y_continuous(limit=c(0,9)) +
  geom_point(col = ifelse(mers_data_3$e>0, 'red', ifelse(mers_data_3$e==0, 'black', 'blue')), size=3) +
  geom_line(data=mers_1, aes(x=mers_d_1, y=mers_p_1), size=2) +
  geom_line(data=mers_1_a, aes(x=mers_d_1_a, y=mers_p_1_a), size=2) +
  coord_fixed(ratio = 1) +
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=34,face="bold"))
# gdp
sp_mers_gdp_3 <- ggplot(data=mers_data_3, aes(x=gdp_rank, y=paper_rank)) +
  geom_abline(size=2) +
  theme_classic() +
  scale_x_continuous(limit=c(0,9)) +
  scale_y_continuous(limit=c(0,9)) +
  geom_point(col = ifelse(mers_data_3$e_gdp>0, 'red', ifelse(mers_data_3$e_gdp==0, 'black', 'blue')), size=3) +
  geom_line(data=mers_1, aes(x=mers_d_1, y=mers_p_1), size=2) +
  geom_line(data=mers_1_a, aes(x=mers_d_1_a, y=mers_p_1_a), size=2) +
  coord_fixed(ratio = 1) +
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=34,face="bold"))
# Plot PDF of e
e_mers_pdf_3 <- ggplot(mers_data_3, aes(x=e_gdp)) +
  geom_density(color='red', size=2) +
  geom_density(aes(x=e), size=2) +
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=34,face="bold"))
# Calculate the quantile
country_mers_quantile_3 <- quantile(mers_data_3$paper, seq(0, 1, 0.05))
# Calculate the mean
country_mers_mean_3 <- mean(mers_data_3$paper)
# Calculate the standard deviation
country_mers_sd_3 <- sd(mers_data_3$paper)
# Calculate the Coefficient of Variation
country_mers_cv_3 <-  country_mers_sd_3/country_mers_mean_3

## Sars
# Covert quantitative data to ranking data
sars_paper_rank <- rank(sars_data$paper, ties.method='average')
sars_dia_rank <- rank(sars_data$diagnose, ties.method='average')
sars_gdp_rank <- rank(sars_data$gdp, ties.method='average')
# Add the ranking data to data frame
sars_data$dia_rank <- sars_dia_rank
sars_data$paper_rank <- sars_paper_rank
sars_data$gdp_rank <- sars_gdp_rank
## 计算spearman系数
# sars 诊断数与论文数
cor.test(sars_data$diagnose, sars_data$paper, method='spearman')
cor.test(sars_data$diagnose, sars_data$paper, method='kendall')
# sars gpd与论文数
cor.test(sars_data$gdp, sars_data$paper, method='spearman')
cor.test(sars_data$gdp, sars_data$paper, method='kendall')
# Calculate the e
e_sars <- rep(1, 30)
for (i in 1:30) {
  if(sars_dia_rank[i] > sars_paper_rank[i]){
    e_sars[i] <- (sars_paper_rank[i]-sars_dia_rank[i])/sqrt(sars_dia_rank[i])
  }else if(sars_dia_rank[i] == sars_paper_rank[i]){
    e_sars[i] <- 0
  }else{
    e_sars[i] <- (sars_paper_rank[i]-sars_dia_rank[i])/sqrt(sars_paper_rank[i])
  }
}
sars_data$e <- e_sars
sum(e_sars)
# gdp 的e
e_sars_gdp <- rep(1, 30)
for (i in 1:30) {
  if(sars_gdp_rank[i] > sars_paper_rank[i]){
    e_sars_gdp[i] <- (sars_paper_rank[i]-sars_gdp_rank[i])/sqrt(sars_gdp_rank[i])
  }else if(sars_gdp_rank[i] == sars_paper_rank[i]){
    e_sars_gdp[i] <- 0
  }else{
    e_sars_gdp[i] <- (sars_paper_rank[i]-sars_gdp_rank[i])/sqrt(sars_paper_rank[i])
  }
}
sars_data$e_gdp <- e_sars_gdp
sum(e_sars_gdp)
# add e=3 curve
sars_d_3 <- seq(1,30,0.1)
sars_p_3 <- (sars_d_3 - 3*sqrt(sars_d_3))
sars_3 <- data.frame(sars_d_3, sars_p_3)
sars_p_3_a <- seq(1,30,0.1)
sars_d_3_a <- (sars_p_3_a - 3*sqrt(sars_p_3_a))
sars_3_a <- data.frame(sars_d_3_a, sars_p_3_a)
# plot scatter plot
sp_sars <- ggplot(data=sars_data, aes(x=dia_rank, y=paper_rank)) +
  geom_point(col = ifelse(sars_data$e>0, 'red', ifelse(sars_data$e==0, 'black', 'blue')), size=3) +
  geom_abline(size=2) +
  theme_classic() +
  scale_x_continuous(limit=c(0,30)) +
  scale_y_continuous(limit=c(0,30)) +
  geom_line(data=sars_3, aes(x=sars_d_3, y=sars_p_3), size=2) +
  geom_line(data=sars_3_a, aes(x=sars_d_3_a, y=sars_p_3_a), size=2) +
  coord_fixed(ratio = 1) +
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=34,face="bold"))
# gdp
sp_sars_gdp <- ggplot(data=sars_data, aes(x=gdp_rank, y=paper_rank)) +
  geom_point(col = ifelse(sars_data$e_gdp>0, 'red', ifelse(sars_data$e_gdp==0, 'black', 'blue')), size=3) +
  geom_abline(size=2) +
  theme_classic() +
  scale_x_continuous(limit=c(0,30)) +
  scale_y_continuous(limit=c(0,30)) +
  geom_line(data=sars_3, aes(x=sars_d_3, y=sars_p_3), size=2) +
  geom_line(data=sars_3_a, aes(x=sars_d_3_a, y=sars_p_3_a), size=2) +
  coord_fixed(ratio = 1) +
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=34,face="bold"))
# Plot PDF of e
e_sars_pdf <- ggplot(sars_data, aes(x=e_gdp)) +
  geom_density(color='red', size=2) +
  geom_density(aes(x=e), size=2) +
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=34,face="bold"))
# Calculate the quantile
country_sars_quantile <- quantile(sars_data$paper, seq(0, 1, 0.05))
# Calculate the mean
country_sars_mean <- mean(sars_data$paper)
# Calculate the standard deviation
country_sars_sd <- sd(sars_data$paper)
# Calculate the Coefficient of Variation
country_sars_cv <-  country_sars_sd/country_sars_mean


# 患病数大于三例的
sars_data_3 <- subset(sars_data, diagnose>3)
sars_paper_rank_3 <- rank(sars_data_3$paper, ties.method='average')
sars_dia_rank_3 <- rank(sars_data_3$diagnose, ties.method='average')
sars_gdp_rank_3 <- rank(sars_data_3$gdp, ties.method='average')
# Add the ranking data to data frame
sars_data_3$dia_rank <- sars_dia_rank_3
sars_data_3$paper_rank <- sars_paper_rank_3
sars_data_3$gdp_rank <- sars_gdp_rank_3
## 计算spearman系数
# sars 诊断数与论文数
cor.test(sars_data_3$diagnose, sars_data_3$paper, method='spearman')
cor.test(sars_data_3$diagnose, sars_data_3$paper, method='kendall')
# sars gdp与论文数
cor.test(sars_data_3$gdp, sars_data_3$paper, method='spearman')
cor.test(sars_data_3$gdp, sars_data_3$paper, method='kendall')
# Calculate the e
e_sars_3 <- rep(1, 15)
for (i in 1:15) {
  if(sars_dia_rank_3[i] > sars_paper_rank_3[i]){
    e_sars_3[i] <- (sars_paper_rank_3[i]-sars_dia_rank_3[i])/sqrt(sars_dia_rank_3[i])
  }else if(sars_dia_rank_3[i] == sars_paper_rank_3[i]){
    e_sars_3[i] <- 0
  }else{
    e_sars_3[i] <- (sars_paper_rank_3[i]-sars_dia_rank_3[i])/sqrt(sars_paper_rank_3[i])
  }
}
sars_data_3$e <- e_sars_3
sum(e_sars_3)
# e_gdp
e_sars_gdp_3 <- rep(1, 15)
for (i in 1:15) {
  if(sars_gdp_rank_3[i] > sars_paper_rank_3[i]){
    e_sars_gdp_3[i] <- (sars_paper_rank_3[i]-sars_gdp_rank_3[i])/sqrt(sars_gdp_rank_3[i])
  }else if(sars_dia_rank_3[i] == sars_paper_rank_3[i]){
    e_sars_gdp_3[i] <- 0
  }else{
    e_sars_gdp_3[i] <- (sars_paper_rank_3[i]-sars_gdp_rank_3[i])/sqrt(sars_paper_rank_3[i])
  }
}
sars_data_3$e_gdp <- e_sars_gdp_3
sum(e_sars_gdp_3)
# add e=1 curve
sars_d_1 <- seq(1,15,0.1)
sars_p_1 <- (sars_d_1 - sqrt(sars_d_1))
sars_1 <- data.frame(sars_d_1, sars_p_1)
sars_p_1_a <- seq(1,15,0.1)
sars_d_1_a <- (sars_p_1_a - sqrt(sars_p_1_a))
sars_1_a <- data.frame(sars_d_1_a, sars_p_1_a)
# plot scatter plot
sp_sars_3 <- ggplot(data=sars_data_3, aes(x=dia_rank, y=paper_rank)) +
  geom_abline(size=2) +
  theme_classic() +
  scale_x_continuous(limit=c(0,9)) +
  scale_y_continuous(limit=c(0,9)) +
  geom_point(col = ifelse(sars_data_3$e>0, 'red', ifelse(sars_data_3$e==0, 'black', 'blue')), size=3) +
  geom_line(data=sars_1, aes(x=sars_d_1, y=sars_p_1), size=2) +
  geom_line(data=sars_1_a, aes(x=sars_d_1_a, y=sars_p_1_a), size=2) +
  coord_fixed(ratio = 1) +
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=34,face="bold"))
# gdp
sp_sars_gdp_3 <- ggplot(data=sars_data_3, aes(x=gdp_rank, y=paper_rank)) +
  geom_abline(size=2) +
  theme_classic() +
  scale_x_continuous(limit=c(0,9)) +
  scale_y_continuous(limit=c(0,9)) +
  geom_point(col = ifelse(sars_data_3$e_gdp>0, 'red', ifelse(sars_data_3$e_gdp==0, 'black', 'blue')), size=3) +
  geom_line(data=sars_1, aes(x=sars_d_1, y=sars_p_1), size=2) +
  geom_line(data=sars_1_a, aes(x=sars_d_1_a, y=sars_p_1_a), size=2) +
  coord_fixed(ratio = 1) +
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=34,face="bold"))
# Plot PDF of e
e_sars_pdf_3 <- ggplot(sars_data_3, aes(x=e_gdp)) +
  geom_density(color='red', size=2) +
  geom_density(aes(x=e), size=2) +
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=34,face="bold"))
# Calculate the quantile
country_sars_quantile_3 <- quantile(sars_data_3$paper, seq(0, 1, 0.05))
# Calculate the mean
country_sars_mean_3 <- mean(sars_data_3$paper)
# Calculate the standard deviation
country_sars_sd_3 <- sd(sars_data_3$paper)
# Calculate the Coefficient of Variation
country_sars_cv_3 <-  country_sars_sd_3/country_sars_mean_3

## Covid-19
# Covert quantitative data to ranking data
cov_paper_rank <- rank(cov19_data$paper, ties.method='average')
cov_dia_rank <- rank(cov19_data$diagnose, ties.method='average')
cov_gdp_rank <- rank(cov19_data$gdp, ties.method='average')
# Add the ranking data to data frame
cov19_data$dia_rank <- cov_dia_rank
cov19_data$paper_rank <- cov_paper_rank
cov19_data$gdp_rank <- cov_gdp_rank
## 计算spearman系数
# cov 诊断数与论文数
cor.test(cov19_data$diagnose, cov19_data$paper, method='spearman')
cor.test(cov19_data$diagnose, cov19_data$paper, method='kendall')
# cov gpd与论文数
cor.test(cov19_data$gdp, cov19_data$paper, method='spearman')
cor.test(cov19_data$gdp, cov19_data$paper, method='kendall')
# Calculate the e
e_cov <- rep(1, 124)
for (i in 1:124) {
  if(cov_dia_rank[i] > cov_paper_rank[i]){
    e_cov[i] <- (cov_paper_rank[i]-cov_dia_rank[i])/sqrt(cov_dia_rank[i])
  }else if(cov_dia_rank[i] == cov_paper_rank[i]){
    e_cov[i] <- 0
  }else{
    e_cov[i] <- (cov_paper_rank[i]-cov_dia_rank[i])/sqrt(cov_paper_rank[i])
  }
}
cov19_data$e <- e_cov
sum(e_cov)
# gdp 的e
e_cov_gdp <- rep(1, 124)
for (i in 1:124) {
  if(cov_gdp_rank[i] > cov_paper_rank[i]){
    e_cov_gdp[i] <- (cov_paper_rank[i]-cov_gdp_rank[i])/sqrt(cov_gdp_rank[i])
  }else if(cov_gdp_rank[i] == cov_paper_rank[i]){
    e_cov_gdp[i] <- 0
  }else{
    e_cov_gdp[i] <- (cov_paper_rank[i]-cov_gdp_rank[i])/sqrt(cov_paper_rank[i])
  }
}
cov19_data$e_gdp <- e_cov_gdp
sum(e_cov_gdp)
# add e=3 curve
cov_d_3 <- seq(1,124,0.1)
cov_p_3 <- (cov_d_3 - 3*sqrt(cov_d_3))
cov_3 <- data.frame(cov_d_3, cov_p_3)
cov_p_3_a <- seq(1,124,0.1)
cov_d_3_a <- (cov_p_3_a - 3*sqrt(cov_p_3_a))
cov_3_a <- data.frame(cov_d_3_a, cov_p_3_a)
# plot scatter plot
sp_cov <- ggplot(data=cov19_data, aes(x=dia_rank, y=paper_rank)) +
  geom_point(col = ifelse(cov19_data$e>0, 'red', ifelse(cov19_data$e==0, 'black', 'blue')), size=3) +
  geom_abline(size=2) +
  theme_classic() +
  scale_x_continuous(limit=c(0,124)) +
  scale_y_continuous(limit=c(0,124)) +
  geom_line(data=cov_3, aes(x=cov_d_3, y=cov_p_3), size=2) +
  geom_line(data=cov_3_a, aes(x=cov_d_3_a, y=cov_p_3_a), size=2) +
  coord_fixed(ratio = 1) +
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=34,face="bold"))
# gdp
sp_cov_gdp <- ggplot(data=cov19_data, aes(x=gdp_rank, y=paper_rank)) +
  geom_point(col = ifelse(cov19_data$e_gdp>0, 'red', ifelse(cov19_data$e_gdp==0, 'black', 'blue')), size=3) +
  geom_abline(size=2) +
  theme_classic() +
  scale_x_continuous(limit=c(0,124)) +
  scale_y_continuous(limit=c(0,124)) +
  geom_line(data=cov_3, aes(x=cov_d_3, y=cov_p_3), size=2) +
  geom_line(data=cov_3_a, aes(x=cov_d_3_a, y=cov_p_3_a), size=2) +
  coord_fixed(ratio = 1) +
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=34,face="bold"))
# Plot PDF of e
e_cov_pdf <- ggplot(cov19_data, aes(x=e_gdp)) +
  geom_density(color='red', size=2) +
  geom_density(aes(x=e), size=2) +
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=34,face="bold"))
# Calculate the quantile
country_cov_quantile <- quantile(cov19_data$paper, seq(0, 1, 0.05))
# Calculate the mean
country_cov_mean <- mean(cov19_data$paper)
# Calculate the standard deviation
country_cov_sd <- sd(cov19_data$paper)
# Calculate the Coefficient of Variation
country_cov_cv <-  country_cov_sd/country_cov_mean


# 论文大于20的
cov19_data_3 <- subset(cov19_data, paper>20)
cov_paper_rank_3 <- rank(cov19_data_3$paper, ties.method='average')
cov_dia_rank_3 <- rank(cov19_data_3$diagnose, ties.method='average')
cov_gdp_rank_3 <- rank(cov19_data_3$gdp, ties.method='average')
# Add the ranking data to data frame
cov19_data_3$dia_rank <- cov_dia_rank_3
cov19_data_3$paper_rank <- cov_paper_rank_3
cov19_data_3$gdp_rank <- cov_gdp_rank_3
## 计算spearman系数
# cov 诊断数与论文数
cor.test(cov19_data_3$diagnose, cov19_data_3$paper, method='spearman')
cor.test(cov19_data_3$diagnose, cov19_data_3$paper, method='kendall')
# cov gdp与论文数
cor.test(cov19_data_3$gdp, cov19_data_3$paper, method='spearman')
cor.test(cov19_data_3$gdp, cov19_data_3$paper, method='kendall')
# Calculate the e
e_cov_3 <- rep(1, 63)
for (i in 1:63) {
  if(cov_dia_rank_3[i] > cov_paper_rank_3[i]){
    e_cov_3[i] <- (cov_paper_rank_3[i]-cov_dia_rank_3[i])/sqrt(cov_dia_rank_3[i])
  }else if(cov_dia_rank_3[i] == cov_paper_rank_3[i]){
    e_cov_3[i] <- 0
  }else{
    e_cov_3[i] <- (cov_paper_rank_3[i]-cov_dia_rank_3[i])/sqrt(cov_paper_rank_3[i])
  }
}
cov19_data_3$e <- e_cov_3
sum(e_cov_3)
# e_gdp
e_cov_gdp_3 <- rep(1, 63)
for (i in 1:63) {
  if(cov_gdp_rank_3[i] > cov_paper_rank_3[i]){
    e_cov_gdp_3[i] <- (cov_paper_rank_3[i]-cov_gdp_rank_3[i])/sqrt(cov_gdp_rank_3[i])
  }else if(cov_dia_rank_3[i] == cov_paper_rank_3[i]){
    e_cov_gdp_3[i] <- 0
  }else{
    e_cov_gdp_3[i] <- (cov_paper_rank_3[i]-cov_gdp_rank_3[i])/sqrt(cov_paper_rank_3[i])
  }
}
cov19_data_3$e_gdp <- e_cov_gdp_3
sum(e_cov_gdp_3)
# add e=1 curve
cov_d_1 <- seq(1,63,0.1)
cov_p_1 <- (cov_d_1 - sqrt(cov_d_1))
cov_1 <- data.frame(cov_d_1, cov_p_1)
cov_p_1_a <- seq(1,63,0.1)
cov_d_1_a <- (cov_p_1_a - sqrt(cov_p_1_a))
cov_1_a <- data.frame(cov_d_1_a, cov_p_1_a)
# plot scatter plot
sp_cov_3 <- ggplot(data=cov19_data_3, aes(x=dia_rank, y=paper_rank)) +
  geom_abline(size=2) +
  theme_classic() +
  scale_x_continuous(limit=c(0,63)) +
  scale_y_continuous(limit=c(0,63)) +
  geom_point(col = ifelse(cov19_data_3$e>0, 'red', ifelse(cov19_data_3$e==0, 'black', 'blue')), size=3) +
  geom_line(data=cov_1, aes(x=cov_d_1, y=cov_p_1), size=2) +
  geom_line(data=cov_1_a, aes(x=cov_d_1_a, y=cov_p_1_a), size=2) +
  coord_fixed(ratio = 1) +
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=34,face="bold"))
# gdp
sp_cov_gdp_3 <- ggplot(data=cov19_data_3, aes(x=gdp_rank, y=paper_rank)) +
  geom_abline(size=2) +
  theme_classic() +
  scale_x_continuous(limit=c(0,63)) +
  scale_y_continuous(limit=c(0,63)) +
  geom_point(col = ifelse(cov19_data_3$e_gdp>0, 'red', ifelse(cov19_data_3$e_gdp==0, 'black', 'blue')), size=3) +
  geom_line(data=cov_1, aes(x=cov_d_1, y=cov_p_1), size=2) +
  geom_line(data=cov_1_a, aes(x=cov_d_1_a, y=cov_p_1_a), size=2) +
  coord_fixed(ratio = 1)  +
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=34,face="bold"))
# Plot PDF of e
e_cov_pdf_3 <- ggplot(cov19_data_3, aes(x=e_gdp)) +
  geom_density(color='red', size=2) +
  geom_density(aes(x=e), size=2) +
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=34,face="bold"))
# Calculate the quantile
country_cov_quantile_3 <- quantile(cov19_data_3$paper, seq(0, 1, 0.05))
# Calculate the mean
country_cov_mean_3 <- mean(cov19_data_3$paper)
# Calculate the standard deviation
country_cov_sd_3 <- sd(cov19_data_3$paper)
# Calculate the Coefficient of Variation
country_cov_cv_3 <-  country_cov_sd_3/country_cov_mean_3

#####把所有的图集中#####
plot_grid(sp_mers, sp_mers_gdp, e_mers_pdf, sp_mers_3, sp_mers_gdp_3, e_mers_pdf_3, 
          sp_sars, sp_sars_gdp, e_sars_pdf, sp_sars_3, sp_sars_gdp_3, e_sars_pdf_3, 
          sp_cov, sp_cov_gdp, e_cov_pdf, sp_cov_3, sp_cov_gdp_3, e_cov_pdf_3,
          nrow=6, ncol=3, labels='auto', label_size=50)


#####绘制国家合作网#####
NetMatrix_ccn_cov <- biblioNetwork(Mccn_cov, analysis = "collaboration", network = "countries", sep = ";")
# Plot the network
net_ccn_cov=networkPlot(NetMatrix_ccn_cov, n = dim(NetMatrix_ccn_cov)[1], 
                            Title = "Country Collaboration", type = "sphere", size=TRUE, 
                            remove.multiple=FALSE,labelsize=2,cluster="none")

############################结束############################









