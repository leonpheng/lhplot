
#' DV vs IPRED
#'
#' Generate GOF1
#' @param data data frame
#' @keywords lh_dv_ipred
#' @export
#' @examples p1<-lh_dv_ipred(dat1,type="lin",IPREDN="IPRED")
#' @examples p2<-lh_dv_ipred(dat1,type="log")


lh_dv_ipred<-function(data,y="DV",
                 x="IPRED",type="log",scale=c(0.1,100),
              IPREDN="Individual Predicted Concentration (ng/mL)",
              PREDN="Population Predicted Concentration (ng/mL)",
              DVN="Observed Concentration (ng/mL)",
              TADN="Time After Dose (h)",
              RTIMEN="Time After First Dose (h)",
              IVARN="Time After First Dose (h)",
              CWRESN="Conditional Weighted Residuals"
){
  r<-data[,c(x,y)]
  names(r)<-c("x","y")
  if("auto"%in%scale){
  limx <- range(r$x, r$y)
if(min(limx)==0){
  limx1 <-c(0.01,10^ceiling(log10(max(limx))))}else{
  limx1 <-c(10^floor(log10(min(limx))),10^ceiling(log10(max(limx))))  
  }}else{
    limx <-scale
    limx1 <-c(10^floor(log10(scale[1])),10^ceiling(log10(scale[2])))
  }

cols <- c("Observed"="#A6CEE3")
cols1 <- c("Identity"="#1F78B4")

p<-ggplot2::ggplot(r,aes(x=x,y=y))+
  geom_point(aes(col="Observed"))+
  xlab(IPREDN)+ylab(DVN)+
  geom_abline(slope=1,size=1,col="blue")+
  geom_line(aes(x=0,y=0,linetype="Identity"))+
  geom_smooth(method="loess", method.args=list(span=2/3, degree=1, family="symmetric"), se=F,aes(linetype="LOESS"))+
  scale_colour_manual(name="",values=cols) +
  scale_linetype_discrete(name = "")+
  theme_bw()
if(type=="lin"){
p=p+scale_x_continuous(limits=limx)+scale_y_continuous(limits=limx)
}else{
  p=p+scale_x_log10(limits=limx1)+scale_y_log10(limits=limx1)
}
p
}

#' DV vs PRED
#'
#' Generate GOF1
#' @param data data frame
#' @keywords lh_dv_pred
#' @export
#' @examples p1<-lh_dv_pred(dat1,type="lin",IPREDN="IPRED")
#' @examples p2<-lh_dv_pred(dat1,type="log",scale="auto")

lh_dv_pred<-function(data=dat1,y="DV",
                      x="PRED",type="log",scale=c(0.1,100),
                      IPREDN="Individual Predicted Concentration (ng/mL)",
                      PREDN="Population Predicted Concentration (ng/mL)",
                      DVN="Observed Concentration (ng/mL)",
                      TADN="Time After Dose (h)",
                      RTIMEN="Time After First Dose (h)",
                      IVARN="Time After First Dose (h)",
                      CWRESN="Conditional Weighted Residuals"
){
  r<-data[,c(x,y)]
  names(r)<-c("x","y")
  if("auto"%in%scale){
    limx <- range(r$x, r$y)
    if(min(limx)==0){
      limx1 <-c(0.01,10^ceiling(log10(max(limx))))}else{
        limx1 <-c(10^floor(log10(min(limx))),10^ceiling(log10(max(limx))))  
      }}else{
        limx <-scale
        limx1 <-c(10^floor(log10(scale[1])),10^ceiling(log10(scale[2])))
      }
  
  cols <- c("Observed"="#A6CEE3")
  cols1 <- c("Identity"="#1F78B4")
  
  p<-ggplot2::ggplot(r,aes(x=x,y=y))+
    geom_point(aes(col="Observed"))+
    xlab(PREDN)+ylab(DVN)+
    geom_abline(slope=1,size=1,col="blue")+
    geom_line(aes(x=0,y=0,linetype="Identity"))+
    geom_smooth(method="loess", method.args=list(span=2/3, degree=1, family="symmetric"), se=F,aes(linetype="LOESS"))+
    scale_colour_manual(name="",values=cols) +
    scale_linetype_discrete(name = "")+
    theme_bw()
  if(type=="lin"){
    p=p+scale_x_continuous(limits=limx)+scale_y_continuous(limits=limx)
  }else{
    p=p+scale_x_log10(limits=limx1)+scale_y_log10(limits=limx1)
  }
  p
}

#' CWRES vs TAD
#'
#' Generate GOF1
#' @param data data frame
#' @keywords lh_cwres_tad
#' @export
#' @examples data$TAD<-data$TAD;p1<-lh_cwres_tad(dat1)
#' @examples p2<-lh_cwres_tad(dat1)


lh_cwres_tad<-function(data=dat1,y="CWRES",
                     x="TAD",type="log",scale=c(0.1,100),
                     IPREDN="Individual Predicted Concentration (ng/mL)",
                     PREDN="Population Predicted Concentration (ng/mL)",
                     DVN="Observed Concentration (ng/mL)",
                     TADN="Time After Dose (h)",
                     RTIMEN="Time After First Dose (h)",
                     IVARN="Time After First Dose (h)",
                     CWRESN="Conditional Weighted Residuals"
){
  r<-data[,c(x,y)]
  names(r)<-c("x","y")
  if("auto"%in%scale){
    limx <- range(r$x, r$y)
    if(min(limx)==0){
      limx1 <-c(0.01,10^ceiling(log10(max(limx))))}else{
        limx1 <-c(10^floor(log10(min(limx))),10^ceiling(log10(max(limx))))  
      }}else{
        limx <-scale
        limx1 <-c(10^floor(log10(scale[1])),10^ceiling(log10(scale[2])))
      }
  
  cols <- c("Observed"="#A6CEE3")
  cols1 <- c("Identity"="#1F78B4")
  
  p<-ggplot2::ggplot(r,aes(x=x,y=y))+
    geom_point(aes(col="Observed"))+
    xlab(TADN)+ylab(CWRESN)+
    geom_hline(aes(yintercept=0),linetype="solid") +
    geom_hline(yintercept=c(-4,4,-6,6),linetype = "dashed",col="grey")+
    geom_smooth(method="loess", method.args=list(span=2/3, degree=1, family="symmetric"), se=F,linetype="dashed")+
    scale_x_continuous(labels = function(x) format(x, scientific =F))+
    scale_y_continuous(limits=c(-8,8),breaks=seq(-8,8,2))+
    scale_colour_manual(name="",values=cols) +
    scale_linetype_discrete(name = "")+
    theme_bw()
p  
}


#' CWRES vs RTIME or TIME
#'
#' Generate GOF1
#' @param data data frame
#' @keywords lh_cwres_time
#' @export
#' @examples data$TIME<-data$RTIME;p1<-lh_cwres_time(dat1)
#' @examples p2<-lh_cwres_time(dat1)


lh_cwres_time<-function(data=dat1,y="CWRES",
                       x="TIME",type="log",scale=c(0.1,100),
                       IPREDN="Individual Predicted Concentration (ng/mL)",
                       PREDN="Population Predicted Concentration (ng/mL)",
                       DVN="Observed Concentration (ng/mL)",
                       TADN="Time After Dose (h)",
                       RTIMEN="Time After First Dose (h)",
                       IVARN="Time After First Dose (h)",
                       CWRESN="Conditional Weighted Residuals"
){
  r<-data[,c(x,y)]
  names(r)<-c("x","y")
  if("auto"%in%scale){
    limx <- range(r$x, r$y)
    if(min(limx)==0){
      limx1 <-c(0.01,10^ceiling(log10(max(limx))))}else{
        limx1 <-c(10^floor(log10(min(limx))),10^ceiling(log10(max(limx))))  
      }}else{
        limx <-scale
        limx1 <-c(10^floor(log10(scale[1])),10^ceiling(log10(scale[2])))
      }
  
  cols <- c("Observed"="#A6CEE3")
  cols1 <- c("Identity"="#1F78B4")
  
  p<-ggplot2::ggplot(r,aes(x=x,y=y))+
    geom_point(aes(col="Observed"))+
    xlab(RTIMEN)+ylab(CWRESN)+
    geom_hline(aes(yintercept=0),linetype="solid") +
    geom_hline(yintercept=c(-4,4,-6,6),linetype = "dashed",col="grey")+
    geom_smooth(method="loess", method.args=list(span=2/3, degree=1, family="symmetric"), se=F,linetype="dashed")+
    scale_x_continuous(labels = function(x) format(x, scientific =F))+
    scale_y_continuous(limits=c(-8,8),breaks=seq(-8,8,2))+
    scale_colour_manual(name="",values=cols) +
    scale_linetype_discrete(name = "")+
    theme_bw()
  p  
}

#' CWRES vs PRED
#'
#' Generate GOF1
#' @param data data frame
#' @keywords lh_cwres_pred
#' @export
#' @examples p1<-lh_cwres_pred(dat1)
#' @examples p2<-lh_cwres_pred(dat1)

lh_cwres_pred<-function(data=dat1,y="CWRES",
                        x="TIME",type="log",scale=c(0.1,100),
                        IPREDN="Individual Predicted Concentration (ng/mL)",
                        PREDN="Population Predicted Concentration (ng/mL)",
                        DVN="Observed Concentration (ng/mL)",
                        TADN="Time After Dose (h)",
                        RTIMEN="Time After First Dose (h)",
                        IVARN="Time After First Dose (h)",
                        CWRESN="Conditional Weighted Residuals"
){
  r<-data[,c(x,y)]
  names(r)<-c("x","y")
  if("auto"%in%scale){
    limx <- range(r$x, r$y)
    if(min(limx)==0){
      limx1 <-c(0.01,10^ceiling(log10(max(limx))))}else{
        limx1 <-c(10^floor(log10(min(limx))),10^ceiling(log10(max(limx))))  
      }}else{
        limx <-scale
        limx1 <-c(10^floor(log10(scale[1])),10^ceiling(log10(scale[2])))
      }
  
  cols <- c("Observed"="#A6CEE3")
  cols1 <- c("Identity"="#1F78B4")
  
  p<-ggplot2::ggplot(r,aes(x=x,y=y))+
    geom_point(aes(col="Observed"))+
    xlab(PREDN)+ylab(CWRESN)+
    geom_hline(aes(yintercept=0),linetype="solid") +
    geom_hline(yintercept=c(-4,4,-6,6),linetype = "dashed",col="grey")+
    geom_smooth(method="loess", method.args=list(span=2/3, degree=1, family="symmetric"), se=F,linetype="dashed")+
    scale_x_continuous(labels = function(x) format(x, scientific =F))+
    scale_y_continuous(limits=c(-8,8),breaks=seq(-8,8,2))+
    scale_colour_manual(name="",values=cols) +
    scale_linetype_discrete(name = "")+
    theme_bw()
  p  
}

#' SAVE BATCH PLOTS 
#'
#' Generate GOF1
#' @param data data frame
#' @keywords lh_gof
#' @export
#' @examples p1<-lh_gof("test.png")
#' @examples p2<-lh_gof()

lh_gof<-function(file.name=NULL){
p1<-lh_dv_ipred(dat1,type="lin",IPREDN="IPRED")
p2<-lh_dv_ipred(dat1,type="log")
p3<-lh_dv_pred(dat1,type="lin")
p4<-lh_dv_pred(dat1,type="log")
p5<-lh_cwres_time(dat1)
p6<-lh_cwres_pred(dat1) 
if(is.null(file.name)){
p<-ggpubr::ggarrange(p1,p3,p5,p2,p4,p6,ncol=3, nrow=2, common.legend = TRUE, legend="bottom")
}else{
  ggsave(file.name,
         ggpubr::ggarrange(p1,p3,p5,p2,p4,p6,ncol=3, nrow=2, common.legend = TRUE, legend="bottom"),dpi = 300, width =12, height = 8,units = c("in"))  
}
p
}

