
#' GOF with Facet
#'
#' Generate GOF with stratification
#' @param data data frame
#' @keywords lh_dv_strat
#' @export
#' @examples lh_dv_strat(...)

lh_dv_strat<-function(data,x,y,id,xlab="Individual Predictions (\U03BCg/mL)",ylab="Observed\n Concentrations (\U03BCg/mL)",by="Age_Group"){
  data$x<-data[,x]
  data$y<-data[,y]
  data$id<-data[,id]
  data$strat<-data[,by]
ggplot(data,aes(x,y,label=id))+
  geom_blank(aes(x,y))  +
  geom_abline(data= data.frame(slope = 1,intercept =0 ),
              aes(slope=slope,intercept=intercept,alpha = "Identity"),
              col = "blue",size=1.5,linetype=1, key_glyph = "path")+
  geom_smooth(aes(linetype = "LOESS"),method = "loess",
              method.args = list(span = 2/3,degree = 1, family = "symmetric"),
              se=F,size=1.5)+
  geom_point(pch=16, size=1.5, aes(col="Observed"),alpha=0.2) +
  facet_grid(~strat,margins = "strat")+
  scale_colour_manual(name=NULL,values="black",guide = guide_legend(order=1,label.hjust = 1))+
  scale_alpha_manual(name=NULL,values = c(0.2), guide=guide_legend(order=2,label.hjust = 1) )+
  scale_linetype_manual(name=NULL,values=2, guide=guide_legend(order=3,label.hjust = 1))+
  labs(x=xlab, y=ylab)+
  theme_bw(base_size = 13)+
  theme(aspect.ratio = 1,legend.position="bottom",
        legend.key.width = unit(0.5,"in"),
        legend.key.height=unit(0.125,"in"))
}

#' CWRES with Facet
#'
#' Generate GOF with stratification
#' @param data data frame
#' @keywords lh_cwres_strat
#' @export
#' @examples lh_cwres_strat(...)
lh_cwres_strat<-function(data,x,y,id,xlab="Population Predictions (\U03BCg/mL)",ylab="Conditional Weighted Residual",by="Age_Group"){
  data$x<-data[,x]
  data$y<-data[,y]
  data$id<-data[,id]
  data$strat<-data[,by]

ggplot(data,aes(x,y,label=id))+
  geom_hline(yintercept=0,size=1.5,color="blue",alpha=0.2) +
  geom_hline(yintercept=4,linetype=2,color="gray") +
  geom_hline(yintercept=-4,linetype=2,color="gray") +
  geom_line(aes(x = 0, y = 0, alpha = "Identity"),size=1.5,color="blue")+
  geom_smooth(aes(linetype = "LOESS"),method = "loess",
              method.args = list(span = 2/3,degree = 1, family = "symmetric"),
              se=F,size=1.5)+
  geom_point(pch=16, size=1.5, aes(col="Observed"),alpha=0.2) +
  facet_grid(~strat,margins = "strat") +
  ylab(ylab)+ xlab(xlab)+
  ylim(-6,6)+
  scale_y_continuous(breaks=seq(-6,6,2),labels=c("-6","-4","-2","   0","2","4","6"))+
  scale_colour_manual(name=NULL,values="black",guide = guide_legend(order=1,label.hjust = 1))+
  scale_alpha_manual(name=NULL,values = c(0.2), guide=guide_legend(order=2,label.hjust = 1) )+
  scale_linetype_manual(name=NULL,values=2, guide=guide_legend(order=3,label.hjust = 1))+
  theme_bw(base_size = 13)+
  theme(aspect.ratio = 1,legend.position="bottom",
        legend.key.width = unit(0.5,"in"),
        legend.key.height=unit(0.125,"in"))
}


#' BOXPLOT MC SIMULATION
#'
#' Generate boxplot with targets and stats
#' @param data data frame
#' @keywords lhboxplot2
#' @export
#' @examples lhboxplot2(data=rall,
#' @examples y="Cmaxss",
#' @examples x="Label",
#' @examples x.title="Group",
#' @examples y.title="AUCss",
#' @examples low.targ.line=rf$Cmaxsslow,
#' @examples high.targ.line=rf$Cmaxsshi,
#' @examples add.target="yes",
#' @examples add.obs.point="yes",
#' @examples add.stats="yes",
#' @examples stat.label.space=c(0.1,0.2),jit=c(0.15,0)) +theme_bw()+ theme(axis.text.x #' @examples = element_text(angle = 45, hjust = 1))+scale_y_log10()

lhboxplot2<-function(data,y,x,x.title,y.title,low.targ.line,high.targ.line,add.target="no",add.obs.point="no",add.stats="no",stat.label.space=c(0.1,0.2),jit=c(0.1,0.1))
{
  #Compute target attainment
  library(ggplot2)
  xti=x.title
  yti=y.title
  lowt=low.targ.line
  hit=high.targ.line
  obs=add.obs.point
  prop=add.stats
  space=stat.label.space
  setdiff(c(x,y,lowt,hit),names(df))
  s1<-data[,c(x,y)]
  s1<-s1[order(s1[,x]),]
  s1$hit<-hit
  s1$lowt<-lowt

  s1$y<-s1[,y]
  s1$x<-s1[,x]
  s1$f<-"dum"##s1[,facet]

  s1$hi<-with(s1,ifelse(y>hit,1,0))
  s1$wi<-with(s1,ifelse(y<=hit&y>=lowt,1,0))
  s1$lo<-with(s1,ifelse(y<lowt,1,0))

  nn<-addvar(s1,x,x,"length(x)","var","n")
  #s1<-join(s1,nn)
  str(s1)

  targ<-addvar(s1,x,"hi","sum(x)","var","hi")
  targ<-dplyr::left_join(targ,addvar(s1,x,"wi","sum(x)","var","wi"))
  targ<-dplyr::left_join(targ,addvar(s1,x,"lo","sum(x)","var","lo"))
  targ<-dplyr::left_join(targ,nn);
  targ$hi<-with(targ,round(hi/n*100,1))
  targ$wi<-with(targ,round(wi/n*100,1))
  targ$lo<-with(targ,round(lo/n*100,1))
  #Highlight greater % by inclreasing font size
  for(i in 1:nrow(targ)){
    targ$si1[i]<-ifelse(targ$wi[i]==max(targ[i,c("hi","wi","lo")]),3,2)
    targ$si2[i]<-ifelse(targ$hi[i]==max(targ[i,c("hi","wi","lo")]),3,2)
    targ$si3[i]<-ifelse(targ$lo[i]==max(targ[i,c("hi","wi","lo")]),3,2)
  }

  var<-y
  targ$hi<-with(targ,paste0(hi,"%"))
  targ$wi<-with(targ,paste0(wi,"%"))
  targ$lo<-with(targ,paste0(lo,"%"))
  targ$allt<-with(targ,paste(hi,wi,lo,sep="\n"))

  head(s1)
  s1[,c("f","lowt")]
  lowtar<-nodup(s1,"f","add",c("lowt","hit","y"))

  var<-y
  coord<-addvar(s1,"f",var,"max(x)","add","max1")
  targ$f<-"dum"
  targ2<-dplyr::left_join(targ,coord)

  targ2$ma1<-targ2$max1+(space[2]*targ2$max1)
  targ2$me1<-targ2$max1+(space[1]*targ2$max1)
  targ2$mi1<-targ2$max1
  targ2$x<-targ2[,x]
  #targ2$f<-targ2[,facet]

  head(s1)

  p0<-ggplot(s1,aes(x=x,y=y))
  if(obs=="yes"){
    p0<-p0+geom_jitter(width =jit[1], height =jit[2],col="gray",size=0.4)
  }

  p0<- p0+ geom_boxplot(outlier.shape ="", alpha = 0.5)+
    xlab(xti)+ylab(yti)
  # if(!is.null(facet)){
  #   p0<-p0+facet_wrap(~f,scale="free")
  # }
  #if(obs=="yes"){
  #    p0<-p0+geom_jitter(width =jit[1], height =jit[2],col="gray",size=0.4)
  #  }



  if(add.target=="yes"){

    p0<-p0+ geom_hline(data=lowtar,aes(yintercept=lowt), linetype="dashed", color = "green4")+
      geom_hline(data=lowtar,aes(yintercept=hit), linetype="dashed", color = "green4")

  }

  if(prop=="yes"){
    if(hit==lowt){
      p0<-p0+geom_text(data=targ2, aes(x=x, y=me1, label=lo), col='red', size=3)+
        #geom_text(data=targ2, aes(x=x, y=me1, label=wi), col='blue', size=targ2$si1)+
        geom_text(data=targ2, aes(x=x, y=ma1, label=hi), col='green4', size=2)
    }else{p0<-p0+geom_text(data=targ2, aes(x=x, y=mi1, label=lo), col='red', size=2)+
      geom_text(data=targ2, aes(x=x, y=me1, label=wi), col='green4', size=2)+
      geom_text(data=targ2, aes(x=x, y=ma1, label=hi), col='red', size=2)
    # +
    # annotate("rect", xmin=hitar$x[1], xmax=hitar$x[2], ymin=lowtar$y[1], ymax=hitar$y[2], alpha=0.2, fill="red")
    }
  }
  p0
}





#' DV vs X with STRAT
#'
#' Generate GOF1
#' @param data data frame
#' @keywords lh_dv_ipred
#' @export
#' @examples p1<-lh_dv_ipred(dat1,type="lin",IPREDN="IPRED")
#' @examples p2<-lh_dv_ipred(dat1,type="log")
lh_dv_x<-function (data, y = "DV", x = "IPRED", type = "log",
                   scale = c(0.1, 100), IPREDN = "Individual Predicted Concentration (ng/mL)",
                   DVN = "Observed Concentration (ng/mL)", col.obs = "#A6CEE3",
                   col.ident = "#1F78B4",col.point=NULL,shape.point=NULL)
{
  cw <- data#[, c(x, y, strat)]
  cw$x <- cw[,x];cw$y <- cw[,y]
  if("auto" %in% scale) {
    limx <- range(cw$x, cw$y)
    if (min(limx) == 0) {
      limx1 <- c(0.01, 10^ceiling(log10(max(limx))))
    }else {
      limx1 <- c(10^floor(log10(min(limx))), 10^ceiling(log10(max(limx))))
    }
  } else {
    limx <- scale
    limx1 <- c(10^floor(log10(scale[1])), 10^ceiling(log10(scale[2])))
  }
  cols <- c(Observed = col.obs)
  cols1 <- c(Identity = col.ident)
  p <- ggplot(cw, aes_string(x = x, y = y))
  p<-p+geom_point(aes_string(col=col.point,shape=shape.point))+
    xlab(IPREDN) + ylab(DVN) + geom_abline(slope = 1, size = 1, col = "blue") + geom_line(aes(x = 0, y = 0, linetype = "Identity")) + geom_smooth(method = "loess",method.args = list(span = 2/3, degree = 1, family = "symmetric"),se = F, aes(linetype = "LOESS")) + scale_linetype_discrete(name = "") +
    guides(col = guide_legend(title = strat)) + theme_bw()
}



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
              CWRESN="Conditional Weighted Residuals",
              col.obs="#A6CEE3",col.ident="#1F78B4"
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

cols <- c("Observed"=col.obs)
cols1 <- c("Identity"=col.ident)

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

lh_dv_pred<-function(data,y="DV",
                      x="PRED",type="log",scale=c(0.1,100),
                      IPREDN="Individual Predicted Concentration (ng/mL)",
                      PREDN="Population Predicted Concentration (ng/mL)",
                      DVN="Observed Concentration (ng/mL)",
                      TADN="Time After Dose (h)",
                      RTIMEN="Time After First Dose (h)",
                      IVARN="Time After First Dose (h)",
                      CWRESN="Conditional Weighted Residuals",
                     col.obs="#A6CEE3",col.ident="#1F78B4"
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

  cols <- c("Observed"=col.obs)
  cols1 <- c("Identity"=col.ident)

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


lh_cwres_tad<-function(data,y="CWRES",
                     x="TAD",type="log",scale=c(0.1,100),
                     IPREDN="Individual Predicted Concentration (ng/mL)",
                     PREDN="Population Predicted Concentration (ng/mL)",
                     DVN="Observed Concentration (ng/mL)",
                     TADN="Time After Dose (h)",
                     RTIMEN="Time After First Dose (h)",
                     IVARN="Time After First Dose (h)",
                     CWRESN="Conditional Weighted Residuals",
                     col.obs="#A6CEE3",col.ident="#1F78B4"
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

  cols <- c("Observed"=col.obs)
  cols1 <- c("Identity"=col.ident)

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


lh_cwres_time<-function(data,y="CWRES",
                       x="TIME",type="log",scale=c(0.1,100),
                       IPREDN="Individual Predicted Concentration (ng/mL)",
                       PREDN="Population Predicted Concentration (ng/mL)",
                       DVN="Observed Concentration (ng/mL)",
                       TADN="Time After Dose (h)",
                       RTIMEN="Time After First Dose (h)",
                       IVARN="Time After First Dose (h)",
                       CWRESN="Conditional Weighted Residuals",
                       col.obs="#A6CEE3",col.ident="#1F78B4",strat=NULL
){
  if(!is.null(strat)){
    cw<-data[,c(x,y,strat)]
    names(cw)<-c("x","y","strat")}else{
      cw<-data[,c(x,y)]
      names(cw)<-c("x","y")}
  if("auto"%in%scale){
    limx <- range(cw$x, cw$y)
    if(min(limx)==0){
      limx1 <-c(0.01,10^ceiling(log10(max(limx))))}else{
        limx1 <-c(10^floor(log10(min(limx))),10^ceiling(log10(max(limx))))
      }}else{
        limx <-scale
        limx1 <-c(10^floor(log10(scale[1])),10^ceiling(log10(scale[2])))
      }

  cols <- c("Observed"=col.obs)
  cols1 <- c("Identity"=col.ident)
  if(is.null(strat)){
    p<-ggplot2::ggplot(cw,aes(x=x,y=y))+
      geom_point(aes(col=factor(z)))+
      xlab(PREDN)+ylab(CWRESN)+
      geom_hline(aes(yintercept=0),linetype="solid") +
      geom_hline(yintercept=c(-4,4,-6,6),linetype = "dashed",col="grey")+
      geom_smooth(method="loess", method.args=list(span=2/3, degree=1, family="symmetric"), se=F,linetype="dashed")+
      scale_x_continuous(labels = function(x) format(x, scientific =F))+
      scale_y_continuous(limits=c(-8,8),breaks=seq(-8,8,2))+
      scale_linetype_discrete(name = "")+
      guides(col=guide_legend(title=strat))+
      theme_bw()
  }else{
    cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
              "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    cw[,"strat"]<-factor(cw[,"strat"])
    p<-ggplot2::ggplot(cw,aes(x=x,y=y))+
      geom_point(aes(col=strat))+
      xlab(PREDN)+ylab(CWRESN)+
      geom_hline(aes(yintercept=0),linetype="solid") +
      geom_hline(yintercept=c(-4,4,-6,6),linetype = "dashed",col="grey")+
      geom_smooth(method="loess", method.args=list(span=2/3, degree=1, family="symmetric"), se=F,linetype="dashed")+
      scale_x_continuous(labels = function(x) format(x, scientific =F))+
      scale_y_continuous(limits=c(-8,8),breaks=seq(-8,8,2))+
      scale_colour_manual(name="Observed",values=cbp1) +
      scale_linetype_discrete(name = "")+
      theme_bw()
  }
  p
}

#' CWRES vs X with Strat
#'
#' Generate GOF1
#' @param data data frame
#' @keywords lh_cwres_pred
#' @export
#' @examples p1<-lh_cwres_pred(dat1)
#' @examples p2<-lh_cwres_pred(dat1)

lh_cwres_x<-function(data,y="CWRES",
                        x="PREDN",type="log",scale=c(0.1,100),
                        PREDN="Population Predicted Concentration (ng/mL)",
                        CWRESN="Conditional Weighted Residuals",
                        col.obs="#A6CEE3",col.ident="#1F78B4",strat=NULL){

  if(!is.null(strat)){
    cw<-data[,c(x,y,strat)]
    names(cw)<-c("x","y","strat")}else{
      cw<-data[,c(x,y)]
      names(cw)<-c("x","y")}
  if("auto"%in%scale){
    limx <- range(cw$x, cw$y)
    if(min(limx)==0){
      limx1 <-c(0.01,10^ceiling(log10(max(limx))))}else{
        limx1 <-c(10^floor(log10(min(limx))),10^ceiling(log10(max(limx))))
      }}else{
        limx <-scale
        limx1 <-c(10^floor(log10(scale[1])),10^ceiling(log10(scale[2])))
      }

  cols <- c("Observed"=col.obs)
  cols1 <- c("Identity"=col.ident)
if(is.null(strat)){
    p<-ggplot2::ggplot(cw,aes(x=x,y=y))+
      geom_point(aes(col=factor(z)))+
      xlab(PREDN)+ylab(CWRESN)+
      geom_hline(aes(yintercept=0),linetype="solid") +
      geom_hline(yintercept=c(-4,4,-6,6),linetype = "dashed",col="grey")+
      geom_smooth(method="loess", method.args=list(span=2/3, degree=1, family="symmetric"), se=F,linetype="dashed")+
      scale_x_continuous(labels = function(x) format(x, scientific =F))+
      scale_y_continuous(limits=c(-8,8),breaks=seq(-8,8,2))+
      scale_linetype_discrete(name = "")+
      guides(col=guide_legend(title=strat))+
      theme_bw()
  }else{
    cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
              "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    cw[,"strat"]<-factor(cw[,"strat"])
    p<-ggplot2::ggplot(cw,aes(x=x,y=y))+
      geom_point(aes(col=strat))+
      xlab(PREDN)+ylab(CWRESN)+
      geom_hline(aes(yintercept=0),linetype="solid") +
      geom_hline(yintercept=c(-4,4,-6,6),linetype = "dashed",col="grey")+
      geom_smooth(method="loess", method.args=list(span=2/3, degree=1, family="symmetric"), se=F,linetype="dashed")+
      scale_x_continuous(labels = function(x) format(x, scientific =F))+
      scale_y_continuous(limits=c(-8,8),breaks=seq(-8,8,2))+
      scale_colour_manual(name="Observed",values=cbp1) +
      scale_linetype_discrete(name = "")+
      theme_bw()
  }
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

lh_cwres_pred<-function(data,y="CWRES",
                        x="PRED",type="log",scale=c(0.1,100),
                        IPREDN="Individual Predicted Concentration (ng/mL)",
                        PREDN="Population Predicted Concentration (ng/mL)",
                        DVN="Observed Concentration (ng/mL)",
                        TADN="Time After Dose (h)",
                        RTIMEN="Time After First Dose (h)",
                        IVARN="Time After First Dose (h)",
                        CWRESN="Conditional Weighted Residuals",
                        col.obs="#A6CEE3",col.ident="#1F78B4",strat=NULL){
  if(!is.null(strat)){
    cw<-data[,c(x,y,strat)]
    names(cw)<-c("x","y","strat")}else{
      cw<-data[,c(x,y)]
      names(cw)<-c("x","y")}
  if("auto"%in%scale){
    limx <- range(cw$x, cw$y)
    if(min(limx)==0){
      limx1 <-c(0.01,10^ceiling(log10(max(limx))))}else{
        limx1 <-c(10^floor(log10(min(limx))),10^ceiling(log10(max(limx))))
      }}else{
        limx <-scale
        limx1 <-c(10^floor(log10(scale[1])),10^ceiling(log10(scale[2])))
      }

  cols <- c("Observed"=col.obs)
  cols1 <- c("Identity"=col.ident)
  if(is.null(strat)){
    p<-ggplot2::ggplot(cw,aes(x=x,y=y))+
      geom_point(aes(col=factor(z)))+
      xlab(PREDN)+ylab(CWRESN)+
      geom_hline(aes(yintercept=0),linetype="solid") +
      geom_hline(yintercept=c(-4,4,-6,6),linetype = "dashed",col="grey")+
      geom_smooth(method="loess", method.args=list(span=2/3, degree=1, family="symmetric"), se=F,linetype="dashed")+
      scale_x_continuous(labels = function(x) format(x, scientific =F))+
      scale_y_continuous(limits=c(-8,8),breaks=seq(-8,8,2))+
      scale_linetype_discrete(name = "")+
      guides(col=guide_legend(title=strat))+
      theme_bw()
  }else{
    cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
              "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    cw[,"strat"]<-factor(cw[,"strat"])
    p<-ggplot2::ggplot(cw,aes(x=x,y=y))+
      geom_point(aes(col=strat))+
      xlab(PREDN)+ylab(CWRESN)+
      geom_hline(aes(yintercept=0),linetype="solid") +
      geom_hline(yintercept=c(-4,4,-6,6),linetype = "dashed",col="grey")+
      geom_smooth(method="loess", method.args=list(span=2/3, degree=1, family="symmetric"), se=F,linetype="dashed")+
      scale_x_continuous(labels = function(x) format(x, scientific =F))+
      scale_y_continuous(limits=c(-8,8),breaks=seq(-8,8,2))+
      scale_colour_manual(name="Observed",values=cbp1) +
      scale_linetype_discrete(name = "")+
      theme_bw()
  }
  p
}


#' SAVE BATCH PLOTS
#'
#' Generate GOF1
#' @param data data frame
#' @keywords lh_gof
#' @export
#' @examples lh_gof("test.png")
#' @examples lh_gof()

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

####BOXPLOT
#' CATEGORICAL BOXPLOT vs ETA
#'
#' Generate COVAR
#' @param data Data frame, merged ETA and COVAR data
#' @param lst.eta List of ETA names
#' @param lst.cov List of covariate names. Plots are generated in loop if more than one covariate
#' @keywords lh_cat_cov
#' @export
#' @examples p1<-lh_cat_cov(data=cateta,lst.eta=keta,lst.cov=cat,save.path=NULL)
#' @examples p2<-lh_gof()
lh_cat_cov<-function(data,lst.eta=c("ETA1"),lst.cov=c("SEX","RACE"),save.path=NULL,fancy="yes"){
  cat1 <- lhlong(data, lst.cov)
  names(cat1)[names(cat1) == "variable"] <- "Covariate"
  names(cat1)[names(cat1) == "value"] <- "Categorical"
  cat1 <- chclass(cat1, c("Covariate", "Categorical"), "char")
  cat1 <- addvar(cat1, c("Covariate", "Categorical"), lst.eta[1],
                 "length(x)", "yes", "count")

  cat1$Cat1 <- paste0(cat1$Categorical, "\n (n=", cat1$count,
                      ")")

  cat1 <- lhlong(cat1, lst.eta)
  # head(cat1)
  cat1 <- chclass(cat1, c("Covariate", "Categorical", "variable"),
                  "char")
  unique(cat1$Categorical)
  head(cat1)
  cat1$variable <- factor(cat1$variable, levels = lst.eta)
  catnum <- addvar(nodup(cat1, c("Covariate", "Categorical"),
                         "var"), "Covariate", "Categorical", "length(x)", "no",
                   "catnumber")
  for (i in lst.cov[lst.cov %in%catnum$Covariate[catnum$catnumber >
                                                 0]]) {
    #def1$VARN <- tolower(def1$VARN)
    #z <- def1$LABEL[def1$VARN == i]

    dcat <- cat1[cat1$Covariate %in% i, ]
    ord <- sort(unique(data[, i]))
    label <- nodup(dcat, c("Categorical", "Cat1"), "var")
    label$Categorical <- factor(label$Categorical, levels = ord)
    lablel <- label$Cat1[order(label$Categorical)]
    dcat$Cat1 <- factor(dcat$Cat1, levels = lablel)
    head(dcat)
    if (!is.null(fancy)) {
      dcat$variable1 <- gsub("ETA", "", dcat$variable)
      dcat$variable1 <- paste0("\U03B7",dcat$variable1)
      dcat<-lhfactor(dcat,"variable","variable1")
    }else{dcat$variable1<-dcat$variable}


    p <- ggplot2::ggplot(dcat, aes(x = Cat1, y = value)) +
      geom_boxplot(outlier.shape = NA) + geom_jitter(position = position_jitter(0.1),
                                                     col = "grey") + geom_hline(yintercept = 0, linetype = 2,
                                                                                color = "red", size = 1) + ylab("Individual Random Effect") +
      xlab("") + facet_wrap(~variable1, scale = "free",
                            ncol = 2) + theme_bw() + theme(axis.text.x = element_text(angle = 45,
                                                                                      hjust = 0.1, vjust = 0.4, size = 10))
    if (!is.null(save.path)) {
      nm <- paste0(save.path, z, "_boxplot.png")
      ggsave(nm, p, width = 12, height = 12)
    }else {
      p
    }
  }
  p
}


####ETA DISTRIBUTION
#' DISTRIBUTION OF ETA
#'
#' Generate GOF
#' @param data Data frame, merged ETA and COVAR data
#' @param lst.eta List of ETA names
#' @keywords lh_con_cov
#' @export
#' @examples p1<-lh_con_cov(data=cateta,lst.eta=keta,lst.cov=cat,save.path=NULL)

lh_eta_dist<-function(data,lst.eta=c("ETACL","ETAVC","ETAV2","ETAQ"),ncol=3,nrow=2,fancy="yes"){
  if(!is.null(fancy)){
    xname<-gsub("ETA","",toupper(lst.eta))
 xname<-paste0("\U03B7",xname)
  }else{xname<-lst.eta}
  head(data)
myplots <- list()
for(i in 1:length(lst.eta)){
p1<-ggplot2::ggplot(data,aes_string(x=lst.eta[i]))+geom_density(fill="royalblue3",col=NA,alpha=0.3,bin=60)+
  geom_histogram(aes_string(x=lst.eta[i],y="..density.."),fill=NA,col="black")+
  xlab(xname[i])+ylab("Density")+
  geom_vline(aes(xintercept=0,col="Zero",linetype = "Zero"),size=1.2)+
  geom_vline(aes(xintercept=mean(data[,lst.eta[i]]),col="Mean",linetype = "Mean"),size=1.2)+
  scale_colour_manual(name="",
                      values=c(Mean="red", Zero="blue"))+
  scale_linetype_manual(name="", values = c(Mean = "dashed", Zero = "dashed"))+
  theme_bw()+
  theme(legend.position = c(0.8, 0.8),legend.background = element_rect(fill=NULL,colour =NULL))
myplots[[i]]<-p1
}
ggarrange(plotlist = myplots,ncol=ncol, nrow=nrow, common.legend = TRUE, legend="bottom")
}

#' DISTRIBUTION OF CWRES internal
#'
#' Generate GOF
#' @param data Data frame, merged ETA and COVAR data
#' @keywords qqplot.cwres
#' @export
#' @examples p1<-lh_con_cov(data=cateta,lst.eta=keta,lst.cov=cat,save.path=NULL)

qqplot.cwres <- function(dat, ...) {
  ylim <-c(-10, 10)
  xlim <-c(-4, 4)
  xlab <- "Quantiles of Standard Normal"
  ylab <- "Conditional Weighted Residuals"
  with(dat, qqnorm(cwres, ylim=ylim, xlim=xlim, xlab=xlab, ylab=ylab, ...))
  with(dat, qqline(cwres))
  abline(a=0,b=1,col="red") #add the standard normal qqplot
}


#' DISTRIBUTION OF CWRES internal
#'
#' Generate GOF
#' @param data Data frame, merged ETA and COVAR data
#' @keywords histogram.cwres
#' @export
#' @examples p1<-lh_con_cov(data=cateta,lst.eta=keta,lst.cov=cat,save.path=NULL)

hist.cwres <- function(dat, ...) {
  ylim <-c(0, 0.55)
  xlab <- "Conditional Weighted Residuals"
  with(dat, hist(cwres, ylim=ylim, xlab=xlab, main="", freq=FALSE, ...))
  abline(v=0, lty=2, lwd=3, col="gray")
  xs <- seq(-10, 140, len=100)
  lines(xs, dnorm(xs), col="gray", lwd=3)
}

#' DISTRIBUTION OF CWRES
#'
#' Generate GOF
#' @param data Data frame, merged ETA and COVAR data
#' @keywords lh_cwres_dist
#' @export
#' @examples p1<-lh_con_cov(data=cateta,lst.eta=keta,lst.cov=cat,save.path=NULL)

lh_cwres_dist<-function(data,cwres="cwres",file.name="dist_QQ_CWRES.png"){
names(data)[names(data)==cwres]<-"cwres"

png(file =file.name , width = 8, height = 6, units = 'in', res = 300)
par(mfrow = c(1, 2))
qqplot.cwres(data)
hist.cwres(data)
dev.off()}


#' Individual plots
#'
#' Generate GOF
#' @param data Data frame
#' @param output.name Listing is generated in word document.
#' @keywords lh_indiv_plot
#' @export
#' @examples plh_indiv_plot(data=dat1,id="usubjid",n.plots.page=9,time="time",dv="dv",ipred="ipred"#' #' @examples ,pred="pred",type="linear",
#' @examples xtit="Time after first dose (h)",
#' @examplesytit="Concentration (ng/mL)",output.name="Individiual.docx")

lh_indiv_plot<-function(data,id="usubjid",n.plots.page=9,time="time",dv="dv",ipred="ipred",pred="pred",type="linear",xtit="Time after first dose (h)",ytit="Concentration (ng/mL)",output.name="./test.docx")
  {
  library(scales)
  library(ggplot2)
  npepage<-n.plots.page
  n<-length(unique(data[,id]))
  page<-1:ceiling(n/npepage)
  nb_pg2<-page*9
  nb_pg1<-c(1,nb_pg2[1:(length(nb_pg2)-1)]+1)

  doc<-officer::read_docx()
  for(i in 1:length(nb_pg2)){
    ddat<-data[data[,id]%in%unique(data[,id])[nb_pg1[i]:nb_pg2[i]],]
    ddat$usubjid<-ddat[,id]
    break2<-c(0.0001,0.0005,0.001,0.005,0.1,0.5,1,5,10,100,10^3,10^4,10^5,10^6)
  p<-ggplot2::ggplot(ddat,aes_string(x=time,y=dv))+
      ggplot2::geom_point(aes(col="Observed"))
    if(!is.null(ipred)){
      p<-p+ggplot2::geom_line(aes(y=ipred,col="IPRED"))
    }
    if(!is.null(pred)){
      p<-p+ggplot2::geom_line(aes(y=pred,col="PRED"))
    }
    p<-p+ggplot2::facet_wrap(~usubjid,scales="free")
    if(type=="log"){
    p=p+ ggplot2::scale_y_log10(breaks = break2)}
    p=p+ ggplot2::scale_x_continuous()+
      ggplot2::xlab(xtit)+ggplot2::ylab(ytit)+
      ggplot2::theme_bw()+
      ggplot2::theme(legend.title=element_blank())

    doc<-officer::body_add_gg(doc,p,width=7.08,height=5.9)
  }
  if(!is.null(output.name)){
  print(doc,output.name)}else{doc}
}


#' EXPLORATORY Individual plots
#'
#' Generate DATA
#' @param data Data frame
#' @param output.name Listing is generated in word document.
#' @keywords lh_indiv_plot
#' @export
#' @examples lh_explor_ind(data=dat1,id="usubjid",n.plots.page=9,time="time",dv="dv",ipred="ipred"#' #' @examples ,pred="pred",type="linear",
#' @examples xtit="Time after first dose (h)",
#' @examplesytit="Concentration (ng/mL)",output.name="Individiual.docx")

lh_explor_ind<-function(data,dose="amt",id="id",n.plots.page=9,time="time",dv="dv",ipred=NULL,pred=NULL,type="linear",xtit="Time after first dose (h)",ytit="Concentration (ng/mL)",output.name="./test.docx")
{
  data[,dose][!is.na(data[,dose])&data[,dose]==0]<-NA
  data[,dv][!is.na(data[,dose])]<-NA
  head(data)
  library(scales)
  library(ggplot2)
  npepage<-n.plots.page
  n<-length(unique(data[,id]))
  page<-1:ceiling(n/npepage)
  nb_pg2<-page*9
  nb_pg1<-c(1,nb_pg2[1:(length(nb_pg2)-1)]+1)
  head(data)

  doc<-officer::read_docx()
  for(i in 1:length(nb_pg2)){
    ddat<-data[data[,id]%in%unique(data[,id])[nb_pg1[i]:nb_pg2[i]],]
    ddat$usubjid<-ddat[,id]
ddat<-dplyr::left_join(ddat,lhtool2::addvar(ddat[!is.na(ddat[,dv])&ddat[,dv]>0,c(id,dv)],id,dv,"min(x)/2","no","dose"))
  ddat$time<-ddat[,time]
    break2<-c(0.0001,0.0005,0.001,0.005,0.1,0.5,1,5,10,100,10^3,10^4,10^5,10^6)
p<-ggplot2::ggplot(ddat[is.na(ddat$amt),],aes_string(x=time,y=dv))+
      ggplot2::geom_point(aes(col="Observed"))+
      ggplot2::geom_line()
    if(!is.null(ipred)){
      p<-p+ggplot2::geom_line(aes(y=ipred,linetype="IPRED"),col="blue")
    }
    if(!is.null(pred)){
      p<-p+ggplot2::geom_line(aes(y=pred,linetype="PRED"),col="red")
    }
    p<-p+ggplot2::geom_point(data=ddat[!is.na(ddat$amt),],aes(x=time,y=dose,col="Dose"),shape=17)
    p<-p+ggplot2::facet_wrap(~usubjid,scales="free")
    if(type=="log"){
      p=p+ ggplot2::scale_y_log10(breaks = break2)}
    p=p+ ggplot2::scale_x_continuous()+
      ggplot2::xlab(xtit)+ggplot2::ylab(ytit)+
      ggplot2::theme_bw()+
      ggplot2::theme(legend.title=element_blank())

    doc<-officer::body_add_gg(doc,p,width=7.08,height=5.9)
  }
  if(!is.null(output.name)){
    print(doc,output.name)}else{doc}
}


#############################################
#' CONTINUOUS SCATTER vs ETA
#'
#' Generate COVAR
#' @param data Data frame, merged ETA and COVAR data
#' @param lst.eta List of ETA names
#' @param lst.cov List of covariate names. Add plot name to save.path
#' @keywords lh_con_cov
#' @export
#' @examples p1<-lh_con_cov(data=cateta,lst.eta=keta,lst.cov=cat,save.path=NULL)

lh_con_cov<-function(data,lst.eta=c("eta1","eta2"),lst.cov=c("AGE","WT"),save.path="./scatter.png",fancy="yes"){
  library(lattice)
  library(grid)
  if(!is.null(fancy)){
    names(data)[names(data)%in%lst.eta]<-gsub("ETA","",names(data)[names(data)%in%lst.eta])
    lst.eta<-gsub("ETA","",lst.eta)
    names(data)[names(data)%in%lst.eta]<-paste0("\U03B7",names(data)[names(data)%in%lst.eta])
    lst.eta<-paste0("\U03B7",lst.eta)
  }

  png(save.path,width=768,height=768,pointsize = 16)
  print(gpairs(x=data[,c(lst.eta,lst.cov)],
               upper.pars = list(conditional = 'boxplot', scatter = 'loess'),
               lower.pars = list(scatter = 'stats',conditional = "barcode"),
               diag.pars = list(fontsize = 9, show.hist = TRUE, hist.color = "gray"),
               stat.pars =list(fontsize = 11, signif =F, verbose =T, use.color = TRUE, missing = 'missing', just = 'centre'),
               scatter.pars = list(pch = 20)))
  dev.off()}



##########################################################################################################################
#' CONTINUOUS GP EXTERNAL FUNCTION FOR INTERNAL USE
#'
#' Generate COVAR
#' @param x not to be used internally. Could be downloaded from Github
#' @keywords gpair
#' @export


gpairs<-function (x, upper.pars = list(scatter = "points", conditional = "barcode",
                                       mosaic = "mosaic"), lower.pars = list(scatter = "points",
                                                                             conditional = "boxplot", mosaic = "mosaic"), diagonal = "default",
                  outer.margins = list(bottom = unit(2, "lines"), left = unit(2,
                                                                              "lines"), top = unit(2, "lines"), right = unit(2, "lines")),
                  xylim = NULL, outer.labels = NULL, outer.rot = c(0, 90),
                  gap = 0.05, buffer = 0.02, reorder = NULL, cluster.pars = NULL,
                  stat.pars = NULL, scatter.pars = NULL, bwplot.pars = NULL,
                  stripplot.pars = NULL, barcode.pars = NULL, mosaic.pars = NULL,
                  axis.pars = NULL, diag.pars = NULL, whatis = FALSE)
{

  if (!is.data.frame(x)) {
    if (is.matrix(x))
      x <- as.data.frame(x)
    else stop("What did you give me? You might want to use Excel. (Only one column in argument to gpairs.\n\n")
  }
  zc <- function(x) length(unique(x)) <= 1
  if (any(sapply(x, zc), na.rm = TRUE)) {
    warning(paste(sum(sapply(x, zc), na.rm = TRUE), "columns with less than two distinct values eliminated"))
    x <- x[, !(sapply(x, zc))]
  }
  if (!is.null(lower.pars) & !is.list(lower.pars)) {
    warning("lower.pars is not a list, proceed with caution.")
  }
  if (!is.null(upper.pars) & !is.list(upper.pars)) {
    warning("upper.pars is not a list, proceed with caution.")
  }
  if (!is.null(reorder)) {
    if (pmatch(reorder, "cluster", nomatch = FALSE)) {
      if (is.null(cluster.pars)) {
        cluster.pars <- list(dist.method = "euclidean",
                             hclust.method = "complete")
      }
      x.num <- as.matrix(as.data.frame(lapply(x, as.numeric)))
      x.clust <- hclust(dist(t(x.num), method = cluster.pars$dist.method),
                        method = cluster.pars$hclust.method)
      x <- x[, x.clust$order]
    }
  }
  if (is.null(lower.pars$scatter.pars)) {
    lower.pars$scatter.pars <- "points"
  }
  if (is.null(lower.pars$conditional)) {
    lower.pars$conditional <- "boxplot"
  }
  if (is.null(lower.pars$mosaic)) {
    lower.pars$mosaic <- "mosaic"
  }
  if (is.null(upper.pars$scatter.pars)) {
    upper.pars$scatter.pars <- "points"
  }
  if (is.null(upper.pars$conditional)) {
    upper.pars$conditional <- "barcode"
  }
  if (is.null(upper.pars$mosaic)) {
    upper.pars$mosaic <- "mosaic"
  }
  if (!is.list(outer.margins)) {
    if (length(outer.margins) == 4) {
      if (is.unit(outer.margins[1])) {
        outer.margins <- list(bottom = outer.margins[1],
                              left = outer.margins[2], top = outer.margins[3],
                              right = outer.margins[4])
      } else {
        outer.margins <- list(bottom = unit(outer.margins[1],
                                            "lines"), left = unit(outer.margins[2], "lines"),
                              top = unit(outer.margins[3], "lines"), right = unit(outer.margins[4],
                                                                                  "lines"))
      }
    } else {
      stop("outer.margins are not valid.")
    }
  }
  if (is.null(outer.labels)) {
    outer.labels$top <- rep(FALSE, ncol(x))
    outer.labels$top[seq(2, ncol(x), by = 2)] <- TRUE
    outer.labels$left <- rep(FALSE, ncol(x))
    outer.labels$left[seq(2, ncol(x), by = 2)] <- TRUE
    outer.labels$right <- !outer.labels$left
    outer.labels$bottom <- !outer.labels$top
  } else {
    if (pmatch(as.character(outer.labels), "all", nomatch = FALSE)) {
      all.labeling <- TRUE
    } else if (pmatch(as.character(outer.labels), "none", nomatch = FALSE)) {
      all.labeling <- FALSE
    } else {
      stop("argument to outer.labels not understood\n")
    }
    outer.labels <- NULL
    outer.labels$top <- rep(all.labeling, ncol(x))
    outer.labels$left <- rep(all.labeling, ncol(x))
    outer.labels$bottom <- rep(all.labeling, ncol(x))
    outer.labels$right <- rep(all.labeling, ncol(x))
  }
  if (is.null(stat.pars$fontsize)) {
    stat.pars$fontsize <- 7
  }
  if (is.null(stat.pars$signif)) {
    stat.pars$signif <- 0.05
  }
  if (is.null(stat.pars$verbose)) {
    stat.pars$verbose <- FALSE
  }
  if (is.null(stat.pars$use.color)) {
    stat.pars$use.color <- TRUE
  }
  if (is.null(stat.pars$missing)) {
    stat.pars$missing <- "missing"
  }
  if (is.null(stat.pars$just)) {
    stat.pars$just <- "centre"
  }
  if (is.null(scatter.pars$pch)) {
    scatter.pars$pch <- 1
  }
  if (is.null(scatter.pars$size)) {
    scatter.pars$size <- unit(0.25, "char")
  }
  if (is.null(scatter.pars$col)) {
    scatter.pars$col <- "black"
  }
  if (is.null(scatter.pars$plotpoints)) {
    scatter.pars$plotpoints <- TRUE
  }
  if (is.null(axis.pars$n.ticks)) {
    axis.pars$n.ticks <- 5
  }
  if (is.null(axis.pars$fontsize)) {
    axis.pars$fontsize <- 9
  }
  if (axis.pars$n.ticks < 3) {
    axis.pars$n.ticks <- 3
    warning("Fewer than 3 axis ticks might cause problems.")
  }
  if (is.null(diag.pars$fontsize)) {
    diag.pars$fontsize <- 9
  }
  if (is.null(diag.pars$show.hist)) {
    diag.pars$show.hist <- TRUE
  }
  if (is.null(diag.pars$hist.color)) {
    diag.pars$hist.color <- "black"
  }
  if (is.null(stripplot.pars$pch)) {
    stripplot.pars$pch <- 1
  }
  if (is.null(stripplot.pars$size)) {
    stripplot.pars$size <- unit(0.5, "char")
  }
  if (is.null(stripplot.pars$col)) {
    stripplot.pars$col <- "black"
  }
  if (is.null(stripplot.pars$jitter)) {
    stripplot.pars$jitter <- FALSE
  }
  if (is.null(barcode.pars$nint)) {
    barcode.pars$nint <- 0
  }
  if (is.null(barcode.pars$ptsize)) {
    barcode.pars$ptsize <- unit(0.25, "char")
  }
  if (is.null(barcode.pars$ptpch)) {
    barcode.pars$ptpch <- 1
  }
  if (is.null(barcode.pars$bcspace)) {
    barcode.pars$bcspace <- NULL
  }
  if (is.null(barcode.pars$use.points)) {
    barcode.pars$use.points <- FALSE
  }
  if (is.null(mosaic.pars$gp_labels)) {
    mosaic.pars$gp_labels <- gpar(fontsize = 9)
  }
  if (is.null(mosaic.pars$gp_args)) {
    mosaic.pars$gp_args <- list()
  }
  draw.axis <- function(x, y, axis.pars, xpos, ypos, cat.labels = NULL,
                        horiz = NULL, xlim = NULL, ylim = NULL) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    if (is.null(xlim)) {
      px <- pretty(x, axis.pars$n.ticks)
      px <- px[px > min(x, na.rm = TRUE) & px < max(x,
                                                    na.rm = TRUE)]
    } else {
      px <- pretty(xlim, axis.pars$n.ticks)
      px <- px[px > min(xlim, na.rm = TRUE) & px < max(xlim,
                                                       na.rm = TRUE)]
    }
    if (is.null(ylim)) {
      py <- pretty(y, axis.pars$n.ticks)
      py <- py[py > min(y, na.rm = TRUE) & py < max(y,
                                                    na.rm = TRUE)]
    } else {
      py <- pretty(ylim, axis.pars$n.ticks)
      py <- py[py > min(ylim, na.rm = TRUE) & py < max(ylim,
                                                       na.rm = TRUE)]
    }
    k <- length(cat.labels)
    if (!is.null(xpos)) {
      if (!is.null(cat.labels) && !horiz) {
        grid.text(cat.labels, x = unit(1:k, "native"),
                  y = unit(rep(1 * (1 - xpos), k), "npc") + unit(rep(-1 *
                                                                       xpos + 1 * (1 - xpos), k), "lines"), rot = outer.rot[1],
                  gp = gpar(fontsize = axis.pars$fontsize))
      } else grid.xaxis(at = px, gp = gpar(fontsize = axis.pars$fontsize),
                        main = xpos)
    }
    if (!is.null(ypos)) {
      if (!is.null(cat.labels) && horiz) {
        grid.text(cat.labels, y = unit(1:k, "native"),
                  x = unit(rep(1 * (1 - ypos), k), "npc") + unit(rep(-1 *
                                                                       ypos + 1 * (1 - ypos), k), "lines"), rot = outer.rot[2],
                  gp = gpar(fontsize = axis.pars$fontsize))
      }else grid.yaxis(at = py, gp = gpar(fontsize = axis.pars$fontsize),
                       main = ypos)
    }
  }
  qq.panel <- function(x, y, scatter.pars, axis.pars, xpos,
                       ypos, xlim, ylim) {
    pushViewport(viewport(xscale = xlim, yscale = ylim))
    draw.axis(x, y, axis.pars, xpos, ypos, NULL, NULL, xlim,
              ylim)
    popViewport(1)
    pushViewport(viewport(xscale = xlim, yscale = ylim, clip = TRUE))
    grid.rect(gp = gpar(fill = scatter.pars$frame.fill, col = scatter.pars$border.col))
    x <- sort(x)
    y <- sort(y)
    grid.lines(unit(x, "native"), unit(y, "native"))
    popViewport(1)
  }
  scatterplot.panel <- function(x, y, type, scatter.pars, axis.pars,
                                xpos, ypos, xylim) {
    if (is.null(xylim)) {
      xlim <- range(x, na.rm = TRUE) + c(-buffer * (max(x,
                                                        na.rm = TRUE) - min(x, na.rm = TRUE)), buffer *
                                           (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
      ylim <- range(y, na.rm = TRUE) + c(-buffer * (max(y,
                                                        na.rm = TRUE) - min(y, na.rm = TRUE)), buffer *
                                           (max(y, na.rm = TRUE) - min(y, na.rm = TRUE)))
    }    else {
      xlim <- xylim
      ylim <- xylim
    }
    pushViewport(viewport(xscale = xlim, yscale = ylim))
    draw.axis(x, y, axis.pars, xpos, ypos, NULL, NULL, xlim,
              ylim)
    popViewport(1)
    pushViewport(viewport(xscale = xlim, yscale = ylim, clip = TRUE))
    grid.rect(gp = gpar(fill = scatter.pars$frame.fill, col = scatter.pars$border.col))
    if (scatter.pars$plotpoints & (type == "points" || type ==
                                   "lm" || type == "ci" || type == "symlm" || type ==
                                   "loess")) {
      grid.points(x, y, pch = scatter.pars$pch, size = scatter.pars$size,
                  gp = gpar(col = scatter.pars$col))
    }
    if (type == "lm") {
      xy.lm <- lm(y ~ x)
      panel.abline(xy.lm$coef[1], xy.lm$coef[2], col = "red",
                   lwd = 2)
    }
    if (type == "ci") {
      xy.lm <- lm(y ~ x)
      xy <- data.frame(x = seq(min(x, na.rm = TRUE), max(x,
                                                         na.rm = TRUE), length.out = 20))
      yhat <- predict(xy.lm, newdata = xy, interval = "confidence")
      ci <- data.frame(lower = yhat[, "lwr"], upper = yhat[,
                                                           "upr"])
      grid.lines(x = c(xy$x), y = c(ci$lower), default.units = "native")
      grid.lines(x = c(xy$x), y = c(ci$upper), default.units = "native")
      grid.polygon(x = c(xy$x, xy$x[length(xy$x):1]), y = c(ci$lower,
                                                            ci$upper[length(ci$upper):1]), gp = gpar(fill = "grey"),
                   default.units = "native")
    }
    if (type == "loess") {
      junk <- try(panel.loess(x, y, color = "red", span = 1))
      if (class(junk) == "try-error")
        warning("An error in loess occurred and was ignored; no line was plotted.")
    }
    if (type == "symlm") {
      pcs <- try(prcomp(cbind(x, y)))
      if (class(pcs) == "try-error")
        warning("An error in symlm occurred and was ignored; no line was plotted.")
      else {
        slope <- abs(pcs$rotation[1, 2]/pcs$rotation[1,
                                                     1])
        if (cor(x, y) < 0)
          slope <- -1 * slope
        panel.abline(pcs$center[2] - slope * pcs$center[1],
                     slope, col = "blue")
      }
    }
    if (type == "corrgram") {
      pear.test <- cor.test(x, y, method = "pearson", alternative = "two.sided")
      corr <- format(pear.test$estimate, digits = 2)
      if (as.numeric(corr) > 0) {
        panel.fill(col = hsv(h = 0.5, s = abs(as.numeric(corr)),
                             v = 1), border = hsv(h = 0.5, s = abs(as.numeric(corr)),
                                                  v = 1))
        grid.lines(x = unit(c(0, 1), "npc"), y = unit(c(0,
                                                        1), "npc"), gp = gpar(col = "white", lwd = 2))
      } else {
        panel.fill(col = hsv(h = 0, s = abs(as.numeric(corr)),
                             v = 1), border = hsv(h = 0, s = abs(as.numeric(corr)),
                                                  v = 1))
        grid.lines(x = unit(c(0, 1), "npc"), y = unit(c(1,
                                                        0), "npc"), gp = gpar(col = "white", lwd = 2))
      }
    }
    if (type == "qqplot") {
      qq.panel(x, y, scatter.pars, axis.pars, xpos, ypos,
               xlim, ylim)
    }
    if (type == "stats") {
      complete.obs <- nrow(na.omit(cbind(x, y)))
      missing <- length(x) - complete.obs
      pear.test <- cor.test(x, y, method = "pearson", alternative = "two.sided")
      corr <- sprintf("%03.2f", pear.test$estimate)
      rho.test <- cor.test(x, y, method = "spearman", alternative = "two.sided")
      tau.test <- cor.test(x, y, method = "kendall", alternative = "two.sided")
      rho <- sprintf("%03.2f", rho.test$estimate)
      tau <- sprintf("%03.2f", tau.test$estimate)
      xy.lm <- lm(y ~ x)
      r2 <- sprintf("%03.2f", summary(xy.lm)$r.squared)
      p <- sprintf("%06.4f", pf(q = as.numeric(summary(xy.lm)$fstatistic)[1],
                                df1 = as.numeric(summary(lm(xy.lm))$fstatistic)[2],
                                df2 = as.numeric(summary(lm(xy.lm))$fstatistic)[3],
                                lower.tail = FALSE))
      bonfp <- stat.pars$signif/(N * (N - 1))/2
      sig <- 1
      sigrho <- NULL
      sigtau <- NULL
      sigcor <- NULL
      sigp <- NULL
      if (pear.test$p.value < bonfp) {
        sig <- sig + 1
        sigcor <- "*"
      }
      if (rho.test$p.value < bonfp) {
        sig <- sig + 1
        sigrho <- "*"
      }
      if (tau.test$p.value < bonfp) {
        sig <- sig + 1
        sigtau <- "*"
      }
      if (as.numeric(p) < bonfp) {
        sig <- sig + 1
        sigp <- "*"
      }
      if (mean(as.numeric(rho), as.numeric(tau), as.numeric(corr)) >
          0) {
        text.color <- "black"
        if (sig == 1)
          box.color <- 0.5
        else if (sig > 1 && sig < 5)
          box.color <- 0.75
        else if (sig == 5)
          box.color <- 1
      } else if (mean(as.numeric(rho), as.numeric(tau), as.numeric(corr)) <
                 0) {
        text.color <- "white"
        if (sig == 1)
          box.color <- 0.5
        else if (sig > 1 && sig < 5)
          box.color <- 0.25
        else if (sig == 5)
          box.color <- 0
      }
      if (!stat.pars$use.color) {
        panel.fill(col = grey(box.color), border = grey(box.color))
      } else {
        text.color <- "black"
        if (as.numeric(corr) > 0) {
          panel.fill(col = hsv(h = 0.5, s = abs(as.numeric(corr)),
                               v = 1), border = hsv(h = 0.5, s = abs(as.numeric(corr)),
                                                    v = 1))
        } else {
          panel.fill(col = hsv(h = 0, s = abs(as.numeric(corr)),
                               v = 1), border = hsv(h = 0, s = abs(as.numeric(corr)),
                                                    v = 1))
        }
      }
      if (!is.na(stat.pars$verbose)) {
        if (stat.pars$verbose == TRUE) {
          # grid.text(bquote(rho == .(rho) * .(sigrho)),
          #          x = 0.5, y = 0.9, just = stat.pars$just,
          #       gp = gpar(fontsize = stat.pars$fontsize,
          #                col = text.color))
          # grid.text(bquote(tau == .(tau) * .(sigtau)),
          #          x = 0.5, y = 0.7, just = stat.pars$just,
          #         gp = gpar(fontsize = stat.pars$fontsize,
          #                  col = text.color))
          grid.text(paste("r=", corr, sigcor, sep = ""),
                    x = 0.5, y = 0.5, just = stat.pars$just,
                    gp = gpar(fontsize = stat.pars$fontsize,
                              col = text.color))
          grid.text(paste("p=", p, sigp, sep = ""), x = 0.5,
                    y = 0.3, just = stat.pars$just, gp = gpar(fontsize = stat.pars$fontsize,
                                                              col = text.color))
          if (missing > 0)
            grid.text(paste(missing, stat.pars$missing),
                      x = 0.5, y = 0.1, just = stat.pars$just,
                      gp = gpar(fontsize = stat.pars$fontsize,
                                col = "red"))
        } else {

          grid.text(paste(corr, sigcor, sep = ""), x = 0.5,
                    y = 0.7, just = stat.pars$just, gp = gpar(fontsize = stat.pars$fontsize,
                                                              col = text.color))
          if (missing > 0)
            grid.text(paste(missing, "missing"), x = 0.5,
                      y = 0.3, just = stat.pars$just, gp = gpar(fontsize = stat.pars$fontsize,
                                                                col = text.color))
        }
      }
    }
    popViewport(1)
  }
  mosaic.panel <- function(x, y, mosaic.pars, axis.pars, xpos,
                           ypos) {
    if (!is.null(xpos) & !is.null(ypos)) {
      strucplot(table(y, x), margins = c(0, 0, 0, 0), newpage = FALSE,
                pop = FALSE, keep_aspect_ratio = FALSE, shade = mosaic.pars$shade,
                legend = FALSE, gp = mosaic.pars$gp, gp_args = mosaic.pars$gp_args,
                labeling_args = list(tl_labels = c(xpos, !ypos),
                                     gp_labels = mosaic.pars$gp_labels, varnames = c(FALSE,
                                                                                     FALSE), rot_labels = c(outer.rot, outer.rot)))
    }  else {
      if (is.null(xpos) & is.null(ypos)) {
        strucplot(table(y, x), margins = c(0, 0, 0, 0),
                  shade = mosaic.pars$shade, legend = FALSE,
                  gp = mosaic.pars$gp, gp_args = mosaic.pars$gp_args,
                  newpage = FALSE, pop = FALSE, keep_aspect_ratio = FALSE,
                  labeling = NULL)
      }  else {
        if (is.null(xpos)) {
          strucplot(table(y, x), margins = c(0, 0, 0,
                                             0), newpage = FALSE, pop = FALSE, keep_aspect_ratio = FALSE,
                    shade = mosaic.pars$shade, legend = FALSE,
                    gp = mosaic.pars$gp, gp_args = mosaic.pars$gp_args,
                    labeling_args = list(labels = c(TRUE, FALSE),
                                         tl_labels = c(ypos, FALSE), gp_labels = mosaic.pars$gp_labels,
                                         varnames = c(FALSE, FALSE), rot_labels = c(outer.rot,
                                                                                    outer.rot)))
        } else {
          strucplot(table(y, x), margins = c(0, 0, 0,
                                             0), newpage = FALSE, pop = FALSE, keep_aspect_ratio = FALSE,
                    shade = mosaic.pars$shade, legend = FALSE,
                    gp = mosaic.pars$gp, gp_args = mosaic.pars$gp_args,
                    labeling_args = list(labels = c(FALSE, TRUE),
                                         tl_labels = c(FALSE, !xpos), gp_labels = mosaic.pars$gp_labels,
                                         varnames = c(FALSE, FALSE), rot_labels = c(outer.rot,
                                                                                    outer.rot)))
        }
      }
    }
  }
  boxplot.panel <- function(x, y, type, axis.pars, xpos, ypos,
                            xylim) {
    xlim <- NULL
    ylim <- NULL
    old.color <- trellis.par.get("box.rectangle")$col
    trellis.par.set(name = "box.rectangle", value = list(col = "black"))
    trellis.par.set(name = "box.umbrella", value = list(col = "black"))
    trellis.par.set(name = "box.dot", value = list(col = "black"))
    trellis.par.set(name = "plot.symbol", value = list(col = "black"))
    if (is.factor(x)) {
      cat.labels <- levels(x)
      k <- length(levels(x))
      cat.var <- as.numeric(x)
      cont.var <- y
      horiz <- FALSE
    } else {
      cat.labels <- levels(y)
      k <- length(levels(y))
      cat.labels <- cat.labels[k:1]
      cat.var <- k + 1 - as.numeric(y)
      cont.var <- x
      horiz <- TRUE
    }
    if (horiz) {
      if (is.null(xylim)) {
        xlim <- range(cont.var, na.rm = TRUE) + c(-buffer *
                                                    (max(cont.var, na.rm = TRUE) - min(cont.var,
                                                                                       na.rm = TRUE)), buffer * (max(cont.var, na.rm = TRUE) -
                                                                                                                   min(cont.var, na.rm = TRUE)))
      } else {
        xlim <- xylim
      }
      pushViewport(viewport(xscale = xlim, yscale = c(0.5,
                                                      max(cat.var, na.rm = TRUE) + 0.5)))
      if (is.null(ypos))
        cat.labels <- NULL
      draw.axis(cont.var, cat.var, axis.pars, xpos, ypos,
                cat.labels, horiz, xlim, ylim)
      popViewport(1)
      pushViewport(viewport(xscale = xlim, yscale = c(0.5,
                                                      max(cat.var, na.rm = TRUE) + 0.5), clip = TRUE))
      if (type == "boxplot")
        panel.bwplot(cont.var, cat.var, horizontal = horiz,
                     col = "black", pch = "|", gp = gpar(box.umbrella = list(col = "black")))
      if (type == "stripplot")
        panel.stripplot(cont.var, cat.var, horizontal = horiz,
                        jitter.data = stripplot.pars$jitter, col = stripplot.pars$col,
                        cex = stripplot.pars$size, pch = stripplot.pars$pch)
    }else {
      if (is.null(xylim)) {
        ylim <- range(cont.var, na.rm = TRUE) + c(-buffer *
                                                    (max(cont.var, na.rm = TRUE) - min(cont.var,
                                                                                       na.rm = TRUE)), buffer * (max(cont.var, na.rm = TRUE) -
                                                                                                                   min(cont.var, na.rm = TRUE)))
      }else {
        ylim <- xylim
      }
      pushViewport(viewport(yscale = ylim, xscale = c(0.5,
                                                      max(cat.var, na.rm = TRUE) + 0.5)))
      if (is.null(xpos))
        cat.labels <- NULL
      draw.axis(cat.var, cont.var, axis.pars, xpos, ypos,
                cat.labels, horiz, xlim, ylim)
      popViewport(1)
      pushViewport(viewport(yscale = ylim, xscale = c(0.5,
                                                      max(cat.var, na.rm = TRUE) + 0.5), clip = TRUE))
      if (type == "boxplot")
        panel.bwplot(cat.var, cont.var, horizontal = horiz,
                     col = "black", pch = "|", gp = gpar(box.umbrella = list(col = "black")))
      if (type == "stripplot")
        panel.stripplot(cat.var, cont.var, horizontal = horiz,
                        jitter.data = stripplot.pars$jitter, col = stripplot.pars$col,
                        cex = stripplot.pars$size, pch = stripplot.pars$pch)
    }
    grid.rect(gp = gpar(fill = NULL))
    popViewport(1)
    trellis.par.set(name = "box.rectangle", value = list(col = old.color))
    trellis.par.set(name = "box.umbrella", value = list(col = old.color))
    trellis.par.set(name = "box.dot", value = list(col = old.color))
    trellis.par.set(name = "plot.symbol", value = list(col = old.color))
  }
  diag.panel <- function(x, varname, diag.pars, axis.pars,
                         xpos, ypos, xylim) {
    x <- x[!is.na(x)]
    if (is.null(xylim)) {
      xlim <- range(as.numeric(x), na.rm = TRUE) + c(-buffer *
                                                       (max(as.numeric(x), na.rm = TRUE) - min(as.numeric(x),
                                                                                               na.rm = TRUE)), buffer * (max(as.numeric(x),
                                                                                                                             na.rm = TRUE) - min(as.numeric(x), na.rm = TRUE)))
    }else {
      xlim <- xylim
    }
    ylim <- xlim
    pushViewport(viewport(xscale = xlim, yscale = ylim))
    draw.axis(as.numeric(x), as.numeric(x), axis.pars, xpos,
              ypos, NULL, NULL, xlim, ylim)
    popViewport(1)
    pushViewport(viewport(xscale = xlim, yscale = ylim, clip = TRUE))
    if (!diag.pars$show.hist) {
      grid.rect()
      grid.text(varname, 0.5, 0.5, gp = gpar(fontsize = diag.pars$fontsize,
                                             fontface = 2))
    }
    popViewport(1)
    if (diag.pars$show.hist) {
      if (!is.factor(x)) {
        pushViewport(viewport(xscale = xlim, yscale = c(0,
                                                        100), clip = TRUE))
        panel.histogram(as.numeric(x), breaks = NULL,
                        type = "percent", col = diag.pars$hist.color)
      }else {
        pushViewport(viewport(xscale = c(min(as.numeric(x),
                                             na.rm = TRUE) - 1, max(as.numeric(x), na.rm = TRUE) +
                                           1), yscale = c(0, 100), clip = TRUE))
        panel.barchart(1:length(table(x)), 100 * table(x)/sum(table(x)),
                       horizontal = FALSE, col = diag.pars$hist.color)
      }
      grid.text(varname, 0.5, 0.85, gp = gpar(fontsize = diag.pars$fontsize))
      popViewport(1)
    }
  }
  grid.newpage()
  N <- ncol(x)
  vp.main <- viewport(x = outer.margins$bottom, y = outer.margins$left,
                      width = unit(1, "npc") - outer.margins$right - outer.margins$left,
                      height = unit(1, "npc") - outer.margins$top - outer.margins$bottom,
                      just = c("left", "bottom"), name = "main", clip = "off")
  pushViewport(vp.main)
  for (i in 1:N) {
    for (j in 1:N) {
      if (diagonal == "default")
        labelj <- j
      else labelj <- N - j + 1
      x[is.infinite(x[, i]), i] <- NA
      x[is.infinite(x[, j]), j] <- NA
      vp <- viewport(x = (labelj - 1)/N, y = 1 - i/N, width = 1/N,
                     height = 1/N, just = c("left", "bottom"), name = as.character(i *
                                                                                     N + j))
      pushViewport(vp)
      vp.in <- viewport(x = 0.5, y = 0.5, width = 1 - gap,
                        height = 1 - gap, just = c("center", "center"),
                        name = paste("IN", as.character(i * N + j)))
      pushViewport(vp.in)
      xpos <- NULL
      if (i == 1 && outer.labels$top[j]) {
        xpos <- FALSE
      }
      if (i == N && outer.labels$bottom[j]) {
        xpos <- TRUE
      }
      ypos <- NULL
      if (j == N && outer.labels$right[i]) {
        ypos <- FALSE
      }
      if (j == 1 && outer.labels$left[i]) {
        ypos <- TRUE
      }
      if (!is.null(ypos) & diagonal != "default") {
        ypos <- !ypos
      }
      if (i == j) {
        diag.panel(x[, i], names(x)[i], diag.pars, axis.pars,
                   xpos, ypos, xylim)
      }else {
        if (is.factor(x[, i]) + is.factor(x[, j]) ==
            1) {
          if (i < j & upper.pars$conditional != "barcode")
            boxplot.panel(x[, j], x[, i], upper.pars$conditional,
                          axis.pars, xpos, ypos, xylim)
          if (i > j & lower.pars$conditional != "barcode")
            boxplot.panel(x[, j], x[, i], lower.pars$conditional,
                          axis.pars, xpos, ypos, xylim)
          if (i < j & upper.pars$conditional == "barcode") {
            if (is.factor(x[, i])) {
              barcode(split(x[, j], x[, i])[length(levels(x[,
                                                            i])):1], horizontal = TRUE, xlim = xylim,
                      labelloc = ypos, axisloc = xpos, labelouter = TRUE,
                      newpage = FALSE, fontsize = axis.pars$fontsize,
                      buffer = buffer, nint = barcode.pars$nint,
                      ptsize = barcode.pars$ptsize, ptpch = barcode.pars$ptpch,
                      bcspace = barcode.pars$bcspace, use.points = barcode.pars$use.points)
            } else {
              if (!is.null(ypos))
                ypos <- !ypos
              barcode(split(x[, i], x[, j])[length(levels(x[,
                                                            j])):1], horizontal = FALSE, xlim = xylim,
                      labelloc = xpos, axisloc = ypos, labelouter = TRUE,
                      newpage = FALSE, fontsize = axis.pars$fontsize,
                      buffer = buffer, nint = barcode.pars$nint,
                      ptsize = barcode.pars$ptsize, ptpch = barcode.pars$ptpch,
                      bcspace = barcode.pars$bcspace, use.points = barcode.pars$use.points)
            }
          }
          if (i > j & lower.pars$conditional == "barcode") {
            if (is.factor(x[, i])) {
              barcode(split(x[, j], x[, i])[length(levels(x[,
                                                            i])):1], horizontal = TRUE, xlim = xylim,
                      labelloc = ypos, axisloc = xpos, labelouter = TRUE,
                      newpage = FALSE, fontsize = axis.pars$fontsize,
                      buffer = buffer, nint = barcode.pars$nint,
                      ptsize = barcode.pars$ptsize, ptpch = barcode.pars$ptpch,
                      bcspace = barcode.pars$bcspace, use.points = barcode.pars$use.points)
            } else {
              if (!is.null(ypos))
                ypos <- !ypos
              barcode(split(x[, i], x[, j])[length(levels(x[,
                                                            j])):1], horizontal = FALSE, xlim = xylim,
                      labelloc = xpos, axisloc = ypos, labelouter = TRUE,
                      newpage = FALSE, fontsize = axis.pars$fontsize,
                      buffer = buffer, nint = barcode.pars$nint,
                      ptsize = barcode.pars$ptsize, ptpch = barcode.pars$ptpch,
                      bcspace = barcode.pars$bcspace, use.points = barcode.pars$use.points)
            }
          }
        }
        if (is.factor(x[, i]) + is.factor(x[, j]) ==
            0) {
          if (i < j)
            type <- upper.pars$scatter
          else type <- lower.pars$scatter
          scatterplot.panel(x[, j], x[, i], type, scatter.pars,
                            axis.pars, xpos, ypos, xylim)
        }
        if (is.factor(x[, i]) + is.factor(x[, j]) ==
            2) {
          if (i < j)
            mosaic.panel(x[, j], x[, i], mosaic.pars,
                         axis.pars, xpos, ypos)
          else mosaic.panel(x[, j], x[, i], mosaic.pars,
                            axis.pars, xpos, ypos)
        }
      }
      popViewport(1)
      upViewport()
    }
  }
  popViewport()
  if (whatis)
    whatis(x)
}

gpar <- function(...) {
  gp <- validGP(list(...))
  class(gp) <- "gpar"
  gp
}

is.gpar <- function(x) {
  inherits(x, "gpar")
}

print.gpar <- function(x, ...) {
  print(unclass(x), ...)
  invisible(x)
}

validGP <- function(gpars) {
  # Check a (non-NULL) gpar is not of length 0
  check.length <- function(gparname) {
    if (length(gpars[[gparname]]) == 0)
      stop(gettextf("'gpar' element '%s' must not be length 0", gparname),
           domain = NA)
  }
  # Check a gpar is numeric and not NULL
  numnotnull <- function(gparname) {
    if (!is.na(match(gparname, names(gpars)))) {
      if (is.null(gpars[[gparname]]))
        gpars[[gparname]] <<- NULL
      else {
        check.length(gparname)
        gpars[[gparname]] <<- as.numeric(gpars[[gparname]])
      }
    }
  }
  # fontsize, lineheight, cex, lwd should be numeric and not NULL
  numnotnull("fontsize")
  numnotnull("lineheight")
  numnotnull("cex")
  numnotnull("lwd")
  numnotnull("lex")
  # gamma defunct in 2.7.0
  if ("gamma" %in% names(gpars)) {
    warning("'gamma' 'gpar' element is defunct")
    gpars$gamma <- NULL
  }
  numnotnull("alpha")
  # col and fill are converted in C code
  # BUT still want to check length > 0
  if (!is.na(match("col", names(gpars)))) {
    if (is.null(gpars$col))
      gpars$col <- NULL
    else
      check.length("col")
  }
  if (!is.na(match("fill", names(gpars)))) {
    if (is.null(gpars$fill))
      gpars$fill <- NULL
    else
      check.length("fill")
  }
  # lty converted in C code
  # BUT still want to check for NULL and check length > 0
  if (!is.na(match("lty", names(gpars)))) {
    if (is.null(gpars$lty))
      gpars$lty <- NULL
    else
      check.length("lty")
  }
  if (!is.na(match("lineend", names(gpars)))) {
    if (is.null(gpars$lineend))
      gpars$lineend <- NULL
    else
      check.length("lineend")
  }
  if (!is.na(match("linejoin", names(gpars)))) {
    if (is.null(gpars$linejoin))
      gpars$linejoin <- NULL
    else
      check.length("linejoin")
  }
  # linemitre should be larger than 1
  numnotnull("linemitre")
  if (!is.na(match("linemitre", names(gpars)))) {
    if (any(gpars$linemitre < 1))
      stop("invalid 'linemitre' value")
  }
  # alpha should be 0 to 1
  if (!is.na(match("alpha", names(gpars)))) {
    if (any(gpars$alpha < 0 || gpars$alpha > 1))
      stop("invalid 'alpha' value")
  }
  # font should be integer and not NULL
  if (!is.na(match("font", names(gpars)))) {
    if (is.null(gpars$font))
      gpars$font <- NULL
    else {
      check.length("font")
      gpars$font <- as.integer(gpars$font)
    }
  }
  # fontfamily should be character
  if (!is.na(match("fontfamily", names(gpars)))) {
    if (is.null(gpars$fontfamily))
      gpars$fontfamily <- NULL
    else {
      check.length("fontfamily")
      gpars$fontfamily <- as.character(gpars$fontfamily)
    }
  }
  # fontface can be character or integer;  map character to integer
  # store value in font
  # Illegal to specify both font and fontface
  if (!is.na(match("fontface", names(gpars)))) {
    if (!is.na(match("font", names(gpars))))
      stop("must specify only one of 'font' and 'fontface'")
    gpars$font <-
      if (is.null(gpars$fontface)) NULL # remove it
    else {
      check.length("fontface")
      if (is.numeric(gpars$fontface))
        as.integer(gpars$fontface)
      else
        vapply(as.character(gpars$fontface),
               function(ch) # returns integer
                 switch(ch,
                        plain = 1L,
                        bold  = 2L,
                        italic=, oblique = 3L,
                        bold.italic = 4L,
                        symbol= 5L,
                        # These are Hershey variants
                        cyrillic=5L,
                        cyrillic.oblique=6L,
                        EUC   = 7L,
                        stop("invalid fontface ", ch)), 0L)
    }
  }
  gpars
}

# Method for subsetting "gpar" objects
`[.gpar` <- function(x, index, ...) {
  if (length(x) == 0)
    return(gpar())
  maxn <- do.call("max", lapply(x, length))
  newgp <- lapply(x, rep, length.out=maxn)
  newgp <- lapply(X = newgp, FUN = "[", index, ...)
  class(newgp) <- "gpar"
  newgp
}

# possible gpar names
# The order must match the GP_* values in grid.h
.grid.gpar.names <- c("fill", "col", "gamma", "lty", "lwd", "cex",
                      "fontsize", "lineheight", "font", "fontfamily",
                      "alpha", "lineend", "linejoin", "linemitre",
                      "lex",
                      # Keep fontface at the end because it is never
                      # used in C code (it gets mapped to font)
                      "fontface")

set.gpar <- function(gp) {
  if (!is.gpar(gp))
    stop("argument must be a 'gpar' object")
  temp <- grid.Call(L_getGPar)
  # gamma defunct in 2.7.0
  if ("gamma" %in% names(gp)) {
    warning("'gamma' 'gpar' element is defunct")
    gp$gamma <- NULL
  }
  # Special case "cex" (make it cumulative)
  if (match("cex", names(gp), nomatch=0L))
    tempcex <- temp$cex * gp$cex
  else
    tempcex <- temp$cex
  # Special case "alpha" (make it cumulative)
  if (match("alpha", names(gp), nomatch=0L))
    tempalpha <- temp$alpha * gp$alpha
  else
    tempalpha <- temp$alpha
  # Special case "lex" (make it cumulative)
  if (match("lex", names(gp), nomatch=0L))
    templex <- temp$lex * gp$lex
  else
    templex <- temp$lex
  # All other gpars
  temp[names(gp)] <- gp
  temp$cex <- tempcex
  temp$alpha <- tempalpha
  temp$lex <- templex
  # Do this as a .Call.graphics to get it onto the base display list
  grid.Call.graphics(L_setGPar, temp)
}

get.gpar <- function(names=NULL) {
  if (is.null(names)) {
    result <- grid.Call(L_getGPar)
    # drop gamma
    result$gamma <- NULL
  } else {
    if (!is.character(names) ||
        !all(names %in% .grid.gpar.names))
      stop("must specify only valid 'gpar' names")
    # gamma deprecated
    if ("gamma" %in% names) {
      warning("'gamma' 'gpar' element is defunct")
      names <- names[-match("gamma", names)]
    }
    result <- unclass(grid.Call(L_getGPar))[names]
  }
  class(result) <- "gpar"
  result
}

# When editing a gp slot, only update the specified gpars
# Assume gp is NULL or a gpar
# assume newgp is a gpar (and not NULL)
mod.gpar <- function(gp, newgp) {
  if (is.null(gp))
    gp <- newgp
  else
    gp[names(newgp)] <- newgp
  gp
}
