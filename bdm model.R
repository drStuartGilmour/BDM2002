# program to make datasets for DiD modeling, where an observational unit is only in tx or only in control
data.genBl<-function(txUnits,txPeriod,dataRep,baseMean,varLim,txEff,spillEff,rndError,plSlope){
  
  mean.mod<-varLim
  n.treat<-txUnits
  n.unit<-2*n.treat
  n.time<-txPeriod
  n.rep<-dataRep
  n.period<-n.time/2
  slope.val<-plSlope
  # now let's set up teh baseline differences
  # unit 1 has a value of 5; every additional unit varies with a uniform distribution around this
  # random error set here
  # have a treatment effect of treat.eff, and a spillover effect of spill.eff
  base.alpha<-baseMean
  
  mean.diff<-0
  base.err<-rndError
  time.vec<-1:n.time
  treat.eff<-txEff
  spill.eff<-spillEff
  # for complete regularity, let's set the treatment group to be even and the intervention group to be odd
  
  treat.base<-rep(c(rep(0,times=n.time),rep(1,times=n.time)),times=n.treat)
  period.base<-rep(c(rep(0,times=n.period),rep(1,times=n.period)),times=n.unit)
  time.base<-rep(time.vec,times=n.unit)
  unit.id<-rep(1,times=n.time)
  base.mean<-rep(base.alpha,times=n.time)
  for (i in 2:n.unit){
    unit.id<-c(unit.id,rep(i,times=n.time))
    if (mean.mod>0){
      base.mean<-c(base.mean,base.alpha+rep(runif(1,min=-mean.mod,max=mean.mod),times=n.time))
    } else {
      base.mean<-c(base.mean,rep(base.alpha,times=n.time))
    }
  }
  total.eff<-base.mean+treat.eff*treat.base*period.base+spill.eff*(1-treat.base)*period.base+treat.base*mean.diff+plSlope*time.base
  dvar.vec<-treat.base*period.base
  ind.id<-rep(1,times=n.time*n.treat)
  base.mat<-cbind(unit.id,ind.id,time.base,treat.base,period.base,dvar.vec,total.eff)
  up.mat<-base.mat
  if (n.rep>1) {  
    for (j in 2:n.rep){
      up.mat[,2]<-rep(j,times=n.time*n.treat)
      base.mat<-rbind(base.mat,up.mat)
    }
  }
  r.vec<-rnorm(n.rep*n.time*n.unit,mean=0,sd=base.err)
  base.mat[,7]<-base.mat[,7]+r.vec
  
  temp.mat<-data.frame(base.mat)
  names(temp.mat)<-c("unit","individual","time","treatment","period","twfeTerm","response")
  temp.mat<-temp.mat[with(temp.mat, order(unit,individual,time)), ]
  return(temp.mat)
}

# test generate with 4 reps and 4 time points, no variation between units
data.t<-data.genBl(10,4,4,24,0,0,0,1,0)

# test generate with 1 rep and 2 time points, baseline difference in units given by unif(-10,10)
data.t<-data.genBl(10,2,1,24,10,0,0,1,0)

# test generate with 300 reps, 10 treatment units and 2 time points
data.t<-data.genBl(10,2,300,24,10,0,0,1,0)


# test models

lm.true<-lm(response~as.factor(period)+as.factor(treatment)+as.factor(treatment):as.factor(period),data=data.t)
lm.twfe<-lm(response~as.factor(unit)+as.factor(time)+as.factor(twfeTerm),data=data.t)

summary(lm.true)
summary(lm.twfe)

# test generate with 300 reps, 25 treatment units and 20 time points
data.t<-data.genBl(25,20,300,24,10,0,0,1,0)


# test models

lm.true<-lm(response~as.factor(period)+as.factor(treatment)+as.factor(treatment):as.factor(period),data=data.t)
lm.twfe<-lm(response~as.factor(unit)+as.factor(time)+as.factor(twfeTerm),data=data.t)

summary(lm.true)
summary(lm.twfe)

#######################
#
# Simulation 1: Effect of underlying baseline structure
#
####################

# in this simulation we model increasing degrees of separation between baseline units by keeping time and number of units
# and observations fixed, but slowly increasing hte limits of the uniform distribution that determines baseline differences
# between treatment units

xv.vec<-seq(from=0,to=20,by=0.1)
step.num<-length(xv.vec)
sim1.mat<-data.frame(matrix(NA,nrow=step.num,ncol=7))
names(sim1.mat)<-c("baseVar","lmErr","lmErrIQR_l","lmErrIQR_u","twfeErr","twfeErrIQR_l","twfeErrIQR_u")
sim.num<-100
base.val<-20
tx.num<-25
time.num<-20
rep.num<-10

for (i in 1:step.num){
  lmErr.vec<-rep(0,times=sim.num)
  twfeErr.vec<-rep(0,times=sim.num)
  for (j in 1:sim.num){
    data.t<-data.genBl(tx.num,time.num,rep.num,base.val,xv.vec[i],0,0,1,0)
    lm.true<-lm(response~as.factor(period)+as.factor(treatment)+as.factor(treatment):as.factor(period),data=data.t)
    lm.twfe<-lm(response~as.factor(twfeTerm)+as.factor(unit)+as.factor(time),data=data.t)
    lmErr.vec[j]<-summary(lm.true)$coefficients[4,2]
    twfeErr.vec[j]<-summary(lm.twfe)$coefficients[2,2]
  }

  sim1.mat[i,1]<-xv.vec[i]
  sim1.mat[i,2]<-mean(lmErr.vec)
  sim1.mat[i,3]<-quantile(lmErr.vec,0.25)
  sim1.mat[i,4]<-quantile(lmErr.vec,0.75)
  sim1.mat[i,5]<-mean(twfeErr.vec)
  sim1.mat[i,6]<-quantile(twfeErr.vec,0.25)
  sim1.mat[i,7]<-quantile(twfeErr.vec,0.75)
}
#####################
#
# Simulation 2: checking rejection rates
##
#####################

### we can simulate this by varying: the uniform distribution; the time points; the reps

## start with 1 rep and go up, with no uniform difference in treatment units
sim.num<-100
unit.num<-25
studyLength<-20
tx.var<-0
rep.prog<-c(1,2,4,8,16,32,64,128,256,512)

out.mat<-data.frame(matrix(NA,nrow=length(rep.prog),ncol=4))
names(out.mat)<-c("TreatmentDifference","Observations","LMSig","TWFESig")

for (i in 1:length(rep.prog)){
  lm.sig<-rep(0,times=sim.num)
  twfe.sig<-rep(0,times=sim.num)
  for (j in 1:sim.num){
    data.t<-data.genBl(unit.num,studyLength,rep.prog[i],24,tx.var,0,0,1,0)
    lm.true<-lm(response~as.factor(period)+as.factor(treatment)+as.factor(treatment):as.factor(period),data=data.t)
    lm.twfe<-lm(response~as.factor(twfeTerm)+as.factor(unit)+as.factor(time),data=data.t)
    lm.sig[j]<-summary(lm.true)$coefficients[4,4]<0.05
    twfe.sig[j]<-summary(lm.twfe)$coefficients[2,4]<0.05
  }
  out.mat[i,1]<-tx.var
  out.mat[i,2]<-rep.prog[i]
  out.mat[i,3]<-sum(lm.sig)/sim.num
  out.mat[i,4]<-sum(twfe.sig)/sim.num
}
#exactly the same ...
# now let's build some alternative models, with structural differences between the treatment units at baseline. 
# The CPS data had differences of up to 10% at baseline between states. With 50 states and a baseline of 24, that would
# mean a uniform distribution between -5 and 5 (10/50=0.2) would be sufficient
# so let's expand the treatment value to 10 in steps of 2

for (i in c(2,4,6,8,10)){
  tx.var<-i
  out.matT<-data.frame(matrix(NA,nrow=length(rep.prog),ncol=4))
  names(out.matT)<-c("TreatmentDifference","Observations","LMSig","TWFESig")
  
  for (i in 1:length(rep.prog)){
    lm.sig<-rep(0,times=sim.num)
    twfe.sig<-rep(0,times=sim.num)
    for (j in 1:sim.num){
      data.t<-data.genBl(unit.num,studyLength,rep.prog[i],24,tx.var,0,0,1,0)
      lm.true<-lm(response~as.factor(period)+as.factor(treatment)+as.factor(treatment):as.factor(period),data=data.t)
      lm.twfe<-lm(response~as.factor(twfeTerm)+as.factor(unit)+as.factor(time),data=data.t)
      lm.sig[j]<-summary(lm.true)$coefficients[4,4]<0.05
      twfe.sig[j]<-summary(lm.twfe)$coefficients[2,4]<0.05
    }
    out.matT[i,1]<-tx.var
    out.matT[i,2]<-rep.prog[i]
    out.matT[i,3]<-sum(lm.sig)/sim.num
    out.matT[i,4]<-sum(twfe.sig)/sim.num
  }  
  out.mat<-rbind(out.mat,out.matT)
}
names(out.mat)<-c("TreatmentDifference","Observations","LMSig","TWFESig")

## now check to see whether a slope makes a difference

data.t<-data.genBl(50,20,10,24,10,0,0,1,2)
lm.true<-lm(response~as.factor(period)+as.factor(treatment)+as.factor(treatment):as.factor(period),data=data.t)
lm.twfe<-lm(response~as.factor(twfeTerm)+as.factor(unit)+as.factor(time),data=data.t)

summary(lm.true)
summary(lm.twfe)

# nah! So there is nothing going on here!

### plot model results for the twfe
plot(out.mat$Observations[out.mat$TreatmentDifference==0],out.mat$TWFESig[out.mat$TreatmentDifference==0])
lines(out.mat$Observations[out.mat$TreatmentDifference==4],out.mat$TWFESig[out.mat$TreatmentDifference==4])
lines(out.mat$Observations[out.mat$TreatmentDifference==8],out.mat$TWFESig[out.mat$TreatmentDifference==8])

#
# Simulation 3: a brief assessment of power
##
#####################

# now we examine the power of the two models. Let's run through simulations with 1, 2, 4 or 16 individuals per observational unit
# and test them with a treatment effect that is bigger than the control effect by 0.5, 1 and 2. We will set teh spillover 
# effect to 0 for simplicity
sim.num<-100
unit.num<-25
studyLength<-20

rep.prog<-c(1,2,4,8,16)
rep.tx<-c(0.5,1,2)
sim.length<-length(rep.prog)*length(rep.tx)
sim3.mat<-data.frame(matrix(NA,nrow=length(sim.length),ncol=6))
names(sim3.mat)<-c("Intervention effect","Observations","LMBias","LMSig","TWFEBias","TWFESig")
step.set<-0
for (i in 1:length(rep.prog)){
  for (j in 1:length(rep.tx)){
    lm.sig<-rep(0,times=sim.num)
    twfe.sig<-rep(0,times=sim.num)
    lm.est<-rep(0,times=sim.num)
    twfe.est<-rep(0,times=sim.num)
    step.set<-step.set+1
    tx.var<-rep.tx[j]
    rep.num<-rep.prog[i]
    for (k in 1:sim.num){
      data.t<-data.genBl(unit.num,studyLength,rep.num,24,10,tx.var,0,1,0)
      lm.true<-lm(response~as.factor(period)+as.factor(treatment)+as.factor(treatment):as.factor(period),data=data.t)
      lm.twfe<-lm(response~as.factor(twfeTerm)+as.factor(unit)+as.factor(time),data=data.t)
      lm.est[k]<-summary(lm.true)$coefficients[4,1]
      lm.sig[k]<-summary(lm.true)$coefficients[4,4]<0.05
      twfe.sig[k]<-summary(lm.twfe)$coefficients[2,4]<0.05
      twfe.est[k]<-summary(lm.twfe)$coefficients[2,1]
    }
  
  sim3.mat[step.set,1]<-tx.var
  sim3.mat[step.set,2]<-rep.num
  sim3.mat[step.set,3]<-mean(lm.est)-tx.var
  sim3.mat[step.set,4]<-sum(lm.sig)/sim.num
  sim3.mat[step.set,5]<-mean(twfe.est)-tx.var
  sim3.mat[step.set,6]<-sum(twfe.sig)/sim.num
  }
}
