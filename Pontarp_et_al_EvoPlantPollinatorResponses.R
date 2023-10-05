# This code is made by Lingzi Wang at Lund University.

# pamameters
parms<- list(
  rr<-c(gr<-1, gb<-1),  # growth rate of plant and pollinator populations
  rm<-c(rm1<-0, rm2<-0), # optimal traits of plants
  bm<-c(bm1<-0, bm2<-0), # optimal traits of pollinators
  resource_var<-c(rm_sig<-1, bm_sig<-1), # environmental niche widths 
  cp<-c(cp1<-1,cp2<-0.5), # competition coefficient
  KR<-c(KR1<-3,KR2<-3), # plant carrying capacity
  KB<-c(KB1<-1,KB2<-1), # pollinator carry capacity
  gamma<-c(gamma1<-1,gamma2<-3), # interactive trait values
  inter_trait<-c(r_sig<-1, b_sig<-1), # interactive niche widths
  cR<-c(cRH<-0.1,cRL<-0.05),  # convergent coefficient from pollinators to plants
  cB<-c(cBH<-0.1,cBL<-0.05), # convergent coefficient from plants to pollinators
  mu1<-0.5, # mutation rate
  mu2<-1, # mutation rate
  mu_sig<-1, # mutation range
  mu_step<-0.1 # mutation steps
)

# the interaction between the environment and plant and pollinator populations
ra_fun<-function(r,b,R_pop,B_pop){
  pop<-R_pop+B_pop
  ra<-rep(NA,pop)
  for (j in 1:R_pop) {
    ra[j]<-gr*exp(-(1/2)*((r[j]-rm[j])/rm_sig)^2)
  }
  for (j in 1:B_pop) {
    ra[R_pop+j]<-gb*exp(-(1/2)*((b[j]-bm[j])/bm_sig)^2)
  }
  return(ra)
}

# the interaction inside and between plants and pollinators
interaction<-function(r,b,beta,R_pop,B_pop){
  pop<-R_pop+B_pop
  inter<-matrix(rep(NA,pop*pop),nrow=pop)
  KR[1]<-KR1
  KR[2]<-KR2
  KB[1]<-KB1
  KB[2]<-KB2
  for (j in 1:R_pop) {
    for (k in 1:R_pop) {if (k==j) {inter[j,k]<--cp1/KR[j]} 
      else {inter[j,k]<--cp2/KR[j]}}
    for (k in 1:B_pop) {if (k==j) {inter[j,k+R_pop]<-cRH*exp(-(1/2)*((gamma[j]-beta[k])/r_sig)^2)} # mutualist interaction
      else {inter[j,k+R_pop]<-cRL*exp(-(1/2)*((gamma[j]-beta[k])/r_sig)^2)}}
  }
  for (j in 1:B_pop) {
    for (k in 1:R_pop) {if (k==j) {inter[j+R_pop,k]<-cBH*exp(-(1/2)*((beta[j]-gamma[k])/r_sig)^2)}  # mutualist interaction
      else {inter[j+R_pop,k]<-cBL*exp(-(1/2)*((beta[j]-gamma[k])/r_sig)^2)}}
    for (k in 1:B_pop) if (k==j) {inter[j+R_pop,k+R_pop]<--cp1/KB[j]}   
    else   {inter[j+R_pop,k+R_pop]<--cp2/KB[j]}
  }
  return(inter)
}

# Jacobian matrix
Jacob<-function (r,b,R,B) {
  dR1<-expression(R1*(gr*exp(-(1/2)*((r1-rm1)/rm_sig)^2)-(cp1*R1+cp2*R2)/KR1+cRH*exp(-(1/2)*((gamma1-beta1)/r_sig)^2)*B1+cRL*exp(-(1/2)*((gamma1-beta2)/r_sig)^2)*B2))
  dR2<-expression(R2*(gr*exp(-(1/2)*((r2-rm2)/rm_sig)^2)-(cp2*R1+cp1*R2)/KR2+cRL*exp(-(1/2)*((gamma2-beta1)/r_sig)^2)*B1+cRH*exp(-(1/2)*((gamma2-beta2)/r_sig)^2)*B2))
  dB1<-expression(B1*(gb*exp(-(1/2)*((b1-bm1)/bm_sig)^2)-(cp1*B1+cp2*B2)/KB1+cBH*exp(-(1/2)*((beta1-gamma1)/b_sig)^2)*R1+cBL*exp(-(1/2)*((beta1-gamma2)/b_sig)^2)*R2))
  dB2<-expression(B2*(gb*exp(-(1/2)*((b2-bm2)/bm_sig)^2)-(cp2*B1+cp1*B2)/KB2+cBL*exp(-(1/2)*((beta2-gamma1)/b_sig)^2)*R1+cBH*exp(-(1/2)*((beta2-gamma2)/b_sig)^2)*R2))
  
  R_pop<-sum(is.na(R)==0)
  B_pop<-sum(is.na(B)==0)
  
  R1<-R[1] 
  R2<-R[2]
  B1<-B[1]
  B2<-B[2]
  r1<-r[1]
  r2<-r[2]
  b1<-b[1]
  b2<-b[2]
  
  Jacob3<-matrix(c(eval(D(dR1,'R1')),eval(D(dR1,'R2')),eval(D(dR1,'B1')),eval(D(dR1,'B2')),
                   eval(D(dR2,'R1')),eval(D(dR2,'R2')),eval(D(dR2,'B1')),eval(D(dR2,'B2')),
                   eval(D(dB1,'R1')),eval(D(dB1,'R2')),eval(D(dB1,'B1')),eval(D(dB1,'B2')),
                   eval(D(dB2,'R1')),eval(D(dB2,'R2')),eval(D(dB2,'B1')),eval(D(dB2,'B2'))),
                 nrow=4,byrow = T)
  Jacob<-Jacob3
  eig<-eigen(Jacob)$values
  return(eig)
}

# invasion fitness
R1_invasion_fitness_for_derivative_fun<-function (r1_t,R1,R2,B1,B2){
  invasion_fitness<-expression(gr*exp(-(1/2)*((r1_t-rm1)/rm_sig)^2)-(cp1*R1+cp2*R2)/KR1+cRH*exp(-(1/2)*((gamma1-beta1_t)/r_sig)^2)*B1+cRL*exp(-(1/2)*((gamma1-beta2_t)/r_sig)^2)*B2)
  return(invasion_fitness)
}

R2_invasion_fitness_for_derivative_fun<-function (r2_t,R1,R2,B1,B2){
  invasion_fitness<-expression(gr*exp(-(1/2)*((r2_t-rm2)/rm_sig)^2)-(cp2*R1+cp1*R2)/KR2+cRL*exp(-(1/2)*((gamma2-beta1_t)/r_sig)^2)*B1+cRH*exp(-(1/2)*((gamma2-beta2_t)/r_sig)^2)*B2)
  return(invasion_fitness)
}

B1_invasion_fitness_for_derivative_fun<-function (b1_t,beta1_t,R1,R2,B1,B2){
  invasion_fitness<-expression(gb*exp(-(1/2)*((b1_t-bm1)/bm_sig)^2)-(cp1*B1+cp2*B2)/KB1+cBH*exp(-(1/2)*((beta1_t-gamma1)/b_sig)^2)*R1+cBL*exp(-(1/2)*((beta1_t-gamma2)/b_sig)^2)*R2)
  return(invasion_fitness)
}

B2_invasion_fitness_for_derivative_fun<-function (b2_t,beta2_t,R1,R2,B1,B2){
  invasion_fitness<-expression(gb*exp(-(1/2)*((b2_t-bm2)/bm_sig)^2)-(cp2*B1+cp1*B2)/KB2+cBL*exp(-(1/2)*((beta2_t-gamma1)/b_sig)^2)*R1+cBH*exp(-(1/2)*((beta2_t-gamma2)/b_sig)^2)*R2)
  return(invasion_fitness)
}

# selection greadient
R1_sele_gradient<-function(r1_t,R,B){
  fitness<-R1_invasion_fitness_for_derivative_fun(r1_t,R1,R2,B1,B2)
  sele_gradient_exp1<-D(fitness,'r1_t')
  R1<-R[1] 
  R2<-R[2]
  B1<-B[1] 
  B2<-B[2]
  sele_gradient<-NA
  sele_gradient<-eval(sele_gradient_exp1)
  return(sele_gradient)
}

R2_sele_gradient<-function(r2_t,R,B){
  fitness<-R2_invasion_fitness_for_derivative_fun(r2_t,R1,R2,B1,B2)
  sele_gradient_exp1<-D(fitness,'r2_t')
  R1<-R[1] 
  R2<-R[2]
  B1<-B[1] 
  B2<-B[2]
  sele_gradient<-NA
  sele_gradient<-eval(sele_gradient_exp1)
  return(sele_gradient)
}

B1_sele_gradient<-function(b1_t,beta1_t,R,B){
  fitness<-B1_invasion_fitness_for_derivative_fun(b1_t,beta1_t,R1,R2,B1,B2)
  sele_gradient_exp1<-D(fitness,'b1_t')
  sele_gradient_exp2<-D(fitness,'beta1_t')
  R1<-R[1] 
  R2<-R[2]
  B1<-B[1] 
  B2<-B[2]
  sele_gradient<-rep(NA,2)
  sele_gradient[1]<-eval(sele_gradient_exp1)
  sele_gradient[2]<-eval(sele_gradient_exp2)
  return(sele_gradient)
}

B2_sele_gradient<-function(b2_t,beta2_t,R,B){
  fitness<-B2_invasion_fitness_for_derivative_fun(b2_t,beta2_t,R1,R2,B1,B2)
  sele_gradient_exp1<-D(fitness,'b2_t')
  sele_gradient_exp2<-D(fitness,'beta2_t')
  R1<-R[1] 
  R2<-R[2]
  B1<-B[1] 
  B2<-B[2]
  sele_gradient<-rep(NA,2)
  sele_gradient[1]<-eval(sele_gradient_exp1)
  sele_gradient[2]<-eval(sele_gradient_exp2)
  return(sele_gradient)
}

#curvatures
R1_curvature<-function(r1_t,R,B){
  fitness<-R1_invasion_fitness_for_derivative_fun(r1_t,R1,R2,B1,B2)
  curve_exp1<-D(D(fitness,'r1_t'),'r1_t')
  R1<-R[1] 
  R2<-R[2]
  B1<-B[1] 
  B2<-B[2]
  curve<-NA
  curve<-eval(curve_exp1)
  return(curve)
}

R2_curvature<-function(r2_t,R,B){
  fitness<-R2_invasion_fitness_for_derivative_fun(r2_t,R1,R2,B1,B2)
  curve_exp1<-D(D(fitness,'r2_t'),'r2_t')
  R1<-R[1] 
  R2<-R[2]
  B1<-B[1] 
  B2<-B[2]
  curve<-NA
  curve<-eval(curve_exp1)
  return(curve)
}

B1_curvature<-function(b1_t,beta1_t,R,B){
  fitness<-B1_invasion_fitness_for_derivative_fun(b1_t,beta1_t,R1,R2,B1,B2)
  curve_exp1<-D(D(fitness,'b1_t'),'b1_t')
  curve_exp2<-D(D(fitness,'beta1_t'),'beta1_t')
  R1<-R[1] 
  R2<-R[2]
  B1<-B[1] 
  B2<-B[2]
  curve<-rep(NA,2)
  curve[1]<-eval(curve_exp1)
  curve[2]<-eval(curve_exp2)
  return(curve)
}

B2_curvature<-function(b2_t,beta2_t,R,B){
  fitness<-B2_invasion_fitness_for_derivative_fun(b2_t,beta2_t,R1,R2,B1,B2)
  curve_exp1<-D(D(fitness,'b2_t'),'b2_t')
  curve_exp2<-D(D(fitness,'beta2_t'),'beta2_t')
  R1<-R[1] 
  R2<-R[2]
  B1<-B[1] 
  B2<-B[2]
  curve<-rep(NA,2)
  curve[1]<-eval(curve_exp1)
  curve[2]<-eval(curve_exp2)
  return(curve)
}

beta_fitness_landscape<-function(b,beta,R,B){
  start<-min(beta[!is.na(beta)])-3
  end<-max(beta[!is.na(beta)])+3
  seq<-seq(start,end,by=0.0001)
  N<-length(seq)
  R1<-R[1]
  R2<-R[2]
  B1<-B[1]
  B2<-B[2]
  b1_t<-b[1]
  b2_t<-b[2]
  beta1_t<-beta[1]
  beta2_t<-beta[2]
  
  beta1_invasion_fitness<-NA
  beta2_invasion_fitness<-NA
  bb1_fit<-eval(B1_invasion_fitness_for_derivative_fun (b1_t,beta1_t,R1,R2,B1,B2))
  bb2_fit<-eval(B2_invasion_fitness_for_derivative_fun (b2_t,beta2_t,R1,R2,B1,B2))
  
  for(i in 1:N){
    beta1_t<-seq[i]
    beta1_invasion_fitness[i]<-eval(B1_invasion_fitness_for_derivative_fun (b1_t,beta1_t,R1,R2,B1,B2))
  }
  for(i in 1:N){
    beta2_t<-seq[i]
    beta2_invasion_fitness[i]<-eval(B2_invasion_fitness_for_derivative_fun (b2_t,beta2_t,R1,R2,B1,B2))
  }
  plot(seq,beta1_invasion_fitness,ylim=c(min(beta1_invasion_fitness,beta2_invasion_fitness),max(beta1_invasion_fitness,beta2_invasion_fitness)+0.01),type = "l",xlab = "trait a",ylab="fitness")
  points(beta[1],bb1_fit,col="red",pch=16)
  lines(seq,beta2_invasion_fitness,type = "l",,xlab = "trait b",ylab="fitness",col="blue")
  points(beta[2],bb2_fit,col="green",pch=16)
}

# initial values

# populations nums
R_pop<-2
B_pop<-2

# population sizes
R<-rep(NA,R_pop)
B<-rep(NA,B_pop)
r<-rep(NA,R_pop)
b<-rep(NA,B_pop)

# trait values
initials<-list(
  r<-c(r1<-0,r2<-0),
  b<-c(b1<-0,b2<-0),
  beta<-c(beta1<-1.5, beta2<-1.5)
)

t<-1; 
r1_t<-r[1]; r2_t<-r[2];b1_t<-b[1]; b2_t<-b[2]
beta1_t<-beta1; beta2_t<-beta2

# the system ESS

repeat { 
  r[1]<-r1_t
  r[2]<-r2_t
  b[1]<-b1_t
  b[2]<-b2_t
  beta[1]<-beta1_t
  beta[2]<-beta2_t
  ra<-ra_fun(r,b,R_pop,B_pop)
  inter<-interaction(r,b,beta,R_pop,B_pop)
  pop<-solve(inter,-ra)
  R<-pop[1:2]
  B<-pop[3:4]
  lambd<-Jacob(r,b,R,B)
  sel_R1<-R1_sele_gradient(r1_t,R,B)
  sel_R2<-R2_sele_gradient(r2_t,R,B)
  sel_B1<-B1_sele_gradient(b1_t,beta1_t,R,B)
  sel_B2<-B2_sele_gradient(b2_t,beta2_t,R,B)
  curve_R1<-R1_curvature(r1_t,R,B)
  curve_R2<-R2_curvature(r2_t,R,B)
  curve_B1<-B1_curvature(b1_t,beta1_t,R,B)
  curve_B2<-B2_curvature(b2_t,beta2_t,R,B)
  if(t%%10==0){cat(t,"traits","r1=",r1_t,"r2=",r2_t,"b1=",b1_t,"b2=",b2_t,"beta1=",beta1_t,"beta2=",beta2_t,"\n")}
  if ((any(R[which(!is.na(B)==1)]<0))|(any(B[which(!is.na(B)==1)]<0))){
    cat("extinction")
    break
  }
  if (any(Re(lambd)>0)){
    cat("unstable")
    break
  }
  if((any(abs(sel_R1)<1e-08)) | (any(abs(sel_R2)<1e-08)) | (any(abs(sel_B1)<1e-08)) | (any(abs(sel_B2)<1e-08)) ){
    
    
    if(all(curve_R1<0.0001)&all(curve_R2<0.0001)&all(curve_B1<0.0001)&all(curve_B2<0.0001)& (all(abs(sel_R1)<1e-08)) & (all(abs(sel_R2)<1e-08)) & (all(abs(sel_B1)<1e-08)) & (all(abs(sel_B2)<1e-08)) ){
      ESS_time<-t
      cat("ESS","at time t=",t,'!!!!!!\n')
      cat("traits","r1=",r1_t,"r2=",r2_t,"b1=",b1_t,"b2=",b2_t,"beta1=",beta1_t,"beta2=",beta2_t)
      break
    } else if((abs(sel_R1[1])<1e-08)&(curve_R1[1]>=0.0001)) {cat("branching");break}
    else if ((abs(sel_R2[1])<1e-08)&(curve_R2[1]>=0.0001)) {cat("branching");break}
    else if ((abs(sel_B1[1])<1e-08)&(curve_B1[1]>=0.0001)) {cat("branching");break}
    else if ((abs(sel_B1[2])<1e-08)&(curve_B1[2]>=0.0001)) {cat("branching");break}
    else if ((abs(sel_B2[1])<1e-08)&(curve_B2[1]>=0.0001)) {cat("branching");break}
    else if ((abs(sel_B2[2])<1e-08)&(curve_B2[2]>=0.0001)) {cat("branching");break}
    
    else {
      r1_t<-r1_t+1/2*mu1*mu_sig^2*mu_step*sel_R1[1]*R[1]
      r2_t<-r2_t+1/2*mu1*mu_sig^2*mu_step*sel_R2[1]*R[2]
      b1_t<-b1_t+1/2*mu1*mu_sig^2*mu_step*sel_B1[1]*B[1]
      beta1_t<-beta1_t+1/2*mu2*mu_sig^2*mu_step*sel_B1[2]*B[1]
      b2_t<-b2_t+1/2*mu1*mu_sig^2*mu_step*sel_B2[1]*B[2]
      beta2_t<-beta2_t+1/2*mu2*mu_sig^2*mu_step*sel_B2[2]*B[2]
      t<-t+1
    }
  } else {
    r1_t<-r1_t+1/2*mu1*mu_sig^2*mu_step*sel_R1[1]*R[1]
    r2_t<-r2_t+1/2*mu1*mu_sig^2*mu_step*sel_R2[1]*R[2]
    b1_t<-b1_t+1/2*mu1*mu_sig^2*mu_step*sel_B1[1]*B[1]
    beta1_t<-beta1_t+1/2*mu2*mu_sig^2*mu_step*sel_B1[2]*B[1]
    b2_t<-b2_t+1/2*mu1*mu_sig^2*mu_step*sel_B2[1]*B[2]
    beta2_t<-beta2_t+1/2*mu2*mu_sig^2*mu_step*sel_B2[2]*B[2]
    t<-t+1
  }
}

# population sizes & fitness landscapes
if ((all(R[which(!is.na(B)==1)]>0))&(all(B[which(!is.na(B)==1)]>0))&(all(lambd<0))) {     
  par(mfrow=c(1,1))
  par(mar = c(2, 3, 2, 2))
  barplot(c(R[1],R[2],B[1],B[2]),names.arg = c("R1","R2","B1","B2"))
  title(main = "population size at ESS",sub = "My subtitle")
  beta_fitness_landscape(b,beta,R,B)
  title(main = "fitness landscape at ESS",sub = "My subtitle")
  legend(4.8, 0, expression(paste(beta[1]),paste(beta[2])), col = c("black","blue"), lty = 1,
         pch = "", ncol = 1, cex = 0.8)
} else {cat ("extinct or unstable")}

pop_ess<-pop


# disturbance

KR1<-4 # carrying capacity changes

if ((all(R[which(!is.na(B)==1)]>0))&(all(B[which(!is.na(B)==1)]>0))&(all(lambd<0))) {
  par(mfrow=c(1,1))
  par(mar = c(2, 3, 2, 2))
  barplot(c(R[1],R[2],B[1],B[2]),names.arg = c("R1","R2","B1","B2"))
  title(main = "population size with disturbance",sub = "My subtitle")
  beta_fitness_landscape(b,beta,R,B)
  title(main = "fitness landscape with disturbance",sub = "My subtitle")
  legend(4.8, 0, expression(paste(beta[1]),paste(beta[2])), col = c("black","blue"), lty = 1,
         pch = "", ncol = 1, cex = 0.8)
  
}else {"extinction or unstable"}

pop_disturb<-pop


# the new system ESS after disturbance

t<-1; 
r1_t<-r[1]; r2_t<-r[2];b1_t<-b[1]; b2_t<-b[2]
beta1_t<-beta1; beta2_t<-beta2

repeat { 
  r[1]<-r1_t
  r[2]<-r2_t
  b[1]<-b1_t
  b[2]<-b2_t
  beta[1]<-beta1_t
  beta[2]<-beta2_t
  ra<-ra_fun(r,b,R_pop,B_pop)
  inter<-interaction(r,b,beta,R_pop,B_pop)
  pop<-solve(inter,-ra)
  R<-pop[1:2]
  B<-pop[3:4]
  lambd<-Jacob(r,b,R,B)
  sel_R1<-R1_sele_gradient(r1_t,R,B)
  sel_R2<-R2_sele_gradient(r2_t,R,B)
  sel_B1<-B1_sele_gradient(b1_t,beta1_t,R,B)
  sel_B2<-B2_sele_gradient(b2_t,beta2_t,R,B)
  curve_R1<-R1_curvature(r1_t,R,B)
  curve_R2<-R2_curvature(r2_t,R,B)
  curve_B1<-B1_curvature(b1_t,beta1_t,R,B)
  curve_B2<-B2_curvature(b2_t,beta2_t,R,B)
  if(t%%10==0){cat(t,"traits","r1=",r1_t,"r2=",r2_t,"b1=",b1_t,"b2=",b2_t,"beta1=",beta1_t,"beta2=",beta2_t,"pop_size=",pop,"\n")}
  if ((any(R[which(!is.na(B)==1)]<0))|(any(B[which(!is.na(B)==1)]<0))){
    cat("extinction")
    break
  }
  if (any(Re(lambd)>0)){
    cat("unstable")
    break
  }
  
  if((any(abs(sel_R1)<1e-08)) | (any(abs(sel_R2)<1e-08)) | (any(abs(sel_B1)<1e-08)) | (any(abs(sel_B2)<1e-08)) ){
    
    if(all(curve_R1<0.0001)&all(curve_R2<0.0001)&all(curve_B1<0.0001)&all(curve_B2<0.0001)& (all(abs(sel_R1)<1e-08)) & (all(abs(sel_R2)<1e-08)) & (all(abs(sel_B1)<1e-08)) & (all(abs(sel_B2)<1e-08)) ){
      ESS_time<-t
      cat("ESS","at time t=",t,'!!!!!!\n')
      cat("traits","r1=",r1_t,"r2=",r2_t,"b1=",b1_t,"b2=",b2_t,"beta1=",beta1_t,"beta2=",beta2_t,"pop_size=",pop)
      break
    } else if((abs(sel_R1[1])<1e-08)&(curve_R1[1]>=0.0001)) {cat("branching");break}
     else if ((abs(sel_R2[1])<1e-08)&(curve_R2[1]>=0.0001)) {cat("branching");break}
    else if ((abs(sel_B1[1])<1e-08)&(curve_B1[1]>=0.0001)) {cat("branching");break}
    else if ((abs(sel_B1[2])<1e-08)&(curve_B1[2]>=0.0001)) {cat("branching");break}
    else if ((abs(sel_B2[1])<1e-08)&(curve_B2[1]>=0.0001)) {cat("branching");break}
    else if ((abs(sel_B2[2])<1e-08)&(curve_B2[2]>=0.0001)) {cat("branching");break}
    
    else {
      r1_t<-r1_t+1/2*mu1*mu_sig^2*mu_step*sel_R1[1]*R[1]
      r2_t<-r2_t+1/2*mu1*mu_sig^2*mu_step*sel_R2[1]*R[2]
      b1_t<-b1_t+1/2*mu1*mu_sig^2*mu_step*sel_B1[1]*B[1]
      beta1_t<-beta1_t+1/2*mu2*mu_sig^2*mu_step*sel_B1[2]*B[1]
      b2_t<-b2_t+1/2*mu1*mu_sig^2*mu_step*sel_B2[1]*B[2]
      beta2_t<-beta2_t+1/2*mu2*mu_sig^2*mu_step*sel_B2[2]*B[2]
      t<-t+1
    }
  } else {
    r1_t<-r1_t+1/2*mu1*mu_sig^2*mu_step*sel_R1[1]*R[1]
    r2_t<-r2_t+1/2*mu1*mu_sig^2*mu_step*sel_R2[1]*R[2]
    b1_t<-b1_t+1/2*mu1*mu_sig^2*mu_step*sel_B1[1]*B[1]
    beta1_t<-beta1_t+1/2*mu2*mu_sig^2*mu_step*sel_B1[2]*B[1]
    b2_t<-b2_t+1/2*mu1*mu_sig^2*mu_step*sel_B2[1]*B[2]
    beta2_t<-beta2_t+1/2*mu2*mu_sig^2*mu_step*sel_B2[2]*B[2]
    t<-t+1
    
  }
}

# population sizes and fitness landscapes
if ((all(R[which(!is.na(B)==1)]>0))&(all(B[which(!is.na(B)==1)]>0))&(all(lambd<0))) {     
  par(mfrow=c(1,1))
  par(mar = c(2, 3, 2, 2))
  barplot(c(R[1],R[2],B[1],B[2]),names.arg = c("R1","R2","B1","B2"))
  title(main = "population size at new ESS",sub = "My subtitle")
  beta_fitness_landscape(b,beta,R,B)
  title(main = "fitness landscape at new ESS",sub = "My subtitle")
  legend(4.8, 0, expression(paste(beta[1]),paste(beta[2])), col = c("black","blue"), lty = 1,
         pch = "", ncol = 1, cex = 0.8)
} else {cat ("extinct or unstable")}

pop_new_ess<-pop

delta_disturb <- pop_disturb - pop_ess
delta_ess <- pop_new_ess - pop_ess

barplot(delta_disturb,names.arg = c("R1","R2","B1","B2"),horiz=F)
title(main = "population change with disturbance",sub = "My subtitle")

barplot(delta_ess,names.arg = c("R1","R2","B1","B2"),horiz=F)
title(main = "population change at new ESS",sub = "My subtitle")





