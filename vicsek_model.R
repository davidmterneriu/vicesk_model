################################################################################
#Vicsek model
#David Terner 
#Spring 2022
################################################################################

rm(list=ls())

library(readr)
library(readxl)
library(tidyverse)
library(ggplot2)
library(tidyverse)
library(latex2exp)
options(dplyr.summarise.inform = FALSE)


#https://www.r-bloggers.com/2013/05/pairwise-distances-in-r/
vectorized_pdist <- function(A,B){ 
  an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
  bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))
  m = nrow(A)
  n = nrow(B)
  tmp = matrix(rep(an, n), nrow=m) 
  tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
  return(sqrt( tmp - 2 * tcrossprod(A,B) ))
}

A1=matrix(runif(100*2),ncol = 2,nrow = 100)
B1=matrix(runif(100*2),ncol = 2,nrow = 100)

sum((A1[1,]-B1[1,])^2)%>%sqrt()



ex1=vectorized_pdist(A=A1,B=B1)
ex1[1,1]




#Starting parameters

periodic_dist=function(pos_vec,L){
  
  n_birds=dim(pos_vec)[1]
  add_L=c(-L,0,L)
  L_grid=expand.grid(xL=add_L,yL=add_L)
  
  org_ind=which(L_grid[,1]==0 &L_grid[,2]==0)
  
  big_array=array(0,dim=c(9,n_birds,2))
  
  for(i in 1:nrow(L_grid)){
    big_array[i,1:n_birds,1]<-pos_vec[,1]+ L_grid[i,1]
    big_array[i,1:n_birds,2]<-pos_vec[,2]+ L_grid[i,2]
  }
  test_indx=expand.grid(org_ind,1:9)
  indx <- !duplicated(t(apply( test_indx, 1, sort)))
  test_indx=test_indx[ indx,]
  test_array_d=nrow(test_indx)
  
  dist_array=array(0,dim=c(test_array_d,n_birds,n_birds))
  
  for(i in 1:test_array_d){
    dist_array[i,1:n_birds,1:n_birds]<-vectorized_pdist(A=big_array[test_indx[i,1],1:n_birds,1:2],
                                                        B=big_array[test_indx[i,2],1:n_birds,1:2])
  }
  r1=apply(dist_array, c(2,3), min)

  return(r1)
}

nbrid=100
A1=matrix(runif(nbrid*2),ncol = 2,nrow = nbrid)
system.time(ex1<-periodic_dist(pos_vec = A1,L=1))
ex1




vicsek_fn=function(N=100,L=1,vel=0.5,eta=0.3,R=1,steps=100){
  pos_0=data.frame(bird=seq(1,N),pos_x=runif(N,0,L),pos_y=runif(N,0,L))
  pos_0=pos_0%>%mutate(z=pos_x+1i*pos_y,
                       angle=2*pi*(runif(N,0,L)-0.5))
  order_par_df=data.frame(steps=1:steps,order_par=0)
  dt=1/1000
  
  bird_df=pos_0%>%mutate(steps=0)
  
  
  for(s in 1:steps){
    order_par_df$order_par[s]= pos_0%>%mutate(v=vel*exp(angle*1i))%>%
      summarise(v1=abs(mean(v))/vel)%>%
      unlist()%>%
      as.numeric()
    
    d=periodic_dist(pos_vec = pos_0[,c("pos_x","pos_y")],L=L)
    d[d>=R]<-0
    y1=1/N*d%*%sin(pos_0$angle)
    x1=1/N*d%*%cos(pos_0$angle)
    avg_angle=atan2( y1,x1)
    
    #avg_angle=1/N*(d%*%pos_0$angle)
    next_angle= avg_angle+runif(N,-eta/2,eta/2)
    #Update
    pos_0=pos_0%>%
      mutate(pos_x=pos_x+vel*dt*cos(angle),
             pos_y=pos_y+vel*dt*sin(angle),
             pos_x=ifelse(pos_x<0,L+pos_x,
                          ifelse(L<pos_x,pos_x-L,pos_x)),
             pos_y=ifelse(pos_y<0,L+pos_y,
                          ifelse(L<pos_y,pos_y-L,pos_y)),
             z=pos_x+1i*pos_y)
    
    bird_df=rbind(bird_df,
                  pos_0%>%select(bird,pos_x,pos_y,z,angle)%>%mutate(steps=s))
    
    
    pos_0$angle=next_angle
    
  }
  
  
  
  return( list(bird_df=bird_df,order=order_par_df))
  
}


################################################################################
# Rho and eta 
################################################################################


param_df2=expand.grid(rho=10^seq(-1,1,length.out=20),eta=seq(0,5,length.out=20))%>%
  as.data.frame()
param_df2$N=40
param_df2$L=sqrt(param_df2$N/param_df2$rho)


fig2_df=data.frame()

tictoc::tic()
#nrow(fig1_par_df)
for(i in 1:nrow(param_df2)){
  N1=param_df2$N[i]
  L1=param_df2$L[i]
  eta1=param_df2$eta[i]
  temp_res=data.frame()
  
  
  for(t in 1:10){
    mes1= paste0("i: ",signif(i/nrow(param_df2),3),"|| t: ",t/50)
    print(mes1)
    temp_res=rbind( temp_res,vicsek_fn(N=N1,L=L1,vel=0.03,R=1,eta= eta1)$order)
  }
  res_vec=temp_res$order_par[temp_res$steps>79]
  mean_order=mean(res_vec)
  se_mean=sd(res_vec)/sqrt(length(res_vec))
  
  fig2_df=rbind(fig2_df,
                data.frame(N=N1,L=L1,
                           order= mean_order,
                           std.err=se_mean))
  
  
  
}

tictoc::toc()


fig2_df$eta=param_df2$eta
fig2_df$rho=param_df2$rho

fig2_df%>%
  ggplot(aes(x=eta,y=rho,z=order))+
  geom_contour_filled(binwidth = 0.05)+
  scale_fill_viridis_d(option="C")+
  theme_minimal()+
  labs(x=TeX("$\\eta$"),
       y=TeX("$\\rho$"),
       fill=TeX("$v_a$"))+
  scale_y_continuous(breaks = scales::pretty_breaks())+
  scale_x_continuous(breaks = scales::pretty_breaks())+
  guides(fill=guide_legend(nrow = 2))+
  theme(legend.position = "bottom",
        axis.title.y = element_text(angle = 0),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
  
  





###############################################################################


library(gganimate)
library(gifski)
library(transformr)


df1=vicsek_fn(R=0.1,steps=50,vel = 10,L=1)
df1$order


bird_temp=df1$bird_df


g1=bird_temp%>%
  ggplot(aes(x=pos_x,y=pos_y,color=angle))+
  geom_point()+
  scale_color_viridis_c()+
  theme_bw()+
  labs(x="",y="",color="Angle")+
  scale_x_continuous(breaks = scales::pretty_breaks())+
  scale_y_continuous(breaks = scales::pretty_breaks())

#g1


g1=g1+
  transition_time(steps) +
  labs(subtitle = "Step: {frame_time}")+
  ease_aes('linear')

animate(g1, width = 400, height = 400, fps = 2, duration = 20,
        end_pause = 5, res = 100,renderer = gifski_renderer())

anim_save("~/Desktop/IU/INFO/I609/Final Project/ex1.gif")

################################################################################
# Figure 2a Replication
################################################################################

rho1=4
fig1_par_df=data.frame(N=c(40,100,400,4000,10000))
fig1_par_df=fig1_par_df%>%
  mutate(L=sqrt(N/rho1))

eta_seq=seq(0,5,length.out=25)

fig1_df=data.frame()

tictoc::tic()
#nrow(fig1_par_df)
for(i in 1:3){
  N1=fig1_par_df$N[i]
  L1=fig1_par_df$L[i]
  for(j in 1:length(eta_seq)){
    temp_res=data.frame()
    
    
    for(t in 1:50){
      mes1= paste0("i: ",signif(i/nrow(fig1_par_df),3),"|| j: ", signif(j/length(eta_seq),3),"|| t:",t/10)
      print(mes1)
      temp_res=rbind( temp_res,vicsek_fn(N=N1,L=L1,vel=0.03,R=1,eta=eta_seq[j])$order)
    }
    res_vec=temp_res$order_par[temp_res$steps>79]
    mean_order=mean(res_vec)
    se_mean=sd(res_vec)/sqrt(length(res_vec))
    res1=t.test(temp_res$order_par)
    data.frame()
    
    fig1_df=rbind(fig1_df,
                  data.frame(N=N1,L=L1,eta=eta_seq[j],
                             order= mean_order,
                             std.err=se_mean))
  }
  
  
}

tictoc::toc()

write_csv(fig1_df,"~/Desktop/IU/INFO/I609/Final Project/fig1_df.csv")

fig1_df=read_csv("~/Desktop/IU/INFO/I609/Final Project/fig1_df.csv")

fig1_df%>%
  mutate(N1=as.character(N)%>%factor(c(40,100,400)))%>%
  ggplot(aes(x=eta,y=order,color=N1,shape=N1))+
  geom_point(size=2)+
  labs(x=TeX("$\\eta$"),y=TeX("$v_a$"),color="N",shape="N")+
  ggthemes::theme_few()+
  ggthemes::scale_color_calc()+
  theme(axis.title.y = element_text(angle = 0))+
  scale_y_continuous(breaks = seq(0,1,0.2),limits = c(0,1))

################################################################################
# Figure 2b Replication
################################################################################



param_df2=data.frame(N=4*10^seq(1,3,length.out=20))
param_df2$L=20
param_df2$rho=param_df2$N/param_df2$L^2

param_df2=data.frame(rho=10^seq(-1,1,length.out=40))
param_df2$N=40
param_df2$L=sqrt(param_df2$N/param_df2$rho)


fig2_df=data.frame()

tictoc::tic()
#nrow(fig1_par_df)
for(i in 1:nrow(param_df2)){
  N1=param_df2$N[i]
  L1=param_df2$L[i]
  temp_res=data.frame()
  
  
  for(t in 1:50){
    mes1= paste0("i: ",signif(i/nrow(param_df2),3),"|| t: ",t/50)
    print(mes1)
    temp_res=rbind( temp_res,vicsek_fn(N=N1,L=L1,vel=0.03,R=1,eta=0.01)$order)
  }
  res_vec=temp_res$order_par[temp_res$steps>79]
  mean_order=mean(res_vec)
  se_mean=sd(res_vec)/sqrt(length(res_vec))
  
  fig2_df=rbind(fig2_df,
                data.frame(N=N1,L=L1,
                           order= mean_order,
                           std.err=se_mean))
  
  
  
}

tictoc::toc()

write_csv(fig2_df,"~/Desktop/IU/INFO/I609/Final Project/fig1b_df.csv")

################################################################################
fig2_df<-read_csv("~/Desktop/IU/INFO/I609/Final Project/fig1b_df.csv")

fig2_df%>%
  mutate(rho=N/(L^2))%>%View()

fig2_df%>%
  mutate(rho=N/(L^2))%>%
  ggplot(aes(x=rho,y=order))+
  geom_point()+
  geom_line()+
  geom_point(size=2)+
  labs(x=TeX("$\\rho$"),y=TeX("$v_a$"))+
  ggthemes::theme_few()+
  ggthemes::scale_color_calc()+
  theme(axis.title.y = element_text(angle = 0))+
  scale_y_continuous(breaks = seq(0.3,1,0.1),limits = c(0.3,1))+
  scale_x_continuous(breaks = seq(0,10,2))
