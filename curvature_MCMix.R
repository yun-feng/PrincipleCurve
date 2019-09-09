t1=(1:100)/100

#t=floor(rnorm(100,mean=0.5,sd=0.2)*100)/100
#t[1:50]=floor(rnorm(50,mean=0.3,sd=0.1)*100)/100
#t[51:100]=floor(rnorm(50,mean=0.6,sd=0.1)*100)/100

#crossing half circle
v=1*pi
x<-cos(v*t1)
y<-sin(v*t1)
x<-c(x,cos(v*t1))
y<-c(y,1-sin(v*t1))


#crissing square
"x=v*t1
y=v*((2*t1)^2)#+rnorm(100,0,1)
x<-c(x,x)
y<-c(y,-y+7)/5"
data<-matrix(0,2,100*2)
data[1,]=x
data[2,]=y

x1_0=runif(2)
x2_0=runif(2)

alpha1_0=runif(2)
alpha1_0=alpha1_0/sum(alpha1_0*alpha1_0)^0.5
alpha2_0=runif(2)
alpha2_0=alpha2_0/sum(alpha2_0*alpha2_0)^0.5

z=runif(200)

v=1*pi
#crossing sqaure
#v=2.2*pi

alpha1_0=alpha1_0/sqrt(sum(alpha1_0^2))*v

x1_0=x1_0-alpha1_0*0.5

alpha2_0=alpha2_0/sqrt(sum(alpha2_0^2))*v

x2_0=x2_0-alpha2_0*0.5

N=250
alpha1=matrix(rep(c(alpha1_0),2*N),2,N)
alpha2=matrix(rep(c(alpha2_0),2*N),2,N)

prob_t1=matrix(runif(200*N),nrow=200,ncol=N)
prob_t1=prob_t1/apply(prob_t1,1,sum)
prob_t2=matrix(runif(200*N),nrow=200,ncol=N)
prob_t2=prob_t2/apply(prob_t2,1,sum)
x1=matrix(nrow=2,ncol=N)
x2=matrix(nrow=2,ncol=N)

#rho=10

xlim_plot=c(-1,1)
ylim_plot=c(0,1)
#crossing half circle
t_act=10
k=1
#crossing half circle
rho_delta=1e-2*N
#crossing square
#rho_delta=9e-3*N
#t_act=6
rho=(1/rho_delta)/2/N

for ( c in 1:1000){
  
  
  
  Er1=0
  temp_x=x1_0
  prob_t1[,1]=-t_act*apply((data-temp_x)^2,2,sum)
  Er1=Er1+exp(-t_act*apply((data-temp_x)^2,2,sum))
  x1[,1]=x1_0
  for( i in 2:N){
    temp_x=temp_x+alpha1[,i-1]/N
    x1[,i]=temp_x
    prob_t1[,i]=-t_act*apply((data-temp_x)^2,2,sum)
    Er1=Er1+exp(-t_act*apply((data-temp_x)^2,2,sum))
  }
  prob_t1=prob_t1-apply(prob_t1,1,max)
  prob_t1=exp(prob_t1)
  prob_t1=prob_t1/apply(prob_t1,1,sum)
  
  
  Er2=0
  temp_x=x2_0
  prob_t2[,1]=-t_act*apply((data-temp_x)^2,2,sum)
  Er2=Er2+exp(-t_act*apply((data-temp_x)^2,2,sum))
  x2[,1]=x2_0
  for( i in 2:N){
    temp_x=temp_x+alpha2[,i-1]/N
    x2[,i]=temp_x
    prob_t2[,i]=-t_act*apply((data-temp_x)^2,2,sum)
    Er2=Er2+exp(-t_act*apply((data-temp_x)^2,2,sum))
  }
  prob_t2=prob_t2-apply(prob_t2,1,max)
  prob_t2=exp(prob_t2)
  prob_t2=prob_t2/apply(prob_t2,1,sum)
  
  Er1=Er1^k
  Er2=Er2^k
  z=Er1/(Er1+Er2)
  Er=sum((cbind(alpha1,alpha1[,N])-cbind(alpha1[,1],alpha1))^2)+
    sum((cbind(alpha2,alpha2[,N])-cbind(alpha2[,1],alpha2))^2)#+
  #4*sum(alpha1^2)+4*sum(alpha2^2)
  
  #Er=Er/(2*rho*N)
  Er=Er*rho_delta
  Er=Er+sum(-log(Er1+Er2))*1/t_act
  print(Er)
  
  
  
  '
  
    K=sum(z)
    h=0
    temp_x=x1_0
    h=h-apply(t(data-temp_x)*prob_t1[,1]*z,2,sum)
    K=K-sum(prob_t1[,1]*z)
    
    alpha_0=alpha1[,1]
    l=h+K*alpha_0*1/N
    alpha1[,1]=0*alpha1[,1]+sum(alpha1[,2]*alpha1[,1])*alpha1[,2]+rho*l
    alpha1[,1]=alpha1[,1]/sum(alpha1[,1]*alpha1[,1])^0.5*v
    h=h-c((alpha1[,1]-alpha_0)/N)*K
    
    
    la<-c()
    for(i in 2:(N-1)){
      temp_x=temp_x+alpha1[,i-1]/N
      h=h-apply(t(data-temp_x)*prob_t1[,i]*z,2,sum)
      K=K-sum(prob_t1[,i]*z)
      
      alpha_0=alpha1[,i]
      l=h+K*alpha_0/N
      alpha1[,i]=-0*alpha1[,i]+
        sum(alpha1[,i+1]*alpha1[,i])*alpha1[,i+1]+
        sum(alpha1[,i-1]*alpha1[,i])*alpha1[,i-1]+rho*l
      #alpha1[,i]=alpha1[,i]/sum(alpha1[,i]*alpha1[,i])^0.5*v
      
      lambda=v+abs(sum(alpha1[,i+1]*alpha1[,i-1]))/v
      if(sum(alpha1[,i]*alpha1[,i])^0.5 >= lambda){
        alpha1[,i]=alpha1[,i]/sum(alpha1[,i]*alpha1[,i])^0.5*v
      }
      else{
        
        alpha1[,i]=alpha_0
        print("wow")
        
      }
      
      
      
      h=h-c((alpha1[,i]-alpha_0)/N)*K
      
      la<-cbind(la,l)
      
    }
    alpha1[,N]=alpha1[,N-1]
    
    
    K=sum((1-z))
    h=0
    temp_x=x2_0
    h=h-apply(t(data-temp_x)*prob_t2[,1]*(1-z),2,sum)
    K=K-sum(prob_t2[,1]*(1-z))
    
    alpha_0=alpha2[,1]
    l=h+K*alpha_0*1/N
    alpha2[,1]=0*alpha2[,1]+sum(alpha2[,2]*alpha2[,1])*alpha2[,2]+rho*l
    alpha2[,1]=alpha2[,1]/sum(alpha2[,1]*alpha2[,1])^0.5*v
    h=h-c((alpha2[,1]-alpha_0)/N)*K
    
    
    la<-c()
    
    for(i in 2:(N-1)){
      temp_x=temp_x+alpha2[,i-1]/N
      h=h-apply(t(data-temp_x)*prob_t2[,i]*(1-z),2,sum)
      K=K-sum(prob_t2[,i]*(1-z))
      
      alpha_0=alpha2[,i]
      l=h+K*alpha_0/N
      alpha2[,i]=-0*alpha2[,i]+
        sum(alpha2[,i+1]*alpha2[,i])*alpha2[,i+1]+
        sum(alpha2[,i-1]*alpha2[,i])*alpha2[,i-1]+rho*l
      #alpha2[,i]=alpha2[,i]/sum(alpha2[,i]*alpha2[,i])^0.5*v
      
      lambda=v+abs(sum(alpha2[,i+1]*alpha2[,i-1]))/v
      if(sum(alpha2[,i]*alpha2[,i])^0.5 >= lambda){
        alpha2[,i]=alpha2[,i]/sum(alpha2[,i]*alpha2[,i])^0.5*v
      }
      else{
        
        alpha2[,i]=alpha_0
        
        
      }
      
      h=h-c((alpha2[,i]-alpha_0)/N)*K
      
      la<-cbind(la,l)
      
    }
    alpha2[,N]=alpha2[,N-1]
  
  
  
    
    Er1=0
    temp_x=x1_0
    prob_t1[,1]=-t_act*apply((data-temp_x)^2,2,sum)
    Er1=Er1+exp(-t_act*apply((data-temp_x)^2,2,sum))
    x1[,1]=x1_0
    for( i in 2:N){
      temp_x=temp_x+alpha1[,i-1]/N
      x1[,i]=temp_x
      prob_t1[,i]=-t_act*apply((data-temp_x)^2,2,sum)
      Er1=Er1+exp(-t_act*apply((data-temp_x)^2,2,sum))
    }
    prob_t1=prob_t1-apply(prob_t1,1,max)
    prob_t1=exp(prob_t1)
    prob_t1=prob_t1/apply(prob_t1,1,sum)
    
    
    Er2=0
    temp_x=x2_0
    prob_t2[,1]=-t_act*apply((data-temp_x)^2,2,sum)
    Er2=Er2+exp(-t_act*apply((data-temp_x)^2,2,sum))
    x2[,1]=x2_0
    for( i in 2:N){
      temp_x=temp_x+alpha2[,i-1]/N
      x2[,i]=temp_x
      prob_t2[,i]=-t_act*apply((data-temp_x)^2,2,sum)
      Er2=Er2+exp(-t_act*apply((data-temp_x)^2,2,sum))
    }
    prob_t2=prob_t2-apply(prob_t2,1,max)
    prob_t2=exp(prob_t2)
    prob_t2=prob_t2/apply(prob_t2,1,sum)
    
  '
  n=runif(2)
  n=n/sqrt(sum(n^2))
  ori=n
  temp_l=x1-x1_0
  temp_h=t(t(data-x1_0)*z)
  S=-temp_l[,2:N]%*%t(prob_t1[,2:N])%*%t(temp_h)
  S=S-2*ori%*%t(ori)%*%S+t(S)-2*t(S)%*%ori%*%t(ori)
  n=sum(abs(S))*n+S%*%n
#  n=2*sum(abs(S))*n
 # n=n+S%*%ori-2*ori*sum(ori*S%*%ori)-t(S)%*%ori
  n=n/sqrt(sum(n^2))
  alpha1=(diag(c(1,1))-2*n%*%t(n))%*%(diag(c(1,1))-2*ori%*%t(ori))%*%alpha1
  
  x1[,2:N]=x1_0+((diag(c(1,1))-2*n%*%t(n))%*%(diag(c(1,1))-2*ori%*%t(ori))%*%temp_l)[,2:N]
  
  
  temp_x=x1_0
  for(i in 2:(N-1)){
    n=runif(2)
    n=n/sqrt(sum(n^2))
    ori=n
    temp_x=temp_x+alpha1[,i-1]/N
    temp_l=x1-temp_x
    temp_h=t(t(data-temp_x)*z)
    S=-temp_l[,(i+1):N]%*%t(prob_t1[,(i+1):N])%*%t(temp_h)
    S=S-rho_delta*alpha1[,i]%*%t(alpha1[,i-1])
    S=S-2*ori%*%t(ori)%*%S+t(S)-2*t(S)%*%ori%*%t(ori)
    n=sum(abs(S))*n+S%*%n
    #n=2*sum(abs(S))*n
    #n=n+S%*%ori-2*ori*sum(ori*(S%*%ori))-t(S)%*%ori
    n=n/sqrt(sum(n^2))
    alpha1[,i:N]=(diag(c(1,1))-2*n%*%t(n))%*%(diag(c(1,1))-2*ori%*%t(ori))%*%alpha1[,i:N]
    x1[,(i+1):N]=temp_x+((diag(c(1,1))-2*n%*%t(n))%*%(diag(c(1,1))-2*ori%*%t(ori))%*%temp_l)[,(i+1):N]
    
  }
  

  
  
  n=runif(2)
  n=n/sqrt(sum(n^2))
  ori=n
  temp_l=x2-x2_0
  temp_h=t(t(data-x2_0)*(1-z))
  S=-temp_l[,2:N]%*%t(prob_t2[,2:N])%*%t(temp_h)
  S=S-2*ori%*%t(ori)%*%S+t(S)-2*t(S)%*%ori%*%t(ori)
  n=sum(abs(S))*n+S%*%n
  #n=2*sum(abs(S))*n
  #n=n+S%*%ori-2*ori*sum(ori*S%*%ori)-t(S)%*%ori
  n=n/sqrt(sum(n^2))
  alpha2=(diag(c(1,1))-2*n%*%t(n))%*%(diag(c(1,1))-2*ori%*%t(ori))%*%alpha2
  
  x2[,2:N]=x2_0+((diag(c(1,1))-2*n%*%t(n))%*%(diag(c(1,1))-2*ori%*%t(ori))%*%temp_l)[,2:N]
  
  temp_x=x2_0
  for(i in 2:(N-1)){
    n=runif(2)
    n=n/sqrt(sum(n^2))
    ori=n
    temp_x=temp_x+alpha2[,i-1]/N
    temp_l=x2-temp_x
    temp_h=t(t(data-temp_x)*(1-z))
    S=-temp_l[,(i+1):N]%*%t(prob_t2[,(i+1):N])%*%t(temp_h)
    S=S-rho_delta*alpha2[,i]%*%t(alpha2[,i-1])
    S=S-2*ori%*%t(ori)%*%S+t(S)-2*t(S)%*%ori%*%t(ori)
    n=sum(abs(S))*n+S%*%n
    #n=2*sum(abs(S))*n
    #n=n+S%*%ori-2*ori*sum(ori*(S%*%ori))-t(S)%*%ori
    n=n/sqrt(sum(n^2))
    alpha2[,i:N]=(diag(c(1,1))-2*n%*%t(n))%*%(diag(c(1,1))-2*ori%*%t(ori))%*%alpha2[,i:N]
    x2[,(i+1):N]=temp_x+((diag(c(1,1))-2*n%*%t(n))%*%(diag(c(1,1))-2*ori%*%t(ori))%*%temp_l)[,(i+1):N]
    
  }
  
  
  
  
  tc1=apply(t(data)*z,2,sum)/sum(z)
  tc2=apply(t(data)*(1-z),2,sum)/sum(1-z)
  xc1=apply(t(x1)*apply(prob_t1*z,2,sum),2,sum)/sum(z)
  x1_0=x1_0+(tc1-xc1)
  xc2=apply(t(x2)*apply(prob_t2*(1-z),2,sum),2,sum)/sum(1-z)
  x2_0=x2_0+(tc2-xc2)
  
  plot(c(x1[1,],x2[1,]),c(x1[2,],x2[2,]),xlim=xlim_plot,ylim=ylim_plot)
  
  
}



plot(data[1,],data[2,],xlim=xlim_plot,ylim=ylim_plot)
plot(apply(t(prob_t1[which(z>0.5),])*((1:N)/N),2,sum))
plot(apply(t(prob_t2[which(z<0.5),])*((1:N)/N),2,sum))
