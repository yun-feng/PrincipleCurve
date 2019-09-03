t=(1:1000)/1000

#t=floor(rnorm(100,mean=0.5,sd=0.2)*100)/100
#t[1:50]=floor(rnorm(50,mean=0.3,sd=0.1)*100)/100
#t[51:100]=floor(rnorm(50,mean=0.6,sd=0.1)*100)/100
v=2*pi
x<-cos(v*t)
y<-sin(v*t)
#x=v*t
#y=v*t
data<-matrix(0,2,1000)
data[1,]=x+rnorm(1000,0,0.1)
data[2,]=y+rnorm(1000,0,0.1)

full_data=data
data=full_data[,sample.int(1000,100)]

x_0=c(0,0)

x_0=apply(data,1,mean)
#alpha_0=eigen((data-x_0)%*%t(data-x_0))$vectors[,1]
alpha_0=runif(2)
alpha_0=alpha_0/sum(alpha_0^2)^0.5
t=(alpha_0%*%(data-x_0))

temp_data=x_0+alpha_0%*%t

v=2*pi
t=t/v
alpha_0=alpha_0*v

t=t+0.5
x_0=x_0-alpha_0*0.5


N=500
alpha=matrix(rep(c(alpha_0),2*N),2,N)
prob_t=matrix(runif(100*N),nrow=100,ncol=N)
prob_t=prob_t/apply(prob_t,1,sum)
x=matrix(nrow=2,ncol=N)
tc=apply(data,1,mean)
rho=1
plot(temp_data[1,],temp_data[2,],xlim=c(-1.5,1.5),ylim=c(-1.5,1.5))
for ( c in 1:1000){
  temp_x=x_0
  prob_t[,1]=exp(-20*apply((data-temp_x)^2,2,sum))
  x[,1]=x_0
  for( i in 2:N){
    temp_x=temp_x+alpha[,i-1]/N
    x[,i]=temp_x
    prob_t[,i]=exp(-20*apply((data-temp_x)^2,2,sum))
  }
  prob_t=prob_t/apply(prob_t,1,sum)
  xc=apply(t(x)*apply(prob_t,2,sum),2,sum)/100
  x_0=x_0+(tc-xc)
  
  for(bc in 1:1){
    K=100
    h=0
    temp_x=x_0
    h=h-apply(t(data-temp_x)*prob_t[,1],2,sum)
    K=K-sum(prob_t[,1])
    
    alpha_0=alpha[,1]
    l=h+K*alpha_0*1/N
    alpha[,1]=0*alpha[,1]+sum(alpha[,2]*alpha[,1])*alpha[,2]+rho*l
    alpha[,1]=alpha[,1]/sum(alpha[,1]*alpha[,1])^0.5*v
    h=h-c((alpha[,1]-alpha_0)/N)*K
    
    
    la<-c()
    start_flag=1
    end_flag=0
    
    for(i in 2:(N-1)){
        temp_x=temp_x+alpha[,i-1]/N
        h=h-apply(t(data-temp_x)*prob_t[,i],2,sum)
        K=K-sum(prob_t[,i])
        
      if(sum(prob_t[,i])<1e-2){
        if(start_flag){next}
        if (end_flag){next}
        else{end_flag=1;
            end_site=i}
      }
      else{end_flag=0}
      if(start_flag){
        start_flag=0
        alpha[,1:(i-1)]=alpha[,i]
      }
        
      alpha_0=alpha[,i]
      l=h+K*alpha_0/N
      alpha[,i]=-0*alpha[,i]+
        sum(alpha[,i+1]*alpha[,i])*alpha[,i+1]+
        sum(alpha[,i-1]*alpha[,i])*alpha[,i-1]+rho*l
      alpha[,i]=alpha[,i]/sum(alpha[,i]*alpha[,i])^0.5*v

      
      h=h-c((alpha[,i]-alpha_0)/N)*K
      
      la<-cbind(la,l)
      
    }
    alpha[,N]=alpha[,N-1]
    if(end_flag){
      alpha[,end_site:N]=alpha[,end_site-1]
    }
    
    
    "
    for(i in 2:(N-1)){
      temp_x=temp_x+alpha[,i-1]/N
      h=h-apply(t(data-temp_x)*prob_t[,i],2,sum)
      K=K-sum(prob_t[,i])
      
      if(sum(prob_t[,i])<1e-2){
        if(start_flag){next}
        if (end_flag){next}
        else{end_flag=1;
        end_site=i}
      }
      else{end_flag=0}
      if(start_flag){
        start_flag=0
        alpha[,1:(i-1)]=alpha[,i]
      }
      
      alpha_0=alpha[,i]
      l=h+K*alpha_0/N
      alpha[,i]=-0*alpha[,i]+
        sum(alpha[,i+1]*alpha[,i])*alpha[,i+1]+
        sum(alpha[,i-1]*alpha[,i])*alpha[,i-1]+rho*l
      lambda=v+abs(sum(alpha[,i+1]*alpha[,i-1]))/v
      if(sum(alpha[,i]*alpha[,i])^0.5 >= lambda){
        alpha[,i]=alpha[,i]/sum(alpha[,i]*alpha[,i])^0.5*v
      }
      else{
        
        alpha[,i]=alpha[,i]/lambda*v
        
        
      }
      
      
      h=h-c((alpha[,i]-alpha_0)/N)*K
      
      la<-cbind(la,l)
      
    }
    alpha[,N]=alpha[,N-1]
    if(end_flag){
      alpha[,end_site:N]=alpha[,end_site-1]
    }
    
    
    "
  }  

  
  
 
  
   n=runif(2)
    n=n/sqrt(sum(n^2))
    ori=n
    temp_l=x-x_0
    temp_h=data-x_0
    S=-temp_l%*%t(prob_t)%*%t(temp_h)
    n=max(abs(S))*n
    n=n+S%*%ori-2*ori*sum(ori*S%*%ori)-t(S)%*%ori
    n=n/sqrt(sum(n^2))
    alpha=(diag(c(1,1))-2*n%*%t(n))%*%(diag(c(1,1))-2*ori%*%t(ori))%*%alpha
  
  
  
  
  
  
  plot(x[1,],x[2,],xlim=c(-1.5,1.5),ylim=c(-1.5,1.5))
  

}

plot(data[1,],data[2,],xlim=c(-1.5,1.5),ylim=c(-1.5,1.5))
