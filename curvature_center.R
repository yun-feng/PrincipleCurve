####
#####
t=(1:100)/100

#t=floor(rnorm(100,mean=0.5,sd=0.2)*100)/100
#t[1:50]=floor(rnorm(50,mean=0.3,sd=0.1)*100)/100
#t[51:100]=floor(rnorm(50,mean=0.6,sd=0.1)*100)/100
v=1.5*pi
x<-cos(v*t)
y<-sin(v*t)
#x=v*t
#y=v*t
data<-matrix(0,2,100)
data[1,]=x
data[2,]=y

x_0=c(0,0)

x_0=apply(data,1,mean)
alpha_0=eigen((data-x_0)%*%t(data-x_0))$vectors[,1]
t=(alpha_0%*%(data-x_0))

temp_data=x_0+alpha_0%*%t

temp_data=temp_data[,order(t)]
data=data[,order(t)]
t=t[order(t)]

v=3*pi
t=t/v
alpha_0=alpha_0*v

#t=t+0.5
#x_0=x_0-alpha_0*0.5


N=10000
alpha_c=alpha_0
alpha1=matrix(rep(c(alpha_0),2*N),2,N)
alpha2=matrix(rep(c(alpha_0),2*N),2,N)
rho=5
plot(temp_data[1,],temp_data[2,],xlim=c(-1.5,1.5),ylim=c(-1.5,1.5))
for ( c in 1:101){
  
  upper=length(which(t>(1/N)))
  K=length(which(t>(1/N)))
  h=apply(data[,which(t>(1/N))]-temp_data[,which(t>(1/N))],1,sum)
  alpha_0=alpha1[,1]
  l=h+K*alpha_0*1/N
  alpha1[,1]=2*alpha1[,1]+sum(alpha1[,2]*alpha1[,1])*alpha1[,2]+sum(alpha_c*alpha1[,1])*alpha_c+rho*l
  alpha1[,1]=alpha1[,1]/sum(alpha1[,1]*alpha1[,1])^0.5*v
  temp_data=temp_data+c((alpha1[,1]-alpha_0)/N)
  h=h-c((alpha1[,1]-alpha_0)/N)*K
  
  
  for(i in (2):(N-1)){
    if(K>0 && i/N>=t[101-K]){
      h=h-(data[,101-K]-temp_data[,101-K])
      K=K-1
    }
    if(K==0){
      
       alpha1[,i:N]=alpha1[,i-1]
      break;
    }
    alpha_0=alpha1[,i]
    l=h+K*alpha_0/N
    alpha1[,i]=2*alpha1[,i]+sum(alpha1[,i+1]*alpha1[,i])*alpha1[,i+1]+
      sum(alpha1[,i-1]*alpha1[,i])*alpha1[,i-1]+rho*l
    alpha1[,i]=alpha1[,i]/sum(alpha1[,i]*alpha1[,i])^0.5*v
    if(K>0){
      temp_data[,(101-K):100]=temp_data[,(101-K):100]+c((alpha1[,i]-alpha_0)/N)
    }
    h=h-c((alpha1[,i]-alpha_0)/N)*K
    
    

    
  }
  

  lower=length(which(t<(-1/N)))
  K=lower
  h=apply(data[,which(t<(-1/N))]-temp_data[,which(t<(-1/N))],1,sum)
  alpha_0=alpha2[,1]
  l=h-K*alpha_0*1/N
  alpha2[,1]=2*alpha2[,1]+sum(alpha2[,2]*alpha2[,1])*alpha2[,2]+sum(alpha_c*alpha2[,1])*alpha_c+rho*l
  alpha2[,1]=alpha2[,1]/sum(alpha2[,1]*alpha2[,1])^0.5*v
  temp_data=temp_data-c((alpha2[,1]-alpha_0)/N)
  h=h+c((alpha2[,1]-alpha_0)/N)*K
  
  
  la<-c()
  
  for(i in (2):(N-1)){
    if(K>0 && i/N>= -t[K]){
      h=h-(data[,K]-temp_data[,K])
      K=K-1
    }
    if(K==0){
      
      alpha2[,i:N]=alpha2[,i-1]
      break;
    }
    alpha_0=alpha2[,i]
    l=h-K*alpha_0/N
    alpha2[,i]=2*alpha2[,i]+sum(alpha2[,i+1]*alpha2[,i])*alpha2[,i+1]+
      sum(alpha2[,i-1]*alpha2[,i])*alpha2[,i-1]-rho*l
    alpha2[,i]=alpha2[,i]/sum(alpha2[,i]*alpha2[,i])^0.5*v
    if(K>0){
      temp_data[,(101-K):100]=temp_data[,(101-K):100]-c((alpha2[,i]-alpha_0)/N)
    }
    h=h+c((alpha2[,i]-alpha_0)/N)*K
    
    la<-cbind(la,l)
    
  }
  
  h=apply(data[,which(t>(1/N))]-temp_data[,which(t>(1/N))],1,sum)
  h=h-apply(data[,which(t<(-1/N))]-temp_data[,which(t<(-1/N))],1,sum)
  
  alpha_0=alpha_c
  l=h+(length(which(t>(1/N))-length(which(t<(-1/N)))))*alpha_0*1/N
  alpha_c=2*alpha_c+sum(alpha1[,1]*alpha_c)*alpha1[,1]+sum(alpha2[,1]*alpha_c)*alpha2[,1]+rho*l
  alpha_c=alpha_c/sum(alpha_c*alpha_c)^0.5*v
  temp_data[,which(t>(1/N))]=temp_data[,which(t>(1/N))]+c((alpha_c-alpha_0)/N)
  temp_data[,which(t<(-1/N))]=temp_data[,which(t<(-1/N))]-c((alpha_c-alpha_0)/N)

  
  x_0=x_0+apply(data-temp_data,1,mean)
  temp_data=temp_data+apply(data-temp_data,1,mean)
  
  for(i in 1:100){
    
    for(K in 1:50){
      if(t[i]>1/N){
        prime=sum(alpha1[,floor(t[i]*N)]*(data[,i]-temp_data[,i]))
      }
      if(t[i]<(-1/N)){
        prime=sum(alpha2[,floor(abs(t[i])*N)]*(data[,i]-temp_data[,i]))
      }
      if(t[i]<=1/N & t[i]>=(-1/N)){
        prime=sum(alpha_c*(data[,i]-temp_data[,i]))
      }
      if(prime > 0){
        if(t[i]>1/N){
          temp_data[,i]=temp_data[,i]+alpha1[,floor(t[i]*N)]/N
         
        }
        if(t[i]<(-1/N)){
          
          temp_data[,i]=temp_data[,i]+alpha2[,floor(abs(t[i])*N)]/N
        }
        if(t[i]<=1/N & t[i]>=(-1/N)){
          
          temp_data[,i]=temp_data[,i]+alpha_c/N
        }
        
        t[i]=t[i]+1/N
      }
      if(prime < 0 ){
        if(t[i]>1/N){
          
          temp_data[,i]=temp_data[,i]-alpha1[,floor(t[i]*N)]/N
        }
        if(t[i]<(-1/N)){
          temp_data[,i]=temp_data[,i]-alpha2[,floor(abs(t[i])*N)]/N
          t
        }
        if(t[i]<=1/N & t[i]>=(-1/N)){
          temp_data[,i]=temp_data[,i]-alpha_c/N
          
        }
        
        t[i]=t[i]-1/N
        
        
      }
    }
  }
  
  "  
  for(i in 1:100){
  
  min=sum((data[,i]-temp_data[,i])^2)
  t[i]=floor(t[i]*N)/N
  old=data[,i]-temp_data[,i]
  temp=old
  for(K in 1:1){
  temp=temp-alpha[,floor(t[i]*N)-K]/N
  temp_err=sum((temp)^2)
  if(temp_err<min){
  min=temp_err
  t[i]=(floor(t[i]*N)-K)/N
  temp_data[,i]=data[,i]-temp
  }
  }
  temp=old
  for(K in 1:1){
  temp=temp+alpha[,floor(t[i]*N)+K-1]/N
  temp_err=sum((temp)^2)
  if(temp_err<min){
  min=temp_err
  t[i]=(floor(t[i]*N)+K)/N
  temp_data[,i]=data[,i]-temp
  }
  }
  }
  "  
  
  temp_data=temp_data[,order(t)]
  data=data[,order(t)]
  t=t[order(t)]
  
  K=length(which(t<=0))+1
  temp_x=x_0
  for( i in 2:N){
    temp_x=temp_x+alpha1[,i-1]/N
    if(i/N>=t[K]){
      temp_data[,K]=temp_x
      K=K+1
      if(K>100){break}
    }
  }
  K=length(which(t<0))
  temp_x=x_0
  for( i in 2:N){
    temp_x=temp_x-alpha2[,i-1]/N
    if(i/N>=abs(t[K])){
      temp_data[,K]=temp_x
      K=K-1
      if(K<1){break}
    }
  }
  x_0=x_0+apply(data-temp_data,1,mean)
  temp_data=temp_data+apply(data-temp_data,1,mean)
  
  
  plot(temp_data[1,],temp_data[2,],xlim=c(-1.5,1.5),ylim=c(-1.5,1.5),col=col)
  "
  x=matrix(0,2,N)
  x[,1]=x_0
  for(i in 2:(N)){
  x[,i]=x[,i-1]+alpha[,i-1]/N
  }
  plot(x[1,],x[2,],xlim=c(-1.5,1.5),ylim=c(-1.5,1.5))
  "
}
#plot(temp_data[1,],temp_data[2,],xlim=c(-1.5,1.5),ylim=c(0,1))

plot(data[1,],data[2,],xlim=c(-1.5,1.5),ylim=c(-1.5,1.5),col=col)

x=matrix(0,2,N)
x[,1]=x_0
for(i in 2:(N)){
  x[,i]=x[,i-1]+alpha[,i-1]/N
}
plot(x[1,],x[2,],xlim=c(-1.5,1.5),ylim=c(-1.5,1.5))
