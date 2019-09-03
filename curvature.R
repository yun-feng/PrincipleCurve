t=(1:100)/100

v=pi
x<-cos(v*t)
y<-sin(v*t)
data<-matrix(0,2,100)
data[1,]=x
data[2,]=y

x_0=c(0,0)
alpha_0=matrix(c(0,1),2,1)*v

temp_data=x_0+alpha_0%*%t



s=temp_data-x_0
s=(data-apply(data,1,mean))%*%t(s-apply(s,1,mean))
p=svd(s)$u%*%diag(-sign(svd(s)$d))%*%t(svd(s)$v)
alpha_0=p%*%alpha_0
alpha_0=-alpha_0/sum(alpha_0^2)^0.5*v

temp_data=p%*%(temp_data-x_0)
x_0=apply(data+temp_data,1,mean)

temp_data=x_0-temp_data
new=1000000
old=new+1
new_array<-c()

min=new
min_data=temp_data


for ( c in 1:101){
N=10000
alpha=matrix(0,2,N)
Force=matrix(0,2,N)
beta=matrix(0,2,N)
x=matrix(0,2,N)
alpha[,1]=alpha_0
x[,1]=x_0
Force[,1]=c(0,0)
Force_x=0
rho=0.8
K=1
for(i in 2:N){
  
  Force_prime=Force_x+1/rho*Force[,i-1]*(sum(Force[,i-1]*alpha[,i-1]))/(v^2)-alpha[,i-1]*(sum(Force[,i-1]*alpha[,i-1]))^2/(v^4)
  #Force_prime=Force_prime/rho
  Force[,i]=Force[,i-1]+Force_prime*1/N
  if(sum(Force[,i]*Force[,i])>1000000){
  Force[,i]=Force[,i]/sum(Force[,i]*Force[,i])^0.5*1000
  }
  alpha_prime=Force[,i-1]-alpha[,i-1]*(sum(Force[,i-1]*alpha[,i-1]))/(v^2)
  alpha_prime=alpha_prime/rho
  beta[,i]=alpha_prime
  alpha[,i]=alpha[,i-1]+alpha_prime*1/N
  alpha[,i]=alpha[,i]/sum(alpha[,i]*alpha[,i])^0.5*v
  x[,i]=x[,i-1]+alpha[,i-1]*1/N
  if(K<101 & i/N>=t[K]){
    Force_x=Force_x+data[,K]-x[,i]
    temp_data[,K]=x[,i]
    K=K+1
  }
  
}

#plot(temp_data[1,],temp_data[2,],xlim=c(-1.5,1.5),ylim=c(-0.5,1.5))

n=runif(2)
n=n/sqrt(sum(n^2))
ori=n
temp_l=temp_data-x_0
temp_h=data-x_0
S=-temp_l%*%t(temp_h)
n=max(abs(S))*n
n=n+S%*%ori-2*ori*sum(ori*S%*%ori)-t(S)%*%ori
n=n/sqrt(sum(n^2))
alpha_0=(diag(c(1,1))-2*n%*%t(n))%*%(diag(c(1,1))-2*ori%*%t(ori))%*%alpha_0



x_0=x_0+apply(data-temp_data,1,mean)
temp_data=temp_data+apply(data-temp_data,1,mean)



plot(temp_data[1,],temp_data[2,],xlim=c(-1.5,1.5),ylim=c(0,1))



}
#plot(temp_data[1,],temp_data[2,],xlim=c(-1.5,1.5),ylim=c(0,1))

plot(data[1,],data[2,],xlim=c(-1.5,1.5),ylim=c(0,1))








#####
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

t=t+0.5
x_0=x_0-alpha_0*0.5


N=1000
alpha=matrix(rep(c(alpha_0),2*N),2,N)
rho=100
plot(temp_data[1,],temp_data[2,],xlim=c(-1.5,1.5),ylim=c(-1.5,1.5))
for ( c in 1:1000){
  
  for(bc in 1:1){
  K=100
  h=apply(data-temp_data,1,sum)
  #l=h+K*alpha_0*1/N
  #alpha_0=alpha[,1]
  #alpha[,1]=2*alpha[,1]+sum(alpha[,2]*alpha[,1])*alpha[,2]+rho*l
  #alpha[,1]=alpha[,1]/sum(alpha[,1]*alpha[,1])^0.5*v
  #temp_data=temp_data+c((alpha[,1]-alpha_0)/N)
  #h=h-c((alpha[,1]-alpha_0)/N)*K
  
  
  la<-c()
  for(i in (floor(N*t[1])+1):(N-1)){
    if(K>0 && i/N>=t[101-K]){
      h=h-(data[,101-K]-temp_data[,101-K])
      K=K-1
    }
    if(K==0){
      alpha[,i:N]=alpha[,i-1]
      break;
    }
    alpha_0=alpha[,i]
    l=h+K*alpha_0/N
    alpha[,i]=-2*alpha[,i]+
              sum(alpha[,i+1]*alpha[,i])*alpha[,i+1]+
              sum(alpha[,i-1]*alpha[,i])*alpha[,i-1]+rho*l
    alpha[,i]=alpha[,i]/sum(alpha[,i]*alpha[,i])^0.5*v
    if(K>0){
      temp_data[,(101-K):100]=temp_data[,(101-K):100]+c((alpha[,i]-alpha_0)/N)
    }
    h=h-c((alpha[,i]-alpha_0)/N)*K
    
    la<-cbind(la,l)
    
  }

  alpha[,N]=alpha[,N-1]
  temp=alpha[,1:floor(N*t[1])]
  alpha[,1:floor(N*t[1])]=alpha[,floor(N*t[1])+1]
  x_0=x_0-sum((alpha[,1:floor(N*t[1])]-temp))/N
  
  
  x_0=x_0+apply(data-temp_data,1,mean)
  temp_data=temp_data+apply(data-temp_data,1,mean)
  }  
  for(i in 1:100){
    
    for(K in 1:50){
      prime=sum(alpha[,floor(t[i]*N)]*(data[,i]-temp_data[,i]))
      if(prime > 0 & t[i]<(1-2/N)){
        temp_data[,i]=temp_data[,i]+alpha[,floor(t[i]*N)]/N
        t[i]=t[i]+1/N
      }
      if(prime < 0 & t[i]>(3/N)){
        t[i]=t[i]-1/N
        temp_data[,i]=temp_data[,i]-alpha[,floor(t[i]*N)-1]/N
      
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
  
  n=runif(2)
  n=n/sqrt(sum(n^2))
  ori=n
  temp_l=temp_data-x_0
  temp_h=data-x_0
  S=-temp_l%*%t(temp_h)
  n=max(abs(S))*n
  n=n+S%*%ori-2*ori*sum(ori*S%*%ori)-t(S)%*%ori
  n=n/sqrt(sum(n^2))
  alpha=(diag(c(1,1))-2*n%*%t(n))%*%(diag(c(1,1))-2*ori%*%t(ori))%*%alpha
  
  
  
  
  K=1
  temp_x=x_0
  for( i in 2:N){
    temp_x=temp_x+alpha[,i-1]/N
    if(i/N>=t[K]){
      temp_data[,K]=temp_x
      K=K+1
      if(K>100){break}
    }
  }
  x_0=x_0+apply(data-temp_data,1,mean)
  temp_data=temp_data+apply(data-temp_data,1,mean)
  
  
  plot(temp_data[1,],temp_data[2,],xlim=c(-1.5,1.5),ylim=c(-1.5,1.5))
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

plot(data[1,],data[2,],xlim=c(-1.5,1.5),ylim=c(-1.5,1.5))

x=matrix(0,2,N)
x[,1]=x_0
for(i in 2:(N)){
  x[,i]=x[,i-1]+alpha[,i-1]/N
}
plot(x[1,],x[2,],xlim=c(-1.5,1.5),ylim=c(-1.5,1.5))

