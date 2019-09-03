N=100
t<-rep((1:N)/N,2)
#t<-1:200/200
x1<-4*t
x2<-exp(10*(t-0.5))/(1+(exp(10*(t-0.5))))+runif(2*N)*0.2
x1[(N+1):(2*N)]=x1[(N+1):(2*N)]
x2[(N+1):(2*N)]=-2*x2[(N+1):(2*N)]+runif(N)*0.2


x1<-4*t
x1[101:200]=x1[101:200]
x2<-2*t^2#+rnorm(N,0,0.5)
x2[101:200]=-x2[101:200]

w<-matrix(rep(0,2*N*2*N),2*N,2*N)
for(i in 1:(2*N)){
  for(j in 1:(2*N)){
    w[i,j]=exp(1.5*(-(x1[i]-x1[j])^2-(x2[i]-x2[j])^2))
  }
}
d<-apply(w,1,sum)

delta_m<-diag(1/d)%*%w
di_L<-diag(delta_m)
diag(delta_m)<-0
l<-apply(delta_m,2,sum)
#l<-di_L

z1<-runif(2*N)
z2=1-z1


d<-apply(w,1,sum)


L=diag(rep(1,2*N))-diag(1/d)%*%w


L1=diag(z1)%*%L
L2=diag(z2)%*%L

res_1<-Re(eigen(L1)[[2]][,199])
res_2<-Re(eigen(L2)[[2]][,199])
res_1<-res_1-(sum(res_1)/2/N)
res_2<-res_2-(sum(res_2)/2/N)
res_1=res_1/(sqrt(sum(res_1^2)))
res_2=res_2/(sqrt(sum(res_2^2)))
plot(z1)
z1_2=z1
z2_2=z2

#h=2000
h=1000


rho=2
c=0
for(s in 1:1000){
  c=c+1
  
  res_1<-diag(z1)%*%L%*%res_1-res_1  
  res_2=diag(z2)%*%L%*%res_2-res_2  
  res_1<-res_1-(sum(res_1)/2/N)
  res_2<-res_2-(sum(res_2)/2/N)
  res_1=res_1/(sqrt(sum(res_1^2)))
  res_2=res_2/(sqrt(sum(res_2^2)))
  
  for(j in 1:(2*N)){
    s1=sum(delta_m[,j]*(h*(res_1-res_1[j])^2-log(z1)))/l[j]
    s2=sum(delta_m[,j]*(h*(res_2-res_2[j])^2-log(z2)))/l[j]
    alpha1=sum(delta_m[,j]*z1)/l[j]
    alpha2=sum(delta_m[,j]*z2)/l[j]
    grad=s1-s2-alpha1/z1[j]+alpha2/z2[j]+log(z1[j])-log(z2[j])
    if(grad<0){
      s=c(1-1e-4,1e-4)
    }else{s=c(1e-4,1-1e-4)}
    z1[j]=(1-2/(c+1))*z1[j]+2/(c+1)*s[1]
    z2[j]=(1-2/(c+1))*z2[j]+2/(c+1)*s[2]
    
  }
  
  
  "
  
  for(j in 1:(2*N)){
    s1=sum(delta_m[,j]*(h*(res_1-res_1[j])^2-log(z1_2)))/l[j]
    s2=sum(delta_m[,j]*(h*(res_2-res_2[j])^2-log(z2_2)))/l[j]
    
    z1[j]=(z1_2[j]^rho*exp(-s1))^(1/(1+rho))
    z2[j]=(z2_2[j]^rho*exp(-s2))^(1/(1+rho))
    temp=z1[j]+z2[j]
    z1[j]=z1[j]/temp
    z2[j]=z2[j]/temp
    
    alpha1=sum(delta_m[,j]*z1)/l[j]
    alpha2=sum(delta_m[,j]*z2)/l[j]
    z1_2[j]=(alpha1+rho*z1[j])/(1+rho)
    z2_2[j]=(alpha2+rho*z2[j])/(1+rho)
  }
"
  loss=sum(res_1*(diag(z1)%*%L%*%res_1))+
       sum(res_2*(diag(z2)%*%L%*%res_2))
  loss=loss*h
  
  for(j in 1:200){
    loss=loss+sum(delta_m[,j]*z1[j]*(log(z1[j])-log(z1)))
    loss=loss+sum(delta_m[,j]*z2[j]*(log(z2[j])-log(z2)))
  }
  
  "
  for(j in 1:200){
    loss=loss+sum(delta_m[,j]*z1[j]*(log(z1[j])-log(z1_2)))
    loss=loss+sum(delta_m[,j]*z2[j]*(log(z2[j])-log(z2_2)))
    loss=loss+z1[j]*(log(z1[j])-log(z1_2[j]))*rho*l[j]
    loss=loss+z2[j]*(log(z2[j])-log(z2_2[j]))*rho*l[j]
  }
  "
  print(loss)
       
  
  plot(z1)
  
}




plot(res_1[which(z1>0.5)])
color=rep("red",2*N)
color[which(z1<0.5)]="blue"
plot(x1,x2,col=color)



