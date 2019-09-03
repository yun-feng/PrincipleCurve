t1=(1:100)/100

#t=floor(rnorm(100,mean=0.5,sd=0.2)*100)/100
#t[1:50]=floor(rnorm(50,mean=0.3,sd=0.1)*100)/100
#t[51:100]=floor(rnorm(50,mean=0.6,sd=0.1)*100)/100

v=1*pi
x<-cos(v*t1)
y<-sin(v*t1)
x<-c(x,0.5+cos(v*t1))
y<-c(y,1.5-sin(v*t1))

#x=v*t1
#y=v*t1
#x<-c(x,x)
#y<-c(y,y+0.5)
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

for (reg in c(10)){
for ( i in 1:100){


x1_0=apply(data%*%diag(z),1,sum)/sum(z)
alpha1_0=(data-x1_0)%*%diag(z)%*%t(data-x1_0)%*%alpha1_0
alpha1_0=alpha1_0/sum(alpha1_0*alpha1_0)^0.5

x2_0=apply(data%*%diag(1-z),1,sum)/sum(1-z)
alpha2_0=(data-x2_0)%*%diag(1-z)%*%t(data-x2_0)%*%alpha2_0
alpha2_0=alpha2_0/sum(alpha2_0*alpha2_0)^0.5

Er1=reg*(apply((data-x1_0)^2,2,sum)-(t(data-x1_0)%*%alpha1_0)^2)
Er1=Er1[,1]
Er2=reg*(apply((data-x2_0)^2,2,sum)-(t(data-x2_0)%*%alpha2_0)^2)
Er2=Er2[,1]

z=1/(1+exp(Er1-Er2))

}
}
