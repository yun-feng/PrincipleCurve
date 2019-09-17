import numpy as np;

#load data into matrix N*P
data=np.loadtxt("")
#simulate data
data=np.zeros((100,2))
data[:,0]=np.cos(np.pi*np.array(range(100))/100.0)
data[:,1]=np.sin(np.pi*np.array(range(100))/100.0)

#linear separable case
#parameters sigma=0.01
#k1=5
#k2=0.5
#reg=40
#K=2 or 3
#or more strigently k2=1 reg=50
data[:,0]=np.array(range(100))/100.0
data[50:100,0]=np.array(range(50))/100.0
data[:,1]=2*np.array(range(100))/100.0
data[50:100,1]=-2*np.array(range(50))/100.0

#half circles
#parameters sigma=0.05
#k1=8
#k2=1
#reg=80
#K=2 or 3
data=np.zeros((200,2))
data[:100,0]=np.cos(np.pi*np.array(range(100))/100.0)
data[:100,1]=np.sin(np.pi*np.array(range(100))/100.0)
data[100:200,0]=np.cos(np.pi*np.array(range(100))/100.0)
data[100:200,1]=1-np.sin(np.pi*np.array(range(100))/100.0)

#number of datapoints
N=data.shape[0]
#dimension of data
P=data.shape[1]


#hyper parameters
#velocity for curve
vel=np.pi
#number of intervals each curve divided into
Interval=250
#curvature regularizer
rho_delta=2e-2*Interval
rho=(1/rho_delta)/2/Interval
#number of clusters
K=3
#Guassian mixture variance
#tradeoff between data fitting errors
sigma=0.1
#regularization parameters for model selection
#tradeoff between uniform distribution and laplacian distribution
k1=1e4
#tradeoff between cluster size differences
k2=2
#tradeoff between data fitting and number of curves
reg=20


#Variables
#origin for each curve
x_ori=np.random.random((P,K))
#slope for each curve
alpha=np.ones((Interval,P,K))*np.random.random((P,K))
alpha/=np.sqrt((np.square(alpha)).sum(1))[:,np.newaxis,:]/vel
#curve points position
x=np.zeros((Interval,P,K))+x_ori+np.cumsum(alpha,axis=0)*1.0/Interval
#clustering probabilities
z=np.random.random((N,K))
z/=z.sum(1)[:,np.newaxis]
#distance between datapoints and curve points
dist=np.zeros((N,Interval,K))
#pseudotime probabilities
pseudotime_prob=np.random.random((N,Interval,K))
pseudotime_prob/=pseudotime_prob.sum(1)[:,np.newaxis,:]

for cycle in range(1000):
	dist=np.square(data[:,np.newaxis,:,np.newaxis]-x[np.newaxis,:,:,:]).sum(2)
	
	softmin_val_t=dist.min(1)
	pseudotime_prob=np.exp((-dist+softmin_val_t[:,np.newaxis,:])/sigma)
	
	
	group_loss=-np.log(pseudotime_prob.sum(1))+softmin_val_t/sigma
	softmin_val_group=group_loss.min(1)
	group_loss-=softmin_val_group[:,np.newaxis]
	
	pseudotime_prob/=pseudotime_prob.sum(1)[:,np.newaxis,:]
	
	group_L2=np.sqrt((np.square(z)).sum(0))
	
	Loss=(z*group_loss).sum()+(z*np.log(z)).sum()+softmin_val_group.sum()-reg*np.log(1+k1*np.exp(-k2*group_L2)).sum()
	Loss*=sigma
	Loss+=rho_delta*np.square(np.diff(alpha,axis=0)).sum()
	print(Loss)
	
	group_reg=reg*(k1*np.exp(-k2*group_L2))/(1.0+k1*np.exp(-k2*group_L2))*k2/(2*group_L2)
	max_weight=group_reg.max()
	dual_varaible=group_reg*z-max_weight-max_weight*np.log(z)
	z=np.exp(-(group_loss+2*dual_varaible)/(1+2*max_weight))
	z/=z.sum(1)[:,np.newaxis]
	
	
	
	#optimize x_ori
	data_center=(data[:,:,np.newaxis]*z[:,np.newaxis,:]).sum(0)
	x_center=(x[np.newaxis,:,:,:]*z[:,np.newaxis,np.newaxis,:]*pseudotime_prob[:,:,np.newaxis,:]).sum((0,1))
	shift=(data_center-x_center)/(z.sum(0))
	x_ori+=shift
	x=x+shift
	
	#optimize alpha with Householder matrices
	for k in range(K):
		#pivot
		pivot_x=x_ori[:,k]
		for i in range(Interval):
			#Linear loss 
			S=-((x[np.newaxis,i:,:,k]-pivot_x)*z[:,np.newaxis,np.newaxis,k]*pseudotime_prob[:,i:,np.newaxis,k]).sum(1)
			S=(S[:,:,np.newaxis]*(data[:,np.newaxis,:]-pivot_x)).sum(0)#np.dot(S.transpose(),data)
			if i>0:
				S-=rho_delta*alpha[i,:,k,np.newaxis]*alpha[np.newaxis,i-1,:,k]
			#random direction
			n=np.random.random(P)
			n/=np.sqrt(np.square(n).sum())
			ori=np.copy(n)
			S=S-2*np.dot(n[:,np.newaxis]*n[np.newaxis,:],S)+S.transpose()-2*np.dot(S.transpose(),n[:,np.newaxis]*n[np.newaxis,:])
			n=np.sum(np.abs(S))*n+np.dot(S,n)
			n/=np.sqrt(np.square(n).sum())
			#Orthonomal matirx
			Rotate_P=-2*ori[:,np.newaxis]*ori[np.newaxis,:]-2*n[:,np.newaxis]*n[np.newaxis,:]+4*np.sum(ori*n)*ori[:,np.newaxis]*n[np.newaxis,:]
			alpha[i:,:,k]=alpha[i:,:,k]+np.dot(alpha[i:,:,k],Rotate_P)
			x[i:,:,k]=x[i:,:,k]+np.dot(x[i:,:,k]-pivot_x,Rotate_P)
			
			pivot_x=pivot_x+alpha[i,:,k]*1.0/Interval

np.savetxt("/cygdrive/c/pseudo/x1.txt",x[:,:,1], fmt='%.18e',delimiter='\t')
np.savetxt("/cygdrive/c/pseudo/x2.txt",x[:,:,2], fmt='%.18e',delimiter='\t')