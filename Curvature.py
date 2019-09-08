import numpy as np;

#load data into matrix N*P
data=np.loadtxt("")
#simulate data
data=np.zeros((100,2))
data[:,0]=np.cos(np.pi*np.array(range(100))/100.0)
data[:,1]=np.sin(np.pi*np.array(range(100))/100.0)


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
rho_delta=1e-2*N
rho=(1/rho_delta)/2/N
#number of clusters
K=2
#Guassian mixture variance
sigma=15

#Variables
#origin for each curve
x_ori=np.random.random((P,K))
#slope for each curve
alpha=np.ones((Interval,P,K))*np.random.random((P,K))
alpha/=np.sqrt((np.square(alpha)).sum(1))[:,np.newaxis,:]/vel
#curve points position
x=np.zeros((Interval,P,K))+x_ori+np.cumsum(alpha,axis=0)*1.0/N
#clustering probabilities
z=np.random.random((N,K))
z/=z.sum(1)[:,np.newaxis]
#distance between datapoints and curve points
dist=np.zeros((N,Interval,K))
#pseudotime probabilities
pseudotime_prob=np.random.random((N,Interval,K))
pseudotime_prob/=pseudotime_prob.sum(1)[:,np.newaxis,:]

for cycle in range(100):
	dist=(np.square(data[:,np.newaxis,:,np.newaxis]-x[np.newaxis,:,:,:])).sum(2)
	
	softmin_val_t=dist.min(1)
	pseudotime_prob=np.exp((-dist+softmin_val_t[:,np.newaxis,:])/sigma)
	z=np.log(pseudotime_prob.sum(1))-softmin_val_t
	softmin_val_z=z.max(1)
	z=np.exp(z-softmin_val_z[:,np.newaxis])
	Loss=-(softmin_val_z+np.log(z.sum(1))).sum()
	
	pseudotime_prob/=pseudotime_prob.sum(1)[:,np.newaxis,:]
	z/=z.sum(1)[:,np.newaxis]
	
	
	Loss+=rho_delta*np.square(np.diff(alpha,axis=0)).sum()
	print(Loss)
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
			#alpha[i:,:,k]=alpha[i:,:,k]+np.dot(alpha[i:,:,k],Rotate_P)
			#x[i:,:,k]=x[i:,:,k]+np.dot(x[i:,:,k]-pivot_x,Rotate_P)
			
			pivot_x=pivot_x+alpha[i,:,k]*1.0/N
	
	