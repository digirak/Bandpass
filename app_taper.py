from scipy import signal as sig
## Read in fits first
hdlist=fits.open('C_fg.fits')
C_fg=hdlist[0].data
##Get the C_fg and then get limits
BW=20e+6
A=2e-6
N=len(C_fg[:,0,0])
base_lims=len(C_fg[0,0,:])
print 'baselines=',base_lims,'channels=',N
kk=np.empty(N*N*base_lims).reshape(N,N,base_lims)
kk_gau=np.empty(N*N*base_lims).reshape(N,N,base_lims)
C_fgbn=np.empty(N*N*base_lims).reshape(N,N,base_lims)
kkp=np.empty(N*base_lims).reshape(N,base_lims)
kkp_gau=np.empty(N*base_lims).reshape(N,base_lims)
kkpl=np.empty(N/2*base_lims).reshape(N/2,base_lims)
kkpl_gau=np.empty(N/2*base_lims).reshape(N/2,base_lims)
##DFT matrix
i,j=np.meshgrid(np.arange(N),np.arange(N))
omega=np.exp(-2*PI*1J/N)
W=np.power(omega,i*j)
BN=sig.blackman(N,sym=True)
Gau=sig.gaussian(N,8.5)
for i in range(base_lims):
    for j in range(N):
        C_fgbn[:,j,i]=C_fg[:,j,i]*BN
        C_fg[:,j,i]=C_fg[:,j,i]*Gau

for i in range(base_lims):
    for j in range(N):
        C_fgbn[j,:,i]=BN*C_fgbn[j,:,i]
        C_fg[j,:,i]=C_fg[j,:,i]*Gau

for i in range(base_lims):
    kk[:,:,i]=BW*abs(W.dot(C_fgbn[:,:,i]).dot(np.conjugate(W)))*(A)
    kk_gau[:,:,i]=BW*abs(W.dot(C_fg[:,:,i]).dot(np.conjugate(W)))*(A)

for i in range(N):
    kkp[i,:]=kk[i,i,:]
    kkp_gau[i,:]=kk_gau[i,i,:]
del kk

for fLoop in range(int(N/2)):
    kkpl[fLoop,:]=((abs(kkp[fLoop,:])+abs(kkp[N-fLoop-1,:]))/2)*1
    kkpl_gau[fLoop,:]=((abs(kkp_gau[fLoop,:])+abs(kkp_gau[N-fLoop-1,:]))/2)*1

del kkp
df=40e+3
u=np.arange(base_lims/2)
z=(1420e+6/170.6e+6)-1
etahalf=np.arange(-0.5/df,0.5/df,1/(BW))
eta=2*etahalf[etahalf.size/2+1:etahalf.size]
Ez=np.sqrt(0.27*(1+z)**3+0.73)
Dz=c*(1+z)**2/(Ez*1420e+6*70e+3)
kpar=2*PI*eta/(Dz)
kperp=2*PI*u[1:base_lims/2+1]/(Dz*170.6e+6)
#fig,ax=plt.subplots()

## generate 21 cm Power spectra
rad=3
k=np.zeros(N/2-rad*2)
Pk=np.zeros(N/2-rad*2)
Pk_bn=np.zeros(N/2-rad*2)
Pk_bh=np.zeros(N/2-rad*2)

for i in range(rad,N/2-rad):
    k[i-rad]=np.sqrt(np.sum(kpar[i-rad:i+rad]**2)/5+np.sum(kperp[i-rad:i+rad]**2)/5)

for i in range(rad,N/2-rad):
    Pk[i-rad]=k[i-rad]**3*np.sqrt(np.sum(kkpl[i-rad:i+rad,i]**2)/5+np.sum(kkpl[i,i-rad:i+rad]**2)/5)
    Pk_bn[i-rad]=k[i-rad]**3*np.sqrt(np.sum(kkpl_gau[i-rad:i+rad,i]**2)/5+np.sum(kkpl_gau[i,i-rad:i+rad]**2)/5)
 #   Pk_bh[i-rad]=k[i-rad]**3*np.sqrt(np.sum(kkpl_bh[i-rad:i+rad,i]**2)/5+np.sum(kkpl_bh[i,i-rad:i+rad]**2)/5)

xhmean=0.5
deltax=k**3*np.exp(-k**3/1.24)/(2*PI**2)
deltarho=(k/np.max(k))
deltaxrho=deltax*deltarho
P21=(559.21)*(xhmean**2*deltarho+(xhmean-xhmean**2)*deltaxrho+(xhmean-xhmean**2)*deltax)

## get into 2d
theta=np.arctan(kpar/kperp[0:kpar.size])
P21kperp=P21*np.cos(theta[0:P21.size])
P21kpar=P21*np.sin(theta[0:P21.size])
P212d=np.zeros(P21.size*P21.size).reshape(P21.size,P21.size)
for i in range(P21kperp.size):
    P212d[:,i]=(2*PI**2)*P21kpar[:]*P21kperp[i]/(k**3)


