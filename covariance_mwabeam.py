import numpy as np
import matplotlib as mpl
from scipy.fftpack import fft2
from scipy import special as sp
from scipy import signal as sig
#mpl.rcParams['text.usetex']=True
#mpl.rcParams['text.latex.unicode']=True
#from matplotlib import pyplot as plt
from mwa_new_beam import mwa_beam
from multiprocessing import Pool
from itertools import repeat

freqs_st=150e+6
#PI=np.pi
#freq resolution in Hz
df=200e+3
N=int((170e+6-freqs_st)/df)
freqs_new=freqs_st+np.arange(-1,N+1)*df
#baseline resolution in m
freqs=freqs_new[1:N+1]
#freqs[1]=freqs[200]=freqs[100]=df
freq_base=150e+6
dx=1
#x=4+np.arange(180)*dx
x=np.arange(1600)
#x=7+sp.jv(0,np.arange(300))*dx
base_lims=x.size
c=3e+8
alpha=4100
beta=1.6
gamma=0.78
S_max=1
S=1
freqsp=freqs
epsilon=0.42
D=4.2
BW=freqs[N-1]-freqs[0]
#S_max=S_max*np.power(21e-2,2)*1e-26/(2*1.38e-23*(2e-6))
#S=S*np.power(21e-2,2)*1e-26/(2*1.38e-23*(2e-6))
u=x*freq_base/c
term1=alpha/(3-beta)
term2=(S_max**(3-beta)/S**(-beta))*np.power((freqs/freq_base),(-1*gamma))
#f_numat=np.zeros(N*N).reshape(N,N)

C_fg=np.zeros(N*N*base_lims/2).reshape(N,N,base_lims/2)
beam=np.empty(N*base_lims*base_lims,dtype=complex).reshape(base_lims,base_lims,N)
#beam_r=np.zeros(base_lims*base_lims).reshape(base_lims,base_lims)
#beam_i=np.zeros(base_lims*base_lims).reshape(base_lims,base_lims)
#for i in range(freqs.size):
#    beam_r[:,:]=np.real(mwa_beam(base_lims,freqs[i]))
 #   beam_i[:,:]=np.imag(mwa_beam(base_lims,freqs[i]))
 #   beam[:,:,i]=beam_r+np.exp(-((beam_i))*1J)
#del beam_r,beam_i
#C_fgbh=np.zeros(N*N*base_lims).reshape(N,N,base_lims)
#for i in range(N):
p=Pool(6)
beam=p.map(mwa_beam,zip(repeat(base_lims),freqs),chunksize=1)
beam=np.asarray(beam)
p.terminate()
beam=np.transpose(beam)
beam=beam*term2*term1
#C_fg4d=np.empty(base_lims*base_lims*N*N).reshape(base_lims,base_lims,N,N)
#for uLoop in range(base_lims):
 #   for vLoop in range(base_lims):
  #      for nuLoop in range(N):
            #4d 'cube' of C_fg, here we have l,m,nu,nu'
   #         C_fg4d[uLoop,vLoop,nuLoop,:]=beam[uLoop,vLoop,nuLoop]*np.conj(beam[uLoop,vLoop,:])
#del beam
C_fgft=np.zeros(base_lims*base_lims).reshape(base_lims,base_lims)
temp=np.zeros(base_lims)
temp1=np.zeros(base_lims/2)
for nuLoop in range(N):
    for nupLoop in range(N):
        #selecting all (l,m) using nu and nu' parameters

 #       C_fgft[:,:]=abs(fft2(C_fg4d[:,:,nupLoop,nuLoop]))/(1)
        C_fgft[:,:]=abs(fft2(beam[:,:,nuLoop]*np.conj(beam[:,:,nupLoop])))/(1)

        for i in range(base_lims):
 #          temp[i]=C_fgft[i,i]
            temp[i]=np.sqrt((np.sum((C_fgft[i,0:i]/(i+1))**2))+(np.sum((C_fgft[0:i,i]/(i+1))**2)))
        for i in range(base_lims/2):
         #   temp1[base_lims/2-i-1]=(temp[i]+temp[base_lims-i-1])/2
            C_fg[nuLoop,nupLoop,:]=temp[1:base_lims/2+1]
#del C_fg4d,C_fgft,temp,temp1,W,omega
del C_fgft,temp,temp1,W,omega
C_fgbn=np.zeros(N*N*base_lims/2).reshape(N,N,base_lims/2)
C_fgbh=np.zeros(N*N*base_lims/2).reshape(N,N,base_lims/2)
i,j=np.meshgrid(np.arange(N),np.arange(N))
omega=np.exp(-2*PI*1J/N)
W=np.power(omega,i*j)

    

kk_bn=np.zeros(N*N*base_lims/2).reshape(N,N,base_lims/2)
kk_bh=np.zeros(N*N*base_lims/2).reshape(N,N,base_lims/2)

kkp=np.zeros(N*base_lims/2).reshape(N,base_lims/2)
kkp_bn=np.zeros(N*base_lims/2).reshape(N,base_lims/2)
kkp_bh=np.zeros(N*base_lims/2).reshape(N,base_lims/2)

kkpl=np.zeros(int(N/2)*base_lims/2).reshape(int(N/2),base_lims/2)
kkpl_bn=np.zeros(int(N/2)*base_lims/2).reshape(int(N/2),base_lims/2)
kkpl_bh=np.zeros(int(N/2)*base_lims/2).reshape(int(N/2),base_lims/2)

#blacknutt=0.3635819-0.4891775*np.cos(2*PI*np.arange(N)/(N-1))+0.1365995*np.cos(4*PI*np.arange(N)/(N-1))-0.0106411*np.cos(6*PI*np.arange(N)/(N-1))
blacknutt=sig.nuttall(N)
#blacknutt=0.3635-0.4891*np.cos(2*PI*np.arange(N)/N-1)+0.1365995*np.cos(4*PI*np.arange(N)/N-1)
#blacknutt=0.3635819-0.4891775*np.cos(2*PI*np.arange(-shift,N-shift)/(N-1))+0.1365995*np.cos(4*PI*np.arange(-shift,N-shift)/(N-1))-0.0106411*np.cos(6*PI*np.arange(-shift,N-shift)/(N-1))
#blacknutt=0.3635819-0.4891775*np.cos(2*PI*f_numat/(N-1))+0.1365995*np.cos(4*PI*f_numat/(N-1))-0.0106411*np.cos(6*PI*f_numat/(N-1))
#BH=0.35875-0.48829*np.cos(2*PI*np.arange(N)/(N-1))+0.14128*np.cos(4*PI*np.arange(N)/(N-1))-0.01168*np.cos(6*PI*np.arange(N)/(N-1))
BH=sig.blackmanharris(N)
Gau=sig.gaussian(N,std=9.5)
#BH=0.375-0.5*np.cos(2*PI*np.arange(N)/(N-1))+0.125*np.cos(4*PI*np.arange(N)/(N-1))#-0.012604*np.cos(6*PI*np.arange(N)/(N-1))


#for iLoop in range(N):
 #   for fLoop in rane(N):
      #  denom[:,iLoop]=4*(freqs**2+freqsp[iLoop]**2)*D**2

#iLoop=0
          #  C_fgbn[:,fLoop,uLoop]=C_fg[:,fLoop,uLoop]*blacknutt

   #         C_fgbh[:,fLoop,uLoop]=C_fg[:,fLoop,uLoop]*BHa
#
## converting C_fg from Jy^2 to K^2
Jy_to_K=1e-26*(21e-2)**2/(2*1.38e-23*2e-6)
#C_fgbn=C_fg
C_fg=C_fg*Jy_to_K**2
kk=np.zeros(N*N*base_lims/2).reshape(N,N,base_lims/2)
for base in range(base_lims/2):
   kk[:,:,base]=BW*W.dot(C_fg[:,:,base]).dot(np.conj(W))*2e-6
 #   kk[:,:,base]=BW*np.fft.fft2(C_fg[:,:,base])*2e-6
#uLoop=fLoop=0
for uLoop in range(base_lims/2):
    for fLoop in range(N):
       # C_fgbn[fLoop,fLoop,uLoop]=C_fgbn[fLoop,fLoop,uLoop]*(blacknutt[fLoop])
        C_fgbn[:,fLoop,uLoop]=(C_fg[:,fLoop,uLoop]*(blacknutt))
      #  temp=C_fg[:,fLoop,uLoop
  #      C_fgbn[:,fLoop,uLoop]=C_fg[:,fLoop,uLoop]*(blacknutt)
        C_fgbh[:,fLoop,uLoop]=C_fg[:,fLoop,uLoop]*(Gau)
#        C_fgbn[:,fLoop,uLoop]=temp*(blacknutt)
#uLoop=fLoop=0
for uLoop in range(base_lims/2):
      for fLoop in range(N):
         C_fgbn[fLoop,:,uLoop]=C_fgbn[fLoop,:,uLoop]*(blacknutt)
         C_fgbh[fLoop,:,uLoop]=C_fgbh[fLoop,:,uLoop]*(Gau)
#uLoop=0
#del C_fgbn1

#i,j=np.meshgrid(np.arange(N),np.arange(N))
#omega=np.exp(-2*PI*1J/N)
#W=np.power(omega,i*j)
#W=BW*W*2e-6

for uLoop in range(base_lims/2):
    #kk[:,:,uLoop]=(BW*(abs(np.fft.fft2(C_fg[:,:,uLoop])))*2e-6)
 #   kk[:,:,uLoop]=abs(BW*W.dot(C_fg[:,:,uLoop]).dot(np.conj(W))*2e-6)
    kk_bn[:,:,uLoop]=abs(BW*W.dot(C_fgbn[:,:,uLoop]).dot(np.conj(W))*2e-6)#remember to multiply by Sr-Hz
    kk_bh[:,:,uLoop]=abs(BW*W.dot(C_fgbh[:,:,uLoop]).dot(np.conj(W))*2e-6)
 
  #  kk_bn[:,:,uLoop]=(BW*(abs(np.fft.fft2(C_fgbn[:,:,uLoop])))*2e-6)
  #  kk_bh[:,:,uLoop]=BW*abs(np.fft.fft2(C_fgbh[:,:,uLoop]))*2e-6
 #   kk[:,:,uLoop]=BW*abs(np.fft.fft2(C_fg[:,:,uLoop])*np.conj(np.fft.fft2(C_fg[:,:,uLoop])))*2e-6
#N=460-40
#del C_fgbn,C_fg,C_fgbh
for fLoop in range(N):
       kkp[fLoop,:]=kk[fLoop,fLoop,:]
       kkp_bn[fLoop,:]=kk_bn[fLoop,fLoop,:]
       kkp_bh[fLoop,:]=kk_bh[fLoop,fLoop,:]

#del kk,kk_bh,kk_bn
#del kk_bh,kk_bn
for fLoop in range(int(N/2)):
    kkpl[fLoop,:]=((abs(kkp[fLoop,:])+abs(kkp[N-fLoop-1,:]))/2)*1
    kkpl_bn[fLoop,:]=((abs(kkp_bn[fLoop,:])+abs(kkp_bn[N-fLoop-1,:]))/2)*1
    kkpl_bh[fLoop,:]=((abs(kkp_bh[fLoop,:])+abs(kkp_bh[N-fLoop-1,:]))/2)*1
del kkp,kkp_bn#,kkp_bh
u_new=np.zeros(u.size/2)
for i in range(1,u.size/2):
    u_new[i]=(u[i]+u[u.size-i])/2
BW=freqs[N-1]-freqs[0]
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
k=np.zeros(kpar.size-rad*2)
Pk=np.zeros(kpar.size-rad*2)
Pk_bn=np.zeros(kpar.size-rad*2)
Pk_bh=np.zeros(kpar.size-rad*2)

for i in range(rad,kpar.size-rad):
    k[i-rad]=np.sqrt(np.sum(kpar[i-rad:i+rad]**2)/5+np.sum(kperp[i-rad:i+rad]**2)/5)

for i in range(rad,kpar.size-rad):
    Pk[i-rad]=k[i-rad]**3*np.sqrt(np.sum(kkpl[i-rad:i+rad,i]**2)/5+np.sum(kkpl[i,i-rad:i+rad]**2)/5)
    Pk_bn[i-rad]=k[i-rad]**3*np.sqrt(np.sum(kkpl_bn[i-rad:i+rad,i]**2)/5+np.sum(kkpl_bn[i,i-rad:i+rad]**2)/5)
    Pk_bh[i-rad]=k[i-rad]**3*np.sqrt(np.sum(kkpl_bh[i-rad:i+rad,i]**2)/5+np.sum(kkpl_bh[i,i-rad:i+rad]**2)/5)

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


