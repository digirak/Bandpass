import numpy as np
import matplotlib as mpl
from scipy.fftpack import fft2
from scipy import special as sp
from scipy import integrate as intg
mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True
from matplotlib import pyplot as plt
import sys

freqs_st=150e+6
PI=np.pi
#freq resolution in Hz
df=40e+3
N=int((160e+6-freqs_st)/df)
freqs=freqs_st+np.arange(N)*df
#baseline resolution in m
freq_base=150e+6/2
dx=5
x=7+np.arange(300)*dx
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
term2=(S_max**(3-beta)/S**(-beta))*np.power((freqs/freq_base),(-2*gamma))
f_numat=np.zeros(N*N).reshape(N,N)
for ouLoop in range(N):
            f_numat[:,ouLoop]=(freqs-freqsp[ouLoop])/freq_base

#plt.imshow(f_numat,cmap=plt.cm.rainbow)
#plt.colorbar()
#plt.show()
theta=np.deg2rad(np.arange(-90,90))
phi=np.deg2rad(90)

kx=np.zeros(theta.size*freqs.size).reshape(theta.size,freqs.size)
ky=np.zeros(theta.size*freqs.size).reshape(theta.size,freqs.size)

xn=np.arange(1,16)
yn=np.arange(17,32)

fxreal=np.zeros(theta.size*freqs.size).reshape(theta.size,freqs.size)
fxim=np.zeros(theta.size*freqs.size).reshape(theta.size,freqs.size)
for i in range(N):
    kx[:,i]=2*PI*freqs[i]*np.sin(theta)*np.cos(phi)/c
    ky[:,i]=2*PI*freqs[i]*np.sin(theta)*np.sin(phi)/c

for i in range(N):
    for j in range(theta.size):
        fxreal[j,i]=(1./16)*np.sum(np.real(np.exp(1J*(yn*ky[j,i]+xn*kx[j,i]))))
        fxim[j,i]=(1./16)*np.sum(np.imag(np.exp(1J*(yn*ky[j,i]+xn*kx[j,i]))))

Beampow=np.zeros(theta.size*freqs.size).reshape(theta.size,freqs.size)

for i in range(N):
    for j in range(theta.size):
        A=complex(fxreal[j,i],fxim[j,i])
        Beampow[j,i]=A*np.conj(A)

C_fg=np.zeros(N*N*base_lims).reshape(N,N,base_lims)
for uLoop in range(base_lims):
    for fLoop in range(N):
        for fpLoop in range(N):
   #         for tLoop in range(theta.size):
            res=intg.quad(lambda l:l*Beampow[:,fpLoop]*Beampow[:,fLoop]*sp.j0(2*PI*(u[uLoop]*l*f_numat[fLoop,fpLoop])),0,1)
            C_fg[fLoop,fpLoop,uLoop]=res[0]

#sys.exit()

#term3=np.zeros(N*N).reshape(N,N)

#for iLoop in range(N):
 #   term3[:,iLoop]=2*PI*c**2*epsilon**2./(D**2*(freqs**2+freqs[iLoop]**2))

#numer=np.zeros(N*N*base_lims).reshape(N,N,base_lims)

#terme=np.zeros(N*N*base_lims).reshape(N,N,base_lims)
#expterm=np.zeros(N*N*base_lims).reshape(N,N,base_lims)
#denom=np.zeros(N*N).reshape(N,N)

#C_fgbh=np.zeros(N*N*base_lims).reshape(N,N,base_lims)
#C_fgbn=np.zeros(N*N*base_lims).reshape(N,N,base_lims)

kk=np.zeros(N*N*base_lims).reshape(N,N,base_lims)
#kk_bn=np.zeros(N*N*base_lims).reshape(N,N,base_lims)
#kk_bh=np.zeros(N*N*base_lims).reshape(N,N,base_lims)

kkp=np.zeros(N*base_lims).reshape(N,base_lims)
#kkp_bn=np.zeros(N*base_lims).reshape(N,base_lims)
#kkp_bh=np.zeros(N*base_lims).reshape(N,base_lims)

kkpl=np.zeros(int(N/2)*base_lims).reshape(int(N/2),base_lims)
#kkpl_bn=np.zeros(int(N/2)*base_lims).reshape(int(N/2),base_lims)
#kkpl_bh=np.zeros(int(N/2)*base_lims).reshape(int(N/2),base_lims)

blacknutt=0.3635819-0.4891775*np.cos(2*PI*np.arange(N)/(N-1))+0.1365995*np.cos(4*PI*np.arange(N)/(N-1))-0.0106411*np.cos(6*PI*np.arange(N)/(N-1))
BH=0.35875-0.48829*np.cos(2*PI*np.arange(N)/(N-1))+0.14128*np.cos(4*PI*np.arange(N)/(N-1))-0.01168*np.cos(6*PI*np.arange(N)/(N-1))


#for iLoop in range(N):
 #   for fLoop in rane(N):
 #       denom[:,iLoop]=(freqs**2+freqsp[iLoop]**2)*D**2

iLoop=0
#for iLoop in range(base_lims):
 #   numer[:,:,iLoop]=-1*c**2*(f_numat*u[iLoop])**2*epsilon**2
 #   terme[:,:,iLoop]=np.exp(numer[:,:,iLoop]/denom[:,:])
#iLoop=uLoop=0
#for uLoop in range(base_lims):
 #     for fLoop in range(N):
  #          expterm[:,fLoop,uLoop]=term3[:,fLoop]*terme[:,fLoop,uLoop]
  #          C_fg[:,fLoop,uLoop]=term1*np.transpose(term2)*expterm[:,fLoop,uLoop]
          #  C_fgbn[:,fLoop,uLoop]=C_fg[:,fLoop,uLoop]*blacknutt

   #         C_fgbh[:,fLoop,uLoop]=C_fg[:,fLoop,uLoop]*BHa
#
## converting C_fg from Jy^2 to K^2
#Jy_to_K=1e-26*(21e-2)**2/(2*1.38e-23*2e-6)
#C_fgbn=C_fg
#C_fg=C_fg*Jy_to_K**2
uLoop=fLoop=0
#del numer, denom,expterm,terme,term3,term2
#for uLoop in range(base_lims):
 #   for fLoop in range(N):
       # C_fgbn[fLoop,fLoop,uLoop]=C_fgbn[fLoop,fLoop,uLoop]*(blacknutt[fLoop])
     #   C_fgbn[:,fLoop,uLoop]=(C_fg[:,fLoop,uLoop]*np.transpose(blacknutt))
      #  temp=C_fg[:,fLoop,uLoop
  #      C_fgbn[:,fLoop,uLoop]=C_fg[:,fLoop,uLoop]*(blacknutt)
   #     C_fgbh[:,fLoop,uLoop]=C_fg[:,fLoop,uLoop]*(BH)
#        C_fgbn[:,fLoop,uLoop]=temp*(blacknutt)
#uLoop=fLoop=0
#for uLoop in range(base_lims):
 #     for fLoop in range(N):
  #       C_fgbn[fLoop,:,uLoop]=C_fgbn[fLoop,:,uLoop]*(blacknutt)
   #      C_fgbh[fLoop,:,uLoop]=C_fgbh[fLoop,:,uLoop]*(BH)
uLoop=0
#del C_fgbn1

i,j=np.meshgrid(np.arange(N),np.arange(N))
omega=np.exp(-2*PI*1J/N)
W=np.power(omega,i*j)
#W=BW*W*2e-6

for uLoop in range(base_lims):
    #kk[:,:,uLoop]=(BW*(abs(np.fft.fft2(C_fg[:,:,uLoop])))*2e-6)
    kk[:,:,uLoop]=BW*abs(W.dot(C_fg[:,:,uLoop]).dot(np.conj(W)))*2e-6
 #   kk_bn[:,:,uLoop]=abs(BW*W.dot(C_fgbn[:,:,uLoop]).dot(np.conj(W))*2e-6)
  #  kk_bh[:,:,uLoop]=abs(BW*W.dot(C_fgbh[:,:,uLoop]).dot(np.conj(W))*2e-6)
 
  #  kk_bn[:,:,uLoop]=(BW*(abs(np.fft.fft2(C_fgbn[:,:,uLoop])))*2e-6)
  #  kk_bh[:,:,uLoop]=BW*abs(np.fft.fft2(C_fgbh[:,:,uLoop]))*2e-6
 #   kk[:,:,uLoop]=BW*abs(np.fft.fft2(C_fg[:,:,uLoop])*np.conj(np.fft.fft2(C_fg[:,:,uLoop])))*2e-6
#N=460-40
#del C_fgbn,C_fg,C_fgbh
for fLoop in range(N):
        kkp[fLoop,:]=kk[fLoop,fLoop,:]
   #     kkp_bn[fLoop,:]=kk_bn[fLoop,fLoop,:]
    #    kkp_bh[fLoop,:]=kk_bh[fLoop,fLoop,:]

#del kk,kk_bh,kk_bn

for fLoop in range(int(N/2)):
    kkpl[fLoop,:]=((abs(kkp[fLoop,:])+abs(kkp[N-fLoop-1,:]))/2)*1e+6
    #kkpl_bn[fLoop,:]=((abs(kkp_bn[fLoop,:])+abs(kkp_bn[N-fLoop-1,:]))/2)*1e+6
    #kkpl_bh[fLoop,:]=((abs(kkp_bn[fLoop,:])+abs(kkp_bh[N-fLoop-1,:]))/2)*1e+6
#del kkp,kkp_bn,kkp_bh

BW=freqs[N-1]-freqs[0]
z=(1420e+6/freq_base)-1
etahalf=np.arange(-0.5/df,0.5/df,1/(BW))
eta=2*etahalf[etahalf.size/2+1:etahalf.size]
Ez=np.sqrt(0.27*(1+z)**3+0.73)
Dz=c*(1+z)**2/(Ez*1420e+6*70e+3)
kpar=2*PI*eta/(Dz)
kperp=2*PI*u/(Dz*freq_base)
#fig,ax=plt.subplots()
plt.figure('BN')
plt.pcolormesh(np.log10(abs(kkpl_bn[0:kkpl_bn.shape[0]-1,:])),cmap='jet')
#h1=ax.pcolor(np.log10(abs(kkpl[1:250,:])),cmap='jet')
plt.axis([0.61,base_lims,0.61,249])
cb=plt.colorbar()
cb.set_label('$\log_{10}(P_{k}) \mbox{ mK}^{2} h^{-3} \mbox{Mpc}^{3}$')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$k_{\perp}\mbox{ }h\mbox{ Mpc}^{-1}$')
plt.ylabel('$k_{\parallel}\mbox{ }h\mbox{ Mpc}^{-1}$')
plt.yticks([1,10,100,249],(0.01,0.35,3.5,7.2))
plt.xticks([1,10,base_lims],(kperp[0],kperp[10],kperp[100]))
#ax.set_xticks([1,10,30])
#ax.set_yticks([1,10,99])
#ax.set_xticklabels([0.001,0.01,0.1])
#ax.set_yticklabels([0.1,0.3,3.5])
plt.title('With Taper')
#plt.savefig('/home/rakesh/Plots/Figure_1.svg',format='svg',dpi=1200)
plt.figure('REC')
plt.pcolor(np.log10(abs(kkpl[0:249,:])),cmap='jet')
plt.axis([0.61,base_lims,0.61,249])
cb1=plt.colorbar()
cb.set_label('$\log_{10}(P_{k}) \mbox{ mK}^{2} h^{-3} \mbox{Mpc}^{3}$')
plt.xscale('log')
plt.yscale('log')
plt.yticks([1,10,249],(0.01,0.3,1.0))
plt.xticks([1,10,base_lims],(0.001,0.01,0.1))
plt.xlabel('$k_{\perp}\mbox{ }h\mbox{ Mpc}^{-1}$')
plt.title('No Taper')
#ax.set_xticklabels([0.001,0.01,0.1])
#ax.set_yticklabels([0.1,0.3,3.5])

plt.show()

