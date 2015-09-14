#import numpy as np
#from matplotlib import pyplot as plt
#import matplotlib as mpl
from scipy import constants
#mpl.rcParams['text.usetex']=True
#mpl.rcParams['text.latex.unicode']=True
#import cmath
##define constants
c=constants.c
#c=3e+8
PI=constants.pi

#theta=np.deg2rad(np.arange(-90,90,1))
#phi=np.deg2rad(np.arange(-90,90,1))
l=np.arange(-1,1,1./90)
m=np.arange(-1,1,1./90)
#theta=PI/2-np.arcsin(m)
mask=np.zeros(l.size*m.size).reshape(l.size,m.size)
for i in range(l.size):
    for j in range(m.size):
        if((l[i]**2+m[j]**2)<=1):
            mask[i,j]=1

theta=np.arccos(l)-np.deg2rad(0)
#theta=(PI/2)-npa.arccos(l)
phi=np.arcsin(m)-np.deg2rad(0)
#kx=np.zeros(theta.size*phi.size).reshape(phi.size,theta.size)
kx=np.zeros(l.size)
ky=np.zeros(m.size)
#f1_real=np.zeros(theta.size*phi.size).reshape(phi.size,theta.size)
#f1_im=np.zeros(theta.size*phi.size).reshape(phi.size,theta.size)
f1_real=np.zeros(l.size*m.size).reshape(l.size,m.size)
f1_im=np.zeros(l.size*m.size).reshape(l.size,m.size)
#f1=np.zeros(theta.size*phi.size).reshape(phi.size,theta.size)
#freqs=np.arange(100e+6,102e+6,1e+6)
Pwr=np.zeros(theta.size*phi.size).reshape(phi.size,theta.size)
#beampwr=np.zeros(l.size*m.size*freqs.size).reshape(l.size,m.size,freqs.size)
#for f in range(freqs.size):
freq=150e+6
    #lam1=c/(2*150e+6)
#vn=np.cos(PI/2*np.cos(((theta+(PI/2)))))/np.sin(((theta+(PI/2))))
vn=np.cos(PI/2*np.cos(theta))/(np.sin(theta))
#vn=np.cos((2*PI*freq*0.74/c)*np.cos(theta))-np.cos(2*PI*freq*0.74/c)/(np.sin(theta))
#vn=np.sin(PI/2+theta)**2
#vn=vn
#vn=(np.cos(PI/2*np.sin(((phi))))/np.cos(phi))
#vn=np.sin(theta+PI/2)*np.cos(phi)
#vg=2*np.sin(2*0.3*PI*np.cos(theta)*np.sin(phi))
#    vg=
#vn=(vg*vn)
#for i in range(phi.size):
 #   for j in range(theta.size):
 #       kx[i,j]=2*PI*freq*np.sin(theta[j])*np.cos(phi[i])/c
kx=2*PI*freq*l*0.72/c
 #       ky[i,j]=2*PI*freq*np.sin(theta[j])*np.sin(phi[i])/c
ky=2*PI*freq*(m)*0.72/c

dx=1.1
yn=np.arange(-2,3,1)*dx
xn=np.arange(-2,3,1)*dx
#yn=np.zeros(5)
#posmat=np.mgrid[-2:2,-2:2]*dx

#for k in range(3):
for linds in range(l.size):
    for i in range(xn.size):
        for j in range(yn.size):
                #f1_real[:,linds]+=np.real(vn*vg*np.exp(1J*((yn[j]+yshift[k])*(ky)+(xn[i]+xshift[k])*kx[linds])))
            #if(xn[i]==0 or yn[j]==0):
             #   continue
            #else:
                #f1_real[linds,:]+=np.real(vn*vg*np.exp(1J*((yn[j])*(ky)+(xn[i])*kx[linds])))
                #f1_real[linds,:]+=np.real(vn*vg*np.exp(1J*((yn[j])*(ky)+(xn[i])*kx[linds])))
            f1_real[linds,:]+=np.real(vn*np.exp(-1J*((yn[j]*ky)+(xn[i]*kx[linds]))))
            f1_im[linds,:]+=np.imag(vn*np.exp(-1J*(yn[j]*ky)+(xn[i]*kx[linds])))
f1=(f1_real/(xn.size*yn.size*vn.max())*(np.cos(f1_im/(1))+np.sin(f1_im/(1))*1J))
#f1=(f1_real/(xn.size*yn.size))
Pwr=f1*(np.conjugate(f1))*mask
#beampwr[:,:,f]=Pwr


#plt.pcolormesh(np.log10(Pwr[:,:]),vmin=-6)
#cb=plt.colorbar()
#cb.set_label('$\log_{10}$(B(l,m,$\\nu$))')
#plt.xticks([0,45,90,135,180],(-1,-0.5,0,0.5,1))
#plt.yticks([0,45,90,135,180],(-1,-0.5,0,0.5,1))
#plt.xlabel('l')
#plt.ylabel('m')
#plt.title('$B(l,m,\\nu=100\mbox{ MHz})$')
#plt.show()





