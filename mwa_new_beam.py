import numpy as np
from scipy import constants
c=constants.c
PI=constants.pi
def mwa_beam(var):
    (n,f)=var
    l=np.arange(-1,1,2./n)
    m=np.arange(-1,1,2./n)
    theta=np.arccos(m)
    phi=np.arcsin(l)
    vn=np.cos(PI/2*np.cos(theta))/(np.sin(theta))
    freq=f
 #   for f in range(freqs.size):
    
    print freq
    kx=np.zeros(l.size)
    ky=np.zeros(l.size)
    f1_real=np.zeros(l.size*m.size).reshape(l.size,m.size)
    f1_im=np.zeros(l.size*m.size).reshape(l.size,m.size)
    kx=2*PI*freq*l*1/c
    ky=2*PI*freq*m*1/c
    dx=1.1
    yn=np.arange(2,-3,-1)*dx
    xn=np.arange(2,-3,-1)*dx
    vg=2*np.sin(2*0.3*PI*np.sin(theta))
        #vg=1
    for i in range(xn.size):
       for j in range(yn.size):
            for linds in range(l.size):
                f1_real[linds,:]+=np.real(vn*vg*np.exp(1J*((yn[j])*(ky)+(xn[i])*kx[linds])))
                f1_im[linds,:]+=np.imag(vn*vg*np.exp(1J*(((yn[j]))*ky+(xn[i])*kx[linds])))
        
        #f1=(f1_real/(len(xn)*len(yn)*(vn*vg).max()))*(np.cos(f1_im/(len(xn)*len(yn))))+np.sin(f1_im/(len(xn)*len(yn))*1J)
    f1=f1_real/(len(xn)*len(yn)*(vn*vg).max())*np.exp(1J*f1_im)
    del f1_real
    del f1_im
    return f1
