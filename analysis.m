clear
format long
load('covariance.mat','BW','freqs','N','freq_base','f21','c','df','u')
BW=freqs(N)-freqs(1)
z=(f21/freq_base)-1
etahalf=[-0.5/df:1/(BW):0.5/df]
eta=2*etahalf((numel(etahalf)./2+1):numel(etahalf))
Ez=sqrt(0.27*(1+z).^3+0.73)
Dz=c*(1+z).^2./(Ez*f21*70e+3)
kpar=2*pi*eta./(Dz)
kperp=2*pi*u./(Dz*freq_base)

k=sqrt(kpar(1:numel(kperp)).^2+kperp.^2)
%xhmean=1./(1+exp(-(z-10)./0.5))
xhmean=0.50
deltax=k.^3.*exp(-k.^3./1.24)./(2.*pi.^2)
deltarho=(k./max(k)).^1
deltaxrho=deltax.*deltarho
%deltaxrho=0
Pk=(559.21).*(xhmean.^2.*deltarho+(xhmean-xhmean.^2).*deltaxrho+(xhmean-xhmean.^2).*deltax)

load('covariance.mat','kkpl')
load('Tapered.mat','kkpl_bn','kkpl_kai','kkpl_ftop','kkpl_bh')
for i=3:numel(k)-3
	Pk_rec(i-2)=sqrt((sum(kkpl(i-2:i:i+2,i)./i)./5).^2+(sum(kkpl(i,i-2:i:i+2)./5)./i).^2)
	Pk_bn(i-2)=sqrt((sum(kkpl_bn(i-2:i:i+2,i)./5)./i).^2+(sum(kkpl_bn(i,i-2:i:i+2)./5)./i).^2)
	Pk_bh(i-2)=sqrt((sum(kkpl_bh(i-2:i:i+2,i)./5)./i).^2+(sum(kkpl_bh(i,i-2:i:i+2)./5)./i).^2)
	Pk_ftop(i-2)=sqrt((sum(kkpl_ftop(i-2:i:i+2,i)./5)./i).^2+(sum(kkpl_ftop(i,i-2:i:i+2)./5)./i).^2)
	Pk_kai(i-2)=sqrt((sum(kkpl_kai(i-2:i:i+2,i)./5)./i).^2+(sum(kkpl_kai(i,i-2:i:i+2)./5)./i).^2)
end

Pa=Pk(3:97)+(k(3:97).^3.*Pk_rec./(2.*pi.^2))
Pb=Pk(3:97)+(k(3:97).^3.*Pk_bn./(2.*pi.^2))
Pc=Pk(3:97)+(k(3:97).^3.*Pk_ftop./(2.*pi.^2))
Pd=Pk(3:97)+(k(3:97).^3.*Pk_bh./(2.*pi.^2))
Pe=Pk(3:97)+(k(3:97).^3.*Pk_kai./(2.*pi.^2))

R1=(Pa-Pk(3:97))./Pk(3:97)
R2=(Pb-Pk(3:97))./Pk(3:97)
R3=(Pc-Pk(3:97))./Pk(3:97)
R4=(Pd-Pk(3:97))./Pk(3:97)
R5=(Pe-Pk(3:97))./Pk(3:97)

dP_rec=abs(diff(Pk_rec(5:2:95)))./((0.0729))
dP_bn=abs(diff(Pk_bn(5:2:95)))./((0.0729))
dP_kai=abs(diff(Pk_kai(5:2:95)))./((0.0729))
dP_bh=abs(diff(Pk_bh(5:2:95)))./((0.0729))
dP_ftop=abs(diff(Pk_ftop(5:2:95)))./((0.0729))

%%second derivative
dP_rec_sq=diff(dP_rec)./((0.036))
dP_bn_sq=diff(dP_bn)./(abs(0.036))
dP_kai_sq=diff(dP_kai)./(abs(0.036))
dP_bh_sq=diff(dP_bh)./(abs(0.036))
dP_ftop_sq=diff(dP_ftop)./(abs(0.036))

