clear;
%77%%%%%%%%%%%%%%%
format long
freqs_st=150e+6
df=40e+3
N=round((170e+6-freqs_st)/df)
freqs=freqs_st+[0:N-1]*df
freq_base=160e+6
%freq_base=170.6e+6
%dx=(600-10.)/N
dx=1
x=4+[0:99]*dx
f21=1420e+6
base_lims=numel(x)
c=3.0e+8%%speed of light
alpha=1;%%fudge factor
beta=0.67;
gamma=0.7;
S_max=35
%S_max=50.*1e+6
S=1
freqs=(freqs)
freqsp=(freqs)
epsilon=0.42;
D=4;%% this was orignially 4m
scaler=1;
epsilon=epsilon*scaler;
BW=freqs(N)-freqs(1)
%D=3.5;
%S_max=606*S_max
%S=606*S
S_max=S_max.*(21e-2).^2.*1e-26./(2*1.38e-23*sqrt(2e-6))
S=S.*(21e-2).^2.*1e-26./(2.*1.38e-23*sqrt(2e-6))
u=((x)).*freq_base./c;%%observing frequency/central freq
%k=x
%k=x*0
term1=alpha/(3-beta);
term2=(S_max/(1)).^(3-beta)/(S*1).^(-1*beta)*(freqs/freq_base).^(-2*gamma); %%use this later

%term1=1
%term2=1
%%term2=(S/S_max).^beta;
%% section on constants and standard parameters done
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Actualloop for covariance matrix
parpool(4)
parfor i=1:N
f_numat(:,i)=(freqs-freqsp(i))./freq_base
%sigma(:,i)=epsilon*c./(sqrt(freqs.^2+freqsp(i).^2)*D)
end
%j=find(f_numat==0)%% this has to then have some sort f/f' dependence with gamma to include the spectral index
%f_numat(j)=1
%sigma=epsilon*c./(freqs*D);
parfor i=1:N
%term3(:,i)=pi.^2*c.^2*epsilon.^2./(D.^2*(freqs.^2+freqsp(i).^2));
term3(:,i)=2*pi.^1*c.^2*epsilon.^2./(D.^2*(freqs.^2+freqsp(i).^2));
end
%parpool(2)
parfor j=1:base_lims
%%terme is the exponential term in the freqdepgaus. f_numat is the channel to frequency covariance matrix
%%Terme 375x375 in nu-nu' space, sampled by each k(i) which is the baseline in units of wavelength.
%numer(:,:,j)=-c^.2.*(f_numat.*k(j))^.2.*epsilon^.2
%numer(:,:,j)=-c^2.*(f_numat.*k(j))^2.*epsilon^2
numer(:,:,j)=-1*c.^2.*(f_numat.*u(j)).^2.*epsilon.^2
end
parfor i=1:N
%denom(:,i)=(4.*(freqs.^2+freqsp(i).^2)./(freq_base*D).^2)
denom(:,i)=(1.*(freqs.^2+freqsp(i).^2).*(1*D).^2)
end
parfor i=1:base_lims%%loop over u
	terme(:,:,i)=exp(numer(:,:,i)./denom(:,:))
end
parfor i=1:base_lims
for j=1:N
expterm(:,j,i)=term3(:,j).*terme(:,j,i)
end
end
for i=1:base_lims
parfor j=1:N
	C_fg(:,j,i)=term1.*transpose(term2).*(expterm(:,j,i))
end
end
extension=0
blacknutt=0.3635819-0.4891775.*cos(2*pi.*[1:(N)]./(N-1))+0.1365995.*cos(4*pi.*[1:(N)]./(N-1))-0.0106411.*cos(6*pi.*[1:(N)]./(N-1))
bn=abs(freqz(blacknutt))
kai=kaiser(N+extension,10.8)
gau=gausswin(N+extension,4)
parfor iLoop=1:base_lims
%Ft_Rect=abs(fft2(C_fg(:,:,iLoop)))%%added abs
for j=1:N
	C_fgbn(:,j,iLoop)=C_fg(:,j,iLoop).*transpose(blacknutt)
	%C_fgbn(:,j,iLoop)=C_fg(:,j,iLoop).*transpose(blacknutt)
C_fgkai(:,j,iLoop)=C_fg(:,j,iLoop).*kaiser(N,10.8)
	C_fgbh(:,j,iLoop)=C_fg(:,j,iLoop).*blackmanharris(N,'symmetric')
%	C_fgkai(:,j,iLoop)=C_fg(:,j,iLoop).*kai(150:349)
%	C_fggau(:,j,iLoop)=C_fg(:,j,iLoop).*gausswin(N,4)
%	C_fggau(:,j,iLoop)=C_fg(:,j,iLoop).*gau(150:349)
	C_fgftop(:,j,iLoop)=C_fg(:,j,iLoop).*flattopwin(N,'symmetric')
end
%Ift_Rect=abs(conj(fft2(C_fg(:,:,iLoop))))./N%%added abs%% removed ifft
%kk(:,:,iLoop)=Ft_Rect*Ift_Rect
end
B=dftmtx(100)
Bdag=conj(dftmtx(100))./100
B_bn=BW*dftmtx(N)
Bdag_bn=1*conj(dftmtx(N))

%FFT loop for filters

parfor i=1:base_lims
	%ft_bn=abs(fft2(C_fgbn(:,:,i)))%%removed 2%%
	%ift_bn=abs(conj(fft2(C_fgbn(:,:,i))))./N%%removed 2%%removed ifft 

%
%	ft_kai=abs(fft2(C_fgkai(:,:,i)))
%	ift_kai=abs(ifft2(C_fgkai(:,:,i)))
%
%	ft_gau=abs(fft2(C_fggau(:,:,i)))
%	ift_gau=abs(ifft2(C_fggau(:,:,i)))
%	ft_ftop=abs(fft2(C_fgftop(:,:,i))).^2
%	ift_ftop=abs(ifft2(C_fgftop(:,:,i))).^2
%
%	kk_bn(:,:,i)=((ft_bn*C_fgbn(:,:,i)*ift_bn))
kk(:,:,i)=B_bn*C_fg(:,:,i)*Bdag_bn
kk_bn(:,:,i)=B_bn*C_fgbn(:,:,i)*Bdag_bn
%	kk_bn(:,:,i)=((ft_bn*ift_bn))
	%kk_kai(:,:,i)=((ft_kai*C_fgkai(:,:,i)*ift_kai))
	kk_kai(:,:,i)=B_bn*C_fgkai(:,:,i)*Bdag_bn
	kk_bh(:,:,i)=B_bn*C_fgbh(:,:,i)*Bdag_bn
%	kk_gau(:,:,i)=B_bn*C_fggau(:,:,i)*Bdag_bn
%	kk_gau(:,:,i)=((ft_gau*C_fggau(:,:,i)*ift_gau))
	kk_ftop(:,:,i)=((B_bn*C_fgftop(:,:,i)*Bdag_bn))
end
%%Choose along variances

parfor i=1:N
kkp(i,:)=kk(i,i,:)
end
parfor i=1:N
kkp_bn(i,:)=kk_bn(i,i,:)
kkp_kai(i,:)=kk_kai(i,i,:)
kkp_bh(i,:)=kk_bh(i,i,:)
%kkp_gau(i,:)=kk_gau(i,i,:)
kkp_ftop(i,:)=kk_ftop(i,i,:)
end

%%fold
parfor i=1:(N)/2
	kkpl(i,:)=((abs(kkp(i,:))+abs(kkp(N-i,:)))./2)./(2e-6)
end
parfor i=1:N./2
%kkpl(i,:)=abs(kkp(i,:))./2
%kkpl_bn(i,:)=abs(kkp_bn(i,:))./2
	kkpl_bn(i,:)=((abs(kkp_bn(i,:))+abs(kkp_bn(N-i,:)))./2)./2e-6
%	kkpl_gau(i,:)=(abs(kkp_gau(i,:))+abs(kkp_gau(N-i,:)))./2
	kkpl_kai(i,:)=(abs(kkp_kai(i,:)+abs(kkp_kai(N-i,:)))./2)./2e-6
	kkpl_bh(i,:)=(abs(kkp_bh(i,:)+abs(kkp_bh(N-i,:)))./2)./2e-6
%	kkpl_gau(i,:)=abs(kkp_gau(i,:))./2
	kkpl_ftop(i,:)=(abs(kkp_ftop(i,:))+abs(kkp_ftop(N-i,:))./2)./2e-6
end
%kkpl=kkpl*1e+6
%kkpl_bn=kkpl_bn*1e+6
%kkpl_kai=kkpl_kai*1e+6
%kkpl_bh=kkpl_bh*1e+6
%kkpl_ftop=kkpl_ftop*1e+6
BW=freqs(N)-freqs(1)
z=(f21/freq_base)-1
etahalf=[-0.5/df:1/(BW):0.5/df]
eta=2*etahalf((numel(etahalf)./2+1):numel(etahalf))
Ez=sqrt(0.27*(1+z).^3+0.73)
Dz=c*(1+z).^2./(Ez*f21*70e+3)
kpar=2*pi*eta./(Dz)
kperp=2*pi*u./(Dz*freq_base)
%%convert to 1d
%for i=1:numel(kpar)
%for i=2:2:numel(kperp()
%	kperpbins(i)=sum(kperp(i-1:i))./2
%	kperpbins(i)=sqrt(kperp(i-1).^2+kperp(i).^2)
%end
%for i=1:50
%kbold(i)=sqrt(kpar(round(i./2)).^2+kperp(i).^2)
%kbold(i)=sqrt(kpar(i).^2+kperp(i).^2)
%end
%Pkmax=Pkmax./base_lims
%Pkregion1=kkpl(2:18,6:28)
%Pkregion2=kkpl(35:249,50:99)

%Pkregion1_bn=kkpl_bn(2:18,6:28)
%Pkregion2_bn=kkpl_bn(35:249,50:99)

%Pkregion3=kkpl(120:250,1:30)
%Pkregion3_bn=kkpl_bn(120:250,1:30)
%for i=1:numel(kpar)
%mu(i)=kpar(i)./kbold(i)
%end
%mu_avg=sum(mu)/50


%parfor i=1:numel(Pkregion1(:,1))
%	Pk_1(i)=Pkregion1(i,1)
%	Pk_1_bn(i)=Pkregion1_bn(i,1)
%	kbold_1(i)=0
%	for j=1:numel(Pkregion1(1,:))
%		Pk_1(i)=sqrt(Pkregion1(i,j).^2+Pkregion1(i,j+1).^2)
%		Pk_1(i)=sqrt(Pk_1(i).^2+Pkregion1(i,j).^2)
%		Pk_1_bn(i)=sqrt(Pk_1_bn(i).^2+Pkregion1_bn(i,j).^2)
%	kbold_1(i)=sqrt(kpar(1+i).^2+kperp(j+5).^2)
%	end
	%kbold_1(i)=kbold_1(i)./j
	%	Pk_1(i)=kbold_1(i).^3.*Pk_1(i)./(2*pi.^2)
%		Pk_1(i)=Pk_1(i)./j
%		Pk_1_bn(i)=Pk_1_bn(i)./j
%	%	Pk_1_bn(i)=kbold_1(i).^3.*Pk_1_bn(i)./(2*pi.^2)
%end
%Pk_1=Pk_1./numel(Pkregion1(1,:))
%Dfac1=(max(Pk_1)-Pk_1)./max(Pk_1)
%Dfac1_bn=(max(Pk_1_bn)-Pk_1_bn)./max(Pk_1_bn)
%Dfac1_bn_rec=(max(Pk_1)-Pk_1_bn)./max(Pk_1)
%parfor i=1:numel(Pkregion2(:,1))
%	Pk_2(i)=Pkregion2(i,1)
%	Pk_2_bn(i)=Pkregion2_bn(i,1)
%	kbold_2(i)=0
%	for j=2:numel(Pkregion2(1,:))-1
%		Pk_2_bn(i)=sqrt(Pk_2_bn(i).^2+Pkregion2_bn(i,j).^2)
%		Pk_2(i)=sqrt(Pk_2(i).^2+Pkregion2(i,j).^2)
%		kbold_2(i)=sqrt(kpar(35+i).^2+kperp(29+j).^2)
%	end
%		kbold_2(i)=kbold_2(i)
	%	Pk_2(i)=kbold_2(i).^3.*Pk_2(i)./(2.*pi.^2)
%		Pk_2_bn(i)=Pk_2_bn(i)./j
%		Pk_2(i)=Pk_2(i)./j
%
	%	Pk_2_bn(i)=kbold_2(i).^3.*Pk_2_bn(i)./(2.*pi.^2)
%end

%Dfac2=(max(Pk_2)-Pk_2)./max(Pk_2)
%Dfac2_bn=(max(Pk_2_bn)-Pk_2_bn)./max(Pk_2_bn)
%Dfac2_bn_rec=(max(Pk_2)-Pk_2_bn)./max(Pk_2)
%		kbold_3(ouloop)=kperp(1)
%for ouloop=2:numel(Pkregion3(1,:))
%	parfor i=1:numel(Pkregion3(:,1))
%	kbold_3(i)=0
%		Pk_3(i)=Pkregion3(i,1)
%		Pk_3_bn(i)=Pkregion3_bn(i,1)
%		for j=2:numel(Pkregion3(1,:))
%			Pk_3_bn(i)=sqrt(Pk_3_bn(i).^2+Pkregion3_bn(i,j).^2)./1
%			Pk_3(i)=sqrt(Pk_3(i).^2+Pkregion3(i,j).^2)
%		       	kbold_3(i)=sqrt(kpar(119+i).^2+kperp(j).^2)
%end
%Pk_3(i)=Pk_3(i)./j
%Pk_3_bn(i)=Pk_3_bn(i)./j
%		kbold_3(i)=kbold_3(i)./numel(Pkregion3(1,:))
%end

	%	Pk_3(ouloop)=Pk_3(ouloop)./ouloop
	%	Pk_3_bn(ouloop)=Pk_3_bn(ouloop)./ouloop
%end

%Dfac3=((max(Pk_3)-Pk_3)./max(Pk_3))
%Dfac3_bn=(max(Pk_3_bn)-Pk_3_bn)./max(Pk_3_bn)
%Dfac3_bn_rec=(max(Pk_3)-Pk_3_bn)./max(Pk_3)

delete(gcp)
bn_smooth=smooth(bn,100)
drop_1=(max(bn_smooth(2:18))-bn_smooth(2:18))./max(bn_smooth(2:18))
drop_2=(max(bn_smooth(50:250))-bn_smooth(50:250))./max(bn_smooth(50:250))
drop_3=(max(bn_smooth(120:250))-bn_smooth(120:250))./max(bn_smooth(120:250))

%%generate cosmology
k=sqrt(kpar(1:numel(kperp)).^2+kperp.^2)
%xhmean=1./(1+exp(-(z-10)./0.5))
xhmean=0.50
deltax=k.^3.*exp(-k.^3./1.24)./(2.*pi.^2)
deltarho=(k./max(k)).^1
deltaxrho=deltax.*deltarho
%deltaxrho=0
Pk=(559.21).*(xhmean.^2.*deltarho+(xhmean-xhmean.^2).*deltaxrho+(xhmean-xhmean.^2).*deltax)
for i=3:numel(k)-3
	Pk_rec(i-2)=sqrt((sum(kkpl(i-2:i:i+2,i)./i)./5).^2+(sum(kkpl(i,i-2:i:i+2)./5)./i).^2)
%max(kkpl_bn)=max(kkpl_bn).*2.75
%	Pk_bn(i-2)=sqrt((sum(kkpl_bn(1:i,i)./5)./i).^2+(sum(kkpl_bn(i,1:i)./5)./i).^2)
	Pk_bn(i-2)=sqrt((sum(kkpl_bn(i-2:i:i+2,i)./5)./i).^2+(sum(kkpl_bn(i,i-2:i:i+2)./5)./i).^2)
%	Pk_gau(i-2)=sqrt((sum(kkpl_gau(i-2:i:i+2,i)./5)./i).^2+(sum(kkpl_gau(i,i-2:i:i+2)./5)./i).^2)
	Pk_bh(i-2)=sqrt((sum(kkpl_bh(i-2:i:i+2,i)./5)./i).^2+(sum(kkpl_bh(i,i-2:i:i+2)./5)./i).^2)
	Pk_ftop(i-2)=sqrt((sum(kkpl_ftop(i-2:i:i+2,i)./5)./i).^2+(sum(kkpl_ftop(i,i-2:i:i+2)./5)./i).^2)
	Pk_kai(i-2)=sqrt((sum(kkpl_kai(i-2:i:i+2,i)./5)./i).^2+(sum(kkpl_kai(i,i-2:i:i+2)./5)./i).^2)


end
Dfac_rec=(max(Pk_rec)-Pk_rec)./(max(Pk_rec))
Dfac_bn=(max(Pk_bn)-Pk_bn)./(max(Pk_bn))
%Dfac_kai=(max(Pk_kai)-Pk_kai)./(max(Pk_kai))
%Dfac_bh=(max(Pk_bh)-Pk_bh)./(max(Pk_bh))
%Dfac_ftop=(max(Pk_ftop)-Pk_ftop)./(max(Pk_ftop))
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

plot(log10(1./R1))
hold on
plot(log10(1./R2),'r')
%plot(log10(1./R3),'g')
%plot(log10(1./R4),'m')
%for i=1:40
%Pk1d(i)=sqrt(1.87^8*sum(kkpl(1:round(i./2),i)).^2+sum(kkpl(round(i./2),1:i)).^2)
%Pk1d(i)=sqrt(sum(kkpl(20:i+20,i).^2)+sum(kkpl(i+20,1:i).^2))
%end

%for i=1:40
%Pk1d_bn(i)=sqrt(1.87.^8.*sum(kkpl_bn(1:round(i./2),i)).^2+sum(kkpl_bn(round(i./2),1:i)).^2)
%Pk1d_bn(i)=sqrt(sum(kkpl_bn(20:i+20,i).^2)+sum(kkpl_bn(i+20,1:i).^2))
%end

%Dfact_rec=(max(Pk1d)-Pk1d)./max(Pk1d)
%Dfact_bn=(max(Pk1d_bn)-Pk1d_bn)./max(Pk1d_bn)
%save('/home/rakesh/Code/8MHzruns/withfilts_run191to199Mhz.mat')
%save('Gaussbeam/beam1600.mat')
