clear;
%77%%%%%%%%%%%%%%%
freqs_st=167e+6
df=40e+3
N=round((175e+6-freqs_st)/df)
freqs=freqs_st+[0:N-1]*df
freq_obs=150e+6
freq_base=170.6e+6
%dx=(600-10.)/N
dx=1
x=4+[0:100]*dx
base_lims=numel(x)
c=3.0e+8%%speed of light
alpha=1;%%fudge factor
beta=0.67;
gamma=0;
S_max=1
S=1
extension=300
blacknutt=0.3635819-0.4891775.*cos(2*pi.*[1:(N+extension)]./(N+extension))+0.1365995.*cos(4*pi.*[1:(N+extension)]./(N+extension))-0.0106411.*cos(6*pi.*[1:(N+extension)]./(N+extension))
freqs=(freqs).*blacknutt(150:349)
freqsp=(freqs).*blacknutt(150:349)
epsilon=0.42;
D=4;%% this was orignially 4m
scaler=1;
epsilon=epsilon*scaler;
%D=3.5;

k=((x))*freq_obs/c;%%observing frequency/central freq
term1=alpha/(3-beta);
term2=(S_max).^(3-beta)/S.^(-beta)*(freqs/freq_base).^(-gamma); %%use this later
%%term2=(S/S_max).^beta;
%% section on constants and standard parameters done
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Actualloop for covariance matrix
parpool(4)
parfor i=1:N
f_numat(:,i)=(freqs-freqsp(i))/freq_base
%sigma(:,i)=epsilon*c./(sqrt(freqs.^2+freqsp(i).^2)*D)
end
%j=find(f_numat==0)%% this has to then have some sort f/f' dependence with gamma to include the spectral index
%f_numat(j)=1
%sigma=epsilon*c./(freqs*D);
parfor i=1:N
%term3(:,i)=pi*c.^2*epsilon.^2.*freq_base^.2./(D.^2*(freqs.^2+freqsp(i).^2));
term3(:,i)=pi*c.^2*epsilon.^2.*freq_base^0./(D.^2*(freqs.^2+freqsp(i).^2));
end
%parpool(2)
parfor j=1:base_lims
%%terme is the exponential term in the freqdepgaus. f_numat is the channel to frequency covariance matrix
%%Terme 375x375 in nu-nu' space, sampled by each k(i) which is the baseline in units of wavelength.
%numer(:,:,j)=-c^.2.*(f_numat.*k(j))^.2.*epsilon^.2
%numer(:,:,j)=-c^2.*(f_numat.*k(j))^2.*epsilon^2
numer(:,:,j)=-c.^2.*(f_numat.*k(j)).^2.*epsilon^2./1
end
parfor i=1:N
%denom(:,i)=(4.*(freqs.^2+freqsp(i).^2)./(freq_base*D).^2)
denom(:,i)=(4.*(freqs.^2+freqsp(i).^2).*(1*D).^2)
end
parfor i=1:base_lims%%loop over u
	terme(:,:,i)=exp(numer(:,:,i)./denom(:,:))
end
parfor i=1:base_lims
expterm(:,:,i)=term3(:,:)*terme(:,:,i)
end
for i=1:base_lims
parfor j=1:N
	C_fg(:,j,i)=term1.*transpose(term2).*(expterm(:,j,i))
end
end
kai=kaiser(N+extension,10.8)
gau=gausswin(N+extension,4)
parfor iLoop=1:base_lims
%Ft_Rect=abs(fft2(C_fg(:,:,iLoop)))%%added abs
%end
end

B=dftmtx(N)
Bdag=conj(dftmtx(N))./N
%FFT loop for filters
parfor i=1:base_lims
%	ft_bn=abs(fft2(C_fgbn(:,:,i)))%%removed 2%%
%	ift_bn=abs(ifft2(C_fgbn(:,:,i)))%%removed 2%%removed ifft 


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
kk(:,:,i)=B*C_fg(:,:,i)*Bdag
%kk_bn(:,:,i)=B*C_fgbn(:,:,i)*Bdag
%	kk_bn(:,:,i)=((ft_bn*C_fgbn(:,:,i)*ift_bn))
%%Choose along variances
end
for i=1:N
kkp(i,:)=kk(i,i,:)
end
%%fold
parfor i=1:(N)/2
	kkpl(i,:)=(abs(kkp(i,:))+abs(kkp(N-i,:)))./2
end
BW=freqs(N)-freqs(1)
z=7.6
etahalf=[-0.5/df:1/BW:0.5/df]
eta=2*etahalf((N/2)+1:N)
Ez=sqrt(0.27*(1+z).^2+0.73)
Dz=c*(1+z).^1./(Ez*1420e+6*70e+3)
kpar=2*pi*eta./(Dz*(1+z))
kperp=2*pi*x./(Dz*BW)
delete(gcp)
%save('/home/rakesh/Code/8MHzruns/withfilts_run191to199Mhz.mat')
%save('Gaussbeam/beam1600.mat')

