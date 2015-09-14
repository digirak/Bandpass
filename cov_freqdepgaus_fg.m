clear;
%%%%%%%%%%%%%%%%
freqs_st=167e+6
df=80e+3
N=(197e+6-freqs_st)/df
freqs=freqs_st+[0:N-1]*df
freq_obs=150e+6
freq_base=170.6e+6
%dx=(600-10.)/N
dx=1
x=10+[0:N-1]*dx
c=3.0e+8%%speed of light
alpha=1;%%fudge factor
beta=0.67;
gamma=0;
S_max=1
S=1
freqsp=freqs

epsilon=0.42;
D=4;

k=((x))*freq_base/c;%%observing frequency/central freq
term1=alpha/(3-beta);
term2=(S_max).^(3-beta)/S.^(-beta); %%use this later
%%term2=(S/S_max).^beta;
%% section on constants and standard parameters done
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Actualloop for covariance matrix
parpool(5)
parfor i=1:N
f_numat(:,i)=(freqs-freqsp(i))/freq_base
end
j=find(f_numat==0)%% this has to then have some sort f/f' dependence with gamma to include the spectral index
f_numat(j)=1
sigma=epsilon*c./(freqs*D);
parfor i=1:N
term3(:,i)=pi*c.^2*epsilon.^2./(D*(freqs.^2+freqsp(i).^2));
end

blacknutt=0.3635819-0.4891775.*cos(2*pi.*[1:N]/(N-1))+0.1365995.*cos(4*pi.*[1:N]/(N-1))-0.0106411.*cos(6*pi.*[1:N]/(N-1))
parfor j=1:N
%%terme is the exponential term in the freqdepgaus. f_numat is the channel to frequency covariance matrix
%%Terme 375x375 in nu-nu' space, sampled by each k(i) which is the baseline in units of wavelength.
numer(:,:,j)=-c^.2.*(f_numat.*k(j))^.2.*epsilon^.2
denom(:,j)=(4.*(freqs.^2+freqsp(j).^2).*D^.2./(freq_base).^2)
end
parfor i=1:N%%loop over u
	terme(:,:,i)=exp(numer(:,:,i)./denom(:,:))
end
parfor i=1:N
	C_fg(:,:,i)=term1*term2*(terme(:,:,i)*term3(:,:))
end


parfor iLoop=1:N
%Ft_Rect=fft2(C_fg(:,:,iLoop))
for j=1:N
	C_fgbn(:,j,iLoop)=C_fg(:,j,iLoop).*transpose(blacknutt)
	C_fgkai(:,j,iLoop)=C_fg(:,j,iLoop).*kaiser(N,9.5)
	C_fggau(:,j,iLoop)=C_fg(:,j,iLoop).*gausswin(N,3.5)
end
%Ift_Rect=ifft2(C_fg(:,:,iLoop))
%kk(:,:,iLoop)=Ft_Rect*C_fg(:,:,iLoop)*Ift_Rect
end
%FFT loop for filters
parfor i=1:N
	ft_bn=abs(fft2(C_fgbn(:,:,i))).^2
	ift_bn=abs(ifft2(C_fgbn(:,:,i))).^2

	ft_kai=abs(fft2(C_fgkai(:,:,i))).^2
	ift_kai=abs(ifft2(C_fgkai(:,:,i))).^2

	ft_gau=abs(fft2(C_fggau(:,:,i))).^2
	ift_gau=abs(ifft2(C_fggau(:,:,i))).^2

	kk_bn(:,:,i)=((ft_bn*C_fgbn(:,:,i)*ift_bn))
	kk_kai(:,:,i)=((ft_kai*C_fgkai(:,:,i)*ift_kai))
	kk_gau(:,:,i)=((ft_gau*C_fggau(:,:,i)*ift_gau))
end
%%Choose along variances

for i=1:N
%kkp(i,:)=kk(i,i,:)
kkp_bn(i,:)=kk_bn(i,i,:)
kkp_gau(i,:)=kk_gau(i,i,:)
kkp_kai(i,:)=kk_gau(i,i,:)
end
%%fold
parfor i=1:(N-1)/2
%	kkpl(i,:)=abs(kkp(i,:))+abs(kkp(N-i,:))
	kkpl_bn(i,:)=abs(kkp_bn(i,:))+abs(kkp_bn(N-i,:))
	kkpl_kai(i,:)=abs(kkp_kai(i,:))+abs(kkp_kai(N-i,:))
	kkpl_gau(i,:)=abs(kkp_gau(i,:))+abs(kkp_gau(N-i,:))

end


delete(gcp)
save('Gaussbeam/run2_kaiser6gau3.mat')
