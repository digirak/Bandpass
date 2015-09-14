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
dx=3
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
S_max=S_max.*(21e-2).^2.*1e-26./(2*1.38e-23*(2e-6))
S=S.*(21e-2).^2.*1e-26./(2.*1.38e-23*(2e-6))
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
%parpool(4)
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
B_bn=BW*dftmtx(N)*2e-6
Bdag_bn=1*conj(dftmtx(N))

parfor i=1:base_lims
kk(:,:,i)=B_bn*C_fg(:,:,i)*Bdag_bn
end

parfor i=1:N
kkp(i,:)=kk(i,i,:)
end

parfor i=1:(N)/2
	kkpl(i,:)=((abs(kkp(i,:))+abs(kkp(N-i,:)))./2)
end
delete(gcp)
%save('covariance.mat')