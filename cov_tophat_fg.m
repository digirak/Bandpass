clear;%%clear all
freqs_st=167e+6
df=80e+3
%N=100
N=(197e+6-freqs_st)/df
%N=(300e+6-freqs_st)/df
freqs=freqs_st+[0:N-1]*df
freq_obs=150e+6
%%freq_obs=fliplr(freq_obs)
%freq_base=[100:500]
%freq_base=(167+(197-167.)/2)*1e+6
%freq_base=150e+6
%dx=1*(600-10.)./N
%dx=1*(600-10.)./N
dx=1
x=10+[0:N-1]*dx
c=3.0e+8
freq_base=176.2e+6
k=(x.*(freq_base/c));%%observing frequency/central freq
alpha=1;%%fudge factor
beta=0.67;
gamma=0;
S_max=1
S=1;%%fudge factor 2
%for i=1:N
%freqsp(i)=freqs((i))
%end
freqsp=(freqs)
freqs=(freqs)
term1=alpha/(2-beta);
term2=(S_max).^(3-beta)/S.^(-beta); %%use this later
parpool(4)
parfor i=1:N
f_numat(:,i)=(freqs-freqsp(i))/freq_base
end
j=find(f_numat==0)
%f_numat(j)=df/freq_base
%f_numat
f_numat(j)=1
%for i=1:375
%f_numat(j(i))=term1*term2*(freqsp(i)/(freq_obs)).^(-gamma)
%f_numat(j(i))=(freqsp(i)/(freq_obs)).^(gamma)
%f_numat(j)=df/freq_base
%end
W=375
blacknutt=0.3635819-0.4891775.*cos(2*pi.*[1:N]/(N-1))+0.1365995.*cos(4*pi.*[1:N]/(N-1))-0.0106411.*cos(6*pi.*[1:N]/(N-1))
kai=kaiser(N,5)
gau=gausswin(N,2)
parfor i=1:N
u=k(i)
term3=besselj(1,(2*pi*(f_numat).*k(i)))./(f_numat.*k(i))%%sJ_1(2pif(nu)k)/kf(nu)
C_fg(:,:,i)=term1*term2*term3*(freqsp(i)/freq_obs)^(gamma)
end
%C_fg(:,:,i)=term3
parfor j=1:N
for iLoop=1:N
cfg_part=reshape(C_fg(:,iLoop,j),1,N)
C_fgbn(:,iLoop,j)=(cfg_part).*(blacknutt)

C_fgkai(:,iLoop,j)=cfg_part.*transpose(kai)
C_fggau(:,iLoop,j)=cfg_part.*transpose(gau)
end
end
%parfor i=1:N
%C_fgbn(W+1:N,:,:)=C_fg(W+1:N
parfor i=1:N
	Ft=fft2(C_fg(:,:,i))
	ft_bn=fft2(C_fgbn(:,:,i))
ft_kai=fft2(C_fgkai(:,:,i))
ft_gau=fft2(C_fggau(:,:,i))

%Ft=fft(C_fg(:,:,i),[],1)*fft(C_fg(:,:,i),[],2)
	Ftdag=ifft2(C_fg(:,:,i))
	ift_bn=ifft2(C_fgbn(:,:,i))
	ift_kai=ifft2(C_fgkai(:,:,i))
	ift_gau=ifft2(C_fggau(:,:,i))
%Ftdag=ifft(temp)
%Ftdag=ifft(C_fg(:,:,i),[],1)*ifft(C_fg(:,:,i),[],2)
	kk(:,:,i)=Ftdag*C_fg(:,:,i)*Ft
	kkbn(:,:,i)=ift_bn*C_fgbn(:,:,i)*ft_bn
kkkai(:,:,i)=ift_kai*C_fgkai(:,:,i)*ft_kai
kkgau(:,:,i)=ift_gau*C_fggau(:,:,i)*ft_gau
end
%kk=F*Cfgk*f
for i=1:N
kkp(i,:)=kk(i,i,:)
kkp_bn(i,:)=kkbn(i,i,:)
kkp_kai(i,:)=kkkai(i,i,:)
kkp_gau(i,:)=kkgau(i,i,:)
end
%kkpl=kkp(20:180,:)+kkp(181:341,:)
%kkpl=abs(kkp(1:(N-1)/2,:))+abs(kkp(((N-1)/2)+1:N-1,:))
%kkpl_bn=kkp_bn(1:(N-1)/2,:)+kkp_bn(((N-1)/2)+1:N-1,:)
%kkpl_kai=kkp_kai(1:(N-1)/2,:)+kkp_kai(((N-1)/2)+1:N-1,:)
%kkpl_gau=kkp_gau(1:(N-1)/2,:)+kkp_gau(((N-1)/2)+1:N-1,:)
parfor i=1:(N-1)/2
	kkpl(i,:)=abs(kkp(i,:))+abs(kkp(N-i,:))
	kkpl_bn(i,:)=abs(kkp_bn(i,:))+abs(kkp_bn(N-i,:))

	kkpl_kai(i,:)=abs(kkp_kai(i,:))+abs(kkp_kai(N-i,:))
	kkpl_gau(i,:)=abs(kkp_gau(i,:))+abs(kkp_gau(N-i,:))
end
klab={k(30),k(160),k(330)}

%win_rect=(fir1(N-1,[0.0001,0.12],'bandpass',rectwin(N)))
%win_black=(fir1(N-1,[0.01,0.12],'bandpass',blackman(N)));    % FIR filter
%win_bn=(fir1(N-1,[0.0004,0.12],'bandpass',blacknutt));    % FIR filter
%win_kaiser=(fir1(N-1,[0.00010,0.12],'bandpass',kaiser(N,10)));    % FIR filter
%win_gauss=(fir1(N-1,[0.0004,0.10],'bandpass',gausswin(N,20)))
%rect=freqz(win_rect,512,750)
%black=freqz(win_black,512,750)
%bn=freqz(win_bn,512,750)
%kai=freqz(win_kaiser,512,750)
%gau=freqz(win_gauss,512,750)
%'Created Windows'
%parfor i=1:N
%	rect_out(:,i)=kkp(:,i).*abs(rect)
%	black_out(:,i)=kkp(:,i).*abs(black)
%	bn_out(:,i)=kkp(:,i).*abs(bn)
%kaiser_out(:,i)=kkp(:,i).*abs(kai)

%gauss_out(:,i)=kkp(:,i).*abs(gau)
%end

%delete(gcp)

%parpool(4)

%parfor j=1:N
%cfg_black(:,j)=C_fg(:,j).*blackman(N)
%cfg_bn(:,j)=C_fg(:,j).*(blacknutt(N))
%cfg_kai(:,j)=C_fg(:,j).*kaiser(N,5)
%cfg_gau(:,j)=C_fg(:,j).*gausswin(N,2)
%end
delete(gcp)
%kk_black=fft(cfg_black,[],2)*cfg_black*ifft(cfg_black,[],2)
%kk_bn=fft(cfg_bn,[],2)*cfg_bn*ifft(cfg_bn,[],2)
%kk_kai=fft(cfg_kai,[],2)*cfg_kai*ifft(cfg_kai,[],2)
%kk_gau=fft(cfg_gau,[],2)*cfg_gau*ifft(cfg_gau,[],2)
%eta=abs(fft(f_nu))
z=(1420.0e+6./freqs)-1
omega_m=0.73
omega_k=0.27
c_mpc=c*3.24077929e-23
H_0=70e+3*3.24077929e-23
f21=1420e+6
dum=fft(freqs)
eta=dum(1:(N-1)/2)+dum(((N-1)/2)+1:N-1)
etalab={abs(eta(10)),abs(eta(95)),abs(eta(145))}
%save('gamma0.mat')
%Ez=(omega_m*(1+z).^3+omega_k).^(0.5)
%denom=(c_mpc.*(1+z).^2)./(H_0*f21*Ez)
%Dz=denom.*f_nu.*freq_base
%kperp=k.*2*pi./(Dz)
%kpar=eta.*2*pi./denom
%save('tophat.mat')
%imagesc(log(abs(kkp(60:150,150:320))))
%out_rect=abs(ft).*(win_rect);
%out_black=abs(ft).*(win_black);
%out_bn=abs(ft).*(win_bn)
%loglog(abs(out_black))
%ax=gca
%hold(ax,'on')
%loglog(abs(out_rect),'g')
%loglog(abs(out_bn),'r')
%legend('Blackman-Harr6is','rectangular','Blackman-nutall','location','SouthWest')
