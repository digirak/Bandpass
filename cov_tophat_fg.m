clear;%%clear all
freqs_st=167e+6
df=80e+3
%N=100
N=(197e+6-freqs_st)/df
freq_obs=freqs_st+[0:N-1]*df
%freq_obs=freq_obs+80e+3
%freqs=[1:N]*1e+6
%freqs=[100:355]
%freqs=[1:600]
%freq_obs=[100:355]*1e+6%MHza
freqs=freq_obs
%%freq_obs=fliplr(freq_obs)
%freq_base=[100:500]
%freq_base=(167+(197-167.)/2)*1e+6
%freq_base=150e+6
dx=1*(600-10.)./N
x=10+[0:N-1]*dx
for j=1:1%numel(freq_obs)
c=3.0e+8
freq_obs(j)=192.7e+6
%freq_base=freq_obs(j)
freq_base=freqs_st
k=x*freq_base/c;%%observing frequency/central freq
alpha=1;%%fudge factor
beta=0.67;
gamma=0;
S_max=1
S=1;%%fudge factor 2
%%if(freq_obs(j)==freqs(j))
%%	f_nu=1e-10
%%else
	f_nu=(freqs-freq_obs(j))/freq_base;
%	for i=1:N
%		if(freqs(i)==freq_obs(j))
%			f_nu(i)=80e+3/freq_base;
%		else
%			f_nu(i)=(freqs(i)-freq_obs(j))/freq_base
%		end
%	end
%end
f_nu(find(f_nu==0))=df/freq_base
term1=alpha/(3-beta);
term2=(S_max).^(3-beta)/S.^(-beta); %%use this later
term3=besselj(1,2.*pi.*transpose(f_nu)*k)./(transpose(f_nu)*k)%%sJ_1(2pif(nu)k)/kf(nu)

end
C_fg=term1*term2*term3
%C_fg(isnan(C_fg))=0.0

%for j=1:N 
%	F(j,:)=fft(C_fg(j,:));
%	f(:,j)=fft(C_fg(:,j))
%	Fdag(j,:)=(ifft((C_fg(j,:))));
%	fdag(:,j)=ifft(C_fg(:,j))
%end
%%F=fft(C_fg,[],1)
%%Fdag=ifft(C_fg,[],1)

f=fft(C_fg,[],2)
fdag=ifft(C_fg,[],2)

%%kk=F*C_fg*Fdag
kkp=f*C_fg*fdag
%imagesc(log(abs(kkp(50:375,200:350))))
%set(gca,'YDir','normal')
%colorbar
blacknutt=0.3635819-0.4891775.*cos(2*pi.*[1:N]/(N-1))+0.1365995.*cos(4*pi.*[1:N]/(N-1))+0.0106411.*cos(6*pi.*[1:N]/(N-1))
win_rect=(fir1(N-1,[0.093,0.27],'bandpass',rectwin(N)))
win_black=(fir1(N-1,[0.093,0.27],'bandpass',blackman(N)));    % FIR filter
win_bn=(fir1(N-1,[0.093,0.27],'bandpass',blacknutt));    % FIR filter
win_kaiser=(fir1(N-1,[0.10,0.23],'bandpass',kaiser(N,10)));    % FIR filter
win_gauss=(fir1(N-1,[0.1,0.23],'bandpass',gausswin(N,20)))
rect=freqz(win_rect)
black=freqz(win_black)
bn=freqz(win_bn)
kai=freqz(win_kaiser)
gau=freqz(win_gauss)
'Created Windows'
for i=1:N
%	rect_out(:,i)=kkp(:,i).*rect(1:375)
%	black_out(:,i)=kkp(:,i).*black(1:375)
%	bn_out(:,i)=kkp(:,i).*bn(1:375)
%kaiser_out(:,i)=kkp(:,i).*abs(kai(1:375))
	gauss_out(:,i)=kkp(:,i).*abs(gau(1:375))
end

save('tophat.mat')
imagesc(log(abs(kaiser_out(78:116,150:320))))
%out_rect=abs(ft).*(win_rect);
%out_black=abs(ft).*(win_black);
%out_bn=abs(ft).*(win_bn)
%loglog(abs(out_black))
%ax=gca
%hold(ax,'on')
%loglog(abs(out_rect),'g')
%loglog(abs(out_bn),'r')
%legend('Blackman-Harr6is','rectangular','Blackman-nutall','location','SouthWest')
