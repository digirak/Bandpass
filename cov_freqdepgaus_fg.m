clear;
freqs_st=167e+6
df=80e+3
N=(197e+6-freqs_st)/df
freqs=freqs_st+[0:N-1]*df
freq_obs=freqs
freq_base=192.7e+6
dx=(600-10.)/N
x=10+[0:N-1]*dx
c=3.0e+8%%speed of light
alpha=1;%%fudge factor
beta=0.67;
gamma=0;
S_max=3
S=1
for j=1:1
	freq_obs(j)=freq_base
k=((x))*freq_base/c;%%observing frequency/central freq
f_nu=(freqs-freq_obs(j))/freq_base;
f_nu(find(f_nu==0))=df/freq_base
term1=alpha/(3-beta);
term2=(S_max).^(3-beta)/S.^(-beta); %%use this later
%%term2=(S/S_max).^beta;
epsilon=0.42;
D=4;

sigma=epsilon*c./(freqs*D);
term3=pi*c.^2*epsilon.^2./(D*(freqs.^2+freq_obs(j).^2));

denom=(4*(freqs.^2+freq_obs(j).^2)*D.^2)
terme=exp((-c.^2.*(transpose(f_nu)*(k)).^2*epsilon.^2)/(4*(freqs.^2+freq_obs(j).^2)*D.^2))

%C_fg(j,:)=term1.*term2.*term3.*terme
C_fg=term1*term2*(terme*term3)
end
%F=fft(C_fg,[],1)
%Fdag=ifft(C_fg,[],1)

f=fft(C_fg,[],2)
fdag=ifft(C_fg,[],2)

%kk=F*C_fg*Fdag
kkp=f*C_fg*fdag
imagesc(log(abs(kkp(:,220:350))))
set(gca,'YDir','normal')
colorbar
%for i=1:N
%	F(i,:)=fft(C_fg(i,:));
%	f(:,i)=fft(C_fg(:,i))
%	Fdag(i,:)=(ifft(C_fg(i,:)));
%	fdag(:,i)=ifft(C_fg(:,i));
%end
%%kk=Fdag*(C_fg)*F
blacknutt=0.3635819-0.4891775.*cos(2*pi.*[1:N]/(N-1))+0.1365995.*cos(4*pi.*[1:N]/(N-1))+0.0106411.*cos(6*pi.*[1:N]/(N-1))
win_rect=(fir1(N-1,0.03,'low',rectwin(N)))
win_black=(fir1(N-1,0.03,'low',blackman(N)));    % FIR filter
win_bn=(fir1(N-1,0.03,'low',blacknutt));    % FIR filter
win_kaiser=(fir1(N-1,0.020,'low',kaiser(N,10)));
%win_gauss=(fir1(N-1,0.020,'low',gausswin(N,20)))

win_gauss=(fir1(N-1,[0.012,0.020],'bandpass',gausswin(N,20)))
%out_rect=abs(ft).*(win_rect);
%out_black=abs(ft).*(win_black);
%out_bn=abs(ft).*(win_bn)
%%loglog(abs(out_black))
rect=freqz(win_rect)
black=freqz(win_black)
bn=freqz(win_bn)
kai=freqz(win_kaiser)
gau=freqz(win_gauss)

for i=1:N
	%rect_out(:,i)=kkp(:,i).*abs(rect(1:375))
	%black_out(:,i)=kkp(:,i).*abs(black(1:375))
	%bn_out(:,i)=kkp(:,i).*abs(bn(1:375))
	%kaiser_out(:,i)=kkp(:,i).*abs(kai(1:375))
	gauss_out(:,i)=kkp(:,i).*abs(gau(1:375))
	
end
save('gauss-kaiser.mat')

%ax=gca
%hold(ax,'on')
%loglog(abs(out_rect),'g')
%loglog(abs(out_bn),'r')
%legend('Blackman-Harris','rectangular','Blackman-nutall','location','SouthWest')

