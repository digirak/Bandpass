clear;%%clear all
N=10000
freqs=[1:N]
freq_obs=150.0%%MHz
df=1e+6;
freqs=df*freqs;
x=(randn(1,N))%random distribution of the x
%%l=acos(1:360)
%%m=sqrt(abs(abs(x).^2-abs(l).^2))
c=3.0e+8%%speed of light
k=abs((x))*c/freq_obs;%%observing frequency/central freq
alpha=10;%%fudge factor
beta=0.67;
gamma=0;
S_max=1
S=1;%%fudge factor 2
f_nu=(freqs-freq_obs)/freq_obs;
term1=alpha/(3-beta);
%%term2=(S_max).^(3-beta)/S.^(-beta); %%use this later
term2=(S/S_max).^beta;
term3=besselj(1,2*pi*k.*f_nu)./(f_nu.*k)%%sJ_1(2pif(nu)k)/kf(nu)
C_fg=term1.*term2.*term3
ft=fft(C_fg)

blacknutt=0.3635819-0.4891775.*cos(2*pi.*[1:N]/(N-1))+0.1365995.*cos(4*pi.*[1:N]/(N-1))+0.0106411.*cos(6*pi.*[1:N]/(N-1))
win_rect=(fir1(N-1,[0.25,0.75],'bandpass',rectwin(N)))
win_black=(fir1(N-1,[0.25,0.75],'bandpass',blackman(N)));    % FIR filter
win_bn=(fir1(N-1,[0.25,0.75],'bandpass',blacknutt));    % FIR filter
out_rect=C_fg.*win_rect;
out_black=C_fg.*win_black;
out_bn=C_fg.*win_bn

