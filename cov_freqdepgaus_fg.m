clear;

freqs=[1:100]
freq_obs=150.0%%MHz
df=10;
freqs=df*freqs;
x=randn(1,100)%random distribution of the x
c=3.0e+8%%speed of light
k=abs(fft(x))*c/freq_obs;%%observing frequency/central freq
alpha=10;%%fudge factor
beta=0.67;
gamma=0;
S_max=1
S=1
f_nu=(freqs-freq_obs)/freq_obs;
term1=alpha/(3-beta);
%%term2=(S_max).^(3-beta)/S.^(-beta); %%use this later
term2=(S/S_max).^beta;
epsilon=0.42;
D=4;
sigma=epsilon*c./(freqs*D);
term3=pi*c.^2*epsilon.^2./(D*(freqs+freqs.^2+freq_obs.^2));
terme=exp(-k.^2*c.^2.*f_nu.^2*epsilon.^2)./(4*(freqs.^2+freq_obs.^2)*D.^2)
C_fg=term1.*term2.*term3.*terme

