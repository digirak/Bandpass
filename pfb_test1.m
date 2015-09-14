clear; 
N=512; 
 D = 8;
 cof_bit = 12;

 b = fir1(N*D-1,1/N,'low',kaiser(N*D,5));    % FIR filter

 b = (1-2^-(cof_bit-1))*b/max(b);            % scale to less than 1 (signed)
 ROM_cof12b_512x8 = round(b*2^(cof_bit-1));  % 2009 Coeff
tot_freq=10000
df=512*1.28/tot_freq;%%df in Mhz
dt=1./655.35;%%dt in ns
freqs=[1:tot_freq]*df;%%Total frequency range of all sine waves
t=[1:4096]*dt;%%time duration of each sine wave
%coeff=ROM_cof12b_512x8/8192.;%%Coefficient of Window scaled
coeff=b;
parpool(4)
parfor i=1:tot_freq
	wave=sin(2*pi*freqs(i)*t);
	arr=wave.*coeff;
	arr_reshape=reshape(arr,512,8);%%reshape the array to get 512x8
	wsum=sum(arr_reshape,2);%%sum each of the 8 taps
	ft_temp=fft(wsum);%%take an IFFT which is 512 pts
	ft(i,:)=ft_temp(1:256)
	res(i,:)=pow2db(abs(ft_temp(1:256).^2));%%return only 256 unique points
	;%%print i
end
ft_conj=conj(ft)
ft_new=[ft,ft_conj]
%res_new=fft(res,512)
res_new=abs(fft(ft_new));
%%hold on
%%plot(freqs(300:500),res(300:500,10));
%%plot(freqs(300:500),res(300:500,11),'g');
%%plot(freqs(300:500),res(300:500,12),'r');
%%legend('Channel 10', 'channel 11','channel 12','location','SouthWest');
%%xlabel('Frequency in MHz');(
%%ylabel('Power in dB');
%%title('Plot of power in 3 channels');
%%hold off
%%filter section
%win_rect=fft(fir1(511,[0.25,0.75],'bandpass',rectwin(512)));
%win_black=((blackman(256)));
%win_rect=((rectwin(256)));
%win_black=fft(fir1(511,[0.25,0.75],'bandpass',blackman(512)));    % FIR filter
%win_nutt=(fir1(255,[0.25,0.75],'bandpass',nuttallwin(256)));    % FIR filter
%win_flat=fft(fir1(511,[0.25,0.75],'bandpass',flattopwin(512)));    % FIR filter
%win_hann=(fir1(255,[0.25,0.75],'bandpass',hanning(256)));    % FIR filter
%win_nutt=(nuttallwin(256))
%win_kai=fft(fir1(511,[0.25,0.75],'bandpass',kaiser(N,10)))
%win_gauss=fft(fir1(511,[0.25,0.75],'bandpass',gausswin(N,20)))
blacknutt=0.3635819-0.4891775.*cos(2*pi.*[1:512]/511)+0.1365995.*cos(4*pi.*[1:512]/511)+0.0106411.*cos(6*pi.*[1:512]/511)
%%win_bn=fft(fir1(511,[0.2,0.8],'bandpass',blacknutt))
parfor i=1:tot_freq
	%res_black(i,:)=(res_new(i,:)).*(win_black);
%	res_black(i,:)=res_new(i,:).*transpose(abs(freqz(blackman(512),512,512,'whole')))
%	res_rect(i,:)=(res_new(i,:)).*(win_rect);
res_rect(i,:)=res_new(i,:).*(abs(freqz(rectwin(512),512,[-pi:pi/255:pi+pi/255])))
	%res_nutt(i,:)=(res_new(i,:)).*win_nutt;
	%res_flattop(i,:)=(res_new(i,:)).*(win_flat);
	%res_hann(i,:)=res_new(i,:).*(win_hann)
	%res_bn(i,:)=(res_new(i,:)).*(win_bn);
	res_bn(i,:)=res_new(i,:).*(abs(freqz(blacknutt,512,[-pi:pi/255:pi+pi/255])))
	res_kai(i,:)=(res_new(i,:)).*(win_kai);
%res_kai(i,:)=res_new(i,:).*transpose(abs(fft(kaiser(N,12))))
	%res_gauss(i,:)=(res_new(i,:)).*(win_gauss);
%	res_gauss(i,:)=res_new(i,:).*transpose(abs(fft(gausswin(N,5))))
end
nu=abs(fft(freqs))
delete(gcp)


