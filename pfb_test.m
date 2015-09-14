
tot_f=100;
df=256.*1.28/tot_f;
dt=1./655.36;
f=[1:100];
f=f*df
t=[1:4096]*dt;
for i=1:100
    wave(i,:)=sin(2*pi*f(i)*t);
end;
coeff=MwaPfbProtoFilterCoeff2009512x8/8192./32.
for i=1:100
	arr(i,:)=wave(i,:).*coeff;
end
arr_reshape=reshape(arr,100,8,512);
wsum=sum(arr_reshape,2)
    wsum_split=squeeze(wsum(:,1,:))
            for i=1:100
                ft(i,:)=ifft(wsum_split(i,:))
		ftp(i,:)=pow2db(abs(ft(i,:).^2))
            end
	    
            freqs=[1:256]*df
