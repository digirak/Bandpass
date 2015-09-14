clear
format long
load('covariance.mat','C_fg','B_bn','Bdag_bn','base_lims','N','kkpl','freqs','freq_base','f_numat')
t=freqs/freq_base

blacknutt=0.3635819-0.4891775.*cos(2*pi.*[1:(N)]./(N-1))+0.1365995.*cos(4*pi.*[1:(N)]./(N-1))-0.0106411.*cos(6*pi.*[1:(N)]./(N-1))
bn3t=0.375-0.5.*cos(2*pi.*[1:(N)]./(N-1))+0.125.*cos(4*pi.*[1:(N)]./(N-1))%-0.0106411.*cos(6*pi.*[1:(N)]./(N-1))

parpool(4)
parfor iLoop=1:base_lims
for j=1:N
%blacknutt=0.3635819-0.4891775.*cos(2*pi.*f_numat(:,j)/(N-1))+0.1365995.*cos(4*pi.*f_numat(:,j)/(N-1))-0.0106411.*cos(6*pi.*f_numat(:,j)./(N-1))a
%BN=transpose(blacknutt).*f_numat(:,j)
%	C_fgbn(:,j,iLoop)=((C_fg(:,j,iLoop))).*transpose(blacknutt)
	%C_fgbn(:,:,iLoop)=transpose(blacknutt)*(C_fg(:,:,iLoop)*(blacknutt))
%C_fgkai(:,j,iLoop)=C_fg(:,j,iLoop).*kaiser(N,10.8)
	%C_fgbh(:,j,iLoop)=C_fg(:,j,iLoop).*blackmanharris(N,'symmetric')
	%C_fg3t(:,j,iLoop)=C_fg(:,j,iLoop).*flattopwin(N,'symmetric')
	C_fg3t(:,j,iLoop)=C_fg(:,j,iLoop).*transpose(bn3t)
end
end
parfor iLoop=1:base_lims
for j=1:N
%	C_fgbn(j,:,iLoop)=blacknutt.*C_fgbn(j,:,iLoop)
	C_fg3t(j,:,iLoop)=bn3t.*C_fg3t(j,:,iLoop)
	%C_fgbh(j,:,iLoop)=transpose(blackmanharris(N,'symmetric')).*C_fgbh(j,:,iLoop)
end
end
%clearvars C_fg
parfor i=1:base_lims
%kk_bn(:,:,i)=fft2(C_fg(:,:,i))
%kk_bn(:,:,i)=B_bn*C_fgbn(:,:,i)*Bdag_bn
%	kk_kai(:,:,i)=B_bn*C_fgkai(:,:,i)*Bdag_bn
%	kk_bh(:,:,i)=B_bn*C_fgbh(:,:,i)*Bdag_bn
%	kk_ftop(:,:,i)=((B_bn*C_fgftop(:,:,i)*Bdag_bn))
	kk_3t(:,:,i)=((B_bn*C_fg3t(:,:,i)*Bdag_bn))
end


parfor i=1:N
%kkp_bn(i,:)=(kk_bn(i,i,:))
%kkp_kai(i,:)=kk_kai(i,i,:)
%kkp_bh(i,:)=kk_bh(i,i,:)
%kkp_gau(i,:)=kk_gau(i,i,:)
%kkp_ftop(i,:)=kk_ftop(i,i,:)
kkp_3t(i,:)=kk_3t(i,i,:)
end

parfor i=1:N./2
%kkpl(i,:)=abs(kkp(i,:))./2
%kkpl_bn(i,:)=abs(kkp_bn(i,:))./2
%	kkpl_bn(i,:)=((abs(kkp_bn(i,:))+abs(kkp_bn(N-i,:)))./2).*1
%	kkpl_gau(i,:)=(abs(kkp_gau(i,:))+abs(kkp_gau(N-i,:)))./2
%	kkpl_kai(i,:)=(abs(kkp_kai(i,:)+abs(kkp_kai(N-i,:)))./2).*1
%	kkpl_bh(i,:)=(abs(kkp_bh(i,:)+abs(kkp_bh(N-i,:)))./2).*1
%	kkpl_gau(i,:)=abs(kkp_gau(i,:))./2
%	kkpl_ftop(i,:)=(abs(kkp_ftop(i,:))+abs(kkp_ftop(N-i,:))./2).*1
	kkpl_3t(i,:)=(abs(kkp_3t(i,:))+abs(kkp_3t(N-i,:))./2).*1
end
%save('Tapered.mat')
%clear

%kkpl=kkpl*1e+6
