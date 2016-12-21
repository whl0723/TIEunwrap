clc
clear
close all

%%% 参数设置
um=1e-6;
pixelsize=10*um;
N=512;
M=512;
Lx=M*pixelsize;
Ly=N*pixelsize;
lamda=0.63*um;
mx=-N/2:N/2-1;
my= -M/2:M/2-1;
x=mx*pixelsize;
y=my*pixelsize;
% [mx,my]=meshgrid(mx,my);    % 整数坐标
[x,y]=meshgrid(x,y);    

%%% 模拟相位分布
kk=1;
% for i=0:0.03:pi/2-0.03

pset1=double(ones(M,N));
pset1(250:300,:)=0;
% pset1=pset1./max(pset1(:))*2;

b=1; %
alp=1e6;
beta=3e10;
% pset=2*peaks(M);%+10*exp(-0.05*(1*x.^2+1*y.^2));
% pset=alp.*(x.^2+1*y.^2);
% pset=beta.*(x.^3+1*y.^3);
pset=1e3*(x.^1+1*y.^1);
pset=pset.*pset1;
% pset=pset+pset1;
pset=pset-min(pset(:));

% figure;
% subplot 121
% imagesc(pset);colorbar

%%% 包裹相位
uphase=mod(pset,2*pi)-pi; %
% subplot 122
% imagesc(uphase);colorbar

% uphase1=[uphase fliplr(uphase)
%     flipud(uphase) rot90(uphase,2)];
% figure;
% imagesc(uphase1)

fai=TIEunwrap(uphase,pixelsize,lamda);


% figure;
% imagesc(pGS3)
% title('行列')

% %%% 辅助场
% uphase=[uphase fliplr(uphase)
%     flipud(uphase) rot90(uphase,2)];
% 
% u0=b.*exp(1i.*uphase);
% % u00=u0./abs(u0);
% u00=u0;
% % figure;
% % imagesc(imag(u00-u0));
% 
% %%% 衍射
% dz=1*um;
% I0=abs(u00);
% 
% [M,N]=size(I0);
% mx=-N/2:N/2-1;
% my= -M/2:M/2-1;
% [mx,my]=meshgrid(mx,my);    % 整数坐标
% Lx=M*pixelsize;
% Ly=N*pixelsize;
% % figure;imagesc(I0);
% k=2*pi./lamda;
% fx=mx./Lx;
% fy=my./Ly;
% % H=exp(1i*dz.*sqrt(k.^2-4*pi*pi*(fx.^2+fy.^2)));
% H=exp(-1i*pi*lamda*dz*(((mx/Lx).^2)+(my/Ly).^2)); 
% U0=fftshift(fft2(u00));
% u1=ifft2(ifftshift(U0.*H));
% I1=abs(u1);
% 
% % figure;
% % subplot 121
% % imagesc(abs(u1));
% % subplot 122
% % imagesc(angle(u1))
% 
% %%% TIE
% % fx=mx./Lx;
% % fy=my./Ly;
% temp=(I1-I0)./dz;
% figure;imagesc(temp)
% D=fx.^2+fy.^2;
% gamma=eps;
% % gamma=1./max(Lx,Ly);
% fai=ifft2(ifftshift(k*D.*fftshift(fft2(temp))./(4*pi*pi*(D.^2)+gamma)));
% fai=2*fai(1:M/2,1:N/2);
% fai=fai-min(fai(:));

% figure;
% imagesc(fai);
% colorbar

jb = jieb(uphase);
pGS3=djunwrap(uphase);

err=fai-pset;
err2=jb-pset;
err3=pGS3-pset;
% aa(kk)=std2(err);
% bb(kk)=std2(err2);
% cc(kk)=std2(err3);
% kk=kk+1;

% end
% i=0:0.03:pi/2-0.03;
% figure;
% plot(i,aa,'r.')
% hold on
% plot(i,bb,'g.')
% plot(i,cc,'b.')




% figure
% imagesc(jb)
% colorbar
% title('最小二乘')

figure;
subplot 131
imagesc(err);title('TIE的误差')
colorbar
subplot 132
imagesc(err2)
colorbar
title('最小二乘的误差')

subplot 133
imagesc(err3)
colorbar
title('行列的误差')

figure
imagesc(fai)
hold on
plot(240,194:352,'r')

figure;
plot(fai(194:352,240),'r')
hold on
plot(jb(194:352,240),'g')
plot(pGS3(194:352,240),'b')



