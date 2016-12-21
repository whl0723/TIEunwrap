%%% TIE?
clc
clear
close all

%%% 单位初始化
mm=1e-3;
um=1e-6;

%%% 离焦量
dz=500*um;
amplitude=1;

%%% 图像尺寸
Lx=2.56*mm;
Ly=2.56*mm;

%%% 波长
lamda=0.5*um;

%%% 图像采样点数
M=256;
N=256;

%%% 网格化
mx=-M/2:M/2-1;
my= -N/2:N/2-1;
x=linspace(-Lx/2,Lx/2,M+1);x=x(1:M); 
y=linspace(-Ly/2,Ly/2,N+1);y=y(1:N); 
[x1,y1]=meshgrid(x,y);      % 真实大小
[mx,my]=meshgrid(mx,my);    % 整数坐标

%%% 振幅强度
a=2*peaks(M);

%%% 模拟相位分布
b=a; %
alp=1e7;
beta=1e10;
% pset=2*peaks(M);%+10*exp(-0.05*(1*x1.^2+1*y1.^2));
pset=alp.*(x1.^2+1*y1.^2);
% pset=beta.*(x1.^3+1*y1.^3);
pset=pset-min(pset(:));
figure;imagesc(pset);colorbar

%%% 包裹相位
uphase=mod(pset,2*pi)-pi; %
figure;imagesc(uphase);colorbar

%%% 正焦图 I0
f0=b.*exp(1i*pset); %
f0=f0./abs(f0);
I0=abs(f0).^2; %

%%% 欠焦图 I1
H1=exp(-1i*dz*2*pi./lamda*sqrt(1-(lamda*mx/Lx).^2+(lamda*my/Ly).^2));
Ui=fft2(f0);
% H1=exp(-1i*pi*lamda*dz*(((mx/Lx).^2)+(my/Ly).^2)); %
Ud=Ui.*H1;
fd=ifft2((Ud)); %
I1=abs(fd).^2; %

%%%% 过焦图 I2
H2=exp(1i*dz*2*pi./lamda*sqrt(1-(lamda*mx/Lx).^2+(lamda*my/Ly).^2));
Ui=fft2(f0);
% H2=exp(1i*pi*lamda*dz*(((mx/Lx).^2)+(my/Ly).^2)); %
Ud0=Ui.*H2;
fd0=ifft2(Ud0); %
I2=abs(fd0).^2; %

%%% 图像显示
figure;
subplot 131
imshow(I0,[]);
subplot 132
imshow(I1,[]);
subplot 133
imshow(I2,[]);

%% STIE 求解相位
%%% 求解TIE方程等号右边项 TEMP 
TEMP=(I2-I1)./2./dz;
figure;
imagesc(TEMP);

%%% 求解 空间频率 kx， ky
kx=(2*pi/Lx)*mx; ky=(2*pi/Ly)*my;
kx(kx==0)=eps;ky(ky==0)=eps;

%%% 求解相位分布
Nabla=ifft2(fft2(TEMP)./((kx.^2+ky.^2)));
output=-2*pi/lamda.*Nabla./I0;
output=output-min(output(:));

%%% 对比显示
figure;
subplot 121
imagesc(pset)
colorbar
subplot 122
imagesc(output)
colorbar

%%% 误差对比
err=output-pset;
err=err-min(err(:));
figure;
imagesc(err)
colorbar


% 
% phase1=djunwrap(uphase);
% figure;
% imagesc(phase1);title('行列')
% colorbar
% 
% err=phase1-pset;
% figure;
% imagesc(err);title('err')
% colorbar