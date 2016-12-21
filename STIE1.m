function fai=STIE1(I1,I2,pixelsize,dz,lamda)
% I1=[I1 fliplr(I1)
%     flipud(I1) rot90(I1,2)];
% I2=[I2 fliplr(I2)
%     flipud(I2) rot90(I2,2)];
temp=(I2-I1)./I1./dz.*2*pi./lamda;
% temp=medfilt2(temp);
[M,N]=size(I1);
mx=-N/2:N/2-1;
my= -M/2:M/2-1;
[mx,my]=meshgrid(mx,my);
Lx=M*pixelsize;
Ly=N*pixelsize;
% k=2*pi./lamda;
fx=mx./Lx;
fy=my./Ly;
D=fx.^2+fy.^2;
% gamma=eps;
gamma=1./max(Lx,Ly);
D=D+gamma;
%%
d0=10; 
Rh=3; 
Rl=0.0; 
c=1; 
n=1; 
d=sqrt((my-0).^2+(mx-0).^2);
h4=(Rh-Rl).*exp(-c.*(d0./d)).^n+Rl;

XS=h4./D;

%%
% XS=fftshift(XS);
% k2z=1./XS;
% 
% XS=fftshift(XS);
% k2=1./XS;

% figure;
% subplot 121
% imshow(k2,[])
% title('k2')
% 
% subplot 122
% imshow(k2z,[])
% title('k2z')

% fid=fopen('k2.bin','wb');
% fwrite(fid,k2','double');
% fclose(fid)
% 
% fid=fopen('k2z.bin','wb');
% fwrite(fid,k2z','double');
% fclose(fid)
%%


%%

fai=ifft2(ifftshift(XS.*fftshift(fft2(temp))))./25;
% fai=2*fai(1:M/2,1:N/2)./I1;
% fai=fai;
fai=fai-min(fai(:));

%%
% f=fft2(fai);                            %计算图像I的傅里叶变换 
% g=fftshift(f);                        %将图像频谱中心从矩阵原点移到矩阵的中心 
% 
% 
% my= -M/2:M/2-1;mx= -N/2:N/2-1; 
% [J,I]=meshgrid(mx,my);
% 
% g4=h4.*g;
% % g4=h4.*g;
% g4=ifftshift(g4); 
% fai=real(ifft2(g4)); 
