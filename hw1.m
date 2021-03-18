% Clean workspace
clear all; close all; clc

load subdata.mat % Imports the data as the 262144x49 (space by time) matrix called subdata

L = 10; % spatial domain
n = 64; % Fourier modes
x2 = linspace(-L,L,n+1); x = x2(1:n); y = x; z = x;
k = (2*pi/(2*L))*[0:(n/2 - 1) -n/2:-1]; ks = fftshift(k);

[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

% problem1
ave = zeros(n,n,n);
for j=1:49
    Un(:,:,:)=reshape(subdata(:,j),n,n,n);
    Utn=fftn(Un);
    ave = ave+Utn;
end
ave = abs(fftshift(ave))/49;
[A,I] = max(ave(:));
center = [Kx(I),Ky(I),Kz(I)];
%isosurface(X,Y,Z,abs(ave)/max(abs(ave), [], 'all'),0.7)

% problem2
tau = 0.5;
filter = exp(-tau*((Kx-center(1)).^2+(Ky-center(2)).^2+(Kz-center(3)).^2));
indice = zeros(1,49);
for k = 1:49
    Un(:,:,:)=reshape(subdata(:,k),n,n,n);
    Utn = fftn(Un);
    Unft = filter.*fftshift(Utn);
    Unf = ifftn(ifftshift(Unft));
    maximum = max(abs(Unf(:)));
    indice(k) = find(maximum == abs(Unf));
end
sz = [64,64,64];
[Xvec,Yvec,Zvec] = ind2sub(sz,indice);
xp = Kx(Xvec,Yvec,Zvec);
yp = Ky(Xvec,Yvec,Zvec);
zp = Kz(Xvec,Yvec,Zvec);
xxp = squeeze(xp(49,:,49));
yyp = squeeze(yp(:,49,49));
zzp = squeeze(zp(49,49,:));
plot3(xxp,yyp,zzp,"LineWidth",2)
xlabel("x","fontsize",15)
ylabel("y","fontsize",15)
zlabel("z","fontsize",15)
title("Path of the submarine","fontsize",20)
%annotation('textbox','String','The submarine follows a generally smooth spiral path.')

% problem3
T = table(xxp.',yyp);
