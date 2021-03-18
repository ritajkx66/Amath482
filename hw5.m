clear all; close all; clc

video = VideoReader('ski_drop.mov');
numFrames = video.NumberOfFrames;
vidFrames = read(video);
vidFrames_resize = imresize(vidFrames,0.2);
for j = 1:numFrames
    X = vidFrames_resize(:,:,:,j);
    X_gray = rgb2gray(X);
    X_double = im2double(X_gray);
    X_reshape(:,j) = reshape(X_double,[],1);
%     imshow(X); drawnow
end

duration = video.Duration;
% dt = duration/numFrames;
dt = 1/video.Framerate;
t = 0:dt:duration;

X1 = X_reshape(:,1:end-1);
X2 = X_reshape(:,2:end);
[U, Sigma, V] = svd(X1,'econ');
figure(1)
plot(diag(Sigma)/sum(diag(Sigma)),'mo','Linewidth',2)
ylabel("Energy captured");
xlabel("Singular values");
title("Energy of the singular values");

r = 1;
U = U(:,1:r);
Sigma = Sigma(1:r,1:r);
V = V(:,1:r);
S = U'*X2*V*diag(1./diag(Sigma));
[eV, D] = eig(S); % compute eigenvalues + eigenvectors
mu = diag(D); % extract eigenvalues
omega = log(mu)/dt;
Phi = U*eV;

y0 = Phi\X1(:,1); % pseudoinverse to get initial conditions
u_modes = zeros(length(y0),length(t)-1);
for iter = 1:(length(t)-1)
   u_modes(:,iter) = y0.*exp(omega*t(iter)); 
end
u_dmd = Phi*u_modes;

Xspar = X1 - abs(u_dmd);
trueorfalse = Xspar < 0;
R = Xspar.*trueorfalse;
fore = Xspar - R;
back = R + abs(u_dmd);
Xnew = fore + back;

figure(2)
sizevid = size(X_double);
subplot(2,2,1)
pic_spar = reshape(Xspar,sizevid(1),sizevid(2),454);
imshow(pic_spar(:,:,320))
title("Foreground video");
subplot(2,2,2)
pic_dmd = reshape(u_dmd,sizevid(1),sizevid(2),454);
imshow(pic_dmd(:,:,320))
title("background video");
subplot(2,2,3)
pic_fore = reshape(fore,sizevid(1),sizevid(2),454);
imshow(pic_fore(:,:,320))
title("Reconstructed foreground video");
subplot(2,2,4)
pic_back = reshape(back,sizevid(1),sizevid(2),454);
imshow(pic_back(:,:,320))
title("Reconstructed background video");




clear all; close all; clc

video = VideoReader('monte_carlo.mov');
numFrames = video.NumberOfFrames;
vidFrames = read(video);
vidFrames_resize = imresize(vidFrames,0.2);
for j = 1:numFrames
    X = vidFrames_resize(:,:,:,j);
    X_gray = rgb2gray(X);
    X_double = im2double(X_gray);
    X_reshape(:,j) = reshape(X_double,[],1);
    imshow(X); drawnow
end

duration = video.Duration;
dt = duration/numFrames;
dt = 1/video.Framerate;
t = 0:dt:duration;

X1 = X_reshape(:,1:end-1);
X2 = X_reshape(:,2:end);
[U, Sigma, V] = svd(X1,'econ');
figure(1)
plot(diag(Sigma)/sum(diag(Sigma)),'mo','Linewidth',2)
ylabel("Energy captured");
xlabel("Singular values");
title("Energy of the singular values");

lambda = diag(Sigma);
threshold = 0.8;
energy = 0;
r = 0;
while energy <= threshold
    r = r+1;
    energy = energy+lambda(r)/sum(lambda);
end

% U = U(:,1:r);
% Sigma = Sigma(1:r,1:r);
% V = V(:,1:r);
% S = U'*X2*V*(1./diag(Sigma));
S = U(:,1:r)'*X2*V(:,1:r)/Sigma(1:r,1:r);
[eV, D] = eig(S); % compute eigenvalues + eigenvectors
mu = diag(D); % extract eigenvalues
omega = log(mu)/dt;
Phi = U(:,1:r)*eV;

y0 = Phi\X1(:,1); % pseudoinverse to get initial conditions
u_modes = zeros(length(y0),length(t)-2);
for iter = 1:(length(t)-2)
   u_modes(:,iter) = y0.*exp(omega*t(iter)); 
end
u_dmd = Phi*u_modes;

Xspar = X1 - abs(u_dmd);
trueorfalse = Xspar < 0;
R = Xspar.*trueorfalse;
fore = Xspar - R;
back = R + abs(u_dmd);
Xnew = fore + back;

figure(2)
sizevid = size(X_double);
subplot(2,2,1)
pic_spar = reshape(Xspar,sizevid(1),sizevid(2),380);
imshow(pic_spar(:,:,370))
title("Foreground video");
subplot(2,2,2)
pic_dmd = reshape(u_dmd,sizevid(1),sizevid(2),380);
imshow(pic_dmd(:,:,370))
title("Background video");
subplot(2,2,3)
pic_fore = reshape(fore,sizevid(1),sizevid(2),380);
imshow(pic_fore(:,:,370))
title("Reconstructed foreground video");
subplot(2,2,4)
pic_back = reshape(back,sizevid(1),sizevid(2),380);
imshow(pic_back(:,:,370))
title("Reconstructed background video");



