clear all; close all; clc

% Ideal case
load('cam1_1.mat')
numFrames1_1 = size(vidFrames1_1,4);
load('cam2_1.mat')
numFrames2_1 = size(vidFrames2_1,4);
load('cam3_1.mat')
numFrames3_1 = size(vidFrames3_1,4);

% wx = 200;
% wy = 200;
% Nx = 480;
% Ny = 640;

x1 = zeros(numFrames1_1);
y1 = zeros(numFrames1_1);
% filter = zeros(Nx, Ny);
% filter(Nx/4 + wx/4:Nx/2 + wx, Ny/2 - wy/10:Ny/2 + wy/2) = 1;
for j = 1:numFrames1_1
    X = vidFrames1_1(:,:,:,j);
    X = rgb2gray(X);
%     X = im2double(X);
%     [Nx,Ny] = size(X);
%     X = X.*filter;
%     X = im2uint8(X);
    [placex, placey] = find(X > 230);
    x1(j) = mean(placex);
    y1(j) = mean(placey);

end

x2 = zeros(numFrames2_1);
y2 = zeros(numFrames2_1);
% filter = zeros(Nx,Ny);
% filter(Nx/4 - wx/10:Nx/2 + wx, Ny/2 - wy/2:Ny/2 + wy/5) = 1;
for j = 1:numFrames2_1
    X = vidFrames2_1(:,:,:,j);
    X = rgb2gray(X);
%     X = im2double(X);
%     X = X.*filter;
%     X = im2uint8(X);
    [placex, placey] = find(X > 230);
    x2(j) = mean(placex);
    y2(j) = mean(placey);
end

x3 = zeros(numFrames3_1);
y3 = zeros(numFrames3_1);
% filter = zeros(Nx,Ny);
% filter(Nx/2 - wx/10:Nx/2 + wx/2, Ny/4 + wy/2 - wy/10:Ny/2 + wy - wy/10) = 1;
for j = 1:numFrames3_1
    X = vidFrames3_1(:,:,:,j);
    X = rgb2gray(X);
%     X = im2double(X);
%     X = X.*filter;
%     X = im2uint8(X);
    [placex, placey] = find(X > 230);
    x3(j) = mean(placex);
    y3(j) = mean(placey);
end

min = min([numFrames1_1,numFrames2_1,numFrames3_1]);
xs = [x1(1:min);y1(1:min);x2(1:min);y2(1:min);x3(1:min);y3(1:min)];
position = xs - mean(xs,2);
[U,S,V] = svd(position/sqrt(min-1),'econ');
projection = U' * position;
sigma1 = diag(S);

figure(1)
subplot(3,1,1)
plot(1:6,sigma1,'mo','Linewidth',2);
xlabel('Measurements'); ylabel('Singular values');
title("Ideal case: singular values");
subplot(3,1,2)
plot(1:min, position(1,:), 1:min, position(2,:))
xlabel("Time(frames)"); ylabel("Displacement(pixels)"); 
title("Ideal case: original displacement");
legend("z", "x-y plane")
subplot(3,1,3)
plot(1:min,projection(1,:),1:min,projection(2,:),1:min,projection(3,:))
xlabel("Time(frames)"); ylabel("Displacement(pixels)"); 
title("Ideal case: displacement after PCA");
legend('PC1', 'PC2', 'PC3');

% Noisy case
load('cam1_2.mat')
numFrames1_2 = size(vidFrames1_2,4);
load('cam2_2.mat')
numFrames2_2 = size(vidFrames2_2,4);
load('cam3_2.mat')
numFrames3_2 = size(vidFrames3_2,4);


x1 = zeros(numFrames1_2);
y1 = zeros(numFrames1_2);
for j = 1:numFrames1_2
    X = vidFrames1_2(:,:,:,j);
    X = rgb2gray(X);
    [placex, placey] = find(X > 230);
    x1(j) = mean(placex);
    y1(j) = mean(placey);

end

x2 = zeros(numFrames2_2);
y2 = zeros(numFrames2_2);
for j = 1:numFrames2_2
    X = vidFrames2_2(:,:,:,j);
    X = rgb2gray(X);
    [placex, placey] = find(X > 230);
    x2(j) = mean(placex);
    y2(j) = mean(placey);
end

x3 = zeros(numFrames3_2);
y3 = zeros(numFrames3_2);
for j = 1:numFrames3_2
    X = vidFrames3_2(:,:,:,j);
    X = rgb2gray(X);
    [placex, placey] = find(X > 230);
    x3(j) = mean(placex);
    y3(j) = mean(placey);
end

min = min([numFrames1_2,numFrames2_2,numFrames3_2]);
xs = [x1(1:min);y1(1:min);x2(1:min);y2(1:min);x3(1:min);y3(1:min)];
position = xs - mean(xs,2);
[U,S,V] = svd(position/sqrt(min-1),'econ');
projection = U' * position;
sigma2 = diag(S);

figure(2)
subplot(3,1,1)
plot(1:6,sigma2,'mo','Linewidth',2);
xlabel('Measurements'); ylabel('Singular values');
title("Noisy case: singular values");
subplot(3,1,2)
plot(1:min, position(1,:), 1:min, position(2,:))
xlabel("Time(frames)"); ylabel("Displacement(pixels)"); 
title("Noisy case: original displacement");
legend("z", "x-y plane")
subplot(3,1,3)
plot(1:min,projection(1,:),1:min,projection(2,:),1:min,projection(3,:),1:min,projection(4,:))
xlabel("Time(frames)"); ylabel("Displacement(pixels)"); 
title("Noisy case: displacement after PCA");
legend('PC1', 'PC2', 'PC3', 'PC4');

% Horizontal displacement
load('cam1_3.mat')
numFrames1_3 = size(vidFrames1_3,4);
load('cam2_3.mat')
numFrames2_3 = size(vidFrames2_3,4);
load('cam3_3.mat')
numFrames3_3 = size(vidFrames3_3,4);

x1 = zeros(numFrames1_3);
y1 = zeros(numFrames1_3);
for j = 1:numFrames1_3
    X = vidFrames1_3(:,:,:,j);
    X = rgb2gray(X);
    [placex, placey] = find(X > 230);
    x1(j) = mean(placex);
    y1(j) = mean(placey);

end


x2 = zeros(numFrames2_3);
y2 = zeros(numFrames2_3);
for j = 1:numFrames2_3
    X = vidFrames2_3(:,:,:,j);
    X = rgb2gray(X);
    [placex, placey] = find(X > 230);
    x2(j) = mean(placex);
    y2(j) = mean(placey);
end

x3 = zeros(numFrames3_3);
y3 = zeros(numFrames3_3);
for j = 1:numFrames3_3
    X = vidFrames3_3(:,:,:,j);
    X = rgb2gray(X);
    [placex, placey] = find(X > 230);
    x3(j) = mean(placex);
    y3(j) = mean(placey);
end

min = min([numFrames1_3,numFrames2_3,numFrames3_3]);
xs = [x1(1:min);y1(1:min);x2(1:min);y2(1:min);x3(1:min);y3(1:min)];
position = xs - mean(xs,2);
[U,S,V] = svd(position/sqrt(min-1),'econ');
projection = U' * position;
sigma3 = diag(S);

figure(3)
subplot(3,1,1)
plot(1:6,sigma3,'mo','Linewidth',2);
xlabel('Measurements'); ylabel('Singular values');
title("Horizontal displacement: singular values");
subplot(3,1,2)
plot(1:min, position(1,:), 1:min, position(2,:))
xlabel("Time(frames)"); ylabel("Displacement(pixels)"); 
title("Horizontal displacement: original displacement");
legend("z", "x-y plane")
subplot(3,1,3)
plot(1:min,projection(1,:),1:min,projection(2,:),1:min,projection(3,:),1:min,projection(4,:))
xlabel("Time(frames)"); ylabel("Displacement(pixels)"); 
title("Horizontal displacement: displacement after PCA");
legend('PC1', 'PC2', 'PC3', 'PC4');

% Horizontal displacement and rotation
load('cam1_4.mat')
numFrames1_4 = size(vidFrames1_4,4);
load('cam2_4.mat')
numFrames2_4 = size(vidFrames2_4,4);
load('cam3_4.mat')
numFrames3_4 = size(vidFrames3_4,4);

x1 = zeros(numFrames1_4);
y1 = zeros(numFrames1_4);
for j = 1:numFrames1_4
    X = vidFrames1_4(:,:,:,j);
    X = rgb2gray(X);
    [placex, placey] = find(X > 230);
    x1(j) = mean(placex);
    y1(j) = mean(placey);

end

x2 = zeros(numFrames2_4);
y2 = zeros(numFrames2_4);
for j = 1:numFrames2_4
    X = vidFrames2_4(:,:,:,j);
    X = rgb2gray(X);
    [placex, placey] = find(X > 230);
    x2(j) = mean(placex);
    y2(j) = mean(placey);
end

x3 = zeros(numFrames3_4);
y3 = zeros(numFrames3_4);
for j = 1:numFrames3_4
    X = vidFrames3_4(:,:,:,j);
    X = rgb2gray(X);
    [placex, placey] = find(X > 230);
    x3(j) = mean(placex);
    y3(j) = mean(placey);
end

min = min([numFrames1_4,numFrames2_4,numFrames3_4]);
xs = [x1(1:min);y1(1:min);x2(1:min);y2(1:min);x3(1:min);y3(1:min)];
position = xs - mean(xs,2);
[U,S,V] = svd(position/sqrt(min-1),'econ');
projection = U' * position;
sigma4 = diag(S);

figure(4)
subplot(3,1,1)
plot(1:6,sigma4,'mo','Linewidth',2);
xlabel('Measurements'); ylabel('Singular values');
title("Horizontal displacement and rotation: singular values");
subplot(3,1,2)
plot(1:min, position(1,:), 1:min, position(2,:))
xlabel("Time(frames)"); ylabel("Displacement(pixels)"); 
title("Horizontal displacement and rotation: original displacement");
legend("z", "x-y plane")
subplot(3,1,3)
plot(1:min,projection(1,:),1:min,projection(2,:),1:min,projection(3,:),1:min,projection(4,:))
xlabel("Time(frames)"); ylabel("Displacement(pixels)"); 
title("Horizontal displacement and rotation: displacement after PCA");
legend('PC1', 'PC2', 'PC3', 'PC4');

