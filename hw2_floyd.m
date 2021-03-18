clear all; close all; clc

figure(1)
[y, Fs] = audioread('Floyd.m4a');
plot((1:length(y))/Fs,y);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Comfortably Numb');
p8 = audioplayer(y,Fs); playblocking(p8);

S = y';
trgnr = length(y)/Fs; % record time in seconds
L = trgnr;
n = length(y);
t2 = linspace(0,L,n+1); t = t2(1:n);
k = (1/L)*[0:n/2-1 -n/2:-1];
ks = fftshift(k);
Sgt_spec = [];

a = 5500;
tau = 0:1:L;
for i = 1:length(tau)
   g = exp(-a*(t - tau(i)).^2); % Window function
   Sg = g.*S;
   Sgt = fft(Sg);
   Sgt_spec(:,i) = fftshift(abs(Sgt));
end

figure(2)
pcolor(tau,ks,Sgt_spec(1:end-1,:));
shading interp
set(gca,'ylim',[0 1000],'Fontsize',12)
colormap(hot)
colorbar
yticks([78 104 131 185 247]);
yticklabels({'D^#','G^#','C','F^#','B'});
xlabel('time (t)'), ylabel('frequency (k)')
title('Comfortably Numb','Fontsize',16)

figure(3)
pcolor(tau,ks,Sgt_spec(1:end-1,:));
shading interp
set(gca,'ylim',[0 250],'Fontsize',16)
colormap(hot)
colorbar
yticks([78 104 131 185 247]);
yticklabels({'D^#','G^#','C','F^#','B'});
xlabel('time (t)'), ylabel('frequency (k)')
title('Bass in the Comfortably Numb','Fontsize',16)

figure(4)
pcolor(tau,ks,Sgt_spec(1:end-1,:));
shading interp
set(gca,'ylim',[250 1000],'Fontsize',16)
colormap(hot)
colorbar
% yticks([78 98 110 117]);
% yticklabels({'Em', 'G', 'A', 'Bm'});
yticks([262 392 440 587 740]);
yticklabels({'C','G','A','D','F^#'});
xlabel('time (t)'), ylabel('frequency (k)')
title('Guitar in the Comfortably Numb','Fontsize',16)
