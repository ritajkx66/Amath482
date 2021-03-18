clear all; close all; clc

figure(1)
[y, Fs] = audioread('GNR.m4a');
plot((1:length(y))/Fs,y);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Sweet Child O'' Mine');
p8 = audioplayer(y,Fs); playblocking(p8);

S = y';
trgnr = length(y)/Fs; % record time in seconds
L = trgnr;
n = length(y);
t2 = linspace(0,L,n+1); t = t2(1:n);
k = (1/L)*[0:n/2-1 -n/2:-1];
ks = fftshift(k);

a = 1500;
tau = 0:0.1:L;
for i = 1:length(tau)
   g = exp(-a*(t - tau(i)).^2); % Window function
   Sg = g.*S;
   Sgt = fft(Sg);
   Sgt_spec(:,i) = fftshift(abs(Sgt));
end

figure(2)
pcolor(tau,ks,Sgt_spec(1:end,:));
shading interp
set(gca,'ylim',[200 800],'Fontsize',16)
colormap(hot)
colorbar
yticks([277 370 415 554 698 740])
yticklabels({'C^#','F^#','G^#','{C^#}^+','F','{F^#}^+'})
xlabel('time (t)'), ylabel('frequency (k)')
title('Sweet Child O'' Mine','Fontsize',16);



