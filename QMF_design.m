clc;
clear;
close all;
n=0:1023;
Ts=1/1000;
%% Band Pass Data
%  xi = 1344*cos(0.06*pi*n)+864*cos(0.18*n)+8543*cos(0.38*pi*n)-43*cos(0.8*pi*n);
%  xq = 1344*sin(0.06*pi*n)+864*sin(0.18*pi*n)+8543*sin(0.38*pi*n)-43*sin(0.8*pi*n);
%% Random Data
xi=randn(1,1024);
xq=randn(1,1024);
data=xi+j*xq;
y = fft(data);
fs = 1/Ts;
nl = length(data);
fshift = (-nl/2:nl/2-1)*(fs/nl);
yshift = fftshift(y);
figure;
plot(fshift,abs(yshift))
f = (0:length(y)-1)*fs/length(y);
%% Analysis Part
h0_coeff = [-1 0 3 0 -8 0 21 0 -45 0 91 0 -191 0 643 1024 643 0 -191 0 91 0 -45 0 21 0 -8 0 3 0 -1];
n=0:[length(h0_coeff)-1];
h1_coeff = (-1).^n.*h0_coeff;
len=length(h0_coeff);
vl=conv(data,h0_coeff);
vh=conv(data,h1_coeff);
[lp_filter,wl]=freqz(h0_coeff,1,256);
[hp_filter,wh]=freqz(h1_coeff,1,256);
%lp_db_value=20*log10(abs(lp_filter/max(lp_filter)));
[lp_freq]=fft(h0_coeff);
[hp_freq]=fft(h1_coeff);
leng=length(lp_freq);

T1=0.5*((abs(lp_freq)).^2-(abs(hp_freq)).^2);%% Transfer Function
hp_db_value=20*log10(abs(hp_filter/max(hp_filter)));
figure;
plot(2*wl/pi,20*log10(abs(lp_filter/max(lp_filter))),'g');
hold on;
plot(2*wh/pi,20*log10(abs(hp_filter/max(hp_filter))),'r');
down_vl=downsample(vl,2);
down_vh=downsample(vh,2);
up_ag_vl=upsample(down_vl,2);
up_ag_vh=upsample(down_vh,2);
%% Synthesis Part
g0_coeff= h0_coeff;
g1_coeff=-h1_coeff;
g1_coeff_2=h1_coeff;
syn_low=conv(up_ag_vl,g0_coeff);
syn_high=conv(up_ag_vh,g1_coeff);
syn_high_2=conv(up_ag_vh,g1_coeff_2);
xbar=syn_low+syn_high;
xbar_2=syn_low+syn_high_2;
figure;
subplot(211);plot(data)
subplot(212);plot(xbar);
%% MSE
sum=0;
for k=1:1024
    MSEr(k)=abs(data(k)-xbar(k))^2/abs(data(k))^2;
    sum=sum+MSEr(k);
end
MSError1=sqrt(sum)/1024
%% Another Transfer Function
nl1=length(lp_freq);
fl = (-nl1/2:nl1/2-1)*(fs/nl1);
T2=0.5*((abs(lp_freq)).^2+(abs(hp_freq)).^2);
A2= abs(lp_freq.*hp_freq);
figure;
figure;
plot(2*fl/pi, 20*log10(abs(T1)/max(T1)),'y');
hold on;
plot(2*fl/pi, 20*log10(abs(T2)/max(T2)),'b');
hold on;
plot(2*fl/pi, 20*log10(abs(A2)/max(A2)),'g');
sum=0;
for k=1:1024
    MSEr2(k)=abs(data(k)-xbar_2(k))^2/abs(data(k))^2;
    sum=sum+MSEr2(k);
end
MSError2=sqrt(sum)/1024
