clc;
clear all;
close all;

Transmit_data=importdata("C:\Users\DELL\Tn_GroupID_1A_FreqID_60_SNR_10.000000_Foff_25.000000_Nframe_1_Nstart_1.txt");
str = importdata("C:\Users\DELL\Input_GroupID_1A_FreqID_60_SNR_10.000000_Foff_25.000000_Nframe_1_Nstart_1.txt");

for i =1:1024
    data(i)=str(i,1)+i*str(i,2);
end
frame=data;
for i =1:length(Transmit_data)
    T(i)=Transmit_data(i,1)+i*Transmit_data(i,2);
end
 freq_offset=25;
 np=0:(length(T)-1);
offset_signal=cos(2*pi*freq_offset*np)+j*sin(2*pi*freq_offset*np);
corrupted=T.*offset_signal;
noise=wgn(1,length(T),10);
Re=corrupted+noise;

Receive_data=importdata("C:\Users\DELL\Rn_GroupID_1A_FreqID_60_SNR_10.000000_Foff_25.000000_Nframe_1_Nstart_1.txt");
for i =1:length(Receive_data)
    R(i)=Receive_data(i,1)+i*Receive_data(i,2);
end
 subplot(211);plot(corrupted,'r');
% hold on;
subplot(212);plot(R,'g');
