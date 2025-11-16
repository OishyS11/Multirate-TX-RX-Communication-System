clc;
clear all;
close all;

str = importdata("C:\Users\DELL\Input_GroupID_1A_FreqID_60_SNR_10.000000_Foff_25.000000_Nframe_1_Nstart_1.txt");

for i =1:1024
    data(i)=str(i,1)+i*str(i,2);
end
frame=data;

V0_af1=importdata("C:\Users\DELL\V0_GroupID_1A_FreqID_60_SNR_10.000000_Foff_25.000000_Nframe_1_Nstart_1.txt");
V0_a1=V0_af1.rowheaders;
V0_a1=str2double(V0_a1);
V0_a2=[V0_af1.data];
v0_a=V0_a1+j*V0_a2;
V1_af2=importdata("C:\Users\DELL\V1_GroupID_1A_FreqID_60_SNR_10.000000_Foff_25.000000_Nframe_1_Nstart_1.txt");
V1_a1=V1_af2.rowheaders;
V1_a1=str2double(V1_a1);
V1_a2=[V1_af2.data];
%V0_scdq=[V0_scd2{:}]; 
v1_a=V1_a1+j*V1_a2;
V2_af3=importdata("C:\Users\DELL\V2_GroupID_1A_FreqID_60_SNR_10.000000_Foff_25.000000_Nframe_1_Nstart_1.txt");
V3_af4=importdata("C:\Users\DELL\V3_GroupID_1A_FreqID_60_SNR_10.000000_Foff_25.000000_Nframe_1_Nstart_1.txt");

for i =1:length(V2_af3)
    v2_a(i)=V2_af3(i,1)+i*V2_af3(i,2);
end
for i =1:length(V3_af4)
    v3_a(i)=V3_af4(i,1)+i*V3_af4(i,2);
end


h0_coeff_g = [-1 0 3 0 -8 0 21 0 -45 0 91 0 -191 0 643 1024 643 0 -191 0 91 0 -45 0 21 0 -8 0 3 0 -1];
h0_coeff = h0_coeff_g/2050;
n=0:[length(h0_coeff)-1];
h1_coeff = (-1).^n.*h0_coeff;
len=length(h0_coeff);
MSError_frame = 0;
 
%% Analysis Part
for k=1:1
    % Channel 1
    h0_1=upsample(h0_coeff,2);
    h0_2=upsample(h0_coeff,4);
    h0_f= conv(conv(h0_coeff,h0_1),h0_2);
    h0 =downsample(h0_f,8);
    v0_a1=conv(h0_f,frame(k,:));
    v0_ap = downsample(v0_a1,8);
    
    % Channel_2
    h1_1=upsample(h0_coeff,2);
    h1_2=upsample(h1_coeff,4);
    h1_f= conv(conv(h0_coeff,h1_1),h1_2);
    h1 =downsample(h1_f,8);
    v1_a1=conv(h1_f,frame(k,:));
    v1_ap=downsample(v1_a1,8);
    
    % Channel_3
    h2_1=upsample(h1_coeff,2);
    h2_f= conv(h0_coeff,h2_1);
    h2 =downsample(h2_f,4);
    v2_a1=conv(h2_f,frame(k,:));
    v2_ap=downsample(v2_a1,4);
    
    % Channel_4
    h3 =h1_coeff;
    v3_a1=conv(h3,frame(k,:));
    v3_ap=downsample(v3_a1,2);
end
  %% Synthesis Part

    f0_coeff = h0_coeff;
    f1_coeff = -h1_coeff;
    % Channel 1
    f0_u1=upsample(f0_coeff,4);
    f0_u2=upsample(f0_coeff,2);
    f0 = conv(conv(f0_u1,f0_u2),f0_coeff);
    in_0=upsample(v0_ap,8);
    x0_rcon=conv(in_0,f0);
    
    % Channel 2
    f1_u1 = upsample(f1_coeff,4);
    f1_u2 = upsample(f0_coeff,2);
    f1 = conv(conv(f1_u1,f1_u2),f0_coeff);
    in_1=upsample(v1_ap,8);
    x1_rcon=conv(in_1,f1);
    
    % Channel 3
    f2_u1=upsample(f1_coeff,2);
    f2 = conv(f2_u1,f0_coeff);
    in_2=upsample(v2_ap,4);
    x2_rcon=conv(in_2,f2);
    x2_rcon_zero=[x2_rcon, zeros(1,(length(x1_rcon)-length(x2_rcon)))];
    
    % Channel 4
    in_3=upsample(v3_ap,2);
    x3_rcon=conv(in_3,f1_coeff);
    x3_rcon_zero=[x3_rcon, zeros(1,(length(x1_rcon)-length(x3_rcon)))];
    
    xr =x0_rcon+x1_rcon+x2_rcon_zero+x3_rcon_zero;

     %% MSE calculation
    sum=0;
    for k=1:1024
        MSEr2(k)=abs(data(k)-xr(k))^2/abs(data(k))^2;
        sum=sum+MSEr2(k);
    end
    MSError2=sqrt(sum)/1024

subplot(211);plot(xr,'g');
subplot(212);plot(data,'r');