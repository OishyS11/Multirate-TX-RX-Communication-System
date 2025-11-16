clc;
clear all;
close all;
n=0:1023;
Ts=1/1000;
N=1;
frame=zeros(N,length(n));

% Transmitter Side

%Band Pass Data
% for k=1:N
%  xi = 1344*cos(0.06*pi*n)+864*cos(0.18*pi*n)+8543*cos(0.38*pi*n)-43*cos(0.8*pi*n);
% xq = 1344*sin(0.06*pi*n)+864*sin(0.18*pi*n)+8543*sin(0.38*pi*n)-43*sin(0.8*pi*n);
%     data=xi+j*xq;
%     frame(k,:)=data;
% end
str = importdata("C:\Users\DELL\Input_GroupID_1A_FreqID_60_SNR_10.000000_Foff_25.000000_Nframe_1_Nstart_1.txt");

for i =1:1024
    data(i)=str(i,1)+i*str(i,2);
end

str1 = importdata("C:\Users\DELL\Rn_GroupID_1A_FreqID_60_SNR_10.000000_Foff_25.000000_Nframe_1_Nstart_0.txt");

for i =1:1152
    Re(i)=str1(i,1)+i*str1(i,2);
end
 frame=data;
V=zeros(1,4);
% Random Data
% for k=1:N
%     xi=randn(1,1024);
%     xq=randn(1,1024);
%     data=xi+j*xq;
%     frame(k,:)=data;
% end

h0_coeff_g = [-1 0 3 0 -8 0 21 0 -45 0 91 0 -191 0 643 1024 643 0 -191 0 91 0 -45 0 21 0 -8 0 3 0 -1];
h0_coeff = h0_coeff_g/2050;
n=0:[length(h0_coeff)-1];
h1_coeff = (-1).^n.*h0_coeff;
len=length(h0_coeff);
MSError_frame = 0;
SNR=0:20;
% % Analysis Part
%   
% for k1=1:length(SNR)
for k=1:N
%     Channel 1
    h0_1=upsample(h0_coeff,2);
    h0_2=upsample(h0_coeff,4);
    h0_f= conv(conv(h0_coeff,h0_1),h0_2);
    h0 =downsample(h0_f,8);
    v0_a1=conv(h0_f,frame(k,:));
    v0_a = downsample(v0_a1,8);
    
%     Channel_2
    h1_1=upsample(h0_coeff,2);
    h1_2=upsample(h1_coeff,4);
    h1_f= conv(conv(h0_coeff,h1_1),h1_2);
    h1 =downsample(h1_f,8);
    v1_a1=conv(h1_f,frame(k,:));
    v1_a=downsample(v1_a1,8);
    
%     Channel_3
    h2_1=upsample(h1_coeff,2);
    h2_f= conv(h0_coeff,h2_1);
    h2 =downsample(h2_f,4);
    v2_a1=conv(h2_f,frame(k,:));
    v2_a=downsample(v2_a1,4);
    
%     Channel_4
    h3 =h1_coeff;
    v3_a1=conv(h3,frame(k,:));
    v3_a=downsample(v3_a1,2);
end
% 
    V = [bitsll(k,10)+26 bitsll(k,11)+(3*26) bitsll(k,12)+(5*26) bitsll(k,13)+(7*26)];

    for kc=1:4
        PN1(kc)=xor(xor(bitget(V(kc),1),bitget(V(kc),4)),bitget(V(kc),26));
        PN2(kc)=xor(xor(xor(xor(bitget(V(kc),1),bitget(V(kc),2)),bitget(V(kc),3)),bitget(V(kc),4)),bitget(V(kc),26));

        PN1_key(kc)=1-2*PN1(kc);
        PN2_key(kc)=1-2*PN2(kc);
    end
        V0_sc=real(v0_a)*PN1_key(1)+j*imag(v0_a)*PN2_key(1);
        V1_sc=real(v1_a)*PN1_key(2)+j*imag(v1_a)*PN2_key(2);
        V2_sc=real(v2_a)*PN1_key(3)+j*imag(v2_a)*PN2_key(3);
        V3_sc=real(v3_a)*PN1_key(4)+j*imag(v3_a)*PN2_key(4);
%  
%         %% Interleaving
%         S = zeros(1,length(data));
%         for count=1:length(data)
%             if (rem(count,8)==1)
%               S(count)=V0_sc(((count-1)/8)+1);
%             elseif (rem(count,8)==2)
%               S(count)=V3_sc(((count-2)/8)+1);
%             elseif (rem(count,8)==3)
%               S(count)=V1_sc(((count-3)/8)+1);
%             elseif (rem(count,8)==4)
%               S(count)=V3_sc(((count-4)/8)+1);
%             elseif (rem(count,8)==5)
%               S(count)=V2_sc(((count-5)/8)+1);
%             elseif (rem(count,8)==6)
%               S(count)=V3_sc(((count-6)/8)+1);
%             elseif (rem(count,8)==7)
%               S(count)=V2_sc(((count-7)/8)+1);
%             else
%               S(count)=V3_sc(((count-8)/8)+1);
%             end
%             end
%       
%         % Burst Formation
        freq_ID=60;
        freq_omega=(2*pi/128)*freq_ID;
%         n=0:127;
%         burst= cos(freq_omega*n)+j*sin(freq_omega*n);
%         T=[burst S];
        np=0:(length(Re)-1);
%         freq_offset=-1.832+(1.832*2)*rand();
%         offset_signal=cos(2*pi*ft*np)+j*sin(2*pi*ft*np);
%         freq_offset=25;
%         offset_signal=cos(2*pi*freq_offset*np)+j*sin(2*pi*freq_offset*np);
%         corrupted=T.*offset_signal;
%         noise=wgn(1,length(corrupted),20);
%         R=corrupted+noise;
%         
%         R=corrupted+noise;
        
   %% Receiever Data Side
        dft_value=fft(Re,128);
        [M,i]=max(dft_value);
        detct_freq_offset=(i-1)*(2*pi/128)-freq_omega;
        detect_freq_Hz=(detct_freq_offset*128)/(2*pi);
        offset_comp=cos(detct_freq_offset*np)-j*sin(detct_freq_offset*np);
        Sp=Re.*offset_comp;
   %% De-Interleaver
        for cout=1:length(Sp)
            if(rem(cout,8)==1)
                V0_scp(((cout-1)/8)+1)=Sp(cout);
            end
            if(rem(cout,8)==3)
                V1_scp(((cout-3)/8)+1)=Sp(cout);
            end
            if(rem(cout,8)==5) 
                V2_scp(((cout-5)/8)+1)=Sp(cout);
            end
            if(rem(cout,8)==7) 
                V2_scp(((cout-7)/8)+1)=Sp(cout);
            end
            if(rem(cout,2)==0)
                V3_scp(((cout-2)/2)+1)=Sp(cout);
            end
        end
        % Receiver Side_Descrambling and Systhesis
        V0_dsc=real(V0_scp)*PN1_key(1)+j*imag(V0_scp)*PN2_key(1);
        V1_dsc=real(V1_scp)*PN1_key(2)+j*imag(V1_scp)*PN2_key(2);
        V2_dsc=real(V2_scp)*PN1_key(3)+j*imag(V2_scp)*PN2_key(3);
        V3_dsc=real(V3_scp)*PN1_key(4)+j*imag(V3_scp)*PN2_key(4);

    %% Synthesis Part

    f0_coeff = h0_coeff;
    f1_coeff = -h1_coeff;
    % Channel 1
    f0_u1=upsample(f0_coeff,4);
    f0_u2=upsample(f0_coeff,2);
    f0 = conv(conv(f0_u1,f0_u2),f0_coeff);
    in_0=upsample(V0_dsc,8);
    x0_rcon=conv(in_0,f0);
    
    % Channel 2
    f1_u1 = upsample(f1_coeff,4);
    f1_u2 = upsample(f0_coeff,2);
    f1 = conv(conv(f1_u1,f1_u2),f0_coeff);
    in_1=upsample(V1_dsc,8);
    x1_rcon=conv(in_1,f1);
    
    % Channel 3
    f2_u1=upsample(f1_coeff,2);
    f2 = conv(f2_u1,f0_coeff);
    in_2=upsample(V2_dsc,4);
    x2_rcon=conv(in_2,f2);
    x2_rcon_zero=[x2_rcon, zeros(1,(length(x1_rcon)-length(x2_rcon)))];
    
    % Channel 4
    in_3=upsample(V3_dsc,2);
    x3_rcon=conv(in_3,f1_coeff);
    x3_rcon_zero=[x3_rcon, zeros(1,(length(x1_rcon)-length(x3_rcon)))];
    
    xr =x0_rcon+x1_rcon+x2_rcon_zero+x3_rcon_zero;
    %% MSE calculation
    sum=0;
    for k=1:1024
        MSEr2(k)=abs(data(k)-xr(k))^2/abs(data(k))^2;
        sum=sum+MSEr2(k);
    end
    MSError2=sqrt(sum)/1024;
    MSError_frame=MSError_frame+MSError2 ;
     
    MSError_avg = MSError_frame/N
 
%     subplot(211);plot(data);
%     subplot(212);plot(xr);
%     figure;
    plot(Re)
    
 


