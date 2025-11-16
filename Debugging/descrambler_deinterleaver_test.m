clc;
clear all;
close all;

str = importdata("C:\Users\DELL\Rn_GroupID_1A_FreqID_60_SNR_10.000000_Foff_25.000000_Nframe_1_Nstart_1.txt");
for i =1:length(str)
    R(i)= str(i,1)+i*str(i,2);
end
np=0:(length(R)-1);

V0_af1=importdata("C:\Users\DELL\V0_GroupID_1A_FreqID_60_SNR_10.000000_Foff_25.000000_Nframe_1_Nstart_1.txt");
V0_a1=V0_af1.rowheaders;
V0_a1=str2double(V0_a1);
V0_a2=[V0_af1.data];
%V0_scdq=[V0_scd2{:}]; 
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
V = [bitsll(1,10)+26 bitsll(1,11)+(3*26) bitsll(1,12)+(5*26) bitsll(1,13)+(7*26)];

    for kc=1:4
        PN1(kc)=xor(xor(bitget(V(kc),1),bitget(V(kc),4)),bitget(V(kc),26));
        PN2(kc)=xor(xor(xor(xor(bitget(V(kc),1),bitget(V(kc),2)),bitget(V(kc),3)),bitget(V(kc),4)),bitget(V(kc),26));

        PN1_key(kc)=1-2*PN1(kc);
        PN2_key(kc)=1-2*PN2(kc);
    end
%% Burst Formation
        freq_ID=60;
        freq_omega=(2*pi/128)*freq_ID;

 %% Receiever Data Side
        dft_value=fft(R,128);
        [M,i]=max(dft_value);
        detct_freq_offset=(i-1)*(2*pi/128)-freq_omega;
        detect_freq_Hz=(detct_freq_offset*128)/(2*pi);
        offset_comp=cos(detct_freq_offset*np)-j*sin(detct_freq_offset*np);
        Sp=R.*offset_comp;
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
        % Receiver Side_Descrambling
        V0_dsc=real(V0_scp)*PN1_key(1)+j*imag(V0_scp)*PN2_key(1);
        V1_dsc=real(V1_scp)*PN1_key(2)+j*imag(V1_scp)*PN2_key(2);
        V2_dsc=real(V2_scp)*PN1_key(3)+j*imag(V2_scp)*PN2_key(3);
        V3_dsc=real(V3_scp)*PN1_key(4)+j*imag(V3_scp)*PN2_key(4);

        subplot(211);plot(v3_a);
        subplot(212);plot(V3_dsc);
 
