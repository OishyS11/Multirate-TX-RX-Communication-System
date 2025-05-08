clc;
clear all;
close all;

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


V0_scdt=importdata("C:\Users\DELL\ScrambledV0_GroupID_1A_FreqID_60_SNR_10.000000_Foff_25.000000_Nframe_1_Nstart_1.txt");
V0_scd1=V0_scdt.rowheaders;
V0_scd1=str2double(V0_scd1);
V0_scd2=[V0_scdt.data];
%V0_scdq=[V0_scd2{:}]; 
V0_scf=V0_scd1+j*V0_scd2;
V1_scdt=importdata("C:\Users\DELL\ScrambledV1_GroupID_1A_FreqID_60_SNR_10.000000_Foff_25.000000_Nframe_1_Nstart_1.txt");
V1_scd1=V1_scdt.rowheaders;
V1_scd1=str2double(V1_scd1);
V1_scd2=[V1_scdt.data];
%V0_scdq=[V0_scd2{:}]; 
V1_scf=V1_scd1+j*V1_scd2;
V2_scd3=importdata("C:\Users\DELL\ScrambledV2_GroupID_1A_FreqID_60_SNR_10.000000_Foff_25.000000_Nframe_1_Nstart_1.txt");
V3_scd4=importdata("C:\Users\DELL\ScrambledV3_GroupID_1A_FreqID_60_SNR_10.000000_Foff_25.000000_Nframe_1_Nstart_1.txt");

for i =1:length(V2_af3)
    V2_scf(i)=V2_af3(i,1)+i*V2_af3(i,2);
end
for i =1:length(V3_af4)
    V3_scf(i)=V3_af4(i,1)+i*V3_af4(i,2);
end


V = [bitsll(1,10)+26 bitsll(1,11)+(3*26) bitsll(1,12)+(5*26) bitsll(1,13)+(7*26)];

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
main_str=importdata("C:\Users\DELL\ScrambledV2_GroupID_1A_FreqID_60_SNR_10.000000_Foff_25.000000_Nframe_1_Nstart_1.txt")
for i =1:length(main_str)
    main_v2_sc1(i)=main_str(i,1)+i*main_str(i,2);
end
%Receiver Side_Descrambling and Systhesis
        V0_dsc=real(V0_sc)*PN1_key(1)+j*imag(V0_sc)*PN2_key(1);
        V1_dsc=real(V1_sc)*PN1_key(2)+j*imag(V1_sc)*PN2_key(2);
        V2_dsc=real(V2_sc)*PN1_key(3)+j*imag(V2_sc)*PN2_key(3);
        V3_dsc=real(V3_sc)*PN1_key(4)+j*imag(V3_sc)*PN2_key(4);

subplot(211);plot(V3_sc,'g');
subplot(212);plot(V3_scf,'r');
