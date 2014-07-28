Qchannel=5:0.1:12;

Results=zeros(1,length(Qchannel));

Tin=267.2;

Tmod=60; %C

Pout1=10.5;

for resind=1:length(Qchannel)

    

    Results(1:12,resind)=single_channel_Lockhart_Martinelli_Qloss_LPXSTEAM(Qchannel(resind),Tin, Tmod, Pout1);
end


Results2=zeros(10,length(Qchannel));

Pout2=10;

for resind=1:length(Qchannel)

    

    Results2(1:10,resind)=single_channel_Lockhart_Martinelli_Qloss_LP2(Qchannel(resind),Tin,Tmod,Pout1);
end



plot(Results(1,1:length(Qchannel)),Results(2,1:length(Qchannel)),Results2(1,1:length(Qchannel)),Results2(2,1:length(Qchannel)));
