Qchannel=5:0.1:12;

Results=zeros(1,length(Qchannel));

Tin=267.2;

Tmod=60; %C

Pout1=10;

for resind=1:length(Qchannel)

    

    Results(1:12,resind)=single_channel_Lockhart_Martinelli_Qloss_LPXSTEAM4(Qchannel(resind),Tin, Tmod, Pout1);
end


% Results2=zeros(1,length(Qchannel));
% 
% Pout2=10;
% 
% for resind=1:length(Qchannel)
% 
%     
% 
%     Results2(1:12,resind)=single_channel_Lockhart_Martinelli_Qloss_LPXSTEAM3(Qchannel(resind),Tin,Tmod,Pout2);
% end

% Results3=zeros(1,length(Qchannel));
% 
% Pout3=11.38; 
% 
% 
% for resind=1:length(Qchannel)
% 
%     
% 
%     Results3(1:12,resind)=single_channel_Lockhart_Martinelli_Qloss_LPXSTEAM2(Qchannel(resind),Tin,Tmod,Pout3);
% end


plot(Results(1,1:length(Qchannel)),Results(2,1:length(Qchannel)));
