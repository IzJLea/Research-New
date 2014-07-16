Qchannel=5.5:0.1:12;

Results=zeros(6,length(Qchannel));

Tin=267.2;

for ind=1:length(Qchannel)

    

    Results(1:6,ind)=single_channel(Qchannel(ind),Tin);
    
end



plot(Results(1,1:length(Results)),Results(2,1:length(Results)));