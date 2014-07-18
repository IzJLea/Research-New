Qchannel=5:0.1:12;

Results=zeros(9,length(Qchannel));

Tin=267.2;

for resind=1:length(Qchannel)

    

    Results(1:9,resind)=single_channel(Qchannel(resind),Tin);
end

% plot(Results(1,1:length(Qchannel)),Results(2,1:length(Qchannel)));

Resind=1:10:length(Qchannel);

Rchart=Results(1:9,Resind);


Results2=zeros(9,length(Qchannel));

Tin=267.2;

for resind=1:length(Qchannel)

    

    Results2(1:9,resind)=single_channel_alt(Qchannel(resind),Tin);
end

plot(Results(1,1:length(Qchannel)),Results(2,1:length(Qchannel)),Results2(1,1:length(Qchannel)),Results2(2,1:length(Qchannel)));
