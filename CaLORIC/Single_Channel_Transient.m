% Single Channel temperature Profile calculation

Qchannel=0.15; %MW

Hchannel=11; %m

Tenter=100;

PSH=10000;

PVH=10050;

time=3000;           

div=0.01;

Lchannel=6;

Res=Single_Channel(Qchannel,Hchannel,Tenter,PSH,PVH,time,div,Lchannel);

t=zeros(1,length(Res));

for i=2:length(t)
    
    t(1,i)=t(1,i-1)+div;
end

plot(t,Res(3,1:length(Res)),t,Res(1,1:length(Res)));
