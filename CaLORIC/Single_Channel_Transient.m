% Single Channel temperature Profile calculation
tic
Qchannel=0.15; %MW

Hchannel=2.1; %m

Tenter=100;

PSH=114;

PVH=134;

time=3000;

div=0.1;

Lchannel=6;

Res=Single_Channel(Qchannel,Hchannel,Tenter,PSH,PVH,time,div,Lchannel);

t=zeros(1,length(Res));

for i=2:length(t)
    
    t(1,i)=t(1,i-1)+div;
end

plot(t,Res(3,1:length(Res)),t,Res(1,1:length(Res)));
toc