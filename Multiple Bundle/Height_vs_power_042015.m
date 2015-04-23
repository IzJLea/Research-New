%% This script is to calculate the flow, 
tic
time=3000; %time is seconds

div=0.1; %time division (fraction of seconds)

PVH=134; %kPa Vent header pressure

PSH=114; %kPa Supply Vent pressure

Tenter=100; 

Qchannel=linspace(0.05,0.25,5);

Hchannel=linspace(1,11,24);
% 
% Qchannel=0.05;
% 
% Hchannel=1;

Tfuel=zeros(length(Qchannel),length(Hchannel));

Tclad=zeros(length(Qchannel),length(Hchannel));

Tvap=zeros(length(Qchannel),length(Hchannel));

TPT=zeros(length(Qchannel),length(Hchannel));

TCT=zeros(length(Qchannel),length(Hchannel));

mflow=zeros(length(Qchannel),length(Hchannel));

for i=1:length(Qchannel)
    for j=1:length(Hchannel)
        Iteration=Max_Temp_Single_Channels(Qchannel(1,i),Hchannel(1,j),Tenter,PSH,PVH,time,div);
        
        Tfuel(i,j)=Iteration(1,1);
        
        Tclad(i,j)=Iteration(2,1);
        
        Tvap(i,j)=Iteration(3,1);
        
        TPT(i,j)=Iteration(4,1);
        
        TCT(i,j)=Iteration(5,1);
        
        mflow(i,j)=Iteration(6,1);
    end
end

filename = 'Calandria_Temperatures.xlsx';
xlswrite('Calandria_Temperatures.xlsx',Tfuel,1,'B1:Y5');
xlswrite('Calandria_Temperatures.xlsx',Tclad,2,'B1:Y5');
xlswrite('Calandria_Temperatures.xlsx',Tvap,3,'B1:Y5');
xlswrite('Calandria_Temperatures.xlsx',TPT,4,'B1:Y5');
xlswrite('Calandria_Temperatures.xlsx',TCT,5,'B1:Y5');


toc
        

