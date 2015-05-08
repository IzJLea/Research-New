%% This script calculates the Maximum Temperature of the channel elements for liquid heights from 1 to 12 meters and decay power of 0.05 to 0.25 MW

tic
time=3000; %time is seconds

div=0.1; %time division (fraction of seconds)

PVH=134; %kPa Vent header pressure

PSH=114; %kPa Supply Vent pressure

Tenter=100; 

Qchannel=linspace(0.1,0.5,9);

Hchannel=linspace(2,3,11);

Lchannel=6;
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
        Iteration=Max_Temp_Single_Channels(Qchannel(1,i),Hchannel(1,j),Tenter,PSH,PVH,time,div,Lchannel);
        
        Tfuel(i,j)=Iteration(1,1);
        
        Tclad(i,j)=Iteration(2,1);
        
        Tvap(i,j)=Iteration(3,1);
        
        TPT(i,j)=Iteration(4,1);
        
        TCT(i,j)=Iteration(5,1);
        
        mflow(i,j)=Iteration(6,1);
        
        clear Iteration
    end
end

filename = 'Calandria_Temperatures_2to3m.xlsx';
xlswrite('Calandria_Temperatures_2to3m.xlsx',Tfuel,1,'B1:L9');
xlswrite('Calandria_Temperatures_2to3m.xlsx',Tclad,2,'B1:L9');
xlswrite('Calandria_Temperatures_2to3m.xlsx',Tvap,3,'B1:L9');
xlswrite('Calandria_Temperatures_2to3m.xlsx',TPT,4,'B1:L9');
xlswrite('Calandria_Temperatures_2to3m.xlsx',TCT,5,'B1:L9');
xlswrite('Calandria_Temperatures_2to3m.xlsx',mflow,6,'B1:L9');
xlswrite('Calandria_Temperatures_2to3m.xlsx',Qchannel.',7,'C4:C12');

xlswrite('Calandria_Temperatures_2to3m.xlsx',Hchannel,7,'D3:N3');

toc
        

