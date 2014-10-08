function Qloss=Qloss_single_channel2(Tbulk,Tmod,hbulk,hmod,Lbund)
%% Initial dimension Definition

Lchannel=Lbund; % m

DPT=0.10338;

tPT=0.00424;

DCT=0.12869;

tCT=0.0014;

ript=DPT/2;

ropt=ript+tPT;

rict=DCT/2;

roct=rict+tCT;

Tsystem=[Tbulk Tmod]; %C

TevalCO2=mean(Tsystem); %C

Tbulk=Tbulk+273.15;



Tmod=Tmod+273.15;

%% System properties
% pressure tube thermal conductivity

kzircPT=12.767-(5.4348e-4*Tbulk)+(8.9818e-6*Tbulk^2); %W/m.K

% calandria tube thermal conductivity

kzircCT=12.767-(5.4348e-4*Tbulk)+(8.9818e-6*Tbulk^2);  %W/m.K

% CO2 thermal conductivity

kCO2=[14.60e-3 16.23e-3 17.87e-3 19.52e-3 21.18e-3 22.84e-3 27.00e-3 31.12e-3 35.20e-3 39.23e-3];

% CO2 thermal conductivity Temperatures

kCO2Temp=[0 20 40 60 80 100 150 200 250 300];

% CO2 thermal conductivity

kCO2sys=interp1(kCO2Temp,kCO2,TevalCO2);

%% Total Termal Resistance

R1=1/(hbulk*DPT*pi()*Lchannel);

R2=log(ropt/ript)/(2*pi()*kzircPT*Lchannel);

R3=log(rict/ropt)/(2*pi()*kCO2sys*Lchannel);

R4=log(roct/rict)/(2*pi()*kzircCT*Lchannel);

R5=1/(hmod*DCT*pi()*Lchannel);

Rtotal=((1/R1)+(1/R2)+(1/R3)+(1/R4)+(1/R5))^-1;

Qlossval=(Tbulk-Tmod)/Rtotal;

Qloss=[Qlossval R1 R2 R3 R4 R5 Rtotal];

















