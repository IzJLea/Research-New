function Qloss=Qloss_SS(TPT,Tmod,hmod)

Lbund=0.5

Tmean=(TPT+Tmod)/2;

% calandria tube thermal conductivity

kzircCT=12.767-(5.4348e-4*Tmod)+(8.9818e-6*Tmod^2);  %W/m.K

% CO2 thermal conductivity

kCO2=[14.60e-3 16.23e-3 17.87e-3 19.52e-3 21.18e-3 22.84e-3 27.00e-3 31.12e-3 35.20e-3 39.23e-3];

% CO2 thermal conductivity Temperatures

kCO2Temp=[0 20 40 60 80 100 150 200 250 300];

% CO2 thermal conductivity
if Tmean<300
    kCO2sys=interp1(kCO2Temp,kCO2,Tmean);
else
    Tmean=299;
    
    kCO2sys=interp1(kCO2Temp,kCO2,Tmean);
end


PTo=0.1118;

CTi=0.1287;

CTo=CTi+(2*0.0014);

R1=log(CTi/PTo)/(2*pi()*Lbund*kCO2sys);

R2=log(CTo/CTi)/(2*pi()*Lbund*kzircCT);

R3=1/(hmod*Lbund*CTo*pi());

Rtotal=((1/R1)+(1/R2)+(1/R3))^-1;

Qloss=(TPT-Tmod)/Rtotal/1000000;




