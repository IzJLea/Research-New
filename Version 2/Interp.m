Tvap=120;

Tmod=60;

kCO2=[14.60e-3 16.23e-3 17.87e-3 19.52e-3 21.18e-3 22.84e-3 27.00e-3 31.12e-3 35.20e-3 39.23e-3];

    % CO2 thermal conductivity Temperatures

    kCO2Temp=[0 20 40 60 80 100 150 200 250 300];

    % CO2 thermal conductivity
    
    TevalCO2=((Tvap(1,1)+Tmod))/2;    

    kCO2sys=interp1(kCO2Temp,kCO2,TevalCO2);
