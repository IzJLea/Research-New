function k=kCO2(Temp)
% CO2 thermal conductivity for use in interpolation

    kCO2=[0.01051 0.01456 0.01858 0.02257 0.02652 0.03044 0.03814 0.04565 0.05293 0.08491 0.10688 0.11522];

% CO2 thermal conductivity Temperatures for use in interpolation

    kCO2Temp=[-50 0 50 100 150 200 300 400 500 1000 1500 2000];
    
    k=interp1(kCO2Temp,kCO2,Temp);
end
