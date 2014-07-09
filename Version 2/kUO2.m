function k=kUO2(Temp)

t=(Temp+273.15)/1000;

k=((115.8/(7.5408+(17.692*t)+(3.6142*t^2)))+(7410.5*t^(-5/2))*exp(-16.35/t))/1000;
    




%% Will give results in kW/m.K
% Temp is in degrees Celsius

