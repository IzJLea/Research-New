function k=kUO2(Temp)

T=200:10:3000;

t=T./1000;

klist=zeros(1,length(T));

for n=(1:length(T))
    
    klist(n)=(115.8/(7.5408+(17.692*t(n))+(3.6142*t(n)^2)))+(7410.5*t(n)^(-5/2))*exp(-16.35/t(n));
    
end

k=interp1(T,klist,(Temp+273.15))/1000;

%% Will give results in kW/m.K
% Temp is in degrees Celsius

