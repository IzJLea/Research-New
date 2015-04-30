function Cp=CpUO2(Tfuel)

tau=(Tfuel+273.15)/1000;

Cp=(52.1743+(87.951*tau)-(85.2411*tau^2)+(31.542*tau^3)-(2.6334*tau^4)-(0.71391*tau^-2))/270.03*1000; %J/kg.K