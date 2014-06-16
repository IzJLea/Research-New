Qchannel=5.580; % channel power in MW

Psat=10; %channel pressure in Mpa

hout=1355.1; % exit enthalpy in kJ/kg

Msp=25; % mass flowrate 

alpha=LevyCandu(Psat,hout,Qchannel,Msp);

display(alpha,'Void Fraction');
