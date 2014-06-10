%% Initial single pass k calculation- 
% pressure differences in each pass section
%
[~, ~, raw] = xlsread('C:\Users\Izaak\Documents\Research\Generic C9 HTS Data.xlsx','FLOW-DP','I13:I17');


DP = reshape([raw{:}],size(raw));

%%calculate values for k for each section

% normal mass flowrate

M=25.8; % kg/s

Rho=780.6; %kg/m^3


keff=DP/(M^2);

reff=keff*Rho;

display(reff);

%% system information input

Minput=input('Flowrate');

Rhoinl=input('Density of liquid');

Rhoinv=input('Density of gas');

Hlin=input('Liquid enthalpy');

Hvin=input('Gas enthalpy');

Hin=input('Enthalpy of 2 phase flow');

%% Single Phase Pressure drop

DPm=keff*Minput^2;

DPd=DPm*Rho/Rhoinl;

DPt=sum(DPd);

display(DPd, 'Pressure drop per section');

display(DPt, 'Total single phase Pressure drop:');

%% Two Phase Pressure drop

x=(Hin-Hvin)/(Hlin-Hvin);

a0=[0;x;Rhoinl;Rhoinv];

options = optimoptions('fsolve','Display','iter');

[a,fval]= fsolve(Levy, a0, options);


display(a,'Void Fraction= ');

