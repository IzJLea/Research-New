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

Minput=25.8;    %kg/s coolant flowrate

Rhoinl=780.6;   %kg/m^3 liquid water density at saturation conditions

Rhoinv=52.6;    %kg/m^3 water vapour density at saturation conditions

Hlin=1450;  %kJ/kg enthalpy of liquid water at saturation conditions

Hvin=2706;  %kJ/kg enthalpy of gaseous water at saturation conditions

Hin=1550;   %kJ/kg enthalpy of fluid entering channel

%% Single Phase Pressure drop

DPm=keff*Minput^2;

DPd=DPm*Rho/Rhoinl;

DPt=sum(DPd);

display(DPd, 'Pressure drop per section');

display(DPt, 'Total single phase Pressure drop:');

%% Two Phase Pressure drop

Dens=Rhoinl/Rhoinv;

qual=(Hin-Hlin)/(Hvin-Hlin);   % vapour quality (mass)

p=(1-qual)^2;

t=qual^2*Dens;

m=((1-qual)^2)*0.5;

syms x

void=solve((p*(1-x))-(t*x)-(m/((1-x)^2))-0.5);

display(void);

%a=Levyvoid(x,Rhoinv,Rhoinl);

%[result, fval, exitflag, output]=fsolve(@Levy, guess);

%result;

%fval;

%Levy(guess);

%output;
