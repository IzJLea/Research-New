Qin=7.5;

Tenter=267.2;

Tbulk=310.950736486010;

Tmod=60;

hbulk=78000;

%% Initial dimension Definition
Lchannel=5.94; % m

DPT=103.38e-3; % m

tPT=4.24e-3; %m

DCT=128.69e-3; %m

tCT=1.4e-3; %m


%% Property info 
% CO2 thermal conductivity

kCO2=[14.60e-3 16.23e-3 17.87e-3 19.52e-3 21.18e-3 22.84e-3 27.00e-3 31.12e-3 35.20e-3 39.23e-3];

% CO2 thermal conductivity Temperatures

kCO2Temp=[0 20 40 60 80 100 150 200 250 300];

% moderator themperatures

Tmodref=[273.15 275 280 285 290 295 300 305 310 315 320 325 330 335 340 345 350 355 360 365 370 373.15]-273.15;

% moderator dynamic viscosity

mumod=[1750e-6 1652e-6 1422e-6 1225e-6 1080e-6 959e-6 855e-6 769e-6 695e-6 631e-6 577e-6 528e-6 489e-6 453e-6 420e-6 389e-6 365e-6 343e-6 324e-6 306e-6 289e-6 279e-6];

% moderator thermal conductivity

kmod=[569e-3 574e-3 582e-3 590e-3 598e-3 606e-3 613e-3 620e-3 628e-3 634e-3 640e-3 645e-3 650e-3 656e-3 660e-3 668e-3 668e-3 671e-3 674e-3 677e-3 679e-3 680e-3];

% moderator Prandtl number

Prlmod=[12.99 12.22 10.26 8.81 7.56 6.62 5.83 5.20 4.62 4.16 3.77 3.42 3.15 2.88 2.66 2.45 2.29 2.14 2.02 1.91 1.80 1.76];



% moderator expansion coefficient

betamod=[-68.05e-6 -32.74e-6 46.04e-6 114.1e-6 174.0e-6 227.5e-6 276.1e-6 320.6e-6 361.9e-6 400.4e-6 436.7e-6 471.2e-6 504.0e-6 535.5e-6 566.0e-6 595.4e-6 624.2e-6 652.3e-6 697.9e-6 707.1e-6 728.7e-6 750.1e-6];

% liquid heat capacity

Cpmod=[4.217 4.211 4.198 4.189 4.184 4.181 4.179 4.178 4.178 4.179 4.180 4.182 4.184 4.186 4.188 4.191 4.195 4.199 4.203 4.209 4.214 4.217]/1000;

%% System properties
% pressure tube thermal conductivity

kzircPT=(7.51+(0.362e-3*Tbulk)-(0.618e-7*Tbulk^2)+(0.718e-11*Tbulk^3));

% calandria tube thermal conductivity

kzircCT=(7.51+(0.362e-3*Tmod)-(0.618e-7*Tmod^2)+(0.718e-11*Tmod^3));

% CO2 thermal conductivity

Tsystem=[Tbulk Tmod];

TevalCO2=mean(Tsystem);

kCO2sys=interp1(kCO2Temp,kCO2,TevalCO2);

% moderator Prandtl number

Prlmodsys=interp1(Tmodref,Prlmod, Tmod);

% moderator dynamic viscosity

mumodsys=interp1(Tmodref,mumod,Tmod);

% Thermal expansion coefficient

betamodsys=interp1(Tmodref,betamod,Tmod);

% moderator heat capacity

Cpmodsys=interp1(Tmodref,Cpmod,Tmod);

% moderator thermal conductivity

kmodsys=interp1(Tmodref, kmod, Tmod);




%% thermal resistance calculation

ript=DPT/2;

ropt=ript+tPT;

rict=DCT/2;

roct=rict+tCT;

Rtotal=(1/(hbulk*DPT*pi()))+(log(ropt/ript)/(2*pi()*kzircPT*Lchannel))+(log(rict/ropt)/(2*pi()*kCO2sys*Lchannel))+(log(roct/rict)/(2*pi()*kzircCT*Lchannel));

A=kmodsys*pi()*Rtotal*Lchannel;

B=0.60;

C=0.387/(1+((0.559/Prlmodsys)^(9/16)))^(8/27);

D=9.81*betamodsys*Lchannel^3*Cpmodsys/mumodsys/kmodsys;

a=(A*B^2)+1;

b=2*A*C*D^(1/6);

c=A*C^2*D^(1/3);

n1=7/6;

n2=4/3;

dT=Tmod-Tbulk;

fx=@(x) (a*x)+(b*x^(n1))+(c*x^(n2))+dT;

Tdiff=fzero( fx,50);

Ts=Tdiff+Tmod;

Qloss=(Tbulk-Ts)/Rtotal/1000000;