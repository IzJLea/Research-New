function xprime = Levyxprime(Psat,hout,Qchannel,Msp)

%System values

Achannel=0.0034; %m^2

Achannelimp=Achannel*3.2808^2; %ft^2

Aht=9.1224; %m^2 heat transfer area

Ahtimp=Aht*3.2808^2; %ft^2

Dh=0.0074; %m hydraulic diameter

Dhimp=Dh*3.2808; %ft 

g=32.174*3600^2; %gravitational acceleration

gc=1/32.174; % conversion from lbf to lbm

relrough=10^-4; % relative roughness (drawn tubing)

% conversions

Mspimp=Msp*2.2046226*3600; % single phase mass flow

Psatimp=Psat*1000/6.894757; %psia

houtimp=hout/2.326; %Btu/lbm

Qchannelimp=Qchannel*1000*3412.14;

%% Data for Calculated system properties (surface tension liquid and vapour enthalpies, saturation temperature, liquid density, specific heat
% dynamic viscosity, prandtl number)
%Pressure (psia)

Pimp=[24.97 29.82 35.42 41.85 49.18 57.53 66.98 89.60 117.93 152.92 195.60 241.1 422.1 680.0 1046.7 1541 2210 3090];

% Temperature F

Temperatureimp=[240 250 260 270 280 290 300 320 340 360 380 400 450 500 550 600 650 700];

% Liquid Density (lbm/ft^3)

rholimp=[59.09 58.82 58.53 58.24 57.94 57.63 57.31 56.65 55.95 55.22 54.46 53.65 51.46 48.95 45.96 42.32 37.31 27.28];

% vapour Density (lbm/ft^3)

rhovimp=[0.0612 0.0723 0.0850 0.0993 0.1156 0.1390 0.1545 0.2033 0.2637 0.3377 0.4275 0.5359 0.9082 1.479 2.268 3.736 6.152 13.44];

% Specific Heat (liquid-Btu/lbm.R)

Cpimp=[1.013 1.015 1.018 1.020 1.023 1.026 1.029 1.036 1.044 1.054 1.065 1.078 1.121 1.188 1.298 1.509 2.086 13.80];

%dynamic viscosity (lmb/ft.s)

mulimp=[1.625e-4 1.544e-4 1.472e-4 1.406e-4 1.344e-4 1.289e-4 1.236e-4 1.144e-4 1.063e-4 9.972e-5 9.361e-5 8.833e-5 7.722e-5 6.833e-5 6.083e-5 5.389e-5 4.639e-5 3.417e-5];

% Prandtl number

Prlimp=[1.50 1.43 1.37 1.31 1.25 1.21 1.16 1.09 1.02 0.973 0.932 0.893 0.842 0.830 0.864 0.979 1.30 6.68];

% liquid enthalpy (Btu/lbm)

hlimp=[208.49 218.63 228.79 238.98 249.20 259.45 269.73 290.40 311.24 332.28 353.53 375.04 430.18 487.89 549.39 616.92 696.08 822.76];

% vapour enthalpy (Btu/lbm)

hvimp=[1160.5 1164.0 1167.4 1170.7 1173.9 1177.0 1180.0 1185.5 1190.5 1194.8 1198.5 1201.4 1205.1 1202.3 1190.9 1166.6 1119.7 991.1];

% liquid thermal conductivity (Btu/hr.ft.R)

kl=[0.394 0.395 0.395 0.395 0.395 0.395 0.394 0.393 0.391 0.389 0.385 0.382 0.370 0.352 0.329 0.299 0.267 0.254];

% surface tension (requires conversion)

    %temperature
    
Tsurf=[190
195
200
205
210
215
220
225
230
235
240
245
250
255
260
265
270
275
280
285
290
295
300
305
310
315
320
325
330
335
340
345
350
355
360
365
370
374];

    %surface tension
sigma=[39.9500
38.8300
37.6900
36.5500
35.4100
34.2500
33.1000
31.9300
30.7700
29.6000
28.4200
27.2400
26.0600
24.8700
23.6700
22.4800
21.3000
20.1100
18.9400
17.7700
16.6100
15.4500
14.3000
13.1700
12.0400
10.9200
9.8100
8.7300
7.6600
6.6100
5.5900
4.6000
3.6500
2.7500
1.9000
1.1300
0.4500
0];

%% Calculated system properties

%temperature
Tsatimp=LinLook(Pimp,Temperatureimp,Psatimp);

%liquid density

rholsys=LinLook(Pimp,rholimp,Psatimp);

% Gas density

rhogsys=LinLook(Pimp, rhovimp, Psatimp);

% liquid enthalpy

hlsys=LinLook(Pimp, hlimp, Psatimp);

% vapour enthalpy

hvsys=LinLook(Pimp, hvimp,Psatimp);

% heat capacity

Cpsys=LinLook(Pimp, Cpimp,Psatimp);

% liquid dynamic viscosity

mulsys=LinLook(Pimp,mulimp,Psatimp)*3600;

% surface tension

Tsatmet=(Tsatimp-32)/1.8;

sigmasys=LinLook(Tsurf,sigma,Tsatmet);

% liquid thermal conductivity 

klsys=LinLook(Pimp,kl,Psatimp);

% Prandtl Number

Prlsys=LinLook(Pimp,Prlimp,Psatimp);

%% Td calculation
% mass flux

G=Mspimp/Achannelimp;


% friction factor

f=0.0055*(1+((20000*relrough)+(10^6/(G*Dhimp/mulsys)))^(1/3));


% shear stress

tauw=f*G^2/(8*rholsys*gc); 

% heat transfer coefficient

h=0.023*klsys/Dhimp*((G*Dhimp/mulsys)^0.8)*(Prlsys^0.4);



% heat flux

Qf=Qchannelimp/Ahtimp;

% Dimensionless heat flux 

Qlevy=Qf/(rholsys*Cpsys*sqrt(tauw*gc/rholsys));

% Ybplus calculation

Ybplus=0.015*sqrt(sigmasys*gc*Dhimp*rholsys)/mulsys;

% Td calculation

if Ybplus>=0&&Ybplus<5
    
    Td=(Qf/h)-(Qlevy*Prlsys*Ybplus);
    
else if Ybplus>=5&&Ybplus<30;
        
        Td=(Qf/h)-(5*Qlevy*(Prlsys+log(1+(Prlsys*((Ybplus/5)-1)))));
        
    else if Ybplus>=30
            
            Td=(Qf/h)-(5*Qlevy*(Prlsys+log(1+(5*Prlsys))+(0.5*log(Ybplus./30))));
            
        end
    end
end

% mass fraction calculation

x=(houtimp-hlsys)/(hvsys-hlsys);

% x prime calculation

xd=-Cpsys*Td/(hvsys-hlsys);

xprime=x-(xd*exp((x/xd)-1));