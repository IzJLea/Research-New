%%Inputs and property calculation 

Qchannel=5.52; %MW

Qchannelimp=Qchannel*1000*3412.14; %Btu/h

Tin=267.2; %C

Achannel=0.0034;%m^2

Achannelimp=Achannel*3.2808^2; %ft^2

%% Initial single pass k calculation- 
% pressure differences in each pass section
%



DP = [165;
262.5000;
520.5000;
258;
174];


%%calculate values for k for each section

% reference mass flowrate

M=25.8; % kg/s

Rho=780.6; %kg/m^3

keff=DP./1000./(M^2);

reff=keff.*Rho;

reffT=sum(reff);

Pout=8.3; %MPa

Pin=((Pout*1000)-sum(DP))/1000;  %MPa


%% Physical data import
% Temperature for water properties (C)


Tsat=[100
110
120
130
140
150
160
170
180
190
200
210
220
230
240
250
260
270
280
290
300
310
320
330
340
350
360
370
373.9500];

Tsat=reshape(Tsat,29,1);

% Saturation Pressure (MPa)

Psat=[0.1014
0.1434
0.1987
0.2703
0.3615
0.4762
0.6182
0.7922
1.0028
1.2552
1.5549
1.9077
2.3196
2.7971
3.3469
3.9762
4.6923
5.5030
6.4166
7.4418
8.5879
9.8651
11.2840
12.8580
14.6010
16.5290
18.6660
21.0440
22.0640];

Psat=reshape(Psat,29,1);

% Fluid specific volume

svfluid=[0.0010
0.0011
0.0011
0.0011
0.0011
0.0011
0.0011
0.0011
0.0011
0.0011
0.0012
0.0012
0.0012
0.0012
0.0012
0.0013
0.0013
0.0013
0.0013
0.0014
0.0014
0.0014
0.0015
0.0016
0.0016
0.0017
0.0019
0.0022
0.0031];

svfluid=reshape(svfluid,29,1);

% density conversion (kg/m^3)

rhof=1./svfluid;

% vapour specific volume

svvapour=[1.6718
1.2093
0.8912
0.6680
0.5085
0.3925
0.3068
0.2426
0.1938
0.1564
0.1272
0.1043
0.0861
0.0715
0.0597
0.0501
0.0422
0.0356
0.0302
0.0256
0.0217
0.0183
0.0155
0.0130
0.0108
0.0088
0.0069
0.0050
0.0031];

svvapour=reshape(svvapour,29,1);

% density conversion (kg/m^3)

rhov=1./svvapour;

% fluid enthalpy (kJ/kg)

hf=[419.1700
461.4200
503.8100
546.3800
589.1600
632.1800
675.4700
719.0800
763.0500
807.4300
852.2700
897.6300
943.5800
990.1900
1.0376e+03
1.0858e+03
1135
1.1853e+03
1.2369e+03
1290
1345
1.4022e+03
1.4622e+03
1.5259e+03
1.5945e+03
1.6709e+03
1.7617e+03
1.8907e+03
2.0843e+03];

hf=reshape(hf,29,1);

% vapour enthalpy (kJ/kg)

hv=[2.6756e+03
2.6911e+03
2.7059e+03
2.7201e+03
2.7334e+03
2.7459e+03
2.7574e+03
2.7679e+03
2.7772e+03
2.7853e+03
2792
2.7973e+03
2.8009e+03
2.8029e+03
2803
2.8009e+03
2.7966e+03
2.7897e+03
2.7799e+03
2.7667e+03
2.7496e+03
2.7279e+03
2.7006e+03
2666
2.6218e+03
2.5636e+03
2.4815e+03
2.3345e+03
2.0843e+03];

hv=reshape(hv,29,1);

% latent heat (kJ/kg)

hfv=[2.2564e+03
2.2297e+03
2.2021e+03
2.1737e+03
2.1442e+03
2.1137e+03
2.0819e+03
2.0488e+03
2.0141e+03
1.9779e+03
1.9397e+03
1.8997e+03
1.8573e+03
1.8127e+03
1.7654e+03
1.7151e+03
1.6616e+03
1.6044e+03
1543
1.4767e+03
1.4046e+03
1.3257e+03
1.2384e+03
1.1401e+03
1.0273e+03
892.7000
719.8000
443.8000
0];

hfv=reshape(hfv,29,1);
% Compressed water Pressures

Pcomp=[5;10]; % MPa

% Compressed water enthalpy (10 MPa, 260 C)

hcomp=[1.1349e+03;1.1343e+03]; %kJ/kg

% Surface tension Temperatures

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

Tsurf=reshape(Tsurf,38,1);

% Surface tension mN/m

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

sigma=reshape(sigma,38,1);

% Steam heat capacity temperatures

TCp=[160;170;180;190;200;220;240;260;280;300;320;340;360];

% liquid Heat capacity

Cp=[4.340;4.370;4.410;4.460;4.500;4.610;4.760;4.970;5.280;5.750;6.540;8.240;14.96];

% liquid thermal conductivity

kl=[0.680;0.677;0.673;0.669;0.663;0.650;0.632;0.609;0.581;0.548;0.509;0.469;0.427];

% dynamic viscosity

mul=[0.170e-3;0.160e-3;0.150e-3;0.142e-3;0.134e-3;0.122e-3;0.111e-3;0.102e-3;0.094e-3;0.086e-3;0.078e-3;0.070e-3;0.060e-3];

% liquid Prandtl Number

Prl=[1.09;1.03;0.983;0.947;0.910;0.865;0.836;0.832;0.854;0.902;1.00;1.23;2.06];

%% System Properties
% Compressed water enthalpy 

hin=LinLook(Pcomp,hcomp,Pin); %kJ/kg

% Compressed water enthalpy imperial 

hinimp=hin*0.94782/2.2046; %Btu/lb

% saturation temperature 

Tsys=LinLook(Psat,Tsat,Pout); %C

% saturation temp imperial

TsysimpF=(Tsys*1.8)+32;  %Fahrenheit

TsysimpR=TsysimpF+459.67; %Rankine

% saturation liquid enthalpy

hfsys=LinLook(Psat,hf,Pout); %kJ/kg

% saturation liquid enthalpy imperial

hfsysimp=hfsys*0.94782/2.2046; %Btu/lb

% saturation vapour enthalpy

hvsys=LinLook(Psat,hv,Pout); %kJ/kg

% saturation vapour enthalpy imperial

hvsysimp=hvsys*0.94782/2.2046; %Btu/lb

% saturation latent heat

hfvsys=LinLook(Psat,hfv,Pout); %kJ/kg

% saturation latent heat imperial

hfvsysimp=hfvsys*0.94782/2.2046; %Btu/lb

% saturation liquid density

rhofsys=LinLook(Psat,rhof,Pout);  %kg/m^3

% saturation liquid density imperial

rhofsysimp=rhofsys*2.2046/3.2808^3; %lb/ft^3

% saturation gas density

rhovsys=LinLook(Psat,rhov,Pout); %kg/m^3

%saturation gas density imperial

rhovsysimp=rhovsys*2.2046/3.2808^3; %lb/ft^3

% Cp Liquid

Cpl=LinLook(TCp,Cp,Tsys); %kJ/kg.K

% Cp liquid imperial

Cplimp=Cpl*0.97482/2.2046*1.8; %Btu/lb.R

% Thermal Conductivity

ksys=LinLook(TCp,kl,Tsys); %W/m.K

% Thermal Conductivity imperial

ksysimp=ksys*3.41214/3.2808*1.8; %Btu/ft.R

% Dynamic Viscosity

mulsys=LinLook(TCp,mul,Tsys); %kg/m.s

% Dynamic Viscosity imperial

mulsysimp=mulsys*2.2046/3.2808*3600; %lb/ft.h

% Prandtl Number

Prlsys=LinLook(TCp,Prl,Tsys); 

% hout calculation

hout=hin+(Qchannel.*1000./M); %kJ/kg

% hout imperial

houtimp=hout*0.94782/2.2046; %Btu/lb

% Hydraulic Diameter

Dh=0.0074; %m

% Hydraulic Diameter imperial

Dhimp=Dh*3.2808; % ft

% Heat transfer area

Aht=9.1224; %m^2

% Heat transfer area imperial

Ahtimp=Aht*3.2808^2;

% saturation surface tension

sigmasys=LinLook(Tsurf,sigma,Tsys); %N/m

%saturation surface tension imperial

sigmasysimp=sigmasys*0.22481/3.2808; %lb/ft

% lb force to lbmass conversion

gc=1/32.174; %s^2/ft

% imperial gravitational constant

g=32.174*3600^2; %ft/h^2


%% Single Phase mass flow calculation

Mchannel=sqrt((Pout-Pin)*rhofsys/reffT);




%% Alpha calculation 

% delta T departure calculation

Ybplus=0.015*sqrt(sigmasysimp*gc*Dhimp*rhofsysimp)/mulsysimp;

G=Mchannel/Dhimp*2.2046*3600;

f=0.0055*(1+(((20000*5e-6/Dhimp)+(10^6/(G*Dhimp/mulsysimp)))^(1/3)));

tauw=f*G^2/(8*rhofsysimp*gc);

h=0.023*ksysimp/Dhimp*((G*Dhimp/mulsysimp)^0.8)*(Prlsys^0.4);

Qlevy=(Qchannelimp/Achannelimp)/(rhofsysimp*Cplimp*sqrt(tauw*gc/rhofsysimp));

if Ybplus>=0&&Ybplus<5
    Td=(Qchannelimp/Achannelimp/h)-(Qlevy*Prlsys*Ybplus);
else 
    if Ybplus>=5&&Ybplus<30
        Td=(Qchannelimp/Achannelimp/h)-(5*Qlevy*ln(1+(Prlsys*((Ybplus/5)-1))));
    else
        if Ybplus>=30
            Td=(Qchannelimp/Achannelimp/h)-(5*Qlevy*(1+ln(1+(5*Prlsys))+(0.5*ln(Ybplus/30))));
        end
    end
end
display(Ybplus,'Ybplus: ')

display(Td,'Td: ');


%% Single phase mass flow calculation

if hout>hfsys
    display('Two Phase Flow')
end

    
display(Mchannel,'Mass Flow Rate (kg/s): ');


    


