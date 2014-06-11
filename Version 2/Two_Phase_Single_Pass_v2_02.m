%%Inputs and property calculation 

Qchannel=5.52; %MW

Tin=267.2; %C

Dh=0.0074; %m^2

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

keff=DP/(M^2);

reff=keff*Rho;

Pout=10000; %kPa

Pin=Pout-sum(DP);


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

% Surface tension

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

hin=Lagint(Pcomp,hcomp,Pin);

% saturation temperature 

Tsys=Lagint(Psat,Tsat,Pout);

% saturation liquid enthalpy

hfsys=Lagint(Psat,hf,Pout);

% saturation vapour enthalpy

hvsys=Lagint(Psat,hv,Pout);

% saturation latent heat

hfvsys=Lagint(Psat,hfv,Pout);

% saturation liquid density

rhofsys=Lagint(Psat,rhof,Pout);

% saturation gas density

rhovsys=Lagint(Psat,rhov,Pout);

% Cp Liquid

Cpl=Lagint(TCp,Cp,Tsys);

% Thermal Conductivity

ksys=Lagint(TCp,kl,Tsys);

% Dynamic Viscosity

mulsys=Lagint(TCp,mul,Tsys);

% Prandtl Number

Prlsys=Lagint(TCp,Prl,Tsys);

% hout calculation

hout=hin+(Qchannel./M);

%% Two phase property calculation (use Levy from previous version to get averaged properties. If single phase )


%% Conversions  MAKE THIS A PART OF A SEPARATE LEVY ALPHA CALCULATION FUNCTION.
%liquid density

rhofsys=rhofsys/0.45329/3.2808^3; %lb/ft^3

%gaseous density

rhogsys=rhogsys/0.45329/3.2808^3; %lb/ft^3

%surface tension

sigma=sigma/0.22481*3.2808;

% Hydraulic diameter

Dh=Dh/3.2808; %ft

% dynamic viscosity

mulsys=mulsys/0.45329/3.2808;

% thermal conductivity

ksys=ksys*3.4121/3.2808/1.8;

%gravitational acceleration (ft/hr^2)

accimp=9.80665*0.22481/(3600^2);

%Mass flux (lb/h-ft^2)

G=G*0.22481*32.174/5.32/3600/(3.2808^2);

%lbf to lbm conversion

gc=32.174/5.32;




    
    


%% Single Phase Flow (re-do to calculate flowrate with given temperature of coolant DP will be constant. Use averaged properties if two phase)

Mchannel=sqrt(sum(DP)*rhosys/reff);





DPm=keff.*Minput.^2;

DPd=DPm*Rho./rhofsys;

DPt=sum(DPd);

display(DPd, 'Pressure drop per section');

display(DPt, 'Total single phase Pressure drop:');





