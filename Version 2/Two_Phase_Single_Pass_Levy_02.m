%%Inputs and property calculation 

Qchannel=5:0.1:12; %MW

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

Pout=10; %MPa

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

% Compressed water Pressures

Pcomp=[5;10]; % MPa

% Compressed water enthalpy (10 MPa, 260 C)

hcomp=[1.1349e+03;1.1343e+03]; %kJ/kg

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

hin=interp1(Pcomp,hcomp,Pin); %kJ/kg

% saturation temperature 

Tsys=interp1(Psat,Tsat,Pout); %C

% saturation liquid enthalpy

hfsys=interp1(Psat,hf,Pout); %kJ/kg

% saturation vapour enthalpy

hvsys=interp1(Psat,hv,Pout); %kJ/kg

% saturation latent heat

hfvsys=interp1(Psat,hfv,Pout); %kJ/kg

% saturation liquid density

rhofsys=interp1(Psat,rhof,Pout);  %kg/m^3

% saturation gas density

rhovsys=interp1(Psat,rhov,Pout); %kg/m^3

% Cp Liquid

Cpl=interp1(TCp,Cp,Tsys); %kJ/kg.K

% Thermal Conductivity

ksys=interp1(TCp,kl,Tsys); %W/m.K

% Dynamic Viscosity

mulsys=interp1(TCp,mul,Tsys); %kg/m.s

% Prandtl Number

Prlsys=interp1(TCp,Prl,Tsys); 

% hout calculation

hout=hin+(Qchannel.*1000./M); %kJ/kg

% Hydraulic Diameter

Dh=0.0074; %m

% Heat transfer area

Aht=9.1224; %m^2

%% Mass flow calculation

x=(hout-hfsys)./(hvsys-hfsys);

n=length(Qchannel);

Mchannel=zeros(1,n);

alphas=zeros(1,n);

for in=1:n
    if x(in)<=0
    
    Mchannel(in)=sqrt((Pout-Pin)*rhofsys/reffT);
    else
    check=1;
    eta=1;
    alpha=x(in);
    
    while eta>=0.001
        
        xLevy=((alpha*(1-(2*alpha)))+(alpha*sqrt((1-(2*(x(in))))+(alpha*((2*rhovsys/rhofsys*((1-alpha)^2))))+(alpha*(1-(2*alpha))))))/((2*rhovsys/rhofsys*(1-alpha))+(alpha*(1+(2*alpha))));
    
        eta=x(in)-xLevy;
        
        alpha=alpha-eta;
        
    end
    clear eta
    LF=((1-x(in))^1.75)/((1-alpha)^2);
    
    Mchannel(in)=sqrt((Pout-Pin)*rhofsys/reffT/LF);
    
    alphas(in)=alpha;
    
    clear LF
    
    clear alpha
    end
    
    
end

plot(Qchannel,Mchannel);


