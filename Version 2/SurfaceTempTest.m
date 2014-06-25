%% Temperature calculation

% channel as one uniform property channel

Qchannel=5.0; %MW

Tin=267.2; %C

M=25.8; % kg/s

Aflow=0.0034102; %m^2

G=M/Aflow; %kg/m^2.s

Lchannel=6; %m

doutclad=0.0138; %m

roc=doutclad/2; %m

dfuel=0.0122; %m

tclad=0.00038; %m

Q=Qchannel*1000/37; % kW

Dptin=0.1034; % pressure tube inner diameter

Dh=7.4885e-3; %m

%%
% Liquid coolant heat capacity temperatures (D20)

TCp=[160;170;180;190;200;220;240;260;280;300;320;340;360];

% liquid Heat capacity

Cp=[4.340;4.370;4.410;4.460;4.500;4.610;4.760;4.970;5.280;5.750;6.540;8.240;14.96];

% liquid thermal conductivity

kl=1e-3.*[0.680;0.677;0.673;0.669;0.663;0.650;0.632;0.609;0.581;0.548;0.509;0.469;0.427];

% UO2 thermal conductivity temperatures

kUO2T=400:100:2800;

% UO2 thermal conductivity (zero burnup)

kUO20B=1e-3.*[4.74 4.28 3.89 3.55 3.26 3.01 2.79 2.61 2.45 2.32 2.22 2.14 2.09 2.06 2.06 2.08 2.12 2.18 2.26 2.35 2.45 2.56 2.68 2.80 2.93];

% coolant dynamic viscosity

mul=[0.170e-3;0.160e-3;0.150e-3;0.142e-3;0.134e-3;0.122e-3;0.111e-3;0.102e-3;0.094e-3;0.086e-3;0.078e-3;0.070e-3;0.060e-3];

% liquid Prandtl Number

Prl=[1.09;1.03;0.983;0.947;0.910;0.865;0.836;0.832;0.854;0.902;1.00;1.23;2.06];
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




%% Bulk temperature calculation

Cp=LinLook(TCp,Cp,Tin);

Tbulk=(Qchannel*1000/M/Cp)+Tin;

%% Cladding surface temperature

kcool=LinLook(TCp,kl,Tbulk);

Prlcool=LinLook(TCp,Prl,Tbulk);

mulcool=LinLook(TCp,mul,Tbulk);

rholcool=1/LinLook(Tsat, svfluid, Tbulk);

Recool=M*rholcool*Dh/mulcool;

Nucool=0.023*(Recool^0.8)*(Prlcool^0.3);

hcool=Nucool*kcool/Dh; %kW/m^2.K

Tclado=(Q/(hcool*doutclad*Lchannel*pi()))+Tbulk;

%% Inner cladding temperature
kzirc=(7.51+(0.362e-3*Tclado)-(0.618e-7*Tclado^2)+(0.718e-11*Tclado^3))*10^-3;

ric=roc-tclad;

Tcladi=(Q*log(roc/ric)/(2*pi()*Lchannel*kzirc))+Tclado;

%% Fuel pin outer temperature

kgap=0.0476+(0.362e-3*Tcladi)-(0.618e-7*Tcladi^2)+(0.718e-11*Tcladi^3)*10^-3; %kW/m.C

dman=ric-(dfuel/2); %m

djump=10e-6; %m

hgap=kgap/(dman+djump);

Tfuelo=(Q/(dfuel/2*Lchannel)/hgap)+Tcladi;

%% Fuel pin centerline temperature

diff=10;

Tc=Tfuelo;

while diff>1
    
    Tave=(Tfuelo+Tc)/2;
    
    if Tave<400;
        kUO2=kUO20B(1);
        
    else
        kUO2=LinLook(kUO2T,kUO20B,Tave);
        
    end
   
    
    Tcnew=(Q*dfuel/2/4/kUO2)+Tfuelo;
    
    diff=Tcnew-Tc;
    
    Tc=Tcnew;
end

%% Temperature matrix creation

Tsurf=[Tbulk; Tclado; Tcladi; Tfuelo; Tc];


display(Tsurf);


