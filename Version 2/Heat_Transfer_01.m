%% Tsurfaces
%Heat transfer model
% First section determines temperature of outer clad surfaces under liquid
% or vapor coolant flow

% calculation of heat transfer coefficient for liquid and gas
% temperature of surfaces

function Tsurf=Heat_Transfer_01(Tin,Q,M,Lchannel)

% Physical data
doutclad=0.0130; %m

dfuel=0.0122; %m

tclad=0.00038; %m

Q=Q*1000; % kW
%%
% Steam heat capacity temperatures

TCp=[160;170;180;190;200;220;240;260;280;300;320;340;360];

% liquid Heat capacity

Cp=[4.340;4.370;4.410;4.460;4.500;4.610;4.760;4.970;5.280;5.750;6.540;8.240;14.96];

% liquid thermal conductivity

kl=[0.680;0.677;0.673;0.669;0.663;0.650;0.632;0.609;0.581;0.548;0.509;0.469;0.427];

% Zirconium emissivity temperature range

%TezT=[100 150 200 300 400];

% Zirconium emissivity

%Tez=[0.424 0.414 0.416 0.434 0.433];

% UO2 thermal conductivity temperatures

kUO2T=400:100:2800;

% UO2 thermal conductivity (zero burnup)

kUO20B=[4.74 4.28 3.89 3.55 3.26 3.01 2.79 2.61 2.45 2.32 2.22 2.14 2.09 2.06 2.06 2.08 2.12 2.18 2.26 2.35 2.45 2.56 2.68 2.80 2.93];

%% Bulk temperature calculation

Cp=LinLook(TCp,Cp,Tin);

Tbulk=(Q/M/Cp)+Tin;

%% Cladding surface temperature

hcool=LinLook(TCp,kl,Tbulk);

Tclado=(Q/(hcool*doutclad*Lchannel))+Tbulk;

%% Inner cladding temperature

kzirc=(8.8527+(7.0820e-3*Tclado)+(2.5329e-6*Tclado^2)+(2.9918e3/Tclado))/1000; %kW/m.C

roc=doutclad/2;

ric=roc-tclad;

Tcladi=(Q*log(roc/ric)/(2*pi()*Lchannel*kzirc))+Tclado;

%% Fuel pin outer temperature

kgap=0.002e-3; %kW/m.C

Tfuelo=(Q*log(ric/dfuel/2)/(2*pi()*Lchannel*kgap))+Tcladi;

%% Fuel pin centerline temperature

diff=10;

Guess=kUO20B(1);

Tc=Tfuelo;

while diff>1
    
    Tave=(Tfuelo+Tc)/2;
    
    if Tave<400;
        kUO2=Guess;
        
    else
        kUO2=LinLook(kUO2T,kOU20B,Tave);
        
    end
    Tprev=Tc;
    
    Tc=(Q*(dfuel/2)^2/4/kUO2)+Tfuelo;
    
    diff=Tc-Tprev;
end

%% Temperature matrix creation

Tsurf=[Tbulk; Tclado; Tcladi; Tfuelo; Tc];






