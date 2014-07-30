%% Single channel function

function Res=single_channel_Lockhart_Martinelli_Qloss_LPXSTEAM2(Qin, Tenter, Tmod, Pout)




Qchannel=Qin; %MW

Tin=Tenter; %C

Achannel=0.0034;%m^2

doutclad=0.0138; %m

roc=doutclad/2;

dfuel=0.0122; %m

tclad=0.00038; %m

Lchannel=5.94; %m

DPT=0.10338;

tPT=0.00424;

DCT=0.12869;

tCT=0.0014;

ript=DPT/2;

ropt=ript+tPT;

rict=DCT/2;

roct=rict+tCT;

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



Pin=((Pout*1000)+sum(DP))/1000;  %MPa


%% Prandtl Number data import

% Steam heat capacity Pressures (MPa)

PCp=[617.8;791.7;1002.1;1254.4;1553.8;2318;3344;4688;6412;8581;11274;14586;18651]./1000;

% liquid Prandtl Number

Prl=[1.09;1.03;0.983;0.947;0.910;0.865;0.836;0.832;0.854;0.902;1.00;1.23;2.06];





%% System Properties
% Compressed water enthalpy 

Tinsat=XSteam('Tsat_p',(Pin*10));

Cpin=XSteam('CpL_p',(Pin*10));

hsatin=XSteam('hL_p',(Pin*10)); %kJ/kg

hin=hsatin-(Cpin*(Tinsat-Tin));

% saturation temperature 

Tsys=XSteam('Tsat_p',(Pout*10)); %C

% saturation liquid enthalpy

hfsys=XSteam('hL_p',(Pout*10)); %kJ/kg

% saturation vapour enthalpy

hvsys=XSteam('hV_p',(Pout*10)); %kJ/kg

% saturation latent heat

hfvsys=hvsys-hfsys; %kJ/kg

% saturation liquid density

rhofsys=XSteam('rhoL_p',(Pout*10));  %kg/m^3

% saturation gas density

rhovsys=XSteam('rhoV_p',(Pout*10)); %kg/m^3

% Cp Liquid  

Cplsys=XSteam('CpL_p',(Pout*10)); %kJ/kg.K

% liquid Thermal Conductivity

klsys=XSteam('tcL_p', (Pout*10)); %W/m.K

% vapour Termal Conductivity 

%kvsys=interp1(PCp,kv,Pout); %W/m.K

% Liquid Dynamic Viscosity

mulsys=XSteam('my_pT',(Pout*10),(Tsys-1)); %kg/m.s

% Vapour Dynamic Viscosity

muvsys=XSteam('my_pT',(Pout*10),(Tsys+1)); %kg/m.s

% Liquid Prandtl Number

Prlsys=interp1(PCp,Prl,(Pout)); 


% Hydraulic Diameter

Dh=0.0074; %m

% CO2 thermal conductivity W/m.K

kCO2=[15.20 16.55 18.05 19.70 21.2 22.75 24.3 28.3 32.5 36.6 40.7 44.5 48.1 51.7 55.1].*0.001;

% CO2 thermal conductivity Temperatures

kCO2Temp=[208 300 320 340 360 380 400 450 500 550 600 650 700 750 800]-273.15;

%% Mass flow calculation

MchannelLO=sqrt((Pin-Pout)*rhofsys/reffT);

% hout calculation

hout=hin+(Qchannel*1000/MchannelLO); %kJ/kg

x=(hout-hfsys)/(hvsys-hfsys);

rhosys=rhofsys;
if x<=0
    Mchannel=MchannelLO;
end
if x>0
    
    PI2=(mulsys/muvsys)^0.2*rhovsys/rhofsys;
    
    xtt=((1-x)/x)^0.9*PI2^0.5;
    
    if xtt<=10
        
        al=(1+xtt^0.8)^-0.378;
    else 
        
        al=0.823-(0.157*log(xtt));
        
    end
    if al<=0
        Mchannel=MchannelLO;
    else
        LF=1/(1-al)^1.8;
    
        Mchannel=sqrt((Pin-Pout)*rhofsys/reffT/LF);
    
        rhosys=(al)*rhovsys+((1-al)*rhofsys);
    end
end

%% Determination of system reynolds number

% dynamic viscosity

musys=((x/muvsys)+((1-x)/mulsys)).^(-1);
    
 


% Reynold's number calculation

Reynolds=Mchannel.*rhosys.*Dh./musys;


%liquid only Reynold's number

if x<=0
    ReynoldsLO=Reynolds;
else
    ReynoldsLO=Mchannel*(1-x)*rhofsys*Dh/mulsys;
end
%% bulk temperature 
% This is actually the highest fluid temperature reached in the channel,
% but as there is negligible axial heat transfer this is where the clad and
% centerline temperatures will be calculated so that the highest possible
% values are found. Tbulk is typically the average temperature but this
% will not give the highest possible value. heat transfer properties will
% be calculated at the maximum possible values for the channel (node)
hl=0.023*ReynoldsLO.^0.8*Prlsys.^0.4*klsys/Dh;
    
if x<=0
        
    hsys=hl;
        
else
        
    Co=((1-x)/x)^0.8*(rhovsys/rhofsys)^0.5;
        
    Fr=(Mchannel/Achannel)^2/(rhofsys^2*9.81*Dh);
        
    if Fr>0.04
            
            C=0;
    else
            C=0.3;
    end
       
    Bo=Qchannel/(pi()*37*doutclad*Lchannel)/(Mchannel/Achannel*1000*hfvsys);
        
    htpconv=hl*((1.1360*Co.^-0.9*(25*Fr).^C)+(667.2*Bo.^0.7));
        
    htpnuc=hl*((0.6683*Co.^-0.2*(25*Fr).^C)+(1058.0*Bo.^0.7));
        
    if htpnuc>htpconv
            
        hsys=htpnuc;
            
    else
            
        hsys=htpconv;
    end
end

if x<=0
        
    Tbulk=(Qchannel*1000/Mchannel/Cplsys)+Tin;
        
else Tbulk=Tsys;
end
    
if Tbulk>Tsys
    Tbulk=Tsys;
end

Reouter=Mchannel*rhosys*DPT/musys;

Nuouter=0.023*Reouter^0.8*Prlsys^0.3;

houter=Nuouter*klsys/Lchannel;

Qloss=Qloss_single_channel(Tbulk,Tmod,houter);



%% Volumetric heat generation in fuel

Qvol=Qchannel./Lchannel/pi()/dfuel^2/4*1000;






%% Outer Clad Temperature



Q=Qchannel*1000/37;

    
Tclado=(Q/(hsys*doutclad*Lchannel*pi()))+Tbulk;
        


%% Inner Clad Temperature




ric=roc-tclad;



kzirc=(7.51+(0.362e-3*Tclado)-(0.618e-7*Tclado^2)+(0.718e-11*Tclado^3))*10^-3;

Tcladi=(Q*log(roc/ric)/(2*pi()*Lchannel*kzirc))+Tclado;
    


%% outer fuel meat temperature



dman=ric-(dfuel/2); %m

djump=10e-6; %m



    
kgap=0.0476+(0.362e-3*Tcladi)-(0.618e-7*Tcladi^2)+(0.718e-11*Tcladi^3)*10^-3; %kW/m.C

hgap=kgap/(dman+djump);
    
Tfuelo=(Q/(dfuel/2*Lchannel)/hgap)+Tcladi;


%% Fuel centerline Temperature


divf=1000;



rfuel=dfuel/2;



reval=linspace(rfuel,0,divf);

Tfuel=zeros(1,divf);

Tfuel(1)=Tfuelo;

for pev=2:divf
    
    Tfuel(pev)=(Qvol/4/kUO2(Tfuel(pev-1))*(reval(pev-1)^2-reval(pev)^2))+Tfuel(pev-1);
    
end

Tc=Tfuel(divf);


if x<=0
    al=0;
end

%% Pressure Tube inner Temperature

Reouter=Mchannel*rhosys*DPT/musys;

Nuouter=0.023*Reouter^0.8*Prlsys^0.3;

houter=Nuouter*klsys/Lchannel;

R1=1/(houter*Lchannel*DPT*pi()*Lchannel);

TIPT=Tbulk-(Qloss*1000000*R1);

kzircPT=(12.767-(5.4348e-4*(TIPT+273.15))+(8.9818e-6*(TIPT+273.15)^2));

R2=log(ropt/ript)/(2*pi()*kzircPT*Lchannel);

TOPT=TIPT-(Qloss*1000000*R2);

TCO2eval=(TOPT+Tmod)/2;

kCO2gap=interp1(kCO2Temp,kCO2,TCO2eval);

R3=log(rict/ropt)/(2*pi()*kCO2gap*Lchannel);

TICT=TOPT-(Qloss*1000000*R3);

kzircCT=(12.767-(5.4348e-4*(Tmod+273.15))+(8.9818e-6*(Tmod+273.15)^2));

R4=log(roct/rict)/(2*pi()*kzircCT*Lchannel);

TOCT=TICT-(Qloss*1000000*R4);

Res=[Qin;Mchannel;Tclado;Tcladi;Tfuelo;Tc;TOPT;TOCT;x;rhosys;al;Qloss];



