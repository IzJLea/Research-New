%% Initial channel properties 
% included in this section are the physical specifications of the channel
% at the beginning of the simulation and any calculations needed to
% determine essential properties such as mass of elements, flow resistances
% etc...


Qchannel=.150; %MW channel power

H=0; %m height difference from inlet/outlet headers

hmod=1000; %W/m^2.K heat transfer coefficient 

Tenter=100; %C Entering coolant temperature

PSH=210; %kPa supply header pressure

PVH=190; %kPa vent header pressure

Peval=(PSH+PVH)/2/100; % system evaluation pressure

Tmod=60; %C Moderator temperature

Aflow=0.0035; %m^2 Flow area in channel

Dh=0.0074; %m hydraulic diameter of bundle

DPT=0.10338; %m inner pressure tube diameter

DCT=0.12869;%m inner calandria tube thickness

doutclad=0.0138; %m cladding outer diameter

dfuel=0.0122; %m fuel outer diameter

tclad=0.00038; %m thickness of cladding

roc=doutclad/2; %m outer cladding radius

ric=roc-tclad; %m inner cladding radius

Lchannel=5.94; %m length of fuel channel

Lbund=Lchannel; %m length of bundle (for this initial simulation, entire channel will be simulated as one bundle)

tPT=0.00424; %m pressure tube thickness

tCT=0.0014; %m calandria tube thickness

Aipt=pi()*DPT*Lbund;%m^2 inner pressure tube area

Aopt=pi()*(DPT+(2*tPT))*Lbund;%m^2 outer pressure tube area

Aict=pi()*DCT*Lbund;%m^2 inner calandria tube area

Aoct=pi()*(DCT+(2*tCT))*Lbund;%m^2 outer calandria tube area

Afuel=pi()*doutclad*Lbund*37;%m^2 fuel outer area

Aoutfuel=9*pi()*doutclad; % Area of outer fuel elements that sees PT

doxide=3.00e-6; % thickness of oxide layer (assumed to be low for initial simulation)

sigma=5.670373e-8; %Stefan-Boltzmann constant

eclad=0.325+(0.1246e6*doxide); % dimensionless emissivity value for zirconium cladding
        
ePT=eclad;% emissivity for pressure tube
        
eCT=ePT;% emissivity for calandria tube

ript=DPT/2; %inner pressure tube radius

ropt=ript+tPT; %outer pressure tube radius

rict=DCT/2; % inner calandria tube radius

roct=rict+tCT;% outer calandria tube radius

DhPT=4*Aflow/(pi()*DPT); % pressure tube hydraulic diameter

%% pressure drop calculation based on density difference as well as header pressure

deltaP=(PSH)-PVH+((9.81*H*(XSteam('rhoL_p',Peval)-XSteam('rhoV_p',Peval)))/1000);

%% Resistance Calculation from liquid flow data

DP = [165;
262.5000;
520.5000;
258;
174];  % measured reference pressure drops across the inlet feeder, end fitting, core, and exit end fitting and feeder respectively

M=25.8; % kg/s reference mass flowrate

Rho=780.6; %kg/m^3 coolant density at reference measurements

keff=DP./(M^2); 

reff=keff.*Rho;

RCH=reff(3); %effective resistance of fuel channel

RF=sum(reff(4:5)); %effective resistance of end fitting and feeder




%% mass of zirc in elements
Tref=310; % C reference temperature

TrefK=Tref+273.15;

rhoref=255.66+(0.1024*TrefK); %kg/m^3 reference density

mclad=37*((doutclad)^2-(doutclad-(2*tclad))^2)/4*Lbund*rhoref; %mass of zirconium in cladding

mPT=((DPT+(2*tPT))^2-(DPT)^2)/4*Lbund*rhoref; %kg mass of zirconium in pressure tube

mCT=((DCT+(2*tCT))^2-(DPT)^2)/4*Lbund*rhoref; %kg mass of zirconium in calandria tube

%% mass of fuel

rhofuel=10970/(1+(2.04e-5*Tref)+(8.7e-9*Tref^2)); %kg/m^3 reference fuel density

mfuel=Lbund*pi()/4*dfuel^2*rhofuel*37;%kg mass of fuel within fuel channel

%% CO2 thermal conductivity
% as the correlations for CO2 properties are generally very large complex
% and precise mathematical correlations which can take quite a bit of time
% to implement and compute, simple interpolation of measured values will be
% used during this simulation for properties of CO2

% CO2 thermal conductivity for use in interpolation

    kCO2=[0.01051 0.01456 0.01858 0.02257 0.02652 0.03044 0.03814 0.04565 0.05293 0.08491 0.10688 0.11522];

% CO2 thermal conductivity Temperatures for use in interpolation

    kCO2Temp=[-50 0 50 100 150 200 300 400 500 1000 1500 2000];
    
%% Time determination and property matrix creation
% This section sets the length of the simulation and time divisions and
% creates property matrices for the major variables used in the simulation


Time=500;  %seconds total run time

div=0.1; %s time step length

ind=Time/div; % index to determine time step number

Tclad=zeros(1,ind); %C cladding temperature

dTclad=zeros(1,ind);%C change in cladding temperature

Tvap=zeros(1,ind);%C vapor coolant temperature

dTvap=zeros(1,ind);%C change in vapor coolant temperature

TPT=zeros(1,ind);%C pressure tube temperature

dTPT=zeros(1,ind);%C change in pressure tube temperature

Tfuel=zeros(1,ind);%C average fuel temperature

dTfuel=zeros(1,ind);%C change in fuel temperature

TCT=zeros(1,ind);%C calandria tube temperature

dTCT=zeros(1,ind);%C change in calandria tube temperature

Reynolds=zeros(1,ind);% reynolds number of coolant

Prandtl=zeros(1,ind);% coolant prandtl number

Nusselt=zeros(1,ind);% coolant nusselt number

hcool=zeros(1,ind);% coolant convective heat transfer coefficient 

hrad=zeros(1,ind);%W/m^2.K radiative heat transfer coefficient for cladding/pressure tube

kuo2=zeros(1,ind);% W/m.K thermal conductivity of UO2 fuel 

kclad=zeros(1,ind);% W/m.K thermal conductivity of cladding

kPT=zeros(1,ind);% W/m.K pressure tube thermal conductivity

kCT=zeros(1,ind);%W/m.K calandria tube thermal conductivity

kCO2sys=zeros(1,ind);%W/m.K CO2 insulator thermal conductivity

Cpclad=zeros(1,ind);%J/kg.K cladding heat capacity

CpPT=zeros(1,ind);%J/kg.K pressure tube heat capacity

CpCT=zeros(1,ind);%J/kg.K calandria tube heat capacity

Cpvap=zeros(1,ind);%J/kg.K vapor coolant heat capacity

Cpfuel=zeros(1,ind);%J/kg.K fuel heat capacity

hout=zeros(1,ind);% channel exit enthalpy

A1=zeros(1,ind);
    
B1=zeros(1,ind);
    
F1=zeros(1,ind);
    
A2=zeros(1,ind);
    
B2=zeros(1,ind);
    
C2=zeros(1,ind);
    
D2=zeros(1,ind);
       
B3=zeros(1,ind);
    
C3=zeros(1,ind);
    
D3=zeros(1,ind);

F3=zeros(1,ind);

B4=zeros(1,ind);

C4=zeros(1,ind);

D4=zeros(1,ind);

E4=zeros(1,ind);

D5=zeros(1,ind);

E5=zeros(1,ind);

F5=zeros(1,ind);

Rgap=zeros(1,ind); % fuel-cladding gap resistance

R1=zeros(1,ind); 

R2=zeros(1,ind);

R3=zeros(1,ind);

R4=zeros(1,ind);

R5=zeros(1,ind);

CH=zeros(5,ind);

Mcool=zeros(1,ind); %kg mass of coolant within the channel

%% vapor fraction calculation
%This section creates a loop to solve for the vapor fraction (alpha) in the
%channel under low flow conditions

alpha=0;

x=(XSteam('h_pT',Peval,Tenter)-XSteam('hL_p',Peval))/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval));

wsmx=Qchannel*1000/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval));

d=0.01; %dampening factor

delta=0.1;

while delta>=0.0001
    
    wv=(1-alpha)*wsmx/(1+x);
    
    
    hg=XSteam('hV_P',Peval);

    hv=hg+(alpha*Qchannel*1000/wv);
    
    rhov=XSteam('rho_ph',Peval,hv);
    
    rhog=XSteam('rhoV_p',Peval);
    
    rhoave=(rhov+rhog)/2;
    
    a=alpha^2*rhoave/rhog;
    
    alphanew=(1+((1+x)/wsmx*sqrt(deltaP*rhoave/(RCH+(a*RF)))))^-1;
    
    err=alphanew-alpha;
    
    delta=abs(err);
    
    alpha=alpha+(d*err);
    
end
 


%% Initial property calculations
%in this section the initial temperature for the elements and the coolant
%are determined. The temperatures were determined using the assumption that
%saturated fluid was previously flowing in the channel. also calculated is
%the mass flow through the channel

Qbundle=Qchannel; %due to treatment of channel as single bundle to be changed with multiple bundle model

Qel=Qchannel/37; % generation in one element

Qvol=Qel/(pi()/4*dfuel^2*Lchannel); %Volumetric heat generation per pin

mflow=Qbundle*1000*(1-alpha)/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval)); % kg/s 

Tvap(1,1)=XSteam('Tsat_p',Peval); %C initial vapor temperature (saturated)

Reynolds(1,1)=mflow*Dh/Aflow/XSteam('my_pT',Peval,Tvap(1,1)-0.1);% initial reynolds number based on liquid in channel

Prandtl(1,1)=XSteam('Cp_pT',Peval,Tvap(1,1)-0.1)*1000*XSteam('my_pT',Peval,Tvap(1,1)-0.1)/XSteam('tc_pT', Peval,Tvap(1,1)-0.1); % initial prandtl number 

if Reynolds(1,1)<=3000 % initial nusselt number includes consideration for laminar flow
        Nusselt(1,1)=4.36;
else
        
        Nusselt(1,1)=0.023*Reynolds(1,1)^(4/5)*Prandtl(1,1)^0.4;
end

hcool(1,1)=Nusselt(1,1)*XSteam('tcL_p',Peval)/Dh; %J/m^2.K coolant heat transfer coefficient

Tclad(1,1)=(Qbundle*1000000/Afuel/hcool(1,1))+Tvap(1,1); %C initial cladding temperature 

kuo2(1,1)=kUO2(Tclad(1,1))*1000;%W/m.K initial fuel heat capacity

kclad(1,1)=12.767-(5.4348e-4*(Tclad(1,1)+273.15))+(8.9818e-6*(Tclad(1,1)+273.15)^2); %W/m.K initial cladding heat capacity

Rfuel(1,1)=1/(4*pi()*kuo2(1,1)*Lchannel); % resistance of entire fuel meat

Rclad(1,1)=log(roc/ric)/(2*pi()*kclad(1,1)*Lchannel); %resistance of cladding

Rgap(1,1)=0;% gap resistance to be ignored for now. 

Rconv(1,1)=1/(Afuel/hcool(1,1)); % covective resistance of one element

R1(1,1)=(Rfuel(1,1)/2)+Rgap(1,1)+(Rclad(1,1)/2);% resistance between Tfuel and Tclad

TPT(1,1)=Tvap(1,1); % simplified initial pressure tube temperature

TCT(1,1)=Tmod; % simplified initial calandria tube temperature

Tfuel(1,1)=Tclad(1,1)+(Qchannel*1000000/37*R1(1,1)); % initial average fuel temperature

hin=XSteam('hV_p',Peval)*1000; % J/kg enthalpy of steam entering voided section

%% transient calculation
%These calculations determine the five temperatures during the duration of
%the simulation. 

for n=2:ind
    % calculation of vapor heat transfer coefficient. The if loops are to
    % get values if conditions are at saturation conditions 
    
    if Tvap(1,n-1)==XSteam('Tsat_p',Peval)
        
       Reynolds(1,n)=mflow*Dh/Aflow/alpha/XSteam('my_pT',Peval,Tvap(1,1)+0.1);
       
       Prandtl(1,n)=XSteam('Cp_pT',Peval,Tvap(1,1)+0.1)*1000*XSteam('my_pT',Peval,Tvap(1,1)+0.1)/XSteam('tc_pT', Peval,Tvap(1,1)+0.1);
    else
        Reynolds(1,n)=mflow*Dh/Aflow/alpha/XSteam('my_pT',Peval,Tvap(1,n-1));
        
        Prandtl(1,n)=XSteam('Cp_pT',Peval,Tvap(1,n-1))*1000*XSteam('my_pT',Peval,Tvap(1,n-1))/XSteam('tc_pT', Peval,Tvap(1,n-1));
    end
    
    if Reynolds(1,n)<=3000
        
        Nusselt(1,n)=4.36;
        
    else
        
        Nusselt(1,n)=0.023*(Reynolds(1,n)^(4/5))*(Prandtl(1,n)^0.4);
    end
    
    if Tvap(1,n-1)>XSteam('Tsat_p',Peval)
    
        hcool(1,n)=Nusselt(1,n)*XSteam('tc_pT',Peval,Tvap(1,n-1))/Dh;
    else
        hcool(1,n)=Nusselt(1,n)*XSteam('tc_pT',Peval,Tvap(1,n-1)+0.1)/Dh;
    end
    
    % heat transfer coefficient for radiation- view area of fuel pins is
    % assumed to be equal to the outer half of the 18 outer fuel pins
    hrad(1,n)=sigma*((Tclad(1,n-1)+273.15)+(TPT(1,n-1)+273.15))*(((Tclad(1,n-1)+273.15)^2)+((TPT(1,n-1)+273.15)^2))/((1/eclad)+((1-eclad)/eclad*Afuel*9/37/Aipt));
    
    % heat capacity calculations for main elements
    Cpclad(1,n)=255.66+(0.1024*(Tclad(1,n-1)+273.15)); %J/kg.K
    
    CpPT(1,n)=255.66+(0.1024*(TPT(1,n-1)+273.15)); %J/kg.K
    
    CpCT(1,n)=255.66+(0.1024*(TCT(1,n-1)+273.15));
    
    tau=(Tfuel(1,n-1)+273.15)/1000;
    
    Cpfuel(1,n)=(52.1743+(87.951*tau)-(85.2411*tau^2)+(31.542*tau^3)-(2.6334*tau^4)-(0.71391*tau^-2))/270.03*1000; %J/kg.K
    
    if Tvap(1,n-1)>XSteam('Tsat_p',Peval)
    
        Cpvap(1,n)=XSteam('Cp_pT',Peval,Tvap(1,n-1))*1000;
    else
        Cpvap(1,n)=XSteam('CP_pT',Peval,Tvap(1,n-1)+0.1)*1000;
    end
    % Exit enthalpy calculation based off of temperature of vapor in
    % previous time step
     if Tvap(1,n-1)>XSteam('Tsat_p',Peval)
    
        hout(1,n)=XSteam('h_pT',Peval,Tvap(1,n-1))*1000; % J/kg
    else
        hout(1,n)=XSteam('h_pT',Peval,Tvap(1,n-1)+0.1)*1000;
     end
    
    
     %thermal conductivity of main elements
     
       kuo2(1,n)=kUO2(Tfuel(1,n-1))*1000; % W/m.K
    
    kclad(1,n)=12.767-(5.4348e-4*(Tclad(1,n-1)+273.15))+(8.9818e-6*(Tclad(1,n-1)+273.15)^2); %W/m.K
    
    kPT(1,n)=12.767-(5.4348e-4*(TPT(1,n-1)+273.15))+(8.9818e-6*(TPT(1,n-1)+273.15)^2);
    
    kCT(1,n)=12.767-(5.4348e-4*(TCT(1,n-1)+273.15))+(8.9818e-6*(TCT(1,n-1)+273.15)^2);
    
    kCO2sys(1,n)=interp1(kCO2Temp,kCO2,(Tvap(1,n-1)+Tmod)/2);
    
    % mass of steam in voided section
    
    if Tvap(n-1)==XSteam('Tsat_p',Peval)
        
        Mcool(1,n)=XSteam('rhoV_p',Peval)*Aflow*Lchannel*alpha;
    else
        
        Mcool(1,n)=XSteam('rho_pT',Peval,Tvap(1,n-1))*Aflow*Lchannel*alpha;
    end
    
    %calculation of resistances. See accompanying document for derivations
    R1(1,n)=(1/8/pi()/kuo2(1,n)/Lchannel)+Rgap(1,n)+(log(roc/ric)/4/pi()/kclad(1,n)/Lchannel);

    R2(1,n)=(log(roc/ric)/4/pi()/kclad(1,n)/Lchannel)+(1/hcool(1,n)/Afuel);
    
    R3(1,n)=(1/hcool(1,n)/Aipt)+(log(ropt/ript)/4/pi()/kPT(1,n)/Lchannel);
    
    R4(1,n)=(log(ropt/ript)/4/pi()/kPT(1,n)/Lchannel)+(log(rict/ropt)/2/pi()/kCO2sys(1,n)/Lchannel)+(log(roct/rict)/4/pi()/kCT(1,n)/Lchannel);
    
    R5(1,n)=(log(roct/rict)/4/pi()/kCT(1,n)/Lchannel)+(1/hmod/Aoct);
    
    %Eq.1 coefficients
    
    A1(1,n)=-37/R1(1,n)/mfuel/Cpfuel(1,n);
    
    B1(1,n)=37*Tclad(1,n-1)/R1(1,n)/mfuel/Cpfuel(1,n);
    
    F1(1,n)=Qchannel*1000000/mfuel/Cpfuel(1,n);
    
    %Eq.2 coefficients
    
    A2(1,n)=37*Tfuel(1,n-1)/mclad/Cpclad(1,n)/R1(1,n);
    
    B2(1,n)=-1/mclad/Cpclad(1,n)*((37/R1(1,n))+(37/R2(1,n))+(hrad(1,n)*Afuel*9/37));
    
    C2(1,n)=37*Tvap(1,n-1)/mclad/Cpclad(1,n)/R2(1,n);
    
    D2(1,n)=hrad(1,n)*Afuel*9/37/mclad/Cpclad(1,n)*TPT(1,n-1);
    
    %Eq.3 coefficients
    
    B3(1,n)=37*Tclad(1,n-1)/Mcool(1,n)/Cpvap(1,n)/R2(1,n);
    
    C3(1,n)=-1/Mcool(1,n)/Cpvap(1,n)*((37/R2(1,n))+(1/R3(1,n)));
    
    D3(1,n)=TPT(1,n-1)/Mcool(1,n)/Cpvap(1,n)/R3(1,n);
    
    F3(1,n)=mflow/alpha/Mcool(1,n)*(hout(1,n)-hin)/Cpvap(1,n);
    
    %Eq.4 coefficients
    
    B4(1,n)=hrad(1,n)*Afuel*9/37*Tclad(1,n-1)/mPT/CpPT(1,n);
    
    C4(1,n)=Tvap(1,n-1)/mPT/CpPT(1,n)/R3(1,n);
    
    D4(1,n)=-1/mPT/CpPT(1,n)*((1/R3(1,n))+(1/R4(1,n))+(Afuel*9/37*hrad(1,n)));
    
    E4(1,n)=TCT(1,n-1)/mPT/CpPT(1,n)/R4(1,n);
    
    %Eq.5 coefficients
    
    D5(1,n)=TPT(1,n-1)/R4(1,n)/mCT/CpCT(1,n);
    
    E5(1,n)=-1/mCT/CpCT(1,n)*((1/R4(1,n))+(1/R5(1,n)));
    
    F5(1,n)=Tmod/mCT/CpCT(1,n)/R5(1,n);
    
    % calculation of Temperature values
    
    Tfuel(1,n)=(Tfuel(1,n-1)*exp(-A1(1,n)*div))+((1-exp(-A1(1,n)*div))*((B1(1,n)+F1(1,n))/-A1(1,n)));
    
    Tclad(1,n)=(Tclad(1,n-1)*exp(-A1(1,n)*div))+((1-exp(-A1(1,n)*div))*((B1(1,n)+F1(1,n))/-A1(1,n)));
    
    Tvap(1,n)=(Tvap(1,n-1)*exp(-A1(1,n)*div))+((1-exp(-A1(1,n)*div))*((B1(1,n)+F1(1,n))/-A1(1,n)));
    
    TPT(1,n)=(TPT(1,n-1)*exp(-A1(1,n)*div))+((1-exp(-A1(1,n)*div))*((B1(1,n)+F1(1,n))/-A1(1,n)));
    
    TCT(1,n)=(TCT(1,n-1)*exp(-A1(1,n)*div))+((1-exp(-A1(1,n)*div))*((B1(1,n)+F1(1,n))/-A1(1,n)));
    
    
end

    