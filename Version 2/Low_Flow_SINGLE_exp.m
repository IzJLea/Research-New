%% Description of System Power (equal linear distribution of power to bundles)


Qchannel=.150; %MW

H=0; %m

hmod=1000; %W/m^2.K

Tenter=100; %C

PSH=210; %kPa

PVH=190; %kPa

Peval=(PSH+PVH)/2/100;

Tmod=60; %C

Aflow=0.0035; %m^2

Dh=0.0074; %m

DPT=0.10338;

DCT=0.12869;

doutclad=0.0138; %m

dfuel=0.0122; %m

tclad=0.00038; 

roc=doutclad/2;

ric=roc-tclad;

Lchannel=5.94; %m

Lbund=Lchannel;

tPT=0.00424;

tCT=0.0014;

Aipt=pi()*DPT*Lbund;

Aopt=pi()*(DPT+(2*tPT))*Lbund;

Aict=pi()*DCT*Lbund;

Aoct=pi()*(DCT+(2*tCT))*Lbund;

Afuel=pi()*doutclad*Lbund*37;

Aoutfuel=9*pi()*doutclad; % Area of outer fuel elements that sees PT

doxide=3.00e-6; % thickness of oxide layer

sigma=5.670373e-8; %Stefan-Boltzmann constant

ript=DPT/2;

ropt=ript+tPT;

rict=DCT/2;

roct=rict+tCT;

DhPT=4*Aflow/(pi()*DPT);
%% delta P calculation

deltaP=(PSH)-PVH+((9.81*H*(XSteam('rhoL_p',Peval)-XSteam('rhoV_p',Peval)))/1000);

%% Resistance Calculation from liquid flow data

DP = [165;
262.5000;
520.5000;
258;
174];





%% mass of zirc in elements
Tref=310; % Celsius

TrefK=Tref+273.15;

rhoref=255.66+(0.1024*TrefK);

mclad=37*((doutclad)^2-(doutclad-(2*tclad))^2)/4*Lbund*rhoref;

mPT=((DPT+(2*tPT))^2-(DPT)^2)/4*Lbund*rhoref;

mCT=((DCT+(2*tCT))^2-(DPT)^2)/4*Lbund*rhoref;

%% mass of fuel (single pin)

rhofuel=10970/(1+(2.04e-5*Tref)+(8.7e-9*Tref^2));

mfuel=Lbund*pi()/4*dfuel^2*rhofuel;

%% calculate values for k for each section

%reference mass flowrate

M=25.8; % kg/s

Rho=780.6; %kg/m^3

keff=DP./(M^2);

reff=keff.*Rho;

RCH=reff(3);

RF=sum(reff(4:5));

%% alpha calculation

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

Qbundle=Qchannel;

Qel=Qchannel/37;

Qvol=Qchannel/(pi()/4*dfuel^2*Lchannel);

%% Bundle Mass flow

mbundle=((1-alpha)*Qbundle*1000/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval)));
    
    
    


Lrun=10000; %length of run (s)

dt=1; %time division

divf=10000;
                    
time=zeros(1,Lrun/dt);

for tt=2:length(time)
    
    time(tt)=time(tt-1)+dt;
end
%% Property Matrix Creation
% Temperature
Tvap=zeros(1,length(time)); % Mean Vapour Temperature

Tclad=zeros(1,length(time)); % Mean clad Temperature

TPT=zeros(1,length(time)); % Pressure tube temperature

TCT=zeros(1,length(time)); % Calandria Tube Temperature

Tvap(1,1)=XSteam('Tsat_p',Peval); % Initial Vapour temperature assignment

Tfuelo=zeros(1,length(time)); % Outer fuel Temperature

Tfuelmid=zeros(1,length(time)); % fuel centerline Temperature

Tfuelmeat=zeros(1,length(time)); % mean fuel temperature

% Thermal Conductivity

kvap=zeros(1,length(time)); % vapour thermal conductivity

kvap(1,1)=XSteam('tcV_p',Peval); % initial vapour thermal conductivity assignment

kCO2sys=zeros(1,length(time));
% Heat Capacity

Cpvap=zeros(1,length(time)); % Vapour Heat capacity

Cpvap(1,1)=XSteam('CpV_p',Peval)*1000; % initial vapour heat capacity assignment


Cpclad=zeros(1,length(time)); % cladding Heat capacity

CpPT=zeros(1,length(time)); % Pressure Tube Heat capacity

CpCT=zeros(1,length(time)); % Pressure Tube Heat capacity

% Density

rhovap=zeros(1,length(time)); % vapour density

rhovap(1,1)=XSteam('rhoV_p',Peval); % vapour density assignment

% HT coefficient

hcool=zeros(1,length(time)); % heat transfer coefficient between fuel and vapour

hpt=zeros(1,length(time)); % heat transfer coefficient between PT and vapour

hrFPT=zeros(1,length(time)); % radiative heat transfer coefficient for inside channel

hCO2rad=zeros(1,length(time)); % radiative heat transfer coefficient for CO2 gap

A1=zeros(1,length(time));

A2=zeros(1,length(time));

A3=zeros(1,length(time));

A4=zeros(1,length(time));

B1=zeros(1,length(time));

B2=zeros(1,length(time));

B3=zeros(1,length(time));

B4=zeros(1,length(time));

C1=zeros(1,length(time));

C2=zeros(1,length(time));

C3=zeros(1,length(time));

C4=zeros(1,length(time));

D1=zeros(1,length(time));

D2=zeros(1,length(time));

D3=zeros(1,length(time));

D4=zeros(1,length(time));

E1=zeros(1,length(time));

E2=zeros(1,length(time));

E3=zeros(1,length(time));

E4=zeros(1,length(time));

dTclad=zeros(1,length(time));

dTvap=zeros(1,length(time));

dTPT=zeros(1,length(time));

dTCT=zeros(1,length(time));

% Thermal Energy

Qremoved=zeros(1,length(time)); % Total Heat removed by vapour from fuel elements

Qremovedconv=zeros(1,length(time)); % heat removed from fuel elements by convection

Qremovedrad=zeros(1,length(time)); % heat removed from fuel elements by radiation

Qret=zeros(1,length(time)); % heat remaining in the fuel elements

QlossPT=zeros(1,length(time)); % heat lost from PT/CT to moderator

QconvPT=zeros(1,length(time)); % heat removed by convection from PT

QPT=zeros(1,length(time)); % heat remaining within the PT

Qvap=zeros(1,length(time)); % heat going into the vapour

% Nusselt Number

Nu=zeros(1,length(time)); % Nusselt number for fuel pins/vapour interaction


% Prandtl Number

Pr=zeros(1,length(time)); % Prandtl number for PT/vapour interaction



% Reynolds Number

Re=zeros(1,length(time)); % Reynolds number for fuel clad/vapour interaction

% Bundle exit enthalpy

hout=zeros(1,length(time));

% RESISTANCE



%% Initial Property calculations


    
    Re(1,1)=mbundle(1)*Dh/XSteam('my_pT',Peval,Tvap(1,1)-0.1);
    
    Pr(1,1)=XSteam('Cp_pT',Peval,Tvap(1,1)-0.1)*1000*XSteam('my_pT',Peval,Tvap(1,1)-0.1)/XSteam('tc_pT', Peval,Tvap(1,1)-0.1);
    
    if Re<=3000
        Nu(1,1)=4.36;
    else
        
        Nu(1,1)=0.023*Re(1,1)^(4/5)*Pr(1,1)^0.4;
    end
    
    hcool(1,1)=Nu(1,1)*XSteam('tc_pT',Peval, Tvap(1,1)-0.1)/Dh; % W/m^2.K
    
    Tclad(1,1)=(Qbundle*1000000/hcool(1,1)/Afuel)+Tvap(1,1);
    
    
    hpt(1,1)=hcool(1,1);%W/m^2.K
        %% System properties
    % pressure tube thermal conductivity

    kzircPT=12.767-(5.4348e-4*(Tvap(1,1)+273.15))+(8.9818e-6*(Tvap(1,1)+273.15)^2); %W/m.K

    % calandria tube thermal conductivity

    kzircCT=12.767-(5.4348e-4*(Tmod+273.15))+(8.9818e-6*(Tmod+273.15)^2);  %W/m.K

    % CO2 thermal conductivity

    kCO2=[0.01051 0.01456 0.01858 0.02257 0.02652 0.03044 0.03814 0.04565 0.05293 0.08491 0.10688 0.11522]; %W/m.K

    % CO2 thermal conductivity Temperatures

    kCO2Temp=[-50 0 50 100 150 200 300 400 500 1000 1500 2000]; % C

    % CO2 thermal conductivity
    
    TevalCO2=((Tvap(1,1)+Tmod))/2;    

    kCO2sys(1,1)=interp1(kCO2Temp,kCO2,TevalCO2);
    
    R1=1/(hpt(1,1)*Aipt);
    
    R2=log(ropt/ript)/(kzircPT*Aopt);
    
    R3=log(rict/ropt)/(kCO2sys(1,1)*Aict);
    
    R4=log(roct/rict)/(kzircCT*Aoct);
    
    R5=1/(hmod*Aoct);
    
    Rtotal=R1+R2+R3+R4+R5;
    
    QlossPT(1,1)=(Tvap(1,1)-Tmod)/Rtotal;
    
    RPT=R1+R2;
    
    TPT(1,1)=Tvap(1,1)-(QlossPT(1,1)*RPT);
    
    RCTi=R1+R2+R3+R4;
    
    TCT(1,1)=Tvap(1,1)-(QlossPT(1,1)*RCTi);
    
    kgap=0.0476+(0.362e-3*Tclad(1,1))-(0.618e-7*Tclad(1,1)^2)+(0.718e-11*Tclad(1,1)^3)*10^-3; %kW/m.C

    dman=ric-(dfuel/2); %m

    djump=10e-6; %m

    hgap=kgap/(dman+djump);

    Tfuelo(1,1)=(Qel(1)*1000/(dfuel/2*Lchannel)/hgap)+Tclad(1,1);

%% Fuel pin centerline temperature

    
    rfuel=dfuel/2;

    reval=linspace(rfuel,0,divf);

    Tfuel=zeros(1,divf);

    Tfuel(1)=Tfuelo(1,1);

    for pev=2:divf
    
        Tfuel(pev)=(Qvol(1)*1000/4/kUO2(Tfuel(pev-1))*(reval(pev-1)^2-reval(pev)^2))+Tfuel(pev-1);
    
    end
    
    Tfuelmid(1,1)=Tfuel(length(Tfuel));
    
    Tfuelmeat(1,1)=mean(Tfuel);
    
    clear Tfuel
    
    
    hin=XSteam('hV_p',Peval)/1000;


%% Transient calculations

for n=2:length(time)
    
        
        %% HT coefficient for fuel 1in cladding
        
        Re(1,n)=mbundle*Dh/XSteam('my_pT',Peval,Tvap(1,n-1)+0.1);
    
        Pr(1,n)=XSteam('Cp_pT',Peval,Tvap(1,n-1)+0.1)*1000*XSteam('my_pT',Peval,Tvap(1,n-1)+0.1)/XSteam('tc_pT', Peval,Tvap(1,n-1)+0.1);
        
        if Re<=3000
            
            Nu(1,n)=4.36;
        else
        
            Nu(1,n)=0.023*Re(1,n)^(4/5)*Pr(1,n)^0.4;
        
        end
    
        hcool(1,n)=Nu(1,n)*XSteam('tc_pT',Peval, Tvap(1,n-1)+0.1)/Dh;
        
        %% HT coefficient for PT
        
            
          
        hpt(1,n)=hcool(1,n);
        
        eclad=0.325+(0.1246e6*doxide);
        
        ePT=eclad;
        
        eCT=ePT;
        
        hrFPT(1,n)=sigma*(Tclad(1,n-1)+TPT(1,n-1))*((Tclad(1,n-1)^2)+(TPT(1,n-1)^2))/((1/eclad)+((1-ePT)/ePT*Aoutfuel/Aipt));
        
        hCO2rad(1,n)=sigma*(TPT(1,n-1)+TCT(1,n-1))*((TPT(1,n-1)^2)+(TCT(1,n-1)^2))/((1/ePT)+((1-eCT)/eCT*Aopt/Aict));
        
        Cpclad(1,n)=255.66+(0.1024*(Tclad(1,n-1)+273.15)); %J/kg.K
        
        if Tvap(1,n-1)==XSteam('Tsat_p',Peval)
            Cpvap(1,n)=XSteam('CpV_p',Peval)*1000;
        else
            Cpvap(1,n)=XSteam('Cp_pT',Peval,Tvap(n-1))*1000;
        end
        
        if Tvap(1,n-1)==XSteam('Tsat_p',Peval)
            
            hout(1,n)=XSteam('hV_p',Peval)*1000;
        else
            hout(1,n)=XSteam('h_pT',Peval,Tvap(n-1))*1000;
        end
        
        CpPT(1,n)=255.66+(0.1024*(TPT(1,n-1)+273.15)); %J/kg.K
        
        kCO2sys(1,n)=interp1(kCO2Temp,kCO2,((TPT(1,n-1)+TCT(1,n-1))/2)); % W/m.K
        
        hCO2rad(1,n)=sigma*(TPT(1,n-1)+TCT(1,n-1))*((TPT(1,n-1)^2)+(TCT(1,n-1)^2))/((1/ePT)+(((1-eCT)/eCT)*Aopt/Aict));
        
        CpCT(1,n)=255.66+(0.1024*(TCT(1,n-1)+273.15)); %J/kg.K
        
        A1(1,n)=(-Afuel*hcool(1,n)/mclad/Cpclad(1,n))-(Aoutfuel*hrFPT(1,n)/mclad/Cpclad(1,n));
        
        B1(1,n)=(Afuel*hcool(1,n)/mclad/Cpclad(1,n))*Tvap(1,n-1);
        
        C1(1,n)=(Aoutfuel*hrFPT(1,n)/mclad/Cpclad(1,n))*TPT(1,n-1);
        
        D1(1,n)=0;
        
        E1(1,n)=Qbundle/1000000/mclad/Cpclad(1,n);
        
        dTclad(1,n)=(Tclad(1,n-1)*exp(-A1(1,n)*dt))+((1-exp(-A1(1,n)*dt))*(-(B1(1,n)+C1(1,n)+D1(1,n)+E1(1,n))/A1(1,n)));
        
        Tclad(1,n)=Tclad(1,n-1)+dTclad(1,n);
        
        A2(1,n)=hcool(1,n)/mbundle/Cpvap(1,n)*Tclad(1,n-1);
        
        B2(1,n)=-(2*hcool(1,n)/mbundle/Cpvap(1,n));
        
        C2(1,n)=hcool(1,n)/mbundle/Cpvap(1,n)*TPT(1,n-1);
        
        D2(1,n)=0;
        
        E2(1,n)=(hout(1,n)-hin)/Cpvap(1,n);
        
        dTvap(1,n)=(Tvap(1,n-1)*exp(-B2(1,n)*dt))+((1-exp(-A2(1,n)*dt))*(-(A2(1,n)+C2(1,n)+D2(1,n)+E2(1,n))/B2(1,n)));
        
        Tvap(1,n)=Tvap(1,n-1)+dTvap(1,n);
        
        A3(1,n)=hrFPT(1,n)*Aoutfuel/mPT/CpPT(1,n)*Tclad(1,n-1);
        
        B3(1,n)=hcool(1,n)*Aipt/mPT/CpPT(1,n)*Tvap(1,n-1);
        
        C3(1,n)=-(hrFPT(1,n)*Aoutfuel/mPT/CpPT(1,n))-(hcool(1,n)*Aipt/mPT/CpPT(1,n))-(2*pi()*Lbund*kCO2sys(1,n)/mPT/CpPT(1,n)/log(rict/ropt))-(hCO2rad(1,n)*Aopt/mPT/CpPT(1,n));
        
        D3(1,n)=((2*pi()*Lbund*kCO2sys(1,n)/mPT/CpPT(1,n)/log(rict/ropt))+(hCO2rad(1,n)*Aopt/mPT/CpPT(1,n)))*TCT(1,n-1);
        
        E3(1,n)=0;
        
        dTPT(1,n)=(TPT(1,n-1)*exp(-C3(1,n)*dt))+((1-exp(-A3(1,n)*dt))*(-(B3(1,n)+A3(1,n)+D3(1,n)+E3(1,n))/C3(1,n)));
        
        TPT(1,n)=TPT(1,n-1)+dTPT(1,n);
        
        A4(1,n)=0;
        
        B4(1,n)=0;
        
        C4(1,n)=((2*pi()*Lbund*kCO2sys(1,n)/mCT/CpCT(1,n)/log(rict/ropt))+(hCO2rad(1,n)*Aopt/mCT/CpCT(1,n)))*TPT(1,n-1);
        
        D4(1,n)=((-((2*pi()*Lbund*kCO2sys(1,n)/mCT/CpCT(1,n)/log(rict/ropt))+(hCO2rad(1,n)*Aopt/mCT/CpCT(1,n))))-(hmod*Aoct/mCT/CpCT(1,n)));
        
        E4(1,n)=hmod*Aoct/mCT/CpCT(1,n)*Tmod;
        
        dTCT(1,n)=(TCT(1,n-1)*exp(-D4(1,n)*dt))+((1-exp(-A4(1,n)*dt))*(-(B4(1,n)+A4(1,n)+C4(1,n)+E4(1,n))/D4(1,n)));
        
        TCT(1,n)=TCT(1,n-1)+dTCT(1,n);   
end

plot(time,Tclad);