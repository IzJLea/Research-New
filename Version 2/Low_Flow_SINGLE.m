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

Lbund=Lchannel/12;

tPT=0.00424;

tCT=0.0014;

Aipt=pi()*DPT*Lbund;

Aopt=pi()*(DPT+(2*tPT))*Lbund;

Aict=pi()*DCT*Lbund;

Aoct=pi()*(DCT+(2*tCT))*Lbund;

Afuel=pi()*doutclad*Lbund;

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
    
    
    


Lrun=100; %length of run (s)

dt=0.01; %time division

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

% Heat Capacity

Cpvap=zeros(1,length(time)); % Vapour Heat capacity

Cpvap(1,1)=XSteam('CpV_p',Peval); % initial vapour heat capacity assignment

% Density

rhovap=zeros(1,length(time)); % vapour density

rhovap(1,1)=XSteam('rhoV_p',Peval); % vapour density assignment

% HT coefficient

hclad=zeros(1,length(time)); % heat transfer coefficient between fuel and vapour

hpt=zeros(1,length(time)); % heat transfer coefficient between PT and vapour

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

NuPT=zeros(1,length(time)); % Nusselt number for PT/vapour interaction

% Prandtl Number

Pr=zeros(1,length(time)); % Prandtl number for PT/vapour interaction

PrPT=zeros(1,length(time)); % Prandtl number for PT/vapour interaction

% Reynolds Number

Re=zeros(1,length(time)); % Reynolds number for PT/vapour interaction

RePT=zeros(1,length(time)); % Reynolds number for PT/vapour interaction

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
    
    hclad(1,1)=Nu(1,1)*XSteam('tc_pT',Peval, Tvap(1,1)-0.1)/Dh; % W/m^2.K
    
    Tclad(1,1)=(Qel*1000/hclad(1,1)/Afuel)+Tvap(1,1);
    
    RePT(1,1)=mbundle(1)*DhPT/XSteam('my_pT',Peval,Tvap(1,1)-0.1); %%% checked to here!!!! OCTOBER 8 2014
    
    PrPT(1,1)=XSteam('Cp_pT',Peval,Tvap(1,1)-0.1)*1000*XSteam('my_pT',Peval,Tvap(1,1)-0.1)/XSteam('tc_pT', Peval,Tvap(1,1)-0.1);
    
    NuPT(1,1)=0.023*RePT(1,1)^(4/5)*PrPT(1,1)^0.4;
    
    hpt(1,1)=Nu(1,1)*XSteam('tc_pT',Peval, Tvap(1,1)-0.1)/DPT;
        %% System properties
    % pressure tube thermal conductivity

    kzircPT=12.767-(5.4348e-4*(Tvap(1,1)+273.15))+(8.9818e-6*(Tvap(1,1)+273.15)^2); %W/m.K

    % calandria tube thermal conductivity

    kzircCT=12.767-(5.4348e-4*(Tmod+273.15))+(8.9818e-6*(Tmod+273.15)^2);  %W/m.K

    % CO2 thermal conductivity

    kCO2=[14.60e-3 16.23e-3 17.87e-3 19.52e-3 21.18e-3 22.84e-3 27.00e-3 31.12e-3 35.20e-3 39.23e-3];

    % CO2 thermal conductivity Temperatures

    kCO2Temp=[0 20 40 60 80 100 150 200 250 300];

    % CO2 thermal conductivity
    
    TevalCO2=((Tvap(1,1)+Tmod)+273.15)/2;    

    kCO2sys=interp1(kCO2Temp,kCO2,TevalCO2);
    
    R1=1/(hpt(1,1)*Aipt);
    
    R2=log(ropt/ript)/(kzircPT*Aopt);
    
    R3=log(rict/ropt)/(kCO2sys*Aict);
    
    R4=log(roct/rict)/(kzircCT*Aoct);
    
    R5=1/(hmod*Aoct);
    
    Rtotalinv=(1/R1)+(1/R2)+(1/R3)+(1/R4)+(1/R5);
    
    Rtotal=Rtotalinv^-1;
    
    QlossPT(1,1)=(Tvap(1,1)-Tmod)/Rtotal;
    
    RPT=1/((1/R1)+(1/R2));
    
    TPT(1,1)=Tvap(1,1)-(QlossPT(1,1)*RPT);
    
    RCTi=1/((1/R1)+(1/R2)+(1/R3)+(1/R4));
    
    TCT(1,1)=Tvap(1,1)+(QlossPT(1,1)*RCTi);
    
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
    


%% Transient calculations

for n=2:length(time)
    
        
        %% HT coefficient for fuel 1in cladding
        
        Re(1,n)=mbundle(1)*Dh/XSteam('my_pT',Peval,Tvap(1,n-1)+0.1);
    
        Pr(1,n)=XSteam('Cp_pT',Peval,Tvap(1,n-1)+0.1)*1000*XSteam('my_pT',Peval,Tvap(1,n-1)+0.1)/XSteam('tc_pT', Peval,Tvap(1,n-1)+0.1);
        
        if Re<=3000
            
            Nu(1,n)=4.36;
        else
        
            Nu(1,n)=0.023*Re(1,n)^(4/5)*Pr(1,n)^0.4;
        
        end
    
        hclad(1,n)=Nu(1,n)*XSteam('tc_pT',Peval, Tvap(1,n-1)+0.1)/Dh;
        
        %% HT coefficient for PT
        
        RePT(1,n)=mbundle(1)*DPT/XSteam('my_pT',Peval,Tvap(1,n-1)+0.1);
        
        PrPT(1,n)=XSteam('Cp_pT',Peval,Tvap(1,n-1)+0.1)*1000*XSteam('my_pT',Peval,Tvap(1,n-1)+0.1)/XSteam('tc_pT',Peval,Tvap(1,n-1)+0.1);
        if RePT<=3000
            
            NuPT(1,n)=4.36;
        else
            NuPT(1,n)=0.023*RePT(1,n)^(4/5)*PrPT(1,n)^0.4;
        end
        
        hpt(1,n)=NuPT(1,n)*XSteam('tc_pT',Peval,Tvap(1,n-1)+0.1)/DhPT;
        
        %% Calculate heat removed, lost, generated in fuel pins
        
        Qremovedconv(1,n)=hclad(1,n)*Afuel*(Tclad(1,n-1)-Tvap(1,n-1))/1000000;
        
        if Qremovedconv(1,n)>=Qel(1)
            
            Qremovedconv(1,n)=Qel(1);
        end
        
        ef=0.325+(0.1246e6*doxide);
        
        ept=ef;
        
        Qremovedrad(1,n)=((sigma*Afuel*((Tclad(1,n-1))^4-TPT(1,n-1)^4)/((1/ef)+(Afuel/(Aipt/9)*((1/ept)-1))))/1000000);
        
        Qremoved(1,n)=Qremovedconv(1,n)+Qremovedrad(1,n);
        
        Qret(1,n)=(Qbundle(1)/37)-Qremoved(1,n);
        
        %% Calculate heat loss/retention in Pressure Tube
        
        QconvPT(1,n)=hpt(1,n)*Afuel*(Tvap(1,n-1)-TPT(1,n-1))/1000000;
        
        TmeanCO2=(Tvap(1,n-1)+Tmod);
        
        if TmeanCO2>=300
            TmeanCO2=300;
        end
        
        kCO2=[14.60e-3 16.23e-3 17.87e-3 19.52e-3 21.18e-3 22.84e-3 27.00e-3 31.12e-3 35.20e-3 39.23e-3];

        % CO2 thermal conductivity Temperatures

        kCO2Temp=[0 20 40 60 80 100 150 200 250 300];

        % CO2 thermal conductivity

        kCO2sys=interp1(kCO2Temp,kCO2,TmeanCO2);
        
        RCO2=log(rict/ropt)/(2*pi()*kCO2sys*Lchannel);
        
        kzircCT=12.767-(5.4348e-4*Tmod)+(8.9818e-6*Tmod^2);

        RCT=log(roct/rict)/(2*pi()*kzircCT*Lchannel);

        Rmod=1/(hmod*DCT*pi()*Lchannel);
        
        UA=(1/RCO2)+(1/RCT)+(1/Rmod);
        
        QlossPT(1,n)=UA*(TPT(1,n-1)-Tmod)/1000000;
        
        QPT(1,n)=(37*Qremovedrad(1,n))-QlossPT(1,n)-QconvPT(1,n);
        
        Qvap(1,n)=(37*Qremovedconv(1,n))+QconvPT(1,n);
       
        %% Temperature calculations
        % fuel/cladding temperatures
        tau=(Tfuelmeat(1,n-1)+273.15)/1000;
        
        Cpuo2=52.1743+(87.951*tau)-(84.2411*tau^2)+(31.542*tau^3)-(2.6334*tau^4)-(0.71391*tau^-2);
        
        Tfuelmeat(1,n)=(Tfuelmeat(1,n-1)+(Qret(1,n)*1000000/Cpuo2/mfuel))/270.03;
        
        Tfuelo(1,n)=Tfuelo(1,n-1)+(Qret(1,n)*1000000/Cpuo2/mfuel);
        
        Tfuelmid(1,n)=Tfuelmid(1,n-1)+(Qret(1,n)*1000000/Cpuo2/mfuel);
        
        Tclad(1,n)=Tfuelo(1,n);
        
        CpPT=255.66+(0.1024*(TPT(1,n-1)+273.15));
        
        TPT(1,n)=TPT(1,n-1)+(QPT(1,n)*1000000/mPT/CpPT);
        
        if 1==1
            
            hin=XSteam('hV_p',Peval);
        else
            
            hin=1/mbundle(1)*((XSteam('hV_p',Peval)*mchange(1))+(hout(1-1,n-1)*mbundle(1-1)));
        end
        
        Cpvap(1,n)=XSteam('Cp_pT',Peval,Tvap(1,n-1));
        
        hout(1,n)=hin+(Qvap(1,n)*1000/mbundle(1)/Cpvap(1,n));
        
        Tvap(1,n)=XSteam('T_ph',Peval,hout(1,n));
        
        TCT(1,n)=(QlossPT(1,n)*RCO2)+TCT(1,n-1);
   
end