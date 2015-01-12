%% Initial Variables
%

Qchannel=.150; %MW channel power

bunds=12; %number of bundles

H=0; %m height difference from inlet/outlet headers

hmod=1000; %W/m^2.K heat transfer coefficient 

Tenter=200; %C Entering coolant temperature

PSH=210; %kPa supply header pressure

PVH=190; %kPa vent header pressure

Tmod=60; %C Moderator temperature

Aflow=0.0035; %m^2 Flow area in channel

Dh=0.0074; %m hydraulic diameter of bundle

DPT=0.10338; %m inner pressure tube diameter

DCT=0.12869;%m inner calandria tube thickness

doutclad=0.0138; %m cladding outer diameter

dfuel=0.0122; %m fuel outer diameter

tclad=0.00038; %m thickness of cladding

Lchannel=5.94; %m length of fuel channel

tPT=0.00424; %m pressure tube thickness

tCT=0.0014; %m calandria tube thickness

sigma=5.670373e-8; %Stefan-Boltzmann constant

eclad=0.6988; % dimensionless emissivity value for zirconium cladding
        
ePT=eclad;% emissivity for pressure tube
        
eCT=ePT;% emissivity for calandria tube


%% Initial Calculated Variables

Peval=110; % system evaluation pressure

roc=doutclad/2; %m outer cladding radius

ric=roc-tclad; %m inner cladding radius

Lbund=Lchannel/bunds; %m length of bundle (for this initial simulation, entire channel will be simulated as one bundle)

Aipt=pi()*DPT*Lbund;%m^2 inner pressure tube area

Aopt=pi()*(DPT+(2*tPT))*Lbund;%m^2 outer pressure tube area

Aict=pi()*DCT*Lbund;%m^2 inner calandria tube area

Aoct=pi()*(DCT+(2*tCT))*Lbund;%m^2 outer calandria tube area

Afuel=pi()*doutclad*Lbund;%m^2 fuel outer area

Arad=9/37*Afuel; % Area of outer fuel elements that sees PT

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

rhoref=rhoZirc(TrefK);

mclad=((doutclad)^2-(doutclad-(2*tclad))^2)/4*Lbund*rhoref; %mass of zirconium in cladding of one element

mPT=((DPT+(2*tPT))^2-(DPT)^2)/4*Lbund*rhoref; %kg mass of zirconium in pressure tube

mCT=((DCT+(2*tCT))^2-(DPT)^2)/4*Lbund*rhoref; %kg mass of zirconium in calandria tube

%% mass of fuel

mfuel=Lbund*pi()/4*dfuel^2*rhoUO2(Tref); %kg mass of fuel within one fuel element

%% Time determination 
% This section sets the length of the simulation and time divisions an

Time=1000;  %seconds total run time

div=0.1; %s time step length

ind=Time/div; % index to determine time step number

time=zeros(1,(Time/div));

for p=2:length(time)
    
    time(1,p)=time(1,p-1)+div;
end

%% Property matrix Creation

% This section creates all property matrices for calculations
%Temperatures
Tclad=zeros(bunds,ind); %C cladding temperature

Tvap=zeros(bunds,ind);%C vapor coolant temperature

TPT=zeros(bunds,ind);%C pressure tube temperature

Tfuel=zeros(bunds,ind);%C average fuel temperature

TCT=zeros(bunds,ind);%C calandria tube temperature

%Material properties

kuo2=zeros(bunds,ind);% W/m.K thermal conductivity of UO2 fuel 

kclad=zeros(bunds,ind);% W/m.K thermal conductivity of cladding

kPT=zeros(bunds,ind);% W/m.K pressure tube thermal conductivity

kCT=zeros(bunds,ind);%W/m.K calandria tube thermal conductivity

kCO2sys=zeros(bunds,ind);%W/m.K CO2 insulator thermal conductivity

Cpclad=zeros(bunds,ind);%J/kg.K cladding heat capacity

CpPT=zeros(bunds,ind);%J/kg.K pressure tube heat capacity

CpCT=zeros(bunds,ind);%J/kg.K calandria tube heat capacity

Cpvap=zeros(bunds,ind);%J/kg.K vapor coolant heat capacity

Cpfuel=zeros(bunds,ind);%J/kg.K fuel heat capacity

% Heat transfer properties

Reynolds=zeros(bunds,ind);% reynolds number of coolant

Prandtl=zeros(bunds,ind);% coolant prandtl number

Nusselt=zeros(bunds,ind);% coolant nusselt number

hcool=zeros(bunds,ind);% coolant convective heat transfer coefficient 

hrad=zeros(bunds,ind);%W/m^2.K radiative heat transfer coefficient for cladding/pressure tube

hrpt=zeros(bunds,ind);% W/m^2.K radiative heat transfer coefficient for calandria tube/pressure tube

hin=zeros(bunds,ind);

hout=zeros(bunds,ind);% channel exit enthalpy

hfluid=zeros(bunds,ind);

hvap=zeros(bunds,ind);

%Heat Transfer Resistances


R1=zeros(bunds,ind); 

R2=zeros(bunds,ind);

R3=zeros(bunds,ind);

R4=zeros(bunds,ind);

R5=zeros(bunds,ind);

%Calculation Constants

A1=zeros(bunds,ind);
    
B1=zeros(bunds,ind);
    
F1=zeros(bunds,ind);
    
A2=zeros(bunds,ind);
    
B2=zeros(bunds,ind);
    
C2=zeros(bunds,ind);
    
D2=zeros(bunds,ind);
       
B3=zeros(bunds,ind);
    
C3=zeros(bunds,ind);
    
D3=zeros(bunds,ind);

F3=zeros(bunds,ind);

B4=zeros(bunds,ind);

C4=zeros(bunds,ind);

D4=zeros(bunds,ind);

E4=zeros(bunds,ind);

D5=zeros(bunds,ind);

E5=zeros(bunds,ind);

F5=zeros(bunds,ind);

%Miscellaneous

CH=zeros(5,ind);

Mcool=zeros(bunds,ind); %kg mass of coolant within the channel

Qel=zeros(1,bunds);

mvap=zeros(1,bunds);

Alphas=zeros(1,bunds);
%% Power profile calculation

Bundpower=cosPower(bunds,Qchannel)*1000000;

for n=1:bunds
    Qel(1,n)=Bundpower(1,n)/37;
end
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
clear delta
%constant mass flux calculation

Mflux=wsmx/Aflow*((1-alpha)/alpha); % average channel mass flux

mflow=Mflux*Aflow*alpha;

%initial enthalpy values

hin(1,1:ind)=XSteam('h_pT',Peval,Tenter)*1000;

Alphas(1,1)=0;

hout(1,1)=hin(1,1)+(Bundpower(1,1)*Alphas(1,1));

if hout(1,1)<=XSteam('hL_p',Peval)
    Alphas(1,1)=0;
    display(BUTTS);
else
    delta=1;
    mvap=Bundpower(1,1)*(1-Alphas(1,1))/1000/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval));
    while delta>=0.0001
        alphap=mvap/Aflow/Mflux;
        err=alphap-Alphas(1,1);
        delta=abs(err);
        Alphas(1,1)=Alphas(1,1)-(d*err);
        mvap=Bundpower(1,1)*Alphas(1,1)/1000/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval));
    end
    
    hout(1,1)=Alphas(1,1)*XSteam('hV_p',Peval)+(1-Alphas(1,1))*XSteam('hL_p',Peval);
end
        


%% Initial property calculations
%in this section the initial temperature for the elements and the coolant
%are determined. The temperatures were determined using the assumption that
%saturated fluid was previously flowing in the channel. Mass flux is
%assumed to be constant at channel average value.


for m=1:bunds
        
    Tvap(m,1)=XSteam('Tsat_p',Peval);
    
    Reynolds(m,1)=mflow*Dh/Aflow/XSteam('my_pT',Peval,Tvap(m,1)-0.1);
    
    Prandtl(m,1)=XSteam('Cp_pT',Peval,Tvap(m,1)-0.1)*1000*XSteam('my_pT',Peval,Tvap(m,1)-0.1)/XSteam('tc_pT', Peval,Tvap(m,1)-0.1);
    
    if Reynolds(m,1)<=3000
        
        Nusselt(m,1)=4.63;
        
    else
        
        Nusselt(m,1)=0.023*Reynolds(m,1)^(4/5)*Prandtl(m,1)^0.4;
    end
    
    hcool(m,1)=Nusselt(m,1)*XSteam('tcL_T',Tvap(m,1))/Dh; %W/m^2/K
    
    Tclad(m,1)=(Qel(1,m)/hcool(m,1)/Afuel)+Tvap(m,1);
 
    R1(m,1)=0.5*Rzirc(Tclad(m,1),ric,roc,Lbund)+Rgap(Tclad(m,1),Tclad(m,1))+0.5*Rfuel(Tclad(m,1),Lbund);
    
    Tfuel(m,1)=(Qel(1,m)*R1(m,1))+Tclad(m,1);
    
    R3=(1/hcool(m,1)/Aipt)+0.5*Rzirc(Tvap(m,1),ript,ropt,Lbund);
    
    TPT(m,1)=Tvap(m,1);
    
    TCT(m,1)=Tmod;
    
    
        
end

% for n=2:ind
%     for p=1:bunds
        
            
