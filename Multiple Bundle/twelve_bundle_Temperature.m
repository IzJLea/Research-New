%% Initial Variables
%

Qchannel=0.150; %MW channel power

bunds=12; %number of bundles

H=3; %m height difference from inlet/outlet headers

hmod=1000; %W/m^2.K heat transfer coefficient 

Tenter=100; %C Entering coolant temperature

PSH=114; %kPa supply header pressure

PVH=134; %kPa vent header pressure

Tmod=60; %C Moderator temperature

Aflow=0.0035; %m^2 Flow area in channel

Dh=0.0074; %m hydraulic diameter of bundle

DPT=0.10338; %m inner pressure tube diameter

DCT=0.12869;%m inner calandria tube thickness

doutclad=0.0138; %m cladding outer diameter

dfuel=0.0122; %m fuel outer diameter

tclad=0.00038; %m thickness of cladding

Lchannel=6; %m length of fuel channel

tPT=0.00424; %m pressure tube thickness

tCT=0.0014; %m calandria tube thickness

sigma=5.670373e-8; %Stefan-Boltzmann constant

eclad=0.6988; % dimensionless emissivity value for zirconium cladding
        
ePT=eclad;% emissivity for pressure tube
        
eCT=ePT;% emissivity for calandria tube


%% Initial Calculated Variables



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



%% Resistance Calculation from liquid flow data

DP = [165;
262.5000;
520.5000;
258;
174];  % kPa measured reference pressure drops across the inlet feeder, end fitting, core, and exit end fitting and feeder respectively

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


%% initial conditions

iter=100;

damp=1;

err=zeros(1,iter);

Peval=zeros(1,iter);

Tvap=zeros(1,iter);

Peval(1,1)=(((PSH+PVH)/2)+(XSteam('rho_pT',PSH/100,Tenter)*9.81*H/1000))/100;

Tvap(1,1)=XSteam('Tsat_p',Peval(1,1))+1;

Peval(1,1)=(((PSH+PVH)/2)+((XSteam('rho_pT',Peval(1,1),Tenter)-XSteam('rho_pT',Peval(1,1),Tvap(1,1)))*9.81*H/1000))/100;
    
Tvap(1,1)=XSteam('Tsat_p',Peval(1,1))+1;
    
Peval(1,1)=(((PSH+PVH)/2)+((XSteam('rho_pT',Peval(1,1),Tenter)-XSteam('rho_pT',Peval(1,1),Tvap(1,1)))*9.81*H/1000))/100;
    
Tvap(1,1)=XSteam('Tsat_p',Peval(1,1))+1;

Res=Alphacalc(PSH,PVH,Tenter,Tvap(1,1),Qchannel,RCH,RF,H);

