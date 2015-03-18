%% Initial Variables
%

Qchannel=.150; %MW channel power

bunds=12; %number of bundles

H=0; %m height difference from inlet/outlet headers

hmod=1000; %W/m^2.K heat transfer coefficient 

Tenter=300; %C Entering coolant temperature

PSH=100; %kPa supply header pressure

PVH=80; %kPa vent header pressure

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

Peval=(PSH+PVH)/2; % system evaluation pressure

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

%% Voidfraction/Exit temperature relation 

Tvap=zeros(1,100);

Tvap(1,1)=XSteam('Tsat_p',Peval);

Tout=zeros(1,100);

Tout(1,1)=XSteam('Tsat_p',Peval);

M1=zeros(1,100);

M2=zeros(1,100);

deltas=zeros(100,100);

for i=2:length(Tvap)
    
    Tmax=3000; % degrees celsius
    
    Tvap(1,i)=Tvap(1,i-1)+((Tmax-Tvap(1,1))/length(Tvap));
    
    Tout(1,i)=Tout(1,i-1)+((Tmax-Tout(1,1))/length(Tout));
    
    M1(1,i)=Qchannel*1000*(XSteam('h_pT',Peval,Tvap(1,i))-XSteam('h_pT',Peval,Tenter));
    
    Res=Alphacalc(PSH,PVH,Tenter,Tout(1,i),Qchannel,RCH,RF);
    
    M2(1,i)=Res(2,1);
end

plotyy(Tvap,M1,Tvap,M2)



            
