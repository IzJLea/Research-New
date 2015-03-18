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


Peval=(PSH+PVH)/2;

%% Flow Resistance calculation

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

RCH=sum(reff(3:4)); %effective resistance of fuel channel

RF=sum(reff(5)); %effective resistance of end fitting and feeder

%% calculation of average channel alpha

Tvap=XSteam('Tsat_p',Peval)+5;

x=zeros(2,1);

sub=(XSteam('hL_p',Peval)-XSteam('h_pT',Peval,Tenter))/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval));

wsmx=Qchannel*1000/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval));

rhoave=(XSteam('rhoV_p',Peval)+XSteam('rho_pT',Peval,Tvap))/2;

alpha=0.1;

for i=1:200

    a=alpha^2*rhoave/XSteam('rhoV_p',Peval);
    
    alpha=1/(1+((1+sub)/wsmx*sqrt((PSH-PVH)*rhoave/(RCH+(a*RF)))));
       
   
end

Mflow=alpha*sqrt((PSH-PVH)*rhoave/(RCH+(a*RF)));

  
