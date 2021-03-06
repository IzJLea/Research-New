%% Initial Variables
%

Qchannel=0.050; %MW channel power

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

Peval=(((PSH+PVH)/2)+(XSteam('rho_pT',PSH/100,Tenter)*9.81*H/1000))/100; % system evaluation pressure

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

%% Voidfraction/Exit temperature/mass flow relation for heat up to final temperature

%% Set exit temperature and properties

Tout=XSteam('Tsat_p',Peval)+200; %Exit vapor temperature

Temps=linspace(Tout,Tout+100,100);

alpha=zeros(1,length(Temps));

Length=linspace(0,Lchannel,100);

Lsection=Lchannel/length(Length);

Heights=zeros(1,length(Length));

Hmatrix=linspace(0,H+6,100);

Heightchange=zeros(1,length(Length));

Avap=zeros(1,length(Heights));

Achange=zeros(1,length(Heights));

AlphaMax=Alphacalc(PSH,PVH,Tenter,Tout,Qchannel,RCH,RF,H);

mvap=zeros(1,length(Heights));

for i=1:length(Length)
    
    RCHL=RCH/Lchannel*Length(1,i);
    
    QchannelL=Qchannel/Lchannel*Length(1,i);
    
    Res=Alphacalc(PSH,PVH,Tenter,Tout,QchannelL,RCHL,((RCH/Lchannel*(Lchannel-Length(1,i)))+RF),H);
    
    alpha(1,i)=Res(1,1);
    
    mvap(1,i)=Res(2,1);
    
    if i==1
        
        Heights(1,i)=DPT;
        
    else  
        
        Heights(1,i)=DPT-(alpha(1,i)*DPT);      
                  
    end
    
end
Power=linspace(Qchannel,Qchannel+0.4,41);

void=zeros(length(Power),length(Heights));

mflow=zeros(length(Power),length(Heights));

dP=zeros(length(Power),length(Heights));

rhoout=zeros(length(Power),length(Heights));

houtrev=zeros(length(Power),length(Heights));

Toutrev=zeros(length(Power),length(Heights));

moutrev=zeros(length(Power),length(Heights));

for k=1:length(Power)
    for j=1:length(Heights)
    
        Res=Alphacalc(PSH,PVH,Tenter,Tout,Power(1,k),RCH,RF,Hmatrix(1,j));

        void(k,j)=Res(1,1);
    
        mflow(k,j)=Res(2,1);
    
        dP(k,j)=Res(4,1);
    
        rhoout(k,j)=Res(5,1);
    
        houtrev(k,j)=Res(6,1);
        
        Toutrev(k,j)=Res(7,1);
        
        moutrev(k,j)=Res(8,1);
    end
end

minmax=[mflow(1,1:length(mflow));mflow(size(mflow,1),1:length(mflow))];

mcheck=[mflow(21,1:length(mflow));moutrev(21,1:length(mflow))];

plot(Hmatrix,mflow);

figure

plot(Hmatrix,mcheck);

% figure
% 
% plot(Hmatrix,Toutrev(21,1:length(mflow)));


% 
% figure
% 
% plot(Length,Heights);
% 
% figure
% 
% plot(Length,mvap);
% 
% figure
% 
% plot(Hmatrix,dP);


    
    
    