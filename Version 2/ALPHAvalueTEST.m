%% Description of System Power (equal linear distribution of power to bundles)


Qchannel=.100; %MW

H=1; %m

Tenter=100; %C

PSH=134; %kPa

PVH=114; %kPa

Peval=(PSH+(XSteam('rhoL_p',PSH/100)*9.81*H/1000))/100;



%% delta P calculation

deltaP=(PSH)-PVH+((9.81*H*(XSteam('rhoL_p',Peval)-XSteam('rhoV_p',Peval)))/1000);

%% Resistance Calculation from liquid flow data

DP = [165;
262.5000;
520.5000;
258;
174];


%% calculate values for k for each section

% reference mass flowrate

% M=25.8; % kg/s
% 
% Rho=780.6; %kg/m^3
% 
% keff=DP./1000./(M^2);
% 
% reff=keff.*Rho;

RCH=2.06e6;

RF=1.5e6;

%% alpha calculation

alpha=0;

x=(XSteam('h_pT',Peval,Tenter)-XSteam('hL_p',Peval))/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval));

wsmx=Qchannel*1000/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval));

d=0.01; %dampening factor

delta=0.1;

 delta>=0.0001
    
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
    
