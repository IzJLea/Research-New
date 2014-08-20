%% Description of System Power (equal linear distribution of power to bundles)


Qchannel=.150; %MW

H=0; %m

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

dman=ric-(dfuel/2); %m

djump=10e-6; %m

Lchannel=5.94; %m

tPT=0.00424;

tCT=0.0014;
%% delta P calculation

deltaP=(PSH)-PVH+((9.81*H*(XSteam('rhoL_p',Peval)-XSteam('rhoV_p',Peval)))/1000);

%% Resistance Calculation from liquid flow data

DP = [165;
262.5000;
520.5000;
258;
174];




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

Qbundle=zeros(1,12);

Ptotal=Qchannel/(2*sin(pi()/2));

Qbundle(1,1)=(1-sin(5*pi()/12))*Ptotal;

Qbundle(1,12)=Qbundle(1,1);

Qbundle(1,2)=(sin(5*pi()/12)-sin(pi()/3))*Ptotal;

Qbundle(1,11)=Qbundle(1,2);

Qbundle(1,3)=(sin(pi()/3)-sin(pi()/4))*Ptotal;

Qbundle(1,10)=Qbundle(1,3);

Qbundle(1,4)=(sin(pi()/4)-sin(pi()/6))*Ptotal;

Qbundle(1,9)=Qbundle(1,4);

Qbundle(1,5)=(sin(pi()/6)-sin(pi()/12))*Ptotal;

Qbundle(1,8)=Qbundle(1,5);

Qbundle(1,6)=sin(pi()/12)*Ptotal;

Qbundle(1,7)=Qbundle(1,6);

Qtotal=zeros(1,12);

for l=1:length(Qtotal)
    
    Qtotal(l)=sum(Qbundle(1:l));
    
end

%% Bundle Mass flow

mbundle=zeros(1,length(Qbundle));

Tvap=zeros(1,12);

hvap=zeros(1,12);

mchange=zeros(1,12);

min=zeros(1,12);

for m=1:length(Tvap)
    
    hvap(m)=(alpha*Qtotal(m)*1000/wv)+hg;
    
    Tvap(m)=XSteam('T_ph',Peval,hvap(m));
    
    mbundle(m)=sum(mbundle(1:m))+((1-alpha)*Qbundle(m)*1000/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval)));
    
    mchange(m)=((1-alpha)*Qbundle(m)*1000/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval)));
    
    min(m)=mbundle(m)-mchange(m);
    
end

Trun=500; %run time

divt=1; %time step

time=0:divt:500;

vbund=zeros(length(Qbundle),length(time));

Tfs=zeros(length(Qbundle),length(time));

Tvapt=zeros(length(Qbundle),length(time));

hvapt=zeros(length(Qbundle),length(time));

hevap=XSteam('hV_p',Peval)-XSteam('hL_p',Peval);

hin=zeros(length(Qbundle),length(time));

for n=1:length(Qbundle)
    
    Tfs(n,1)=XSteam('Tsat_p',Peval);
    
    Tvapt(n,1)=XSteam('Tsat_p',Peval);
    
    hvapt(n,1)=XSteam('hV_p',Peval);
    
       
end

for tind=2:length(time)
    
    for bundind=1:length(Qbundle)
        
        if bundind==1
            
            hvapt(bundind,tind)=(Qbundle(bundind)*1000*alpha/mbundle)+;
        




