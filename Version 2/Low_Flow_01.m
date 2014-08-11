%% Description of System Power (equal linear distribution of power to bundles)


Qchannel=.150; %MW

H=1; %m

Tenter=100; %C

PSH=134; %kPa

PVH=114; %kPa

Peval=(PSH)/100;

Tmod=60; %C

Aflow=0.0035; %m^2

Dh=0.0074; %m

DPT=0.10338;

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

d=0.001; %dampening factor

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

Tvap=zeros(1,12);

hvap=zeros(1,12);

for m=1:length(Tvap)
    
    hvap(m)=(alpha*Qtotal(m)*1000/wv)+hg;
    
    Tvap(m)=XSteam('T_ph',Peval,hvap(m));
    
end

Atr=2;
    
Ar=alpha*2;
    
Ainr=Atr-Ar;
    
xflow=acos(Ainr-1);
    
if xflow>pi()/2
    
        theta=asin(sin(xflow)/(pi()/2));
        
else if xflow==pi()/2;
            
        theta=pi()/2;
            
    else
        theta=(pi()/2)+asin(sin(xflow)/(pi()/2));
    end
end
        

    
    thetapipe=2*theta;
    
    Dflow=8*Aflow*alpha/thetapipe/DPT;

%% vapor heat transfer coefficient for each bundle

mbundle=zeros(1,length(Qbundle));

vsteam=zeros(1,length(Qbundle));

htfe=zeros(1,length(Qbundle));

Reysteamf=zeros(1,length(Qbundle));

Reysteampt=zeros(1,length(Qbundle));

for n=1:length(Qbundle);
    
    mbundle(n)=Qtotal(n)/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval));
    
    vsteam(n)=mbundle(n)/XSteam('rhoV_p',Peval)*Aflow*alpha;
    
    Reysteamf(n)=mbundle(n)*XSteam('rho_pT',Peval,Tvap(n))*Dh/XSteam('my_pT',Peval,Tvap(n));
    
      
      
    
    Reysteampt(n)=mbundle(n)*XSteam('rho_pT',Peval,Tvap(n))*Dflow*alpha/XSteam('my_pT',Peval,Tvap(n));
    
    %htfe(n)=
end




plotyy(1:12,Tvap,1:12,Qtotal,'plot','plot');

Trun=600; %s

Tr=0:0.01:Trun;

Res=single_channel_Lockhart_Martinelli_Qloss_LPXSTEAM5(Qchannel,Tenter,Tmod,PVH/1000);

TFs=zeros(length(Qbundle),length(Tr));

TPt=zeros(length(Qbundle),length(Tr));

TCt=zeros(length(Qbundle),length(Tr));

TFs(1:length(Qbundle),1)=Res(3,1);

TPt(1:length(Qbundle),1)=Res(7,1);

TCt(1:length(Qbundle),1)=Res(8,1);






