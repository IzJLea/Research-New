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

Tvap=zeros(1,12);

hvap=zeros(1,12);

for m=1:length(Tvap)
    
    hvap(m)=(alpha*Qtotal(m)*1000/wv)+hg;
    
    Tvap(m)=XSteam('T_ph',Peval,hvap(m));
    
end

  
Atr=alpha*2;

lx=acos(1-Atr);

Dflow=2*Aflow*alpha/(pi()*sin(lx)*DPT/2);




%% vapor heat transfer coefficient for each bundle

mbundle=zeros(1,length(Qbundle));

vsteam=zeros(1,length(Qbundle));

htfe=zeros(1,length(Qbundle));

htpt=zeros(1,length(Qbundle));

htctpt=zeros(1,length(Qbundle));

htctnat=1000;  %W/m.K natural ht coefficient from CT to PT

Reysteamf=zeros(1,length(Qbundle));

Reysteampt=zeros(1,length(Qbundle));

Prgsys=zeros(1,length(Qbundle));

Nuf=zeros(1,length(Qbundle));

Nupt=zeros(1,length(Qbundle));

Qel=zeros(1,length(Qbundle));

for n=1:length(Qbundle);
    
    mbundle(n)=Qtotal(n)*(1-alpha)/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval));
    
    vsteam(n)=mbundle(n)/XSteam('rhoV_p',Peval)*Aflow*alpha;
    
    Reysteamf(n)=mbundle(n)*XSteam('rho_pT',Peval,Tvap(n))*Dh/XSteam('my_pT',Peval,Tvap(n));     
    
    Reysteampt(n)=mbundle(n)*XSteam('rho_pT',Peval,Tvap(n))*Dflow*alpha/XSteam('my_pT',Peval,Tvap(n));
    
    Prgsys(n)=XSteam('Cp_pT',Peval,Tvap(n))*XSteam('my_pT',Peval,Tvap(n))/XSteam('tc_pT',Peval,Tvap(n))*1000;
    
    Nuf(n)=0.023*Reysteamf(n)^(4/5)*Prgsys(n)^0.4;
    
    Nupt(n)=0.023*Reysteampt(n)^(4/5)*Prgsys(n)^0.3;
    
    htfe(n)=Nuf(n)*XSteam('tc_pT',Peval,Tvap(n))/Dh/1000; %fue-sheath ht coefficient
    
    htpt(n)=Nupt(n)*XSteam('tc_pT',Peval,Tvap(n))/Dflow/1000; %steam-pt ht coefficient
    
    Tave=(Tvap(n)+Tmod)/2;
    
    kgap=0.0476+(0.362e-3*Tave)-(0.618e-7*Tave^2)+(0.718e-11*Tave^3)*10^-3;
    
    htctpt(n)=kgap/(dman+djump); % ht-pt ht coefficient
    
    Qel(n)=Qbundle(n)/37;
    
    
end

%% Mass of cladding, pressure tube calandria tube
rhoCT=6595.2-(0.1477*(Tmod+273.15));
    
rhoPTclad=6595.2-(0.1477*(mean(Tvap)+273.15));

Mfst=pi()/4*((dfuel^2)-((dfuel-(2*tclad))^2))*Lchannel*rhoPTclad*37;

MPT=pi()/4*(((DPT+(2*tPT))^2)-(DPT^2))*Lchannel*rhoPTclad;

MCT=pi()/4*(((DCT+(2*tCT))^2)-(DCT^2))*Lchannel*rhoCT;

CPfst=(255.66+(0.1024*(mean(Tvap)+273.15)))/1000;

CPpt=CPfst;

CPct=(255.66+(0.1024*(Tmod+273.15)))/1000;

%plotyy(1:12,Tvap,1:12,Qtotal,'plot','plot');

Trun=500; %s

Tr=0:0.1:Trun;

Res=single_channel_Lockhart_Martinelli_Qloss_LPXSTEAM5(Qchannel,Tenter,Tmod,PVH/1000,Lbund);

TFs=zeros(length(Qbundle),length(Tr));

TPt=zeros(length(Qbundle),length(Tr));

TCt=zeros(length(Qbundle),length(Tr));

TFs(1:length(Tvap),1)=Tvap';

TPt(1:length(Tvap),1)=Tvap';

for ind=2:length(Tr)
    
    for bundind=1:length(Qbundle)
        
        Tfinf=((Qel(bundind)*1000)+(htfe(bundind)*pi()*dfuel*Lchannel*Tvap(bundind)))/(htfe(bundind)*Lchannel*pi()*dfuel);
        
        gammf=htfe(bundind)*pi()*dfuel*Lchannel/(pi()/4*(dfuel^2-((dfuel-(2*tclad))^2))*Lchannel*CPfst);
        
        C1=Tfinf-Tvap(bundind)-1;
        
        TFs(bundind,ind)=Tfinf-exp(-gammf*Tr(ind))-C1;
        
    end
end
        
plot(time,TFs);


