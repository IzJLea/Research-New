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

Lbund=Lchannel;

tPT=0.00424;

tCT=0.0014;

Aipt=pi()*DPT*Lbund;

Aopt=pi()*(DPT+(2*tPT))*Lbund;

Aict=pi()*DCT*Lbund;

Aoct=pi()*(DCT+(2*tCT))*Lbund;

Afuel=pi()*doutclad*Lbund*37;

Aoutfuel=9*pi()*doutclad; % Area of outer fuel elements that sees PT

doxide=3.00e-6; % thickness of oxide layer

sigma=5.670373e-8; %Stefan-Boltzmann constant

eclad=0.325+(0.1246e6*doxide);
        
ePT=eclad;
        
eCT=ePT;

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

mCT=((DCT+(2*tCT))^2-(DPT)^2)/4*Lbund*rhoref;

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

mflow=Qbundle*1000*(1-alpha)/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval)); % kg/s 

Time=500;  %seconds

div=0.5;

ind=Time/div;

Tclad=zeros(1,ind);

Tvap=zeros(1,ind);

TPT=zeros(1,ind);

Reynolds=zeros(1,ind);

Prandtl=zeros(1,ind);

Nusselt=zeros(1,ind);

hcool=zeros(1,ind);

hrad=zeros(1,ind);

Cpclad=zeros(1,ind);

CpPT=zeros(1,ind);

Cpvap=zeros(1,ind);

hout=zeros(1,ind);

A1=zeros(1,ind);
    
B1=zeros(1,ind);
    
C1=zeros(1,ind);
    
D1=zeros(1,ind);
    
A2=zeros(1,ind);
    
B2=zeros(1,ind);
    
C2=zeros(1,ind);
    
D2=zeros(1,ind);
    
A3=zeros(1,ind);
    
B3=zeros(1,ind);
    
C3=zeros(1,ind);
    
D3=zeros(1,ind);
%% initial Temperature

Tvap(1,1)=XSteam('Tsat_p',Peval);

Reynolds(1,1)=mflow*Dh/Aflow/XSteam('my_pT',Peval,Tvap(1,1)-0.1);

Prandtl(1,1)=XSteam('Cp_pT',Peval,Tvap(1,1)-0.1)*1000*XSteam('my_pT',Peval,Tvap(1,1)-0.1)/XSteam('tc_pT', Peval,Tvap(1,1)-0.1);

if Reynolds(1,1)<=3000
        Nusselt(1,1)=4.36;
else
        
        Nusselt(1,1)=0.023*Reynolds(1,1)^(4/5)*Prandtl(1,1)^0.4;
end

hcool(1,1)=Nusselt(1,1)*XSteam('tcL_p',Peval)/Dh; 

Tclad(1,1)=(Qbundle*1000000/Afuel/hcool(1,1))+Tvap(1,1);

TPT(1,1)=Tvap(1,1); % This is only temporary while Qloss is equal to zero. Once method works this needs to be changed.

hin=XSteam('hV_p',Peval)*1000; % J/kg

for n=2:ind
    
    if Tvap(1,n-1)==XSteam('Tsat_p',Peval)
        
       Reynolds(1,n)=mflow*Dh/Aflow/alpha/XSteam('my_pT',Peval,Tvap(1,1)+0.1);
       
       
       Prandtl(1,n)=XSteam('Cp_pT',Peval,Tvap(1,1)+0.1)*1000*XSteam('my_pT',Peval,Tvap(1,1)+0.1)/XSteam('tc_pT', Peval,Tvap(1,1)+0.1);
    else
        Reynolds(1,n)=mflow*Dh/Aflow/alpha/XSteam('my_pT',Peval,Tvap(1,n-1));
        
        Prandtl(1,n)=XSteam('Cp_pT',Peval,Tvap(1,n-1))*1000*XSteam('my_pT',Peval,Tvap(1,n-1))/XSteam('tc_pT', Peval,Tvap(1,n-1));
    end
    
    if Reynolds(1,n)<=3000
        
        Nusselt(1,n)=4.36;
        
    else
        
        Nusselt(1,n)=0.023*(Reynolds(1,n)^(4/5))*(Prandtl(1,n)^0.4);
    end
    
    hcool(1,n)=Nusselt(1,n)*XSteam('tcL_p',Peval)/Dh; 
    
    hrad(1,n)=sigma*(Tclad(1,n-1)+TPT(1,n-1))*((Tclad(1,n-1)^2)+(TPT(1,n-1)^2))/((1/eclad)+((1-eclad)/eclad*Afuel/2/Aipt));
    
    Cpclad(1,n)=255.66+(0.1024*(Tclad(1,n-1)+273.15)); %J/kg.K
    
    CpPT(1,n)=255.66+(0.1024*(TPT(1,n-1)+273.15)); %J/kg.K
    
    if Tvap(1,n-1)>XSteam('Tsat_p',Peval)
    
        Cpvap(1,n)=XSteam('Cp_pT',Peval,Tvap(1,n-1))*1000;
    else
        Cpvap(1,n)=XSteam('CP_pT',Peval,Tvap(1,n-1)+0.1)*1000;
    end
    
    if Tvap(1,n-1)>XSteam('Tsat_p',Peval)
    
        hout(1,n)=XSteam('h_pT',Peval,Tvap(1,n-1))*1000; % J/kg
    else
        hout(1,n)=XSteam('h_pT',Peval,Tvap(1,n-1)+0.1)*1000;
    end
    
    A1(1,n)=-((hcool(1,n)*Afuel)+(hrad(1,n)*Afuel/2))/(mclad*Cpclad(1,n));
    
    B1(1,n)=hcool(1,n)*Afuel*Tvap(1,n-1)/mclad/Cpclad(1,n);
    
    C1(1,n)=((hrad(1,n)*Afuel/2*mclad/Cpclad(1,n))+(hcool(1,n)*Aipt/mclad/Cpclad(1,n)))*TPT(1,n-1);
    
    D1(1,n)=Qchannel*1000000/mclad/Cpclad(1,n);
    
    A2(1,n)=hcool(1,n)*Afuel*Tclad(1,n-1)/mflow/Cpvap(1,n);
    
    B2(1,n)=-hcool(1,n)*(Afuel+Aipt)/mflow/Cpvap(1,n);
    
    C2(1,n)=hcool(1,n)*Aipt*TPT(1,n)/mflow/Cpvap(1,n);
    
    D2(1,n)=(hout(1,n)-hin)/Cpvap(1,n);
    
    A3(1,n)=hrad(1,n)*Afuel/2/mPT/CpPT(1,n);
    
    B3(1,n)=hcool(1,n)*Aipt*Tvap(1,n)/mPT/CpPT(1,n);
    
    C3(1,n)=-((hrad(1,n)*Afuel/2)+(hcool(1,n)*Aipt))/mPT/CpPT(1,n);
    
    D3(1,n)=0;
    
    dTclad=(Tclad(1,n-1)*exp(-A1(1,n)*div))+(-(B1(1,n)+C1(1,n)+D1(1,n))/A1(1,n)*(1-exp(-A1(1,n)*div)));
    
    dTvap=(Tvap(1,n-1)*exp(-B2(1,n)*div))+(-(A2(1,n)+C2(1,n)+D2(1,n))/B2(1,n)*(1-exp(-B2(1,n)*div)));
    
    dTPT=(TPT(1,n-1)*exp(-C3(1,n)*div))+(-(B3(1,n)+A3(1,n)+D3(1,n))/C3(1,n)*(1-exp(-C3(1,n)*div)));
    
    Tclad(1,n)=Tclad(1,n-1)+dTclad;
    
    Tvap(1,n)=Tvap(1,n-1)+dTvap;
    
    TPT(1,n)=TPT(1,n-1)+dTPT;
    
end


    
    