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

%% mass of fuel

rhofuel=10970/(1+(2.04e-5*Tref)+(8.7e-9*Tref^2));

mfuel=Lbund*pi()/4*dfuel^2*rhofuel*37;

%% calculate values for k for each section

%reference mass flowrate

M=25.8; % kg/s

Rho=780.6; %kg/m^3

keff=DP./(M^2);

reff=keff.*Rho;

RCH=reff(3);

RF=sum(reff(4:5));

%% Data for CO2 thermal conductivity

% CO2 thermal conductivity

    kCO2=[0.01051 0.01456 0.01858 0.02257 0.02652 0.03044 0.03814 0.04565 0.05293 0.08491 0.10688 0.11522];

% CO2 thermal conductivity Temperatures

    kCO2Temp=[-50 0 50 100 150 200 300 400 500 1000 1500 2000];

%% alpha calculation

alpha=0; %initial vapor fraction guess

x=(XSteam('h_pT',Peval,Tenter)-XSteam('hL_p',Peval))/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval)); %thermodynamic quality

wsmx=Qchannel*1000/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval)); %maximum steam mass flow

d=0.01; %dampening factor

delta=0.1; %difference in alpha value per iteration

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

div=0.1;

ind=Time/div;

Tclad=zeros(1,ind);

dTclad=zeros(1,ind);

Tvap=zeros(1,ind);

dTvap=zeros(1,ind);

TPT=zeros(1,ind);

dTPT=zeros(1,ind);

Tfuel=zeros(1,ind);

dTfuel=zeros(1,ind);

TCT=zeros(1,ind);

dTCT=zeros(1,ind);

Reynolds=zeros(1,ind);

Prandtl=zeros(1,ind);

Nusselt=zeros(1,ind);

hcool=zeros(1,ind);

hrad=zeros(1,ind);

kuo2=zeros(1,ind);

kclad=zeros(1,ind);

kPT=zeros(1,ind);

kCT=zeros(1,ind);

kCO2sys=zeros(1,ind);

Cpclad=zeros(1,ind);

CpPT=zeros(1,ind);

CpCT=zeros(1,ind);

Cpvap=zeros(1,ind);

Cpfuel=zeros(1,ind);

hout=zeros(1,ind);

A1=zeros(1,ind);
    
B1=zeros(1,ind);
    
F1=zeros(1,ind);
    
A2=zeros(1,ind);
    
B2=zeros(1,ind);
    
C2=zeros(1,ind);
    
D2=zeros(1,ind);
       
B3=zeros(1,ind);
    
C3=zeros(1,ind);
    
D3=zeros(1,ind);

F3=zeros(1,ind);

B4=zeros(1,ind);

C4=zeros(1,ind);

D4=zeros(1,ind);

E4=zeros(1,ind);

D5=zeros(1,ind);

E5=zeros(1,ind);

F5=zeros(1,ind);

Rgap=zeros(1,ind);

R1=zeros(1,ind);

R2=zeros(1,ind);

R3=zeros(1,ind);

R4=zeros(1,ind);

R5=zeros(1,ind);

CH=zeros(5,ind);

Mcool=zeros(1,ind);

%% initial Temperatures

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

kuo2(1,1)=kUO2(Tclad(1,1))*1000;

kclad(1,1)=12.767-(5.4348e-4*(Tclad(1,1)+273.15))+(8.9818e-6*(Tclad(1,1)+273.15)^2); %W/m.K

Rfuel(1,1)=1/(4*pi()*kuo2(1,1)*Lchannel);

Rclad(1,1)=log(roc/ric)/(2*pi()*kclad(1,1)*Lchannel);

Rgap(1,1)=0;

Rconv(1,1)=1/(Afuel/37*hcool(1,1));

R1(1,1)=(Rfuel(1,1)/2)+Rgap(1,1)+(Rclad(1,1)/2);

R2(1,1)=(Rclad(1,1)/2)+Rconv(1,1);

TPT(1,1)=Tvap(1,1); % This is only temporary while Qloss is equal to zero. Once method works this needs to be changed.

TCT(1,1)=Tmod;

Tfuel(1,1)=Tclad(1,1)+(Qchannel*1000000/37*R1(1,1));

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
    
    if Tvap(1,n-1)>XSteam('Tsat_p',Peval)
    
        hcool(1,n)=Nusselt(1,n)*XSteam('tc_pT',Peval,Tvap(1,n-1))/Dh;
    else
        hcool(1,n)=Nusselt(1,n)*XSteam('tc_pT',Peval,Tvap(1,n-1)+0.1)/Dh;
    end
    
     
    
    hrad(1,n)=sigma*((Tclad(1,n-1)+273.15)+(TPT(1,n-1)+273.15))*(((Tclad(1,n-1)+273.15)^2)+((TPT(1,n-1)+273.15)^2))/((1/eclad)+((1-eclad)/eclad*Afuel*9/37/Aipt));
    
    Cpclad(1,n)=255.66+(0.1024*(Tclad(1,n-1)+273.15)); %J/kg.K
    
    CpPT(1,n)=255.66+(0.1024*(TPT(1,n-1)+273.15)); %J/kg.K
    
    CpCT(1,n)=255.66+(0.1024*(TCT(1,n-1)+273.15));
    
    tau=(Tfuel(1,n-1)+273.15)/1000;
    
    Cpfuel(1,n)=(52.1743+(87.951*tau)-(85.2411*tau^2)+(31.542*tau^3)-(2.6334*tau^4)-(0.71391*tau^-2))/270.03*1000; %J/kg.K
    
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
    
    kuo2(1,n)=kUO2(Tfuel(1,n-1))*1000; % W/m.K
    
    kclad(1,n)=12.767-(5.4348e-4*(Tclad(1,n-1)+273.15))+(8.9818e-6*(Tclad(1,n-1)+273.15)^2); %W/m.K
    
    kPT(1,n)=12.767-(5.4348e-4*(TPT(1,n-1)+273.15))+(8.9818e-6*(TPT(1,n-1)+273.15)^2);
    
    kCT(1,n)=12.767-(5.4348e-4*(TCT(1,n-1)+273.15))+(8.9818e-6*(TCT(1,n-1)+273.15)^2);
    
    kCO2sys(1,n)=interp1(kCO2Temp,kCO2,(Tvap(1,n-1)+Tmod)/2);
    if Tvap(n-1)==XSteam('Tsat_p',Peval)
        
        Mcool(1,n)=XSteam('rho_pT',Peval,Tvap(1,n-1)+0.1)*Aflow*Lchannel;
    else
        
        Mcool(1,n)=XSteam('rho_pT',Peval,Tvap(1,n-1)+0.1)*Aflow*Lchannel;
    end
    
        
    
    R1(1,n)=(1/8/pi()/kuo2(1,n)/Lchannel)+Rgap(1,n)+(log(roc/ric)/4/pi()/kclad(1,n)/Lchannel);

    R2(1,n)=(log(roc/ric)/4/pi()/kclad(1,n)/Lchannel)+(1/hcool(1,n)/Afuel);
    
    R3(1,n)=(1/hcool(1,n)/Aipt)+(log(ropt/ript)/4/pi()/kPT(1,n)/Lchannel);
    
    R4(1,n)=(log(ropt/ript)/4/pi()/kPT(1,n)/Lchannel)+(log(rict/ropt)/2/pi()/kCO2sys(1,n)/Lchannel)+(log(roct/rict)/4/pi()/kCT(1,n)/Lchannel);
    
    R5(1,n)=(log(roct/rict)/4/pi()/kCT(1,n)/Lchannel)+(1/hmod/Aoct);
    
    A1(1,n)=-37/R1(1,n)/mfuel/Cpfuel(1,n);
    
    B1(1,n)=37*Tclad(1,n-1)/R1(1,n)/mfuel/Cpfuel(1,n);
    
    F1(1,n)=Qchannel*1000000/mfuel/Cpfuel(1,n);
    
    A2(1,n)=37*Tfuel(1,n-1)/mclad/Cpclad(1,n)/R1(1,n);
    
    B2(1,n)=-1/mclad/Cpclad(1,n)*((37/R1(1,n))+(37/R2(1,n))+(hrad(1,n)*Afuel/2));
    
    C2(1,n)=37*Tvap(1,n-1)/mclad/Cpclad(1,n)/R2(1,n);
    
    D2(1,n)=hrad(1,n)*Afuel/2/mclad/Cpclad(1,n)*TPT(1,n-1);
    
    B3(1,n)=37*Tclad(1,n-1)/Mcool(1,n)/Cpvap(1,n)/R2(1,n);
    
    C3(1,n)=-1/Mcool(1,n)/Cpvap(1,n)*((37/R2(1,n))+(1/R3(1,n)));
    
    D3(1,n)=TPT(1,n-1)/Mcool(1,n)/Cpvap(1,n)/R3(1,n);
    
    F3(1,n)=mflow/alpha/Mcool(1,n)*(hout(1,n)-hin)/Cpvap(1,n);
    
    B4(1,n)=hrad(1,n)*Afuel/2*Tclad(1,n-1)/mPT/CpPT(1,n);
    
    C4(1,n)=Tvap(1,n-1)/mPT/CpPT(1,n)/R3(1,n);
    
    D4(1,n)=-1/mPT/CpPT(1,n)*((1/R3(1,n))+(1/R4(1,n))+(Afuel/2*hrad(1,n)));
    
    E4(1,n)=TCT(1,n-1)/mPT/CpPT(1,n)/R4(1,n);
    
    D5(1,n)=TPT(1,n-1)/R4(1,n)/mCT/CpCT(1,n);
    
    E5(1,n)=-1/mCT/CpCT(1,n)*((1/R4(1,n))+(1/R5(1,n)));
    
    F5(1,n)=Tmod/mCT/CpCT(1,n)/R5(1,n);
    
    dTfuel(1,n)=(Tfuel(1,n-1)*exp(-A1(1,n)*div))+((1-exp(-A1(1,n)*div))*(-(B1(1,n)+F1(1,n))/A1(1,n)));
    
    ch1=isreal(dTfuel(1,n));
    
    if ch1==0
        display('dTfuel');
        break
    end
    
    dTclad(1,n)=(Tclad(1,n-1)*exp(-B2(1,n)*div))+((1-exp(-B2(1,n)*div))*(-(A2(1,n)+C2(1,n)+D2(1,n))/B2(1,n)));
    
    ch2=isreal(dTclad(1,n));
    
    if ch2==0
        display('dTclad');
        break
    end
    
    
    dTvap(1,n)=(Tvap(1,n-1)*exp(-C3(1,n)*div))+((1-exp(-C3(1,n)*div))*(-(B3(1,n)+D3(1,n)+F3(1,n))/C3(1,n)));
    
    
    ch3=isreal(dTvap(1,n));
    
    if ch3==0
        display('dTvap');
        break
    end
        
    dTPT(1,n)=(TPT(1,n-1)*exp(-D4(1,n)*div))+((1-exp(-D4(1,n)*div))*(-(B4(1,n)+C4(1,n)+E4(1,n))/D4(1,n)));
    
    ch4=isreal(dTPT(1,n));
    
    if ch4==0
        display('dTPT');
        break
    end
    
    dTCT(1,n)=(TCT(1,n-1)*exp(-E5(1,n)*div))+((1-exp(-E5(1,n)*div))*(-(D5(1,n)+F5(1,n))/E5(1,n)));
    
    ch5=isreal(dTCT(1,n));
    
    if ch5==0
        display('dTCT');
        break
    end
    
    Tfuel(1,n)=dTfuel(1,n);
    
    Tclad(1,n)=dTclad(1,n);
    
    Tvap(1,n)=dTvap(1,n);
    
    TPT(1,n)=dTPT(1,n);
    
    TCT(1,n)=dTCT(1,n);
    
    CH(1:5,n)=[Tfuel(1,n); Tclad(1,n); Tvap(1,n); TPT(1,n); TCT(1,n)];
    
end


    
    