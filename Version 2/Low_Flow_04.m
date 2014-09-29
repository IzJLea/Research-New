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

dman=ric-(dfuel/2); %m

djump=10e-6; %m

Lchannel=5.94; %m

Lbund=Lchannel/12;

tPT=0.00424;

tCT=0.0014;

Aipt=pi()*DPT*Lbund;

Aopt=pi()*(DPT+(2*tPT))*Lbund;

Aict=pi()*DCT*Lbund;

Aoct=pi()*(DCT+(2*tCT))*Lbund;

Afuel=pi()*doutclad*Lbund;

doxide=3.00e-6; % thickness of oxide layer

sigma=5.670373e-8; %Stefan-Boltzmann constant

%% delta P calculation

deltaP=(PSH)-PVH+((9.81*H*(XSteam('rhoL_p',Peval)-XSteam('rhoV_p',Peval)))/1000);

%% Resistance Calculation from liquid flow data

DP = [165;
262.5000;
520.5000;
258;
174];



kCO2=[14.60e-3 16.23e-3 17.87e-3 19.52e-3 21.18e-3 22.84e-3 27.00e-3 31.12e-3 35.20e-3 39.23e-3];

kCO2Temp=[0 20 40 60 80 100 150 200 250 300];

%% mass of zirc in elements
Tref=310; % Celsius

TrefK=Tref+273.15;

rhoref=255.66+(0.1024*TrefK);

mclad=37*((doutclad)^2-(doutclad-(2*tclad))^2)/4*Lbund*rhoref;

mPT=((DPT+(2*tPT))^2-(DPT)^2)/4*Lbund*rhoref;
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

Qel=Qbundle./37;

Qvol=Qel./(pi()/4*(dfuel^2)*Lbund);

for l=1:length(Qtotal)
    
    Qtotal(l)=sum(Qbundle(1:l));
    
end

%% Bundle Mass flow

mbundle=zeros(1,length(Qbundle));

Tvapi=zeros(1,12);

hvap=zeros(1,12);

mchange=zeros(1,12);

min=zeros(1,12);

vbund_guess=zeros(1,12);

Vol_bund=pi()/4*DPT^2*Lbund;

k_init=zeros(1,length(Qbundle));

Cp_init=zeros(1,length(Qbundle));

rho_init=zeros(1,length(Qbundle));

for m=1:length(Tvapi)
    
    hvap(m)=XSteam('hV_p',Peval);
    
    Tvapi(m)=XSteam('T_ph',Peval,hvap(m));
         
    mchange(m)=((1-alpha)*Qbundle(m)*1000/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval)));
    
    mbundle(m)=sum(mchange(1:m));
    
    min(m)=mbundle(m)-mchange(m);
    
    vbund_guess(m)=mbundle(m)/XSteam('rhoV_p',Peval)*pi()/4*DPT*alpha;
    
    k_init(m)=XSteam('tcV_T',Tvapi(m))/1000;
    
    Cp_init(m)=XSteam('CpV_T',Tvapi(m));
    
    rho_init(m)=XSteam('rhoV_T' ,Tvapi(m));
    
    
end

Lrun=100; %length of run (s)

dt=0.1; %time division

divf=10000;
                    
time=zeros(1,Lrun/dt);

for tt=2:length(time)
    
    time(tt)=time(tt-1)+dt;
end
%% Property Matrix Creation
% Temperature
Tvap=zeros(length(Tvapi),length(time)); % Mean Vapour Temperature

Tclad=zeros(length(Tvapi),length(time)); % Mean clad Temperature

TPT=zeros(length(Tvapi),length(time)); % Pressure tube temperature

TCT=zeros(length(Tvapi),length(time)); % Calandria Tube Temperature

Tvap(1:12,1)=Tvapi; % Initial Vapour temperature assignment

Tfuelo=zeros(length(Qbundle),length(time)); % Outer fuel Temperature

Tfuelmid=zeros(length(Qbundle),length(time)); % fuel centerline Temperature

Tfuelmeat=zeros(length(Qbundle),length(time)); % mean fuel temperature

% Thermal Conductivity

kvap=zeros(length(Tvapi),length(time)); % vapour thermal conductivity

kvap(1:12,1)=k_init; % initial vapour thermal conductivity assignment

% Heat Capacity

Cpvap=zeros(length(Tvapi),length(time)); % Vapour Heat capacity

Cpvap(1:12,1)=Cp_init; % initial vapour heat capacity assignment

% Density

rhovap=zeros(length(Tvapi),length(time)); % vapour density

rhovap(1:12,1)=rho_init; % vapour density assignment

% HT coefficient

hclad=zeros(length(Tvapi),length(time)); % heat transfer coefficient between fuel and vapour

hpt=zeros(length(Tvapi),length(time)); % heat transfer coefficient between PT and vapour

% Thermal Energy

Qremoved=zeros(length(Tvapi),length(time)); % Total Heat removed by vapour from fuel elements

Qremovedconv=zeros(length(Tvapi),length(time)); % heat removed from fuel elements by convection

Qremovedrad=zeros(length(Tvapi),length(time)); % heat removed from fuel elements by radiation

Qret=zeros(length(Tvapi),length(time)); % heat remaining in the fuel elements

QlossPT=zeros(length(Tvapi),length(time)); % heat lost from PT/CT to moderator

QconvPT=zeros(length(Tvapi),length(time)); % heat removed by convection from PT

QPT=zeros(length(Tvapi),length(time)); % heat remaining within the PT

% Nusselt Number

Nu=zeros(length(Tvapi),length(time)); % Nusselt number for fuel pins/vapour interaction

NuPT=zeros(length(Tvapi),length(time)); % Nusselt number for PT/vapour interaction

% Prandtl Number

Pr=zeros(length(Tvapi),length(time)); % Prandtl number for PT/vapour interaction

PrPT=zeros(length(Tvapi),length(time)); % Prandtl number for PT/vapour interaction

% Reynolds Number

Re=zeros(length(Tvapi),length(time)); % Reynolds number for PT/vapour interaction

RePT=zeros(length(Tvapi),length(time)); % Reynolds number for PT/vapour interaction

%% Initial Property calculations

for ind=1:length(Qbundle)
    
    Re(ind,1)=mbundle(ind)*Dh/XSteam('my_pT',Peval,Tvap(ind,1)-0.1);
    
    Pr(ind,1)=XSteam('Cp_pT',Peval,Tvap(ind,1)-0.1)*1000*XSteam('my_pT',Peval,Tvap(ind,1)-0.1)/XSteam('tc_pT', Peval,Tvap(ind,1)-0.1);
    
    Nu(ind,1)=0.023*Re(ind,1)^(4/5)*Pr(ind,1)^0.4;
    
    hclad(ind,1)=Nu(ind,1)*XSteam('tc_pT',Peval, Tvap(ind,1)-0.1)/Dh;
    
    Tclad(ind,1)=(Qbundle(ind)*1000/37/hclad(ind,1)/Afuel)+Tvap(ind,1);
    
    RePT(ind,1)=mbundle(ind)*DPT/XSteam('my_pT',Peval,Tvap(ind,1)-0.1);
    
    PrPT(ind,1)=XSteam('Cp_pT',Peval,Tvap(ind,1)-0.1)*1000*XSteam('my_pT',Peval,Tvap(ind,1)-0.1)/XSteam('tc_pT', Peval,Tvap(ind,1)-0.1);
    
    NuPT(ind,1)=0.023*RePT(ind,1)^(4/5)*PrPT(ind,1)^0.4;
    
    hpt(ind,1)=Nu(ind,1)*XSteam('tc_pT',Peval, Tvap(ind,1)-0.1)/DPT;
    
    Qloss=Qloss_single_channel2(Tvap(ind,1),Tmod,hpt(ind,1),hmod,Lbund);
    
    TPT(ind,1)=Tvap(ind,1)-((((1/Qloss(2))+(1/Qloss(3))+(1/Qloss(2)))/2)*Qbundle(ind));
    
    TCT(ind,1)=Tvap(ind,1)-((((1/Qloss(2))+(1/Qloss(3))+(1/Qloss(4))+(1/Qloss(5))+(1/Qloss(2))+(1/Qloss(3))+(1/Qloss(4)))/2)*Qbundle(ind));
    
    kgap=0.0476+(0.362e-3*Tclad(ind,1))-(0.618e-7*Tclad(ind,1)^2)+(0.718e-11*Tclad(ind,1)^3)*10^-3; %kW/m.C

    dman=ric-(dfuel/2); %m

    djump=10e-6; %m

    hgap=kgap/(dman+djump);

    Tfuelo(ind,1)=(Qel(ind)*1000/(dfuel/2*Lchannel)/hgap)+Tclad(ind,1);

%% Fuel pin centerline temperature

    
    rfuel=dfuel/2;

    reval=linspace(rfuel,0,divf);

    Tfuel=zeros(1,divf);

    Tfuel(1)=Tfuelo(ind,1);

    for pev=2:divf
    
        Tfuel(pev)=(Qvol(ind)*1000/4/kUO2(Tfuel(pev-1))*(reval(pev-1)^2-reval(pev)^2))+Tfuel(pev-1);
    
    end
    
    Tfuelmid(ind,1)=Tfuel(length(Tfuel));
    
    Tfuelmeat(ind,1)=mean(Tfuel);
    
    clear Tfuel
    
end

%% Transient calculations

for n=2:length(time)
    for p=1:length(Qbundle)
        
        %% HT coefficient for fuel pin cladding
        
        %% Start by computing the different characteristic numbers using data from the previous time step where properties are needed
        %% Calculate the values for heat leaving and retained within the fuel elements.
       
        %% Calculate the values for heat leaving/entering and retained within the PT
        %% calculate new vapour/cladding/PT/CT/fuel temperatures 
        %% Repeat
        
       
        
        