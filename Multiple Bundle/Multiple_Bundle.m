%% Initial Variables
%

Qchannel=.150; %MW channel power

bunds=12; %number of bundles

H=0; %m height difference from inlet/outlet headers

hmod=1000; %W/m^2.K heat transfer coefficient 

Tenter=200; %C Entering coolant temperature

PSH=210; %kPa supply header pressure

PVH=190; %kPa vent header pressure

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

Peval=110; % system evaluation pressure

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

%% Time determination 
% This section sets the length of the simulation and time divisions an

Time=1000;  %seconds total run time

div=0.1; %s time step length

ind=Time/div; % index to determine time step number

time=zeros(1,(Time/div));

for p=2:length(time)
    
    time(1,p)=time(1,p-1)+div;
end

%% Property matrix Creation

% This section creates all property matrices for calculations
%Temperatures
Tclad=zeros(bunds,ind); %C cladding temperature

Tvap=zeros(bunds,ind);%C vapor coolant temperature

TPT=zeros(bunds,ind);%C pressure tube temperature

Tfuel=zeros(bunds,ind);%C average fuel temperature

TCT=zeros(bunds,ind);%C calandria tube temperature

%Enthalpies

hin=zeros(bunds,ind);

hout=zeros(bunds,ind);% channel exit enthalpy

houtvap=zeros(bunds,ind);

hinvap=zeros(bunds,ind);

Mfluxvap=zeros(bunds,ind);

%Calculation Constants

A1=zeros(bunds,ind);
    
B1=zeros(bunds,ind);
    
F1=zeros(bunds,ind);
    
A2=zeros(bunds,ind);
    
B2=zeros(bunds,ind);
    
C2=zeros(bunds,ind);
    
D2=zeros(bunds,ind);
       
B3=zeros(bunds,ind);
    
C3=zeros(bunds,ind);
    
D3=zeros(bunds,ind);

F3=zeros(bunds,ind);

B4=zeros(bunds,ind);

C4=zeros(bunds,ind);

D4=zeros(bunds,ind);

E4=zeros(bunds,ind);

D5=zeros(bunds,ind);

E5=zeros(bunds,ind);

F5=zeros(bunds,ind);

%Miscellaneous

Qel=zeros(1,bunds);

Alphas=zeros(bunds,ind);
%% Power profile calculation

Bundpower=cosPower(bunds,Qchannel)*1000000;

for n=1:bunds
    Qel(1,n)=Bundpower(1,n)/37;
end
%% vapor fraction calculation
%This section creates a loop to solve for the vapor fraction (alpha) in the
%channel under low flow conditions

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
clear delta
%constant mass flux calculation

Mflux=(1-alpha)*Qchannel*1000/(XSteam('hV_p',Peval)-XSteam('h_pT',Peval,Tenter))/(alpha*Aflow); % average channel mass flux

mflow=Mflux*Aflow*alpha;

%% Initial property calculations
%in this section the initial temperature for the elements and the coolant
%are determined. The temperatures were determined using the assumption that
%saturated fluid was previously flowing in the channel. Mass flux is
%assumed to be constant at channel average value.


for m=1:bunds
        
    Tvap(m,1)=XSteam('Tsat_p',Peval);
    
    Tclad(m,1)=(Qel(1,m)*R2(Tvap(m,1),Tvap(m,1),Tvap(m,1),Afuel,Arad,Aipt,Dh,Mflux,Peval,ric,roc,Lbund))+Tvap(m,1);
 
    Tfuel(m,1)=(Qel(1,m)*R1(Tclad(m,1),Tclad(m,1),ric,roc,Lbund))+Tclad(m,1);
    
    TPT(m,1)=Tvap(m,1)-(Qloss(Tvap(m,1),Tvap(m,1),Tmod,Tmod,ript,ropt,rict,roct,hmod,Dh,Lbund,Mflux,Peval)*R3(Tvap(m,1),Tvap(m,1),Aipt,Dh,ript,ropt,Lbund,Mflux,Peval));
    
    TCT(m,1)=(Qloss(Tvap(m,1),Tvap(m,1),Tmod,Tmod,ript,ropt,rict,roct,hmod,Dh,Lbund,Mflux,Peval)*R5(Tmod,rict,roct,hmod,Lbund))+Tmod;       
end

%% initial alpha values. Liquid inventory/enthalpy is main variable

hin(1,1:ind)=XSteam('h_pT',Peval,Tenter)*1000;

hout(1,1)=hin(1,1)+(Bundpower(1,1)/mflow);

if hout(1,1)>=XSteam('hL_p',Peval)*1000
    
    
    hout(1,1)=hin(1,1)+(Bundpower(1,1)*Alphas(1,1)/mflow);
    
    Mfluxvap(1,1)=(1-Alphas(1,1))*Bundpower(1,1)/((XSteam('hV_p',Peval)-XSteam('hL_p',Peval))*1000*Alphas(1,1)*Aflow);
    
    houtvap(1,1)=XSteam('hV_p',Peval)*1000;
    
    if Alphas(1,1)~=0
        
        hinvap(1,1)=XSteam('hV_p',Peval)*1000;
    end
    
end

for p=2:bunds
    
    hin(p,1)=hout(p-1,1);
    
    hinvap(p,1)=houtvap(p-1,1);
    
    hout(p,1)=hin(p,1)+(Bundpower(1,p)/mflow);
    
    if hout(p,1)>=XSteam('hL_p',Peval)*1000
        
        Alphas(p,1)=((Bundpower(1,p)/mflow)+hin(p,1)-(1000*XSteam('hL_p',Peval)))/((1000*(XSteam('hV_p',Peval)-XSteam('hL_p',Peval)))+(Bundpower(1,p)/mflow));
    
        hout(p,1)=1000*((Alphas(p,1)*XSteam('hV_p',Peval))+((1-Alphas(p,1))*XSteam('hL_p',Peval)));
        
        houtvap(p,1)=XSteam('hV_p',Peval)*1000;
        
        if Alphas(p,1)~=0
        
            hinvap(p,1)=XSteam('hV_p',Peval)*1000;
        end
        
        Mfluxvap(p,1)=(1-Alphas(p,1))*Bundpower(1,p)/((XSteam('hV_p',Peval)-XSteam('hL_p',Peval))*1000*Alphas(p,1)*Aflow)+Mfluxvap(p-1,1);
    end
end
        
        


%%

for n=2:ind
    for p=1:bunds        
            if p>1
                
                hin(p,n)=hout(p-1,n);
                
                if Alphas(p,n-1)~=0
                    hinvap(p,n)=XSteam('hV_p',Peval)*1000;
                else
                    hinvap(p,n)=houtvap(p-1,n);
                end
                
                
                
            end
        
            
            hout(p,n)=hin(p,n)+(Bundpower(1,p)/mflow);
            
            Tvap(p,n)=XSteam('T_ph',Peval,(hout(p,n)/1000));
            
            Tclad(p,n)=(Qel(1,p)*R2(Tclad(p,n-1),Tvap(p,n-1),TPT(p,n-1),Afuel,Arad,Aipt,Dh,Mflux,Peval,ric,roc,Lbund))+Tvap(p,n);
            
            Tfuel(p,n)=(Qel(1,p)*R1(Tfuel(p,n-1),Tclad(p,n-1),ric,roc,Lbund))+Tclad(p,n);
            
            TPT(p,n)=Tvap(p,n)-(Qloss(Tvap(p,n-1),TPT(p,n-1),TCT(p,n-1),Tmod,ript,ropt,rict,roct,hmod,Dh,Lbund,Mflux,Peval)*R3(Tvap(p,n-1),TPT(p,n-1),Aipt,Dh,ript,ropt,Lbund,Mflux,Peval));
            
            TCT(p,n)=(Qloss(Tvap(p,n-1),TPT(p,n-1),TCT(p,n-1),Tmod,ript,ropt,rict,roct,hmod,Dh,Lbund,Mflux,Peval)*R5(TCT(p,n-1),rict,roct,hmod,Lbund))+Tmod;
                        
        if hout(p,n)>=XSteam('hL_p',Peval)*1000
            
            if Tvap(p,n-1)==XSteam('Tsat_p',Peval)
                Tvap(p,n-1)=Tvap(p,n-1)+0.1;
            end
            
            %% Temperature calculation coefficients
            % Fuel temperature
            A1(p,n)=-inv(mfuel*CpUO2(Tfuel(p,n-1))*R1(Tfuel(p,n-1),Tclad(p,n-1),ric,roc,Lbund));
            
            B1(p,n)=Tclad(p,n-1)/(mfuel*CpUO2(Tfuel(p,n-1))*R1(Tfuel(p,n-1),Tclad(p,n-1),ric,roc,Lbund));
            
            F1(p,n)=Qel(1,p)/(mfuel*CpUO2(Tfuel(p,n-1)));
            
            % Cladding Temperature
            
            A2(p,n)=Tfuel(p,n-1)/(mclad*CpZirc(Tclad(p,n-1))*R1(Tfuel(p,n-1),Tclad(p,n-1),ric,roc,Lbund));
            
            B2(p,n)=-1/(mclad*CpZirc(Tclad(p,n-1)))*(inv(R1(Tfuel(p,n-1),Tclad(p,n-1),ric,roc,Lbund))+inv(R2(Tclad(p,n-1),Tvap(p,n-1),TPT(p,n-1),Afuel,Arad,Aipt,Dh,Mfluxvap(p,n-1),Peval,ric,roc,Lbund)));
            
            C2(p,n)=Tvap(p,n-1)/(mclad*CpZirc(Tclad(p,n-1))*R2(Tclad(p,n-1),Tvap(p,n-1),TPT(p,n-1),Afuel,Arad,Aipt,Dh,Mfluxvap(p,n-1),Peval,ric,roc,Lbund));
            
            % coolant temperature
            
            B3(p,n)=Tclad(p,n-1)/(Mcool(Tvap(p,n-1),Peval,Aflow,Alphas(p,n-1),Lbund)*XSteam('Cp_pT',Peval,Tvap(p,n-1))*1000*R2conv(Tclad(p,n-1),Tvap(p,n-1),ric,roc,Afuel,Dh,Mfluxvap(p,n-1),Peval,Lbund));
            
            C3(p,n)=-((1/R2(Tclad(p,n-1),Tvap(p,n-1),TPT(p,n-1),Afuel,Arad,Aipt,Dh,Mfluxvap(p,n-1),Peval,ric,roc,Lbund))+(1/R3(Tvap(p,n-1),TPT(p,n-1),Aipt,Dh,ript,ropt,Lbund,Mfluxvap(p,n-1),Peval)))/(Mcool(Tvap(p,n-1),Peval,Aflow,Alphas(p,n-1),Lbund)*XSteam('Cp_pT',Peval,Tvap(p,n-1))*1000);
            
            D3(p,n)=TPT(p,n-1)/(Mcool(Tvap(p,n-1),Peval,Aflow,Alphas(p,n-1),Lbund)*XSteam('Cp_pT',Peval,Tvap(p,n-1))*1000*R3(Tvap(p,n-1),TPT(p,n-1),Aipt,Dh,ript,ropt,Lbund,Mfluxvap(p,n-1),Peval));
            
            F3(p,n)=(hinvap(p,n-1)-houtvap(p,n-1))*Mfluxvap(p,n-1)*Aflow*Alphas(p,n-1)/(Mcool(Tvap(p,n-1),Peval,Aflow,Alphas(p,n-1),Lbund)*XSteam('Cp_pT',Peval,Tvap(p,n-1))*1000);
            
                
                
            % pressure tube temperature
            
            B4(p,n)=Tclad(p,n-1)/(mPT*CpZirc(TPT(p,n-1))*R2rad(Tclad(p,n-1),TPT(p,n-1),ric,roc,Arad,Aipt,Lbund));
            
            C4(p,n)=Tvap(p,n-1)/(mPT*CpZirc(TPT(p,n-1))*R3(Tvap(p,n-1),TPT(p,n-1),Aipt,Dh,ript,ropt,Lbund,Mfluxvap(p,n-1),Peval));
            
            D4(p,n)=-1/(mPT*CpZirc(TPT(p,n-1)))*(inv(R3(Tvap(p,n-1),TPT(p,n-1),Aipt,Dh,ript,ropt,Lbund,Mfluxvap(p,n-1),Peval))+inv(R2rad(Tclad(p,n-1),TPT(p,n-1),ric,roc,Arad,Aipt,Lbund))+inv(R4(TPT(p,n-1),TCT(p,n-1),ript,ropt,rict,roct,Lbund)));
            
            E4(p,n)=TCT(p,n-1)/(mPT*CpZirc(TPT(p,n-1))*R4(TPT(p,n-1),TCT(p,n-1),ript,ropt,rict,roct,Lbund));
            
            %calandria tube temperature
            
            D5(p,n)=TPT(p,n-1)/(mCT*CpZirc(TCT(p,n-1))*R4(TPT(p,n-1),TCT(p,n-1),ript,ropt,rict,roct,Lbund));
            
            E5(p,n)=-1/(mCT*CpZirc(TCT(p,n-1)))*(inv(R4(TPT(p,n-1),TCT(p,n-1),ript,ropt,rict,roct,Lbund))+inv(R5(TCT(p,n-1),rict,roct,hmod,Lbund)));
            
            F5(p,n)=Tmod/(mCT*CpZirc(TCT(p,n-1))*R5(TCT(p,n-1),rict,roct,hmod,Lbund));
            
            %% Temperature calculation
            
            Tfuel(p,n)=(Tfuel(p,n-1)*exp(A1(p,n)*div))+((1-exp(A1(p,n)*div))*((B1(p,n)+F1(p,n))/-A1(p,n)));
    
            Tclad(p,n)=(Tclad(p,n-1)*exp(B2(p,n)*div))+((1-exp(B2(p,n)*div))*((A2(p,n)+C2(p,n)+D2(p,n))/-B2(p,n)));
    
            Tvap(p,n)=(Tvap(p,n-1)*exp(C3(p,n)*div))+((1-exp(C3(p,n)*div))*((B3(p,n)+D3(p,n)+F3(p,n))/-C3(p,n)));
    
            TPT(p,n)=(TPT(p,n-1)*exp(D4(p,n)*div))+((1-exp(D4(p,n)*div))*((B4(p,n)+C4(p,n)+E4(p,n))/-D4(p,n)));
    
            TCT(p,n)=(TCT(p,n-1)*exp(E5(p,n)*div))+((1-exp(E5(p,n)*div))*((D5(p,n)+F5(p,n))/-E5(p,n)));
            
            if Tvap(p,n)<=XSteam('Tsat_p',Peval) && Alphas(p,n-1)~=0
                
                Tvap(p,n)=XSteam('Tsat_p',Peval)+0.1;
            end
                
            %% Main property calculations
            
            C=(XSteam('hL_p',Peval)*1000)-(Bundpower(1,p)/mflow);
            
            Alphas(p,n)=(hin(p,n-1)-C)/((XSteam('h_pT',Peval,Tvap(p,n-1))*1000)-C);
            
            if Tvap(p,n-1)==XSteam('Tsat_p',Peval)
                Alphas(p,n)=(hin(p,n-1)-C)/((XSteam('h_pT',Peval,Tvap(p,n-1)+0.1)*1000)-C);
            end
            
            if p==1
                Mfluxvap(p,n)=(1-Alphas(p,n-1))*Bundpower(1,p)/((XSteam('hV_p',Peval)-XSteam('hL_p',Peval))*1000*Alphas(p,n-1)*Aflow);
            else
                Mfluxvap(p,n)=(1-Alphas(p,n-1))*Bundpower(1,p)/((XSteam('hV_p',Peval)-XSteam('hL_p',Peval))*1000*Alphas(p,n-1)*Aflow)+Mfluxvap(p-1,n-1);
            end
                        
            houtvap(p,n)=XSteam('h_pT',Peval,Tvap(p,n-1))*1000;
            
            hout(p,n)=hin(p,n)+((1-Alphas(p,n-1))*Bundpower(1,p)/(Mfluxvap(p,n-1)*Aflow*Alphas(p,n-1)));
                        
        end
    end
end
            
            
