function res=Channel_Time_Step(Tfuel,Tclad,Tvap,TPT,TCT,Tmod,div,Peval,alpha,mflow,Qch,hin,hmod,mfuel,mclad,mPT,mCT,moxide,toxide,Lchannel)
% requirement for t should be checked. Replace with time of oxidation for Qzirc calculation. 
%% Channel Properties
if mflow<0
    
    mflow=-mflow;
end

Dh=0.0074; % hydraulic diameter of the channel

Aflow=0.0035; %m^2 Flow area in channel

DPT=0.10338; %m inner pressure tube diameter

DCT=0.12869;%m inner calandria tube thickness

doutclad=0.0138; %m cladding outer diameter

tclad=0.00038; %m thickness of cladding

roc=doutclad/2; %m outer cladding radius

ric=roc-tclad; %m inner cladding radius

Lbund=Lchannel; %m length of bundle (for this initial simulation, entire channel will be simulated as one bundle)

tPT=0.00424; %m pressure tube thickness

tCT=0.0014; %m calandria tube thickness

Aipt=pi()*DPT*Lbund;%m^2 inner pressure tube area

Aopt=pi()*(DPT+(2*tPT))*Lbund;%m^2 outer pressure tube area

Aict=pi()*DCT*Lbund;%m^2 inner calandria tube area

Aoct=pi()*(DCT+(2*tCT))*Lbund;%m^2 outer calandria tube area

Afuel=pi()*doutclad*Lbund;%m^2 fuel outer area

Arad=9/37*Afuel; % Area of outer fuel elements that sees PT

doxide=3.00e-6; % thickness of oxide layer (assumed to be low for initial simulation)

sigma=5.670373e-8; %Stefan-Boltzmann constant

eclad=0.325+(0.1246e6*doxide); % dimensionless emissivity value for zirconium cladding
        
ript=DPT/2; %inner pressure tube radius

ropt=ript+tPT; %outer pressure tube radius

rict=DCT/2; % inner calandria tube radius

roct=rict+tCT;% outer calandria tube radius

Qel=Qch/37; % Power in MW
%% Temperature calculations

if Tvap==XSteam('Tsat_p',Peval)
        
       Reynolds=mflow*Dh/Aflow/alpha/XSteam('my_pT',Peval,Tvap+0.1);
       
       Prandtl=XSteam('Cp_pT',Peval,Tvap+0.1)*1000*XSteam('my_pT',Peval,Tvap+0.1)/XSteam('tc_pT', Peval,Tvap+0.1);
else     
       
            
       Reynolds=mflow*Dh/Aflow/alpha/muvap(Tvap,Peval);
       
       Prandtl=XSteam('Cp_pT',Peval,Tvap)*1000*muvap(Tvap,Peval)/XSteam('tc_pT', Peval,Tvap);
            
         
end

if Reynolds<=3000
        
        Nusselt=4.36;
        
    else
        
        Nusselt=0.023*(Reynolds^(4/5))*(Prandtl^0.4);
end

if Tvap>XSteam('Tsat_p',Peval)
    
        hcool=Nusselt*XSteam('tc_pT',Peval,Tvap)/Dh;
    else
        hcool=Nusselt*XSteam('tc_pT',Peval,Tvap+0.1)/Dh;
end

% heat transfer coefficient for radiation- view area of fuel pins is
    % assumed to be equal to the outer half of the 18 outer fuel pins
hrad=sigma*((Tclad+273.15)+(TPT+273.15))*(((Tclad+273.15)^2)+((TPT+273.15)^2))/((1/eclad)+((1-eclad)/eclad*Arad/Aipt));
    
hrpt=sigma*((TPT+273.15)+(TCT+273.15))*(((TPT+273.15)^2)+((TCT+273.15)^2))/((1/eclad)+((1-eclad)/eclad*Aopt/Aict));

% heat capacity calculations for main elements
Cpclad=CpZirc(Tclad);
    
CpPT=CpZirc(TPT);
    
CpCT=CpZirc(TCT);

Cpfuel=CpUO2(Tfuel);

if Tvap>XSteam('Tsat_p',Peval)
    
    Cpvap=XSteam('Cp_pT',Peval,Tvap)*1000;
else
    Cpvap=XSteam('CP_pT',Peval,Tvap+0.1)*1000;
end

if Tvap>XSteam('Tsat_p',Peval)
    
    hout=XSteam('h_pT',Peval,Tvap)*1000; % J/kg
else
    hout=XSteam('h_pT',Peval,Tvap+0.1)*1000;
end

kuo2=kUO2(Tfuel); %W/mK

kclad=kzirc(Tclad);

kPT=kzirc(TPT);

kCT=kzirc(TCT);

kCO2sys=kCO2((Tvap+Tmod)/2);

Mcool=mflow*div;

Rgap=0;

%calculation of resistances. See accompanying document for derivations

R1=(1/8/pi()/kuo2/Lchannel)+Rgap+(log(roc/ric)/4/pi()/kclad/Lchannel);

R2=(log(roc/ric)/4/pi()/kclad/Lchannel)+(1/hcool/Afuel);
    
R3=(1/hcool/Aipt)+(log(ropt/ript)/4/pi()/kPT/Lchannel);
    
R4=(log(ropt/ript)/4/pi()/kPT/Lchannel)+(log(rict/ropt)/2/pi()/kCO2sys/Lchannel)+(log(roct/rict)/4/pi()/kCT/Lchannel);
    
R5=(log(roct/rict)/4/pi()/kCT/Lchannel)+(1/hmod/Aoct);

%Eq.1 coefficients
    
A1=-1/R1/mfuel/Cpfuel;
    
B1=Tclad/R1/mfuel/Cpfuel;
    
F1=Qel*1000000/mfuel/Cpfuel;
    
%Eq.2 coefficients
    
A2=Tfuel/mclad/Cpclad/R1;
    
B2=-1/mclad/Cpclad*((1/R1)+(1/R2)+(hrad*Arad));
    
C2=Tvap/mclad/Cpclad/R2;
    
D2=TPT*hrad*Arad/mclad/Cpclad;

F=QZirc_Steam(moxide,mclad,Tclad,toxide,div,Afuel,alpha);

dm=F(1,1);

F2=F(2,1)/Cpclad/mclad;

% dm=0;
% 
% F2=0;
    
%Eq.3 coefficients
    
B3=37*alpha*Tclad/Mcool/Cpvap/R2;
    
C3=-1/Mcool/Cpvap*((37*alpha/R2)+(1/R3));
    
D3=TPT/Mcool/Cpvap/R3;
    
F3=-mflow/Mcool*(hout-hin)/Cpvap;
    
%Eq.4 coefficients
    
B4=hrad*Arad*Tclad/mPT/CpPT;
    
C4=Tvap/mPT/CpPT/R3;
    
D4=-1/mPT/CpPT*((1/R3)+(1/R4)+(Arad*hrad)+(hrpt*Aopt));
    
E4=TCT/mPT/CpPT*((1/R4)+(hrpt*Aopt));
    
%Eq.5 coefficients
    
D5=TPT/mCT/CpCT*((1/R4)+(hrpt*Aopt));
    
E5=-1/mCT/CpCT*((1/R4)+(1/R5)+(hrpt*Aopt));
   
F5=Tmod/mCT/CpCT/R5;

% calculation of Temperature values
    
Tfuel=(Tfuel*exp(A1*div))+((1-exp(A1*div))*((B1+F1)/-A1));
    
Tclad=(Tclad*exp(B2*div))+((1-exp(B2*div))*((A2+C2+D2+F2)/-B2));
    
Tvap=(Tvap*exp(C3*div))+((1-exp(C3*div))*((B3+D3+F3)/-C3));
    
TPT=(TPT*exp(D4*div))+((1-exp(D4*div))*((B4+C4+E4)/-D4));
    
TCT=(TCT*exp(E5*div))+((1-exp(E5*div))*((D5+F5)/-E5));

res=[Tfuel;Tclad;Tvap;TPT;TCT;dm;F2*Cpclad*mclad];



