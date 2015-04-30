function res=Channel_Temperature_profile(Qchannel,Hchannel,Lchannel,Tenter,PSH,PVH,time,div,sections)

Lsection=Lchannel/sections;

Qsection=Qchannel/sections;

Lsections=zeros(sections,1);

Lsections(1,1)=Lsection;

for sec=2:sections
    Lsections(sec,1)=Lsections(sec-1,1)+Lsection;
end

Tmod=60;

hmod=1000; %W/m^2.K heat transfer coefficient 

doutclad=0.0138; %m cladding outer diameter

tclad=0.00038; %m thickness of cladding

DPT=0.10338; %m inner pressure tube diameter

DCT=0.12869;%m inner calandria tube thickness

tPT=0.00424; %m pressure tube thickness

tCT=0.0014; %m calandria tube thickness

dfuel=0.0122; %m fuel outer diameter

%% mass of zirc in elements
Tref=310; % C reference temperature

rhoref=rhoZirc(Tref); %kg/m^3 reference density

mclad=((doutclad)^2-(doutclad-(2*tclad))^2)/4*Lchannel*rhoref; %mass of zirconium in cladding of one element

mPT=((DPT+(2*tPT))^2-(DPT)^2)/4*Lchannel*rhoref; %kg mass of zirconium in pressure tube

mCT=((DCT+(2*tCT))^2-(DPT)^2)/4*Lchannel*rhoref; %kg mass of zirconium in calandria tube

%% mass of fuel

rhofuel=rhoUO2(Tref); %kg/m^3 reference fuel density

mfuel=Lchannel*pi()/4*dfuel^2*rhofuel; %kg mass of fuel within one fuel element

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

RCHf=reff(3); %effective resistance of fuel channel

RFf=sum(reff(4:5)); %effective resistance of end fitting and feeder

%% Creation of property matrices

alpha=zeros(sections,(time/div)+1);

hin=zeros(sections,(time/div)+1);

mflow=zeros(sections,(time/div)+1);

moxide=zeros(sections,(time/div)+1);

Peval=zeros(sections,(time/div)+1);

Qzirconium=zeros(sections,(time/div)+1);

Tvap=zeros(sections,(time/div)+1);

Tfuel=zeros(sections,(time/div)+1);

Tclad=zeros(sections,(time/div)+1);

TPT=zeros(sections,(time/div)+1);

TCT=zeros(sections,(time/div)+1);

t=zeros(sections,(time/div)+1);

toxide=zeros(sections,(time/div)+1);

for i=1:sections
    
    Qchannel=Qsection*i;
    
    Lchannel=Lsection*i;
    
    RCH=(RCHf/sections)*i;
    
    RF=RFf+(RCHf-RCH);
    
    mfuel=(mfuel/sections)*i;
    
    mclad=(mclad/sections)*i;
    
    mPT=(mPT/sections)*i;
    
    mCT=(mCT/sections)*i;
    
    for j=2:(time/div)+1
    
        t(i,j)=t(1,j-1)+div;
    end
    
    Peval(i,1)=(PSH+(XSteam('rhoL_T',Tenter)*9.81*Hchannel/1000))/100;


    %% Initial flow calculation

    Tvap(i,1)=XSteam('Tsat_p',Peval(i,1));

    Init=Alphacalc(PSH,PVH,Tenter,Tvap(i,1),Qchannel,RCH,RF,Hchannel);

    mflow(i,1)=Init(2,1);

    %% Initial Temperature calculation

    Init_Temp=Initial_Temp(mflow(i,1),Tvap(i,1),Tmod,Qchannel,Peval(i,1),Lchannel);

    Tfuel(i,1)=Init_Temp(1,1);

    Tclad(i,1)=Init_Temp(2,1);

    TPT(i,1)=Init_Temp(3,1);

    TCT(i,1)=Init_Temp(4,1);
    
    hin(i,1)=XSteam('hV_p',Peval(i,1));
    
    for j=2:(time/div)+1
    
        A=Alphacalc(PSH,PVH,Tenter,Tvap(i,j-1),Qchannel,RCH,RF,Hchannel);
    
        alpha(i,j)=A(1,1);
        
        if isnan(alpha(i,j))
        
        alpha(i,j)=max(alpha(i,1:length(alpha)));
        end
    
        mflow(i,j)=A(2,1);
        
        if isnan(mflow(i,j))
        
            mflow(i,j)=mflow(i,j-1);
        end
    
        Peval(1,j)=A(3,1);
        
        if isnan(Peval(1,j))
        
            Peval(i,j)=Peval(i,j-1);
        end
        
        
        hin=XSteam('hV_p',Peval(1,j));   
            
        B=Channel_Time_Step(Tfuel(i,j-1),Tclad(i,j-1),Tvap(i,j-1),TPT(i,j-1),TCT(i,j-1),Tmod,div,Peval(i,j),alpha(i,j),mflow(i,j),Qchannel,hin,hmod,mfuel,mclad,mPT,mCT,moxide(i,j-1),toxide(i,j-1),Lchannel);
    
        moxide(i,j)=moxide(1,j-1)+B(6,1);
        
        if moxide(i,j)~=0
            toxide(i,j)=toxide(i,j-1)+div;
        end
    
        Tfuel(i,j)=B(1,1);
    
        if isnan(Tfuel(i,j))
            Tf=Tfuel(i,1:j);
            Tfuel(i,j)=max(Tf);
        end
    
        Tclad(i,j)=B(2,1);
    
        if isnan(Tclad(i,j))
            Tc=Tclad(i,1:j);
            Tclad(i,j)=max(Tc);
        end
    
        Tvap(i,j)=B(3,1);
    
        if isnan(Tvap(1,i))
            Tv=Tvap(i,1:j);
            Tvap(i,j)=max(Tv);
        end
    
        TCT(i,j)=B(4,1);
    
        if isnan(TCT(i,j))
            TC=TCT(i,1:j);
            TCT(i,j)=max(TC);
        end
    
        TPT(i,j)=B(5,1);
    
        if isnan(TPT(i,j))
            TP=TPT(i,1:j);
            TPT(i,j)=max(TP);
        end
          
        moxide(i,j)=moxide(i,j)+B(6,1);
    
        Qzirconium(i,j)=B(7,1);
    end
end

res=[Tfuel;Tclad;Tvap;TPT;TCT;alpha;mflow;moxide];


    
    
    
