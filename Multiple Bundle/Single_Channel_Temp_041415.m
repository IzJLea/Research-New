%% Initial channel conditions: liquid at Tenter is flowing in channel followed by initial voiding to saturation temp 
% steam. 

Lchannel=6; %m

Qchannel=0.050; %MW

Hchannel=1; %m

Tenter=100; %Celsius

Tmod=60;

hmod=1000; %W/m^2.K heat transfer coefficient 

time=3000; %s

div=0.1;

PSH=114; 

PVH=134;

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

RCH=reff(3); %effective resistance of fuel channel

RF=sum(reff(4:5)); %effective resistance of end fitting and feeder

%% Creation of property matrices

alpha=zeros(1,(time/div)+1);

mflow=zeros(1,(time/div)+1);

moxide=zeros(1,(time/div)+1);

Peval=zeros(1,(time/div)+1);

Qzirconium=zeros(1,(time/div)+1);

Tvap=zeros(1,(time/div)+1);

Tfuel=zeros(1,(time/div)+1);

Tclad=zeros(1,(time/div)+1);

TPT=zeros(1,(time/div)+1);

TCT=zeros(1,(time/div)+1);

t=zeros(1,(time/div)+1);



for j=2:(time/div)+1
    
    t(1,j)=t(1,j-1)+div;
    
end

%% Initial Evaluation Pressure

Peval(1,1)=(PSH+(XSteam('rhoL_T',Tenter)*9.81*Hchannel/1000))/100;


%% Initial flow calculation

Tvap(1,1)=XSteam('Tsat_p',Peval(1,1));

Init=Alphacalc(PSH,PVH,Tenter,Tvap(1,1),Qchannel,RCH,RF,Hchannel);

mflow(1,1)=Init(2,1);

%% Initial Temperature calculation

Init_Temp=Initial_Temp(mflow(1,1),Tvap(1,1),Tmod,Qchannel,Peval(1,1));

Tfuel(1,1)=Init_Temp(1,1);

Tclad(1,1)=Init_Temp(2,1);

TPT(1,1)=Init_Temp(3,1);

TCT(1,1)=Init_Temp(4,1);

%% loop for calculation temperatures and alpha values for transient

for i=2:(time/div)+1
    
    A=Alphacalc(PSH,PVH,Tenter,Tvap(1,i-1),Qchannel,RCH,RF,Hchannel);
    
    alpha(1,i)=A(1,1);
    if isnan(alpha(1,i))
        
        alpha(1,i)=max(alpha);
    end
    
    mflow(1,i)=A(2,1);
    if isnan(mflow(1,i))
        
        mflow(1,i)=mflow(1,i-1);
    end
    
    Peval(1,i)=A(3,1);
    if isnan(Peval(1,i))
        
        Peval(1,i)=Peval(1,i-1);
    end
    
    hin=XSteam('hV_p',Peval(1,i));   
            
    B=Channel_Time_Step(Tfuel(1,i-1),Tclad(1,i-1),Tvap(1,i-1),TPT(1,i-1),TCT(1,i-1),Tmod,div,Peval(1,i),alpha(1,i),mflow(1,i),Qchannel,hin,hmod,mfuel,mclad,mPT,mCT,moxide(1,i-1),t(1,i));
    
    Tfuel(1,i)=B(1,1);
    
    if isnan(Tfuel(1,i))
        
        Tfuel(1,i)=max(Tfuel);
    end
    
    Tclad(1,i)=B(2,1);
    
    if isnan(Tclad(1,i))
        
        Tclad(1,i)=max(Tclad);
    end
    
    Tvap(1,i)=B(3,1);
    
    if isnan(Tvap(1,i))
        
        Tvap(1,i)=max(Tvap);
    end
    
    TCT(1,i)=B(4,1);
    
    if isnan(TCT(1,i))
        
        TCT(1,i)=max(TCT);
    end
    
    TPT(1,i)=B(5,1);
    
    if isnan(TPT(1,i))
        
        TPT(1,i)=max(TPT);
    end
    
    
        
   
        
        
    
    moxide(1,i)=moxide(1,i-1)+B(6,1);
    
    Qzirconium(1,i)=B(7,1);
end

Temps=[Tfuel;Tclad;Tvap;TCT;TPT];

plot(t,Temps);

figure

plot(t,alpha);

figure 

plot(t,mflow);






