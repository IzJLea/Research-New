function res=Treverse(mflow,Tflow,Tmod,Qchannel,Peval)

Lchannel=6; % 

Dh=0.0074; %m hydraulic diameter of bundle

Aflow=0.0035; %m^2 Flow area in channel

DPT=0.10338; %m inner pressure tube diameter

DCT=0.12869;%m inner calandria tube thickness

doutclad=0.0138; %m cladding outer diameter

tclad=0.00038; %m thickness of cladding

roc=doutclad/2; %m outer cladding radius

ric=roc-tclad; %m inner cladding radius

Lchannel=6; %m length of fuel channel

Lbund=Lchannel; %m length of bundle (for this initial simulation, entire channel will be simulated as one bundle)

tPT=0.00424; %m pressure tube thickness

tCT=0.0014; %m calandria tube thickness

Aipt=pi()*DPT*Lbund;%m^2 inner pressure tube area

Aopt=pi()*(DPT+(2*tPT))*Lbund;%m^2 outer pressure tube area

Aict=pi()*DCT*Lbund;%m^2 inner calandria tube area

Aoct=pi()*(DCT+(2*tCT))*Lbund;%m^2 outer calandria tube area

Afuel=pi()*doutclad*Lbund;%m^2 fuel outer area

ript=DPT/2; %inner pressure tube radius

ropt=ript+tPT; %outer pressure tube radius

rict=DCT/2; % inner calandria tube radius

roct=rict+tCT;% outer calandria tube radius

roc=doutclad/2; %m outer cladding radius

ric=roc-tclad; %m inner cladding radius

Qel=Qchannel/37; % generation in one element

Afuel=pi()*doutclad*Lchannel;%m^2 fuel outer area

Reynolds=mflow*Dh/Aflow/XSteam('my_pT',Peval,Tflow); % initial reynolds number based on liquid in channel

Prandtl=XSteam('Cp_pT',Peval,Tflow-0.1)*1000*XSteam('my_pT',Peval,Tflow-0.1)/XSteam('tc_pT', Peval,Tflow-0.1); % initial prandtl number 

if Reynolds<=3000 % initial nusselt number includes consideration for laminar flow
        Nusselt=4.36;
else
        
        Nusselt=0.023*Reynolds^(4/5)*Prandtl^0.4;
end

hcool=Nusselt*XSteam('tcL_p',Peval)/Dh; %J/m^2.K coolant heat transfer coefficient

Tclad=(Qel*1000000/Afuel/hcool)+Tflow; %C initial cladding temperature 

kuo2=kUO2(Tflow);%W/m.K initial fuel heat capacity

kclad=kzirc(Tflow); %W/m.K 

kPT=kzirc(Tflow);

kCT=kzirc(Tmod);

kCO2sys=kCO2((Tflow+Tmod)/2);

Rgap=0;

%calculation of resistances. See accompanying document for derivations

R1=(1/8/pi()/kuo2/Lchannel)+Rgap+(log(roc/ric)/4/pi()/kclad/Lchannel);

R2=(log(roc/ric)/4/pi()/kclad/Lchannel)+(1/hcool/Afuel);
    
R3=(1/hcool/Aipt)+(log(ropt/ript)/4/pi()/kPT/Lchannel);
    
R4=(log(ropt/ript)/4/pi()/kPT/Lchannel)+(log(rict/ropt)/2/pi()/kCO2sys/Lchannel)+(log(roct/rict)/4/pi()/kCT/Lchannel);
    
Tclad=(Qchannel*1000000*R2)+Tflow;

Tfuel=(Qchannel*1000000*R1)+Tclad;

TPT=(Qchannel*1000000*R3)+Tflow;

TCT=(Qchannel*1000000*R4)+TPT;

res=[Tfuel;Tclad;Tflow;TPT;TCT];