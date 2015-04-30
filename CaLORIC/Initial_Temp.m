function Res=Initial_Temp(mflow,Tflow,Tmod,Qchannel,Peval,Lchannel)

Dh=0.0074; %m hydraulic diameter of bundle

Aflow=0.0035; %m^2 Flow area in channel

doutclad=0.0138; %m cladding outer diameter

tclad=0.00038; %m thickness of cladding

roc=doutclad/2; %m outer cladding radius

ric=roc-tclad; %m inner cladding radius

Qel=Qchannel/37; % generation in one element

Afuel=pi()*doutclad*Lchannel;%m^2 fuel outer area

Reynolds=mflow*Dh/Aflow/XSteam('my_pT',Peval,Tflow-0.1); % initial reynolds number based on liquid in channel

Prandtl=XSteam('Cp_pT',Peval,Tflow-0.1)*1000*XSteam('my_pT',Peval,Tflow-0.1)/XSteam('tc_pT', Peval,Tflow-0.1); % initial prandtl number 

if Reynolds<=3000 % initial nusselt number includes consideration for laminar flow
        Nusselt=4.36;
else
        
        Nusselt=0.023*Reynolds^(4/5)*Prandtl^0.4;
end

hcool=Nusselt*XSteam('tcL_p',Peval)/Dh; %J/m^2.K coolant heat transfer coefficient

Tclad=(Qel*1000000/Afuel/hcool)+Tflow; %C initial cladding temperature 

kuo2=kUO2(Tclad);%W/m.K initial fuel heat capacity

kclad=kzirc(Tclad); %W/m.K 

Rfuel=1/(4*pi()*kuo2*Lchannel); % resistance of entire fuel meat

Rclad=log(roc/ric)/(2*pi()*kclad*Lchannel); %resistance of cladding

Rgap=0; 

R1=(Rfuel/2)+Rgap+(Rclad/2);% resistance between Tfuel and Tclad

TPT=Tflow; % simplified initial pressure tube temperature

TCT=Tmod; % simplified initial calandria tube temperature

Tfuel=Tclad+(Qel*1000000*R1); % initial average fuel temperature

Res=[Tfuel;Tclad;Tflow;TPT;TCT];



