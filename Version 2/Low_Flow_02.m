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

Lbund=Lchannel/12;

tPT=0.00424;

tCT=0.0014;

Aipt=pi()*DPT*Lbund;

Aopt=pi()*(DPT+(2*tPT))*Lbund;

Aict=pi()*DCT*Lbund;

Aoct=pi()*(DCT+(2*tCT))*Lbund;
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

%% Bundle Mass flow

mbundle=zeros(1,length(Qbundle));

Tvap=zeros(1,12);

hvap=zeros(1,12);

mchange=zeros(1,12);

min=zeros(1,12);

for m=1:length(Tvap)
    
    hvap(m)=(alpha*Qtotal(m)*1000/wv)+hg;
    
    Tvap(m)=XSteam('T_ph',Peval,hvap(m));
         
    mchange(m)=((1-alpha)*Qbundle(m)*1000/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval)));
    
    mbundle(m)=sum(mchange(1:m));
    
    min(m)=mbundle(m)-mchange(m);
    
end

Trun=500; %run time

divt=1; %time step

time=0:divt:Trun;

vbund=zeros(length(Qbundle),length(time));

Tfs=zeros(length(Qbundle),length(time));

Tvapt=zeros(length(Qbundle),length(time));

hvapt=zeros(length(Qbundle),length(time));

hsatv=XSteam('hV_p',Peval);

h1=zeros(length(Qbundle),length(time));

hevap=XSteam('hV_p',Peval)-XSteam('hL_p',Peval);

hin=zeros(length(Qbundle),length(time));

T1=zeros(length(Qbundle),length(time));

Tfinf=zeros(length(Qbundle),length(time));

Tpinf=zeros(length(Qbundle),length(time));

Tvapt(1:length(Qbundle),1)=XSteam('Tsat_p',Peval);

hin(1,1:length(time))=0;

hin(2:length(Qbundle),1)=XSteam('hV_p',Peval);

Tfc=zeros(length(Qbundle),length(time));

Reysc=zeros(length(Qbundle),length(time));

Nusc=zeros(length(Qbundle),length(time));

hsc=zeros(length(Qbundle),length(time));

Reyps=zeros(length(Qbundle),length(time));

Nups=zeros(length(Qbundle),length(time));

hcp=zeros(length(Qbundle),length(time));

hcc=zeros(length(Qbundle),length(time));

Tpt=zeros(length(Qbundle),length(time));

Tct=zeros(length(Qbundle),length(time));

Afuel=Lbund*pi()*doutclad;

Apt=Lbund*pi()*DPT;
Apto=Lbund*pi()*(DPT+(2*tPT));



for h1start=1:length(Qbundle)
    
    h1(h1start,1)=(min(1)/mbundle(1)*hin(1,1))+(mchange(1)/mbundle(1)*XSteam('hV_p',Peval));
    
    T1(h1start,1)=XSteam('T_ph',Peval,h1(h1start,1));
    
    Tvapt(h1start,1)=T1(h1start,1)+(Qbundle(h1start)*alpha*1000/mbundle(h1start)/XSteam('Cp_ph',Peval,h1(h1start,1)));
    
end

rhoclad=6595.2-(0.1477*(Tvapt(1,1)+273.15));

Mclad=(dfuel^2-(dfuel-(2*tclad))^2)*pi()/4*Lbund*rhoclad;


for tind=2:length(time)
    
    for bundind=1:length(Qbundle)
        
        hin(bundind,tind)=XSteam('h_pT',Peval,Tvapt(bundind,tind-1));
        
        h1(bundind,tind)=(min(bundind)/mbundle(bundind)*hin(bundind,tind))+(mchange(bundind)/mbundle(bundind)*hsatv);
        
        T1(bundind,tind)=XSteam('T_ph',Peval,h1(bundind,tind));
        
        Tvapt(bundind,tind)=T1(bundind,tind)+(Qbundle(bundind)*alpha*1000/mbundle(bundind)/XSteam('Cp_ph',Peval,h1(bundind,tind)));
    end
end

for tc=1:length(time)
    for bc=1:length(Qbundle)
        
        vbund(bc,tc)=mbundle(bc)/XSteam('rho_pT',Peval,Tvapt(bc,tc))/(Aflow*alpha);
        
        Reysc(bc,tc)=vbund(bc,tc)*Dh*XSteam('rho_pT',Peval,Tvapt(bc,tc))/XSteam('my_pT',Peval,Tvapt(bc,tc));
        
        Nusc(bc,tc)=0.023*(Reysc(bc,tc))^(4/5)*((XSteam('Cp_pT',Peval,Tvapt(bc,tc))*XSteam('my_pT',Peval,Tvapt(bc,tc))/XSteam('tc_pT',Peval,Tvapt(bc,tc)))^0.4);
        
        hsc(bc,tc)=Nusc(bc,tc)*XSteam('tc_pT',Peval,Tvapt(bc,tc))/Dh/1000;
        
        Tfinf(bc,tc)=((Qbundle(bc)/37)+(hsc(bc,tc)*Lbund*pi()*doutclad*Tvapt(bc,tc)))/(hsc(bc,tc)*Lbund*pi()*doutclad);
        
        lamf=(hsc(bc,tc)*Lbund*pi()*doutclad)/(Mclad*((255.66+(0.1024*(Tvapt(bc,tc)+273.15)))/1000));
        
        Reyps(bc,tc)=vbund(bc,tc)*DPT*XSteam('rho_pT',Peval,Tvapt(bc,tc))/XSteam('my_pT',Peval,Tvapt(bc,tc));
        
        Nups(bc,tc)=0.023*(Reyps(bc,tc))^(4/5)*((XSteam('Cp_pT',Peval,Tvapt(bc,tc))*XSteam('my_pT',Peval,Tvapt(bc,tc))/XSteam('tc_pT',Peval,Tvapt(bc,tc)))^0.4);
        
        hcp(bc,tc)=Nusc(bc,tc)*XSteam('tc_pT',Peval,Tvapt(bc,tc))/DPT/1000;
        
        
        
        %Tpinf(bc,tc)=()/()
        
    end
end

% for td=1:length(time)
%     for bd=1:length(Qbundle)
%         if td==1
%             Tfc(bd,td)=Tvapt(bd,td);
%         else
%             Tfc(bd,td)=Tfc(bd,td-1)+(lamf*(Tfinf(bd,length(time))-Tfc(bd,td-1))*divt);
%             
%             
%             
%         end
%         Reyps(bd,td)=vbund(bd,td)*DPT*XSteam('rho_pT',Peval,Tvapt(bd,td))/XSteam('my_pT',Peval,Tvapt(bd,td));
%         
%         Nups(bd,td)=0.023*(Reyps(bd,td))^(4/5)*((XSteam('Cp_pT',Peval,Tvapt(bd,td))*XSteam('my_pT',Peval,Tvapt(bd,td))/XSteam('tc_pT',Peval,Tvapt(bd,td)))^0.4);
%         
%         hcp(bd,td)=Nusc(bd,td)*XSteam('tc_pT',Peval,Tvapt(bd,td))/DPT/1000;
%         
%         if td==1
%             
%             ho=1000; %W/mK
%             Qloss=Qloss_single_channel2(Tvapt(bd,td),Tmod,hcp(bd,td)*1000,ho, Lbund);
%             kzrct=(7.51+(0.362e-3*Tmod)-(0.618e-7*Tmod^2)+(0.718e-11*Tmod^3))/1000;  %kW/m.K
%             
%             kco2=interp1(kCO2Temp,kCO2,(mean(Tvapt(bd,td),Tmod)))/1000;
%             kzrpt=(7.51+(0.362e-3*Tvapt(bd,td))-(0.618e-7*Tvapt(bd,td)^2)+(0.718e-11*Tvapt(bd,td)^3))/1000;
%             Tsurfpt=(Qloss*1000/hcp(bd,td)/Aipt)+Tvapt(bd,td);
%             Topt=(Qloss*1000*tPT/kzrpt/Aopt)+Tsurfpt;
%             Tict=(Qloss*1000*(DCT-(DPT+(2*tPT)))/kco2/Aict)+Topt;
%             Toct=(Qloss*1000*tCT/kzrct/Aoct)+Tict;
%             
%             Tpt(bd,td)=mean(Tsurfpt,Topt);
%             
%             Tct(bd,td)=mean(Tict,Toct);
%         else
%             
%         end
%         
%     end
% end

plot(time,Tvapt);

% hin(1,1:length(time))=0;
% 
% hin(2:length(Qbundle),1)=XSteam('hV_p',Peval);

% for tind=2:length(time)
%     
%     for bundind=1:length(Qbundle)
        
      


kCO2Temp=1;


