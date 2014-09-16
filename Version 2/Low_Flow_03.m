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

Lrun=5; %length of run (s)

dt=0.01; %time division

damp=0.01;
                    
time=zeros(1,Lrun/dt);

for tt=2:length(time)
    
    time(tt)=time(tt-1)+dt;
end

Tvap=zeros(length(Tvapi),length(time));

kvap=zeros(length(Tvapi),length(time));

Cpvap=zeros(length(Tvapi),length(time));

rhovap=zeros(length(Tvapi),length(time));

Tvap(1:12,1)=Tvapi;

kvap(1:12,1)=k_init;

Cpvap(1:12,1)=Cp_init;

rhovap(1:12,1)=rho_init;

for n=2:length(time)
    for p=1:length(Qbundle)
        if p==1
            V=alpha*Lbund*Aflow;
            hin=XSteam('hV_p',Peval);
            Tvap(p,n)=Tvap(p,n-1);
            err=1;
            while err>=0.1
                if Tvap(p,n)==Tvap(p,1)
                    hout=hin;
                    rho=XSteam('rhoV_p',Peval);
                    Cp=XSteam('CpV_p',Peval);
                else
                    Teval=(Tvap(p,n)+Tvap(p,n-1))/2;
                    hout=XSteam('h_pT',Peval,Teval);
                    rho=XSteam('rho_pT',Peval,Teval);
                    Cp=XSteam('Cp_pT',Peval,Teval);
                end
            
                Qst=V*rho*Cp*(Tvap(p,n)-Tvap(p,n-1))/dt;
            
                Qflow=mbundle(p)*(hin-hout)+(Qbundle(p)*1000);
            
                Qdiff=Qflow-Qst;
            
                err=abs(Qdiff);
            
                Tch=Qdiff/V/rho/Cp;
            
                Tvap(p,n)=Tvap(p,n)+(damp*Tch);
               
            end
             clear V
                clear hin
                clear err
                clear hout
                clear rho
                clear Cp
                clear Qst
                clear Qflow
                clear Qdiff
                clear Tch
        else
            
            V=alpha*Lbund*Aflow;
            hin=((XSteam('h_pT',Peval,Tvap(p-1,n-1))*min(p))+(mchange(p)*XSteam('hV_p',Peval)))/(mbundle(p));
            Tvap(p,n)=Tvap(p,n-1);
            err=1;
            while err>=0.1
                
                if Tvap(p,n)==Tvap(p,1)
                    hout=hin;
                    rho=XSteam('rhoV_p',Peval);
                    Cp=XSteam('CpV_p',Peval);
                else
                    Teval=(Tvap(p,n)+Tvap(p,n-1))/2;
                    hout=XSteam('h_pT',Peval,Teval);
                    rho=XSteam('rho_pT',Peval,Teval);
                    Cp=XSteam('Cp_pT',Peval,Teval);
                    
                end
                
                Qst=V*rho*Cp*(Tvap(p,n)-Tvap(p,n-1))/dt;
            
                Qflow=mbundle(p)*(hin-hout)+(Qbundle(p)*1000);
            
                Qdiff=Qflow-Qst;
            
                err=abs(Qdiff);
            
                Tch=Qdiff/V/rho/Cp;
            
                Tvap(p,n)=Tvap(p,n)+(damp*Tch);
                
            end
            
                clear V
                clear hin
                clear err
                clear hout
                clear rho
                clear Cp
                clear Qst
                clear Qflow
                clear Qdiff
                clear Tch
            
        end
    end
end
plot(time,Tvap(1,1:length(time)));

                

            

        
        
           