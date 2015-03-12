channelpower=0.150;

Peval=110;

Tenter=300;

bunds=12;

Mflux=27.7839;

Aflow=0.0035;

mflow=Mflux*Aflow;

alpha=0.11546;

alphaguess=0;

Lchannel=6; %meters

DPT=0.10338;

Wchar=Aflow/DPT;

vmax=mflow/XSteam('rhoL_p',Peval)/Aflow;

Power=linPower(bunds,channelpower*1000);

hin=XSteam('h_pT',Peval,Tenter);

hevap=XSteam('hL_p',Peval);

%% calculation of saturation length

qlin=channelpower*1000/Lchannel;

Lsat=mflow*(hevap-hin)/qlin;
    
Leval=zeros(1,bunds);

count=zeros(1,bunds);

H=zeros(1,bunds);

B=zeros(1,bunds);

alphas=zeros(1,bunds);

Leval(1,1)=0.5;

for j=2:bunds
    
    
    Leval(1,j)=Leval(1,j-1)+0.5;
    
end

Ho=((2*DPT*Lchannel*(1-alpha))-(DPT*(Lsat+Lchannel)))/(Lchannel-Lsat);

for k=1:bunds
    
    if Leval(1,k)>Lsat;
        
        alphas(1,k)=alphaguess;
        
        err=1;
        for l=1:5000
           
            
            B(1,k)=qlin*alpha/DPT/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval))/Wchar/vmax/XSteam('rhoL_p',Peval)/(alpha-alphas(1,k));
            
            H(1,k)=(Leval(1,k)*exp(-B(1,k)*Leval(1,k)))-(Lsat*exp(-B(1,k)*Lsat))+Ho;
            
            alphanew=(1-((((-Leval(1,k)/B(1,k))+(1/B(1,k)^2))*exp(-B(1,k)*Leval(1,k)))-(((-Lsat/B(1,k))+(1/B(1,k)^2))*exp(-B(1,k)*Lsat))-(Lsat*exp(-B(1,k)*Lsat)*(Leval(1,k)-Lsat))+(Ho*(Leval(1,k)-Lsat))))/DPT/Leval(1,k);
            
            error=alphanew-alphas(1,k);
            
            err=abs(error);
            
            alphas(1,k)=alphas(1,k)+(0.000001*error);
            
            if err<0.000001
                
                break
            end
            
            if isnan(alphas(1,k))
                break
            end
            count(1,k)=count(1,k)+1;
        end
    end
end

    
