channelpower=0.150;

Peval=110;

Tenter=300;

bunds=12;

Mflux=140.994;

Aflow=0.0035;

Power=cosPower(bunds,channelpower*1000);

hin=XSteam('h_pT',Peval,Tenter);

hevap=XSteam('hL_p',Peval);

for m=1:bunds
    hout=hin+Power(1,m)/(Mflux*Aflow);
    
        if hout>=hevap
            
            break
        else
            hin=hout;
        end
end

eqs=bunds-(m-1);

Functions=zeros(2*eqs,1);

X=sym('x',[eqs,1]);

Y=sym('y',[eqs,1]);

for p=2:2:2*eqs
    
    if p==2
        
       
        
        Functions(p-1,1)=1-((eval('x.' num2str(p/2)))*(eval('y.'num2str(p/2))))-((1-eval('x.'num2str(p/2)))*(1-eval('y.'num2str(p/2)));
        
        Functions(p,1)=((1-eval('x.'num2str(p/2)))*Power(1,p/2+(m-1))/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval)))-(Mflow*Aflow*eval('x.'num2str(p/2))*eval('y.'num2str(p/2));
    elseif p==(2*eqs)
        
        Functions(p-1,1)=(eval('x.'num2str(p/2-1))*eval('y.'num2str(p/2-1)))+((1-eval('x.'num2str(p/2-1)))*(1-eval('y.'num2str(p/2-1))))-eval('x.'num2str(p/2));
        
        Functions(p,1)=Mflow*Aflow*((eval('x.'num2str(p/2-1))*eval('y.'num2str(p/2-1)))-eval('x.'num2str(p/2)))+((1-eval('x.'num2str(p/2)))*Power(1,p/2+(m-1))/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval)));
        
    else
        
        Functions(p-1,1)=(eval('x.'num2str(p/2-1))*eval('y.'num2str(p/2-1)))+((1-eval('x.'num2str(p/2-1)))*(1-eval('y.'num2str(p/2-1))))-eval('x.'num2str(p/2))-(eval('x.'num2str(p/2))*eval('y.'num2str(p/2))-((1-eval('x.'num2str(p/2)))*(1-eval('y.'num2str(p/2)));
        
        Functions(p,1)=Mflow*Aflow*((eval('x.'num2str(p/2-1))*eval('y.'num2str(p/2-1)))-(eval('x.'num2str(p/2))*eval('y.'num2str(p/2)))+((1-eval('x.'num2str(p/2)))*Power(1,p/2+(m-1))/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval)));
    end
end



%Solution=solve(Functions,X,Y);