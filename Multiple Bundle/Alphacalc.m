function res=Alphacalc(PSH,PVH,Tenter,Tout,Qchannel,RCH,RF,Hchannel)


P1=(((PSH+PVH)/2)+((XSteam('rhoL_T',Tenter))*9.81*Hchannel/1000))/100;

Peval=(((PSH+PVH)/2)+((XSteam('rhoL_T',Tenter)-XSteam('rho_pT',P1,Tout))*9.81*Hchannel/1000))/100; % system evaluation pressure

hin=XSteam('hL_T',Tenter);

rhovap=XSteam('rho_pT',Peval,Tout);

rhoin=XSteam('rho_ph',Peval,hin+1);

deltaP=(PSH-PVH)+((rhoin-rhovap)*9.81*Hchannel/1000);

rhoave=(rhovap+rhoin)/2;

a=rhoave/rhovap;

x=(XSteam('hL_p',Peval)-hin)/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval));

wsmx=Qchannel*1000/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval));

for i=1:100
    
    alpha=1/(1+((1+x)/wsmx*sqrt(deltaP*rhoave/(RCH+(a*RF)))));
    
    a=alpha^2*rhoave/rhovap;
end

mflow=wsmx*(1+x)/(1-alpha);

rhoout=0;

houtrev=0;

Toutrev=0;

mflowrev=0;

if deltaP<0
    
    deltaP=(PVH-PSH);
    
    Peval=(PVH+PSH)/200;
    
    rhoave=rhovap;
    
    rhoout=rhovap;
    
    for i=1:100
        
        mflow=sqrt(deltaP*rhoave/(RCH+(RF*rhoave/rhoout)));
        
        houtrev=XSteam('h_pT',Peval,Tout)+(Qchannel*1000/mflow);
        
        rhoout=XSteam('rho_ph',Peval,houtrev);
        
        if isnan(rhoout)
            
            for j=1:10000
            
                houtrev=houtrev-100;
                
                rhoout=XSteam('rho_ph',Peval,houtrev);
                
                check=isnan(rhoout);
                
                if check==0
                    break
                end
            end
        end
                
    
           
        
        deltaP=((PVH-PSH)+((rhovap-rhoout)*9.81*Hchannel/1000));
        
        Peval=(((PSH+PVH)/2)+((rhovap-rhoout)*9.81*Hchannel/1000))/100;
        
        rhoave=(rhoout+rhovap)/2;

    end
        
    alpha=1;
    
    mflow=-mflow;
    
    Toutrev=XSteam('T_ph',Peval,houtrev);
    
    mflowrev=-sqrt(deltaP*rhoave/(RCH+(2*RF)));
    
    
end
    
       
    

res=[alpha;mflow;Peval;deltaP;rhoout;houtrev;Toutrev;mflowrev];
