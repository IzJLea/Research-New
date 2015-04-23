function res=QZirc_Steam(moxide,mclad,Tclad,t,div,Afuel,alpha)

if Tclad<1227 % initiation temperature in celsius for rapid zirc oxidation
    
    res=[0;0]; % no major reaction will take place
elseif moxide==(mclad*alpha) % no reaction will take place if all the mass available has been reacted
    
    res=[0;0];
    
else
    
    Tclad=Tclad+273;
    
    moxidenew=sqrt(33.3e6*t*exp(-45500/8.3145/Tclad))*Afuel*alpha*10;
    
    if moxidenew>mclad*alpha
        
        moxidenew=mclad*alpha;
        
    end
        
    dm=moxidenew-moxide;
    
    hrxn=140; %kcal/mol
    
    hrxn=hrxn*4184; %J/mol
    
    Qrxn=dm*hrxn/91.224/div/alpha;
    
    res=[dm;Qrxn];
end
    
    
    