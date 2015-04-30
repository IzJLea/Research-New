function res=QZirc_Steam(moxide,mclad,Tclad,toxide,div,Afuel,alpha)
% t is time that the oxidation has been taking place

if Tclad<1227 % initiation temperature in celsius for rapid zirc oxidation
    
    res=[0;0]; % no major reaction will take place
elseif moxide==(mclad*alpha) % no reaction will take place if all the mass available has been reacted
    
    res=[0;0];
    
else
    
    Tclad=Tclad+273;
    
    moxidenew=sqrt(9.25e4*toxide*Afuel*alpha*exp(-45000/8.314/Tclad));
    
    if moxidenew>mclad*alpha
        
        moxidenew=mclad*alpha;
        
    end
        
    dm=moxidenew-moxide;
    
    hrxn=140; %kcal/mol
    
    hrxn=hrxn*4184; %J/mol
    
    Qrxn=dm*hrxn/91.224/div/alpha;
    
    res=[dm;Qrxn];
end
    
    
    