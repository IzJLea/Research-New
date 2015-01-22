function res=Rcool(Tcool,Afuel,Dh,Mflux,Peval)

if Tcool==XSteam('Tsat_p',Peval)
    Tcool=Tcool+0.1;
end

Reynolds=Mflux*Dh/XSteam('my_pT',Peval,Tcool);

RRC=isnan(Reynolds);

if RRC==1
    
    Reynolds=Mflux*Dh/XSteam('my_pT',Peval,Tcool+0.1);
    
end

    Prandtl=XSteam('Cp_pT',Peval,Tcool)*1000*XSteam('my_pT',Peval,Tcool)/XSteam('tc_pT', Peval,Tcool);
    
    PRC=isnan(Prandtl);
    
    if PRC==1
        Prandtl=XSteam('Cp_pT',Peval,Tcool+0.1)*1000*XSteam('my_pT',Peval,Tcool+0.1)/XSteam('tc_pT', Peval,Tcool+0.1);
    end
    
    if Reynolds<=3000
        
        Nusselt=4.63;
        
    else
        
        Nusselt=0.023*Reynolds^(4/5)*Prandtl^0.4;
    end
    
    hcool=Nusselt*XSteam('tc_pT',Peval,Tcool)/Dh; %W/m^2/K
    
    HRC=isnan(hcool);
    
    if HRC==1
        
         
    hcool=Nusselt*XSteam('tc_pT',Peval,Tcool+0.1)/Dh; %W/m^2/K
    
    end
        
    
    res=1/(hcool*Afuel);
end
    
    
    