function res=Rcool(Tcool,Afuel,Dh,Mflux,Peval)

Reynolds=Mflux*Dh/XSteam('my_pT',Peval,Tcool);
    
    Prandtl=XSteam('Cp_pT',Peval,Tcool)*1000*XSteam('my_pT',Peval,Tcool)/XSteam('tc_pT', Peval,Tcool);
    
    if Reynolds<=3000
        
        Nusselt=4.63;
        
    else
        
        Nusselt=0.023*Reynolds^(4/5)*Prandtl^0.4;
    end
    
    hcool=Nusselt*XSteam('tcL_T',Tcool)/Dh; %W/m^2/K
    
    res=1/(hcool*Afuel*37);
    
    
    