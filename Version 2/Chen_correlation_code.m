xtt=((1-x(hind))/x(hind))^0.9*(rhovsys/rhofsys)^0.5*(mulsys/muvsys)^0.1;
        
        if (1/xtt)<=0.1
            
            F=1;
        else
            F=2.35*((1/xtt)+0.213)^0.736;
            
        end
        
        S=1/(1+(2.53e-6*ReynoldsLO(hind)^1.17));
        
        
        
        syms Tclad
        
        hLO=(0.023*ReynoldsLO(hind)^0.8*Prlsys^0.4*klsys/Dh*F);
        
        Tclado(hind)=solve((4*(hLO+(1.22*(klsys^0.79*((1000*Cplsys)^0.45)*rhofsys^0.49/sigmasys^0.5/mulsys^0.29/(1000*hfvsys)^0.24/rhovsys^0.24)*(Tclad-Tbulk(hind))^0.24*(1000*(Pin-Pout))^0.75*S))/pi()/doutclad^2*(Tclad-Tbulk(hind)))-(Qchannel(hind)*1000000));
        
        clear Tclad
        
             
        Nu=0.023*((Mchannel(Rind)*rhofsys*Dh/mulsys)^0.8)*(Prlsys^0.4);
        
        hfluid=Nu*klsys/Dh;
        
        rhoave=1/((x(hind)/rhovsys)+((1-x(hind))/rhofsys));
        
        hfluid=hfluid*sqrt(rhofsys/rhoave);
        
        Tclado(hind)=(Q(hind)/(hfluid*doutclad^2*Lchannel*pi()))+Tbulk(hind);