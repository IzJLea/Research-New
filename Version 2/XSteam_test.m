pressure=2;
Temperature=(0:1:2000);

% viscosity=zeros(length(pressure),length(Temperature));

Cp=zeros(length(pressure),length(Temperature));

tc=zeros(length(pressure),length(Temperature));

enthalpy=zeros(length(pressure),length(Temperature));

rhoV=zeros(length(pressure),length(Temperature));

rho=zeros(length(pressure),length(Temperature));


for n=1:length(pressure)
    
    for m=1:length(Temperature)
        
%         viscosity(n,m)=XSteam('my_pT',pressure(1,n),Temperature(1,m));
%         
%         if isnan(viscosity(n,m))
%             
%             viscosity(n,m)=0;
%         end
        
        Cp(n,m)=XSteam('Cp_pT',pressure(1,n),Temperature(1,m));
        
        if isnan(Cp(n,m))
            
            Cp(n,m)=-1;
        end
                
        enthalpy(n,m)=XSteam('h_pT',pressure(1,n),Temperature(1,m));
        
        if isnan(enthalpy(n,m))
            
            enthalpy(n,m)=-1;
        end
        
        rho(n,m)=XSteam('rho_pT',pressure(1,n),Temperature(1,m));
        
        if isnan(rho(n,m))
            
            rho(n,m)=-1;
        end
        
        tc(n,m)=XSteam('tc_pT',pressure(1,n),Temperature(1,m));
        
        if isnan(tc(n,m))
            
            tc(n,m)=-1;
        end
    end
end

% for p=1:length(pressure)
%     
%     rhoV(1,p)=XSteam('rhoV_p',pressure(1,p));
%     
%     if isnan(rhoV(1,p))
%         rhoV=0;
%     end
% end


subplot(2,3,1)
plot(Temperature,Cp)
title('Heat capacity at 2 bar')
subplot(2,3,2)
plot(Temperature,enthalpy)
title('enthalpy at 2 bar')
subplot(2,3,3)
plot(Temperature(130:2001),rho(130:2001))
title('density at 2 bar')
subplot(2,3,4)
plot(Temperature,tc)
title('thermal conductivity at 2 bar')
% subplot(2,3,5)
% plot(Temperature,viscosity)
% title('viscosity at 2 bar')

