function m=Mcool(Tvap,Peval,Aflow,alpha,Lbund)

if Tvap==XSteam('Tsat_p',Peval)
    
    Tvap=Tvap+0.1;
    
end

m=XSteam('rho_pT',Peval,Tvap)*Aflow*Lbund*alpha;

end