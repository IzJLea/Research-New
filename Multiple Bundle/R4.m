function res=R4(TPT,TCT,ript,ropt,rict,roct,length)

res=0.5*Rzirc(TPT,ript,ropt,length)+RCO2(TPT,TCT,ropt,rict,length)+0.5*Rzirc(TCT,rict,roct,length);

end