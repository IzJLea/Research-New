function res=Rrad(Tclad,TPT,Arad,Aipt)

eclad=0.6988;


hrad=sigma*((Tclad+273.15)+(TPT+273.15))*(((Tclad+273.15)^2)+((TPT+273.15)^2))/((1/eclad)+((1-eclad)/eclad*Arad/Aipt));

res=1/hrad/Arad;