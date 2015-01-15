function res=RCO2(TPT,TCT,innerrad,outerrad,length)

Teval=(TPT+TCT)/2;

res=log(outerrad/innerrad)/(2*pi()*kCO2(Teval)*length);

end
