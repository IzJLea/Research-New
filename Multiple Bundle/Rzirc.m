function res=Rzirc(Tclad,innerrad,outerrad,length)
% resistance for any cylinder of zirconium. For use in cladding, pt and ct


res=log(outerrad/innerrad)/(2*pi()*kzirc(Tclad)*length);

end