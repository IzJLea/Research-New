function res=Rfuel(Tfuel,length)

% average resistance for a cylinder of fuel. Uses average fuel temperature
% assumed to be at midpoint of cylinder

res=1/(4*pi()*kUO2(Tfuel)*length);

end
