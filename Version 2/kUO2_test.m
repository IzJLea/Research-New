T=300:10:3000;

Tc=T+273.15;

k=zeros(1,length(T));

for ch=1:length(T)
    
    k(ch)=kUO2(Tc(ch));
end

plot(T,k);



