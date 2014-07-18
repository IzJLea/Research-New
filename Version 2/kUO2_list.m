Trange=100:1:4000;


kUO2list=zeros(1,length(Trange));

for n=1:length(Trange)
    
    kUO2list(n)=kUO2(Trange(n));
end

plot(Trange,kUO2list);