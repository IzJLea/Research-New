Ttotal=500; %s
dt=1;
time=zeros(1,(Ttotal/dt)+1);
for f=2:length(time);
    time(f)=time(f-1)+dt;
end
Tvap=zeros(1,length(time));

alpha=0.3;
V=0.0035*alpha; %m^3    
Peval=2; %bar
damp=0.001;
T1=120;
hin=XSteam('hV_p',Peval);
m=8.037484870817750e-04; %kg/s


Qin=0.002555563028320; %MW
Tvap(1)=T1;
for t=2:length(time)
    syms T2
    
   Tvap(t)=solve( m*(XSteam('h_pT',Peval,((Tvap(t-1)+T2)/2)+0.1)-hin)+(Qin*1000*alpha)-(V*rho*Cp/dt*(T2-Tvap(t-1))),T2,Real);
           
   
end



