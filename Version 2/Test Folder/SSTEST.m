Ttotal=500; %s
dt=1;
time=zeros(1,(Ttotal/dt)+1);
for f=2:length(time);
    time(f)=time(f-1)+dt;
end
Tvap=zeros(1,length(time));


V=0.0035*0.3; %m^3    
Peval=2; %bar
damp=0.001;
T1=120;
hin=XSteam('hV_p',Peval);
m=8.037484870817750e-04; %kg/s
dt=1;

Qin=0.002555563028320; %MW
Tvap(1)=T1;
for t=2:length(time)
Tvap(t)=Tvap(t-1);
err=1;
while err>=0.1
if Tvap(t)==Tvap(t-1)
    hout=XSteam('hV_p',Peval);
    rho=XSteam('rhoV_p',Peval);
    Cp=XSteam('CpV_p',Peval);
    
else
    Tprop=((Tvap(t-1)+Tvap(t))/2);
    hout=XSteam('h_pT',Peval,Tprop);
    rho=XSteam('rho_pT',Peval,Tprop);
    Cp=XSteam('Cp_pT',Peval,Tprop);
end



Qst=V*rho*Cp*(Tvap(t)-Tvap(t-1))/dt;
Qflow=m*(hin-hout)+(Qin*1000);
Qdiff=Qflow-Qst;
err=abs(Qdiff);
Tch=Qdiff/V/rho/Cp;
Tvap(t)=Tvap(t)+(damp*Tch);
end

end


