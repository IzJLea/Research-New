%% Steam generator natural flow transient calculation

PSHstart=114;

PVH=134;

Hsysinit=13;

SGfrac=0.5;

Tvap=200;

time=100;

div=0.1;

Hsys=zeros(1,time/div+1);

Hsys(1,1)=Hsysinit;

PSH=zeros(1,time/div+1);

PSH(1,1)=PSHstart;

convcheck=zeros(1,time/div+1);

mflow=zeros(1,time/div+1);

t=linspace(1,time,time/div+1);

SGres=Steam_Generator(Hsys(1,1),PVH,PSH(1,1),SGfrac,Tvap,div);

convcheck(1,1)=SGres(4,1);

mflow(1,1)=SGres(2,1);

for i=2:length(Hsys)
    
    Hsys(1,i)=SGres(3,1);
    
    SGres=Steam_Generator(Hsys(1,i),PVH,PSH(1,i-1),SGfrac,Tvap,div);
    
    PSH(1,i)=SGres(1,1);
    
    mflow(1,i)=SGres(2,1);
    
    convcheck(1,i)=SGres(4,1);
    
end
    
    





