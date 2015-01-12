function power=cosPower(Sections,Power)

div=pi()/Sections;

count=1;

cosPo=zeros(1,Sections);
H=Power/2;

P1=-pi()/2;

for count=1:Sections
    
    P2=P1+div;
    
    cosPo(1,count)=H*(sin(P2)-sin(P1));
    
    count=count+1;
    
    P1=P2;
end

power=cosPo;

end

    