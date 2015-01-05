Time=(1:1:5);
Power=[3,4,5,6,7,7.3];

Temperatures=cell(length(Power),length(Time));


for j=1:length(Time)
    for i=1:length(Power)
        Temperatures{i,j}=SteamFlow(Time(1,j),Power(1,i));
    end
end




