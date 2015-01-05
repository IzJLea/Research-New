T3=cell2mat(Temperatures{1,1});

T4=cell2mat(Temperatures{2,1});

T5=cell2mat(Temperatures{3,1});

T6=cell2mat(Temperatures{4,1});

T7=cell2mat(Temperatures{5,1});

T73=cell2mat(Temperatures{6,1});

plot(T73(6,1:length(T73)),T73(1:5,1:length(T73)))

%% matrix consisting of first row:power row 2 through 4 temperatures of fuel, cladding and vapor

Results=zeros(length(Power),4);

Results(1:length(Power),1)=Power;

Results(1,2:4)=T3(1:3,length(T3));

Results(2,2:4)=T4(1:3,length(T4));

Results(3,2:4)=T5(1:3,length(T5));

Results(4,2:4)=T6(1:3,length(T6));

Results(5,2:4)=T7(1:3,length(T7));

Results(6,2:4)=T73(1:3,length(T73));