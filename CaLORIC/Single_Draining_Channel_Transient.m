%% Single_Draining_Channel_Transient.m
%This calculates the transient using the single draining channel function
%and plots the temepratures of the five main components


Qchannel=0.05;

Qchannel2=0.05;

Hchannel=11;

Hchannel2=5;

Lfeeder=10;

Tenter=100;

PSH=10000;

PVH=10100;

time=3000;

div=0.1;

Lchannel=6;

%%tic

Results=Single_Draining_Channel(Qchannel,Hchannel,Lfeeder,Tenter,PSH,PVH,time,div,Lchannel);

%%toc

Tfuel=Results(1,1:length(Results));

Tclad=Results(2,1:length(Results));

Tvap=Results(3,1:length(Results));

TPT=Results(4,1:length(Results));

TCT=Results(5,1:length(Results));

alpha=Results(7,1:length(Results));

mflow=Results(6,1:length(Results));

Hcool=Results(8,1:length(Results));

Lwet=Results(9,1:length(Results));

Peval=Results(10,1:length(Results));

t=linspace(0,time,(time/div)+1);



figure

plot(t,Tfuel,t,Tclad,t,Tvap,t,TPT,t,TCT);

figure

plot(t,mflow);



