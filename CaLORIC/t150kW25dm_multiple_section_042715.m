Qchannel=0.150; %MW power of channel

Hchannel=2.5; %m height of liquid in header

Lchannel=6; %length of channel in meters

Tenter=100; %Temperature of entering liquid

PSH=114; %supply header pressure

PVH=134; %Vent header Pressure

time=3000; %time of run in seconds

div=0.1; %time step

sections=12;

Results=Channel_Temperature_profile(Qchannel,Hchannel,Lchannel,Tenter,PSH,PVH,time,div,sections);
