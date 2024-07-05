close all
%clear all
clc
%%%%%%%%%%%%%%%%%% Simulate MDL file
%open_system('TFBR_12_NEW')
%out = sim('TFBR_12_NEW');
%%%%%%%%%%%%%%%%%%%%%%%%%% Upper Bridge Leg 1
xref=-12;
yref=0;
mp=2;
hf = figure('units','normalized','outerposition',[0 0 1 1]);


%writeobj = VideoWriter('T12pulse','Motion JPEG AVI');
writeobj = VideoWriter('T12pulse_new','Uncompressed AVI');
open(writeobj);

[D1x1,D1x2,D1y1,D1y2,D1tx,D1ty,D2x1,D2x2,D2y1,D2y2,D2tx,D2ty]=draw_diode(xref,yref,mp);
D1 = plot(D1x1,D1y1,'k', D1x2,D1y2,'k');
text(D1tx,D1ty,'D1','FontSize',8);
set(D1,'LineWidth',1.5);
ul1mpx=D1x2(:,1);
ul1mpy=D1y2(:,1);
UD1ex=D1x1(:,10);
UD1ey=D1y1(:,10);
title('TWELVE PULSE RECTIFIER','Color','r');
subtitle('OPERATION AND WAVEFORMS','Color','r');
hold on


D4 = plot(D2x1,D2y1,'k', D2x2,D2y2,'k');
text(D2tx,D2ty,'D4','FontSize',8);
set(D4,'LineWidth',1.5);


UD4ex=D2x2(:,1);
UD4ey=D2y2(:,1);
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%% Upper Bridge Leg 2
xref=xref+2.0;
yref=yref+0;
mp=3;
[D1x1,D1x2,D1y1,D1y2,D1tx,D1ty,D2x1,D2x2,D2y1,D2y2,D2tx,D2ty]=draw_diode(xref,yref,mp);

D3 = plot(D1x1,D1y1,'k', D1x2,D1y2,'k');
text(D1tx-0.25,D1ty,'D3','FontSize',8);
set(D3,'LineWidth',1.5);
ul2mpx=D1x2(:,1);
ul2mpy=D1y2(:,1);

UD3ex=D1x1(:,10);
UD3ey=D1y1(:,10);
hold on

D6 = plot(D2x1,D2y1,'k', D2x2,D2y2,'k');
text(D2tx,D2ty,'D6','FontSize',8);
set(D6,'LineWidth',1.5);

UD6ex=D2x2(:,1);
UD6ey=D2y2(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%% Upper Bridge Leg 3
xref=xref+2.0;
yref=yref+0;
mp=4;
[D1x1,D1x2,D1y1,D1y2,D1tx,D1ty,D2x1,D2x2,D2y1,D2y2,D2tx,D2ty]=draw_diode(xref,yref,mp);

D5 = plot(D1x1,D1y1,'k', D1x2,D1y2,'k');
text(D1tx,D1ty,'D5','FontSize',8);
set(D5,'LineWidth',1.5);

ul3mpx=D1x2(:,1);
ul3mpy=D1y2(:,1);
UD5ex=D1x1(:,10);
UD5ey=D1y1(:,10);
hold on

D2 = plot(D2x1,D2y1,'k', D2x2,D2y2,'k');
text(D2tx,D2ty,'D2','FontSize',8);
set(D2,'LineWidth',1.5);

UD2ex=D2x2(:,1);
UD2ey=D2y2(:,1);
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Line Joining top and bottom  rails of Upper Bridge

UD1D3 = plot([UD1ex,UD3ex],[UD1ey,UD3ey],'k');
set(UD1D3,'LineWidth',1.5);

UD3D5 = plot([UD3ex,UD5ex],[UD3ey,UD5ey],'k');
set(UD3D5,'LineWidth',1.5);

UD4D6 = plot([UD4ex,UD6ex],[UD4ey,UD6ey],'k');
set(UD4D6,'LineWidth',1.5);

UD6D2 = plot([UD6ex,UD2ex],[UD6ey,UD2ey],'k');
set(UD6D2,'LineWidth',1.5);

%%%%%%%%%%%%Line Joing Midpoints of Upper Bridge & Transformer part
% sp=4

xref1=-18;
yref1=-10;

[lx1,ly1]=Indc(xref1,yref1);
ly1(:,1)=ul1mpy;
ul1mp=plot([ul1mpx,lx1(:,end)],[ul1mpy,ul1mpy],'k');
set(ul1mp,'LineWidth',1.5);

xref1=-16.5;
yref1=-10;

[lx2,ly2]=Indc(xref1,yref1);
ly2(:,1)=ul2mpy;
ul2mp=plot([ul2mpx,lx2(:,end)],[ul2mpy,ul2mpy],'k');
set(ul2mp,'LineWidth',1.5);


xref1=-15;
yref1=-10;
[lx3,ly3]=Indc(xref1,yref1);
ly3(:,1)=ul3mpy;
ul3mp=plot([ul3mpx,lx3(:,end)],[ul3mpy,ul3mpy],'k');
set(ul3mp,'LineWidth',1.5);

%%%%%%%%%%%%%%%%%%%%%%%%Transformer
% theta=330;
% %xx=(lx1.*cos(theta))+(ly1.*sin(theta));
% xx=(lx1.*cos(theta));
% %yy=(-lx1.*sin(theta))+(ly1.*cos(theta));
% yy=(-lx1.*sin(theta));
% tr1=plot(xx,yy,'k');


tr1=plot(lx1,ly1,'k');
set(tr1,'LineWidth',1.5);
lx1r=lx1(:,end);
ly1r=ly1(:,end);
text(-18.25,-5.5,'R','FontSize',8);
hold on


tr2=plot(lx2,ly2,'k');
set(tr2,'LineWidth',1.5);
lx2r=lx2(:,end);
ly2r=ly2(:,end);
text(-16.75,-5.5,'Y','FontSize',8);
hold on


tr3=plot(lx3,ly3,'k');
set(tr3,'LineWidth',1.5);
lx3r=lx3(:,end);
ly3r=ly3(:,end);
text(-15.25,-5.5,'B','FontSize',8);
hold on
l12= plot([lx1r lx2r],[ly1r ly2r],'k');
set(l12,'LineWidth',1.5);
hold on
l23= plot([lx2r lx3r],[ly2r ly3r],'k');
set(l23,'LineWidth',1.5);
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%% Lower Bridge Leg 1
xref=-12;
yref=-15;
mp=2;

[D1x1,D1x2,D1y1,D1y2,D1tx,D1ty,D2x1,D2x2,D2y1,D2y2,D2tx,D2ty]=draw_diode(xref,yref,mp);

D1L = plot(D1x1,D1y1,'k', D1x2,D1y2,'k');
text(D1tx,D1ty,'D1','FontSize',8);
set(D1L,'LineWidth',1.5);

bl1mpx=D1x2(:,1);
bl1mpy=D1y2(:,1);
BD1ex=D1x1(:,10);
BD1ey=D1y1(:,10);
hold on

D4L = plot(D2x1,D2y1,'k', D2x2,D2y2,'k');
text(D2tx,D2ty,'D4','FontSize',8);
set(D4L,'LineWidth',1.5);
BD4ex=D2x2(:,1);
BD4ey=D2y2(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%% Lower Bridge Leg 2
xref=xref+2.0;
yref=yref+0;
mp=3;
[D1x1,D1x2,D1y1,D1y2,D1tx,D1ty,D2x1,D2x2,D2y1,D2y2,D2tx,D2ty]=draw_diode(xref,yref,mp);

D3L = plot(D1x1,D1y1,'k', D1x2,D1y2,'k');
text(D1tx,D1ty,'D3','FontSize',8);
set(D3L,'LineWidth',1.5);

bl2mpx=D1x2(:,1);
bl2mpy=D1y2(:,1);
BD3ex=D1x1(:,10);
BD3ey=D1y1(:,10);
hold on

D6L = plot(D2x1,D2y1,'k', D2x2,D2y2,'k');
text(D2tx,D2ty,'D6','FontSize',8);
set(D6L,'LineWidth',1.5);
BD6ex=D2x2(:,1);
BD6ey=D2y2(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%% Lower Bridge Leg 3
xref=xref+2.0;
yref=yref+0;
mp=4;
[D1x1,D1x2,D1y1,D1y2,D1tx,D1ty,D2x1,D2x2,D2y1,D2y2,D2tx,D2ty]=draw_diode(xref,yref,mp);

D5L = plot(D1x1,D1y1,'k', D1x2,D1y2,'k');
text(D1tx,D1ty,'D5','FontSize',8);
set(D5L,'LineWidth',1.5);

bl3mpx=D1x2(:,1);
bl3mpy=D1y2(:,1);
BD5ex=D1x1(:,10);
BD5ey=D1y1(:,10);
hold on

D2L = plot(D2x1,D2y1,'k', D2x2,D2y2,'k');
text(D2tx,D2ty,'D2','FontSize',8);
set(D2L,'LineWidth',1.5);

BD2ex=D2x2(:,1);
BD2ey=D2y2(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Line Joining top and bottom rails of Lower Bridge

BD1D3 = plot([BD1ex,BD3ex],[BD1ey,BD3ey],'k');
set(BD1D3,'LineWidth',1.5);

BD3D5 = plot([BD3ex,BD5ex],[BD3ey,BD5ey],'k');
set(BD3D5,'LineWidth',1.5);

BD4D6 = plot([BD4ex,BD6ex],[BD4ey,BD6ey],'k');
set(BD4D6,'LineWidth',1.5);

BD6D2 = plot([BD6ex,BD2ex],[BD6ey,BD2ey],'k');
set(BD6D2,'LineWidth',1.5);
%%%%%%%%%%%%%%%%%%%%%%Line Joining Top and Bottom bridges
U6B3=plot([BD3ex UD6ex],[BD3ey UD6ey],'r');
set(U6B3,'LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]	);
%%%%%%%%%%%%%%%%%%%%%%%Line Joing Midpoints of bottom Bridge & Transformer part
xrefb1=-18;
yrefb1=-17;

[lbx1,lby1]=Indc(xrefb1,yrefb1);


xrefb2=-16.5;
yrefb2=-17;

[lbx2,lby2]=Indc(xrefb2,yrefb2);

xrefb3=-15;
yrefb3=-17;

[lbx3,lby3]=Indc(xrefb3,yrefb3);


lji3i2= plot([lbx3(:,end) lbx3(:,end)-0.5 lbx2(:,1)+0.5 lbx2(:,1)], [lby3(:,end) lby3(:,end) lby2(:,1) lby2(:,1)],'k');
set(lji3i2,'LineWidth',1.5);

lji2i1= plot([lbx2(:,end) lbx2(:,end)-0.5 lbx1(:,1)+0.5 lbx1(:,1)], [lby2(:,end) lby2(:,end) lby1(:,1) lby1(:,1)],'k');
set(lji2i1,'LineWidth',1.5);

lji1i3= plot([lbx1(:,end) lbx1(:,end)-0.75 lbx1(:,end)-0.75 lbx3(:,1) lbx3(:,1)], [lby1(:,end) lby1(:,end) lby1(:,1)+0.5 lby1(:,1)+0.5 lby3(:,1)],'k');
set(lji1i3,'LineWidth',1.5);

%bl1mp=plot([bl1mpx,bl1mpx-sp],[bl1mpy,bl1mpy],'k');
bl1mp=plot([bl1mpx,lbx3(:,end) lbx3(:,end)],[bl1mpy,bl1mpy,lby3(:,end)],'k');
set(bl1mp,'LineWidth',1.5);

bl2mp=plot([bl2mpx,lbx2(:,end) lbx2(:,end)],[bl2mpy,bl2mpy,lby2(:,end)],'k');
set(bl2mp,'LineWidth',1.5);

bl3mp=plot([bl3mpx,lbx1(:,end) lbx1(:,end)],[bl3mpy,bl3mpy,lby1(:,end)],'k');
set(bl3mp,'LineWidth',1.5);


trb1=plot(lbx1,lby1,'k');
set(trb1,'LineWidth',1.5);
text(-17,-18.5,'B','FontSize',8);
hold on
trb2=plot(lbx2,lby2,'k');
set(trb2,'LineWidth',1.5);
text(-15.75,-17.5,'Y','FontSize',8);
hold on
trb3=plot(lbx3,lby3,'k');
set(trb3,'LineWidth',1.5);
text(-14.25,-16.5,'R','FontSize',8);
hold on
%%%%%%%%% Magnetic Coupling Lines
L1=plot([lx1(:,end)-1.25 lbx1(:,end)-1.25],[ul1mpy bl3mpy],'k');
set(L1,'LineWidth',1.5);

L2=plot([lx1(:,end)-1.5 lbx1(:,end)-1.5],[ul1mpy bl3mpy],'k');
set(L2,'LineWidth',1.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load Resistor
rxref=-6;
ryref=-16/3;

[rx1,ry1]=resi(rxref,ryref);
rx1 =[UD5ex rx1(:,1) rx1 rx1(:,end) BD2ex];
ry1 =[UD5ey UD5ey 3*ry1 BD2ey BD2ey];
res=plot(rx1,ry1,'r');
set(res,'LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]	);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  horizontal inductor source side
xrefs1=-22.5;
yrefs1=-10;

[lbsx1,lbsy1]=Indch(xrefs1,yrefs1,3,1.25);
trs1=plot(lbsx1,lbsy1,'k');
set(trs1,'LineWidth',1.5);

xrefs1=-22.5;
yrefs1=-12;

[lbsx2,lbsy2]=Indch(xrefs1,yrefs1,3,1.25);
trs2=plot(lbsx2,lbsy2,'k');
set(trs2,'LineWidth',1.5);

xrefs1=-22.5;
yrefs1=-14;

[lbsx3,lbsy3]=Indch(xrefs1,yrefs1,3,1.25);
trs3=plot(lbsx3,lbsy3,'k');
set(trs3,'LineWidth',1.5);

lj3_ind=plot([lbsx1(:,end) lbsx2(:,end) lbsx3(:,end)], [lbsy1(:,end) lbsy2(:,end) lbsy3(:,end)],'k');
set(lj3_ind,'LineWidth',1.5);

text(lbsx1(:,1)-2.5,lbsy1(:,1),'V_{RN}(p)','FontSize',9.5);
text(lbsx2(:,1)-2.5,lbsy2(:,1),'V_{YN}(p)','FontSize',9.5);
text(lbsx3(:,1)-2.5,lbsy3(:,1),'V_{BN}(p)','FontSize',9.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Voltages from Simulink File

Vrybn_Sec_Star = squeeze(out.Sec_PN_Star.V.data);
Vryb_Sec_Star = squeeze(out.Sec_PP_Star.V.data);
Vryb_Sec_Delta = squeeze(out.Sec_PP_Delta.V.data);
time=squeeze(out.tout);
L = rescale(time,10,30);


Van_s = Vrybn_Sec_Star(:,1)./max(Vrybn_Sec_Star(:,1));
Vbn_s = Vrybn_Sec_Star(:,2)./max(Vrybn_Sec_Star(:,2));
Vcn_s = Vrybn_Sec_Star(:,3)./max(Vrybn_Sec_Star(:,3));


Vab_s = Vryb_Sec_Star(:,1)./max(Vryb_Sec_Star(:,1));
Vbc_s = Vryb_Sec_Star(:,2)./max(Vryb_Sec_Star(:,2));
Vca_s = Vryb_Sec_Star(:,3)./max(Vryb_Sec_Star(:,3));
Vry_star=rescale(Vab_s,-4,4);
Vyb_star=rescale(Vbc_s,-4,4);
Vbr_star=rescale(Vca_s,-4,4);

Vab_d = Vryb_Sec_Delta(:,1)./max(Vryb_Sec_Delta(:,1));
Vbc_d = Vryb_Sec_Delta(:,2)./max(Vryb_Sec_Delta(:,2));
Vca_d = Vryb_Sec_Delta(:,3)./max(Vryb_Sec_Delta(:,3));
Vry_delta=rescale(Vab_d,-29,-21);
Vyb_delta=rescale(Vbc_d,-29,-21);
Vbr_delta=rescale(Vca_d,-29,-21);

Vry_delta_s=-1*Vryb_Sec_Delta(:,1)./max(-1*Vryb_Sec_Delta(:,1));
Vyb_delta_s=-1*Vryb_Sec_Delta(:,2)./max(-1*Vryb_Sec_Delta(:,2));
Vbr_delta_s=-1*Vryb_Sec_Delta(:,3)./max(-1*Vryb_Sec_Delta(:,3));

Vry_delta_shifted=rescale(Vry_delta_s,-29,-21);
Vyb_delta_shifted=rescale(Vyb_delta_s,-29,-21);
Vbr_delta_shifted=rescale(Vbr_delta_s,-29,-21);

%%%%%%%%%%%%%%%%% Current from simulink file

Iabcs = squeeze(out.Sec_PP_Star.I.data);
Iab_Star = Iabcs(:,1)./max(Iabcs(:,1));
Ibc_Star = Iabcs(:,2)./max(Iabcs(:,2));
Ica_Star = Iabcs(:,3)./max(Iabcs(:,3));


Iabcd = squeeze(out.Sec_PP_Delta.I.data);
Iab_Delta = Iabcd(:,1)./max(Iabcd(:,1));
Ibc_Delta = Iabcd(:,2)./max(Iabcd(:,2));
Ica_Delta = Iabcd(:,3)./max(Iabcd(:,3));
%%%%%% Current Plot
plot([-30 -30],[-5 5],'LineWidth',0.750,'Color',"k");
plot([-30 -20],[0 0],'LineWidth',0.750,'Color',"k");
text(-19.5,0,'\omegat','FontSize',10);
T = rescale(time,-30,-20);
Iab_Star_Scaled= rescale(Iab_Star,-4.5,4.5);
Iab_Delta_Scaled= rescale(Iab_Delta,-4.5,4.5);



Iab_Star_Wave=animatedline('LineWidth',1.0,'Color',"r",'LineStyle',"-");
Iab_Delta_Wave=animatedline('LineWidth',1.0,'Color',"b",'LineStyle',"-");
text(-29.5,6,'SECONDARY- R-PHASE CURRENT','FontSize',8);


plot([-30 -30],[-20 -30],'LineWidth',0.750,'Color',"k");
plot([-30 -20],[-25 -25],'LineWidth',0.750,'Color',"k");
text(-29.5,-20,'PRIMARY SIDE- R-PHASE CURRENT','FontSize',8);
text(-19.5,-25,'\omegat','FontSize',10);

text(-14.5,-27,'Developed By','FontSize',8);
text(-15,-28,'Dr.M.Kaliamoorthy','FontSize',8);

Iab_R_Total=Iab_Star+Iab_Delta;
Iab_Total_Scaled= rescale(Iab_R_Total,-29,-21);


Iab_Total_Wave=animatedline('LineWidth',1.5,'Color',[0.6350 0.0780 0.1840]	,'LineStyle',"-");

%%%%%%%%%%%%%%%%%%%%%%%% Angle calculations


angVr = mod((2*pi*50.*time),2*pi)*(180/pi);
angVy = mod((2*pi*50.*time)-(2*pi/3),2*pi)*(180/pi);
angVb = mod((2*pi*50.*time)+(2*pi/3),2*pi)*(180/pi);

angVry_Star = mod((2*pi*50.*time)-(pi/6),2*pi)*(180/pi);
angVyb_Star = mod((2*pi*50.*time)-(2*pi/3)-(pi/6),2*pi)*(180/pi);
angVbr_Star = mod((2*pi*50.*time)+(2*pi/3)-(pi/6),2*pi)*(180/pi);

angVry_Star_Shifted = mod((2*pi*50.*time)-(pi/6)+pi,2*pi)*(180/pi);
angVyb_Star_Shifted = mod((2*pi*50.*time)-(2*pi/3)-(pi/6)+pi,2*pi)*(180/pi);
angVbr_Star_Shifted = mod((2*pi*50.*time)+(2*pi/3)-(pi/6)+pi,2*pi)*(180/pi);


% angVry_Delta = mod((2*pi*50.*time)-(pi/3),2*pi)*(180/pi);  %Red  angVry_Delta
% angVyb_Delta = mod((2*pi*50.*time)-(2*pi/3)-(pi/3),2*pi)*(180/pi); %Yellow  angVyb_Delta
% angVbr_Delta = mod((2*pi*50.*time)+(2*pi/3)-(pi/3),2*pi)*(180/pi); %Blue  angVbr_Delta
% 
% angVry_Delta_Shifted = mod((2*pi*50.*time)-(pi/3)+pi,2*pi)*(180/pi); % Red Dotted  angVry_Delta_Shifted
% angVyb_Delta_Shifted = mod((2*pi*50.*time)-(2*pi/3)-(pi/3)+pi,2*pi)*(180/pi); %Yellow Dotted  angVyb_Delta_Shifted
% angVbr_Delta_Shifted = mod((2*pi*50.*time)+(2*pi/3)-(pi/3)+pi,2*pi)*(180/pi); % Blue Dotted  angVbr_Delta_Shifted

angVry_Delta = mod((2*pi*50.*time)-(2*pi/3)-(pi/3)+pi,2*pi)*(180/pi);
angVyb_Delta = mod((2*pi*50.*time)+(2*pi/3)-(pi/3)+pi,2*pi)*(180/pi);
angVbr_Delta= mod((2*pi*50.*time)-(pi/3)+pi,2*pi)*(180/pi);

angVry_Delta_Shifted = mod((2*pi*50.*time)-(2*pi/3)-(pi/3),2*pi)*(180/pi);
angVyb_Delta_Shifted = mod((2*pi*50.*time)+(2*pi/3)-(pi/3),2*pi)*(180/pi);
angVbr_Delta_Shifted =mod((2*pi*50.*time)-(pi/3),2*pi)*(180/pi);


Vrn_Phasor=animatedline('LineWidth',1.75,'Color',"r",'LineStyle',"-.");
Vrn_Phasor_cursor=animatedline('LineWidth',1.75,'Color',"r", 'Marker', '.', 'MarkerSize', 20);

Vyn_Phasor=animatedline('LineWidth',1.75,'Color',[0.9290 0.6940 0.1250],'LineStyle',"-.");
Vyn_Phasor_cursor=animatedline('LineWidth',1.75,'Color',[0.9290 0.6940 0.1250], 'Marker', '.', 'MarkerSize', 20);

Vbn_Phasor=animatedline('LineWidth',1.75,'Color',"b",'LineStyle',"-.");
Vbn_Phasor_cursor=animatedline('LineWidth',1.75,'Color',"b", 'Marker', '.', 'MarkerSize', 20);

Vry_Phasor=animatedline('LineWidth',1.75,'Color',"r");
Vry_Phasor_cursor=animatedline('LineWidth',1.75,'Color',"r", 'Marker', '.', 'MarkerSize', 20);

Vyb_Phasor=animatedline('LineWidth',1.75,'Color',[0.9290 0.6940 0.1250]);
Vyb_Phasor_cursor=animatedline('LineWidth',1.75,'Color',[0.9290 0.6940 0.1250], 'Marker', '.', 'MarkerSize', 20);

Vbr_Phasor=animatedline('LineWidth',1.75,'Color',"b");
Vbr_Phasor_cursor=animatedline('LineWidth',1.75,'Color',"b", 'Marker', '.', 'MarkerSize', 20);

Vrys_Phasor=animatedline('LineWidth',1.75,'Color',"r",'LineStyle',":");
Vrys_Phasor_cursor=animatedline('LineWidth',1.75,'Color',"r", 'Marker', '.', 'MarkerSize', 20);

Vybs_Phasor=animatedline('LineWidth',1.75,'Color',[0.9290 0.6940 0.1250],'LineStyle',":");
Vybs_Phasor_cursor=animatedline('LineWidth',1.75,'Color',[0.9290 0.6940 0.1250], 'Marker', '.', 'MarkerSize', 20);

Vbrs_Phasor=animatedline('LineWidth',1.75,'Color',"b",'LineStyle',":");
Vbrs_Phasor_cursor=animatedline('LineWidth',1.75,'Color',"b", 'Marker', '.', 'MarkerSize', 20);

%%%%%% Delta Animated Phasor Lines

    Vry_Delta_Phasor=animatedline('LineWidth',2,'Color',"r");
    Vry_Delta_Phasor_cursor=animatedline('LineWidth',2,'Color',"r", 'Marker', '.', 'MarkerSize', 20);

    Vyb_Delta_Phasor=animatedline('LineWidth',2,'Color',[0.9290 0.6940 0.1250]);
    Vyb_Delta_Phasor_cursor=animatedline('LineWidth',2,'Color',[0.9290 0.6940 0.1250], 'Marker', '.', 'MarkerSize', 20);

    Vbr_Delta_Phasor=animatedline('LineWidth',2,'Color',"b");
    Vbr_Delta_Phasor_cursor=animatedline('LineWidth',2,'Color',"b", 'Marker', '.', 'MarkerSize', 20);

    Vrys_Delta_Phasor=animatedline('LineWidth',2,'Color',"r",'LineStyle',":");
    Vrys_Delta_Phasor_cursor=animatedline('LineWidth',2,'Color',"r", 'Marker', '.', 'MarkerSize', 20);

    Vybs_Delta_Phasor=animatedline('LineWidth',2,'Color',[0.9290 0.6940 0.1250],'LineStyle',":");
    Vybs_Delta_Phasor_cursor=animatedline('LineWidth',2,'Color',[0.9290 0.6940 0.1250], 'Marker', '.', 'MarkerSize', 20);

    Vbrs_Delta_Phasor=animatedline('LineWidth',2,'Color',"b",'LineStyle',":");
    Vbrs_Delta_Phasor_cursor=animatedline('LineWidth',2,'Color',"b", 'Marker', '.', 'MarkerSize', 20);

%%%%%%%%%%%%%%%%%%%%%%%%% Waves Animated Lines
%%%% Star Side
plot([10 10],[4.75,-4.75],'LineWidth',1.0,'Color',"k")
plot([10 30],[0,0],'LineWidth',1.0,'Color',"k")
text(29,1,'\omegat','FontSize',10);
text(10,5,'v_{line}','FontSize',10);
text(13,5,'STAR SIDE- LINE VOLTAGE WAVEFORMS','FontSize',8);


Vry_Star_Wave=animatedline('LineWidth',1.25,'Color',[0.8500 0.3250 0.0980],'LineStyle',"-");
Vyb_Star_Wave=animatedline('LineWidth',1.25,'Color',[0.9290 0.6940 0.1250],'LineStyle',"-");
Vbr_Star_Wave=animatedline('LineWidth',1.25,'Color',"b",'LineStyle',"-");

Vry_Star_Wave_Shifted=animatedline('LineWidth',1.5,'Color',[0.8500 0.3250 0.0980],'LineStyle',":");
Vyb_Star_Wave_Shifted=animatedline('LineWidth',1.5,'Color',[0.9290 0.6940 0.1250],'LineStyle',":");
Vbr_Star_Wave_Shifted=animatedline('LineWidth',1.5,'Color',"b",'LineStyle',":");

lgd = legend([Vry_Star_Wave Vyb_Star_Wave Vbr_Star_Wave Vry_Star_Wave_Shifted Vyb_Star_Wave_Shifted Vbr_Star_Wave_Shifted], 'V_{RY}', 'V_{YB}', 'V_{BR}','V_{YR}','V_{BY}','V_{RB}','Location','northeast','Orientation','horizontal','TextColor','blue','FontSize',10,'AutoUpdate','off');
lgd.NumColumns = 1;


%%% Delta Side
plot([10 10],[-29.75,-20.75],'LineWidth',1.0,'Color',"k")
plot([10 30],[-25,-25],'LineWidth',1.0,'Color',"k")




text(29,-24,'\omegat','FontSize',10);
text(10,-20,'v_{line}','FontSize',10);
text(13,-20,'DELTA SIDE- LINE VOLTAGE WAVEFORMS','FontSize',8);

Vry_Delta_Wave=animatedline('LineWidth',1.25,'Color',[0.8500 0.3250 0.0980],'LineStyle',"-");
Vyb_Delta_Wave=animatedline('LineWidth',1.25,'Color',[0.9290 0.6940 0.1250],'LineStyle',"-");
Vbr_Delta_Wave=animatedline('LineWidth',1.25,'Color',"b",'LineStyle',"-");

Vry_Delta_Wave_Shifted=animatedline('LineWidth',1.5,'Color',[0.8500 0.3250 0.0980],'LineStyle',":");
Vyb_Delta_Wave_Shifted=animatedline('LineWidth',1.5,'Color',[0.9290 0.6940 0.1250],'LineStyle',":");
Vbr_Delta_Wave_Shifted=animatedline('LineWidth',1.5,'Color',"b",'LineStyle',":");



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load Voltage Animated Waves

plot([-4 -4 9],[-6 -19 -19],'LineWidth',1.0,'Color','k');   %%Bridge Voltage axis
plot([-4 9],[-12.5 -12.5],'LineWidth',1.0,'Color','k');

text(-3.0,-6,'UPPER BRIDGE-LOAD VOLTAGE','FontSize',8);
text(1.0,-7,'(UB-LV)','FontSize',8);
text(-3.0,-14,'LOWER BRIDGE-LOAD VOLTAGE','FontSize',8);
text(1.0,-15,'(LB-LV)','FontSize',8);
text(8.5,-11.5,'\omegat','FontSize',10);
text(8.5,-18,'\omegat','FontSize',10);

V_Tload = squeeze(out.V_Load.V_Total.data);
V_up = squeeze(out.V_Load.V_Up.data);
V_low = squeeze(out.V_Load.V_Lower.data);

L1 = rescale(time,-4,9);
V_up_scaled = rescale(V_up,-9,-8);
V_low_scaled = rescale(V_low,-17,-16);
V_Tload_scaled = rescale(V_Tload,-12,-10);

Vup_Wave=animatedline('LineWidth',1.25,'Color',[0.6350 0.0780 0.1840],'LineStyle',"-");
Vlow_Wave=animatedline('LineWidth',1.25,'Color',[0 0.4470 0.7410],'LineStyle',"-");
VTotal_Wave=animatedline('LineWidth',1.25,'Color',[1 0 1],'LineStyle',"-");

plot([10 10 30],[-6 -15 -15],'LineWidth',1.0,'Color','k');
text(15.0,-7,'TOTAL-LOAD VOLTAGE','FontSize',8);
text(16.0,-8,'(UB-LV + LB-LV)','FontSize',8);
text(29,-14,'\omegat','FontSize',10);
text(10,-5.25,'v_{dc}','FontSize',10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nf=4001;
tc = linspace(0,2*pi,Nf); % vector of angles to draw circles
circ_radius_star_L = 4.25/sqrt(3);
circ_radius_star_P = (4.25)/2;
Line_Circle_star= plot(((circ_radius_star_L-0)*sqrt(3)*cos(tc))+0,(circ_radius_star_L-0)*sqrt(3)*sin(tc)-0,'Color','k');
Phase_Circle_star= plot(((circ_radius_star_P-0)*cos(tc))+0,(circ_radius_star_P-0)*sin(tc)-0,'Color','k');
pos1 = get(gca, 'Position');

plot([2*circ_radius_star_P, -2*circ_radius_star_P],[0 0],'LineStyle',":",'LineWidth',1.5,'Color','k');
plot([0 0], [2*circ_radius_star_P, -2*circ_radius_star_P],'LineStyle',":",'LineWidth',1.5,'Color','k');

annotation('textarrow',[0.57 0.54],[0.77 0.77],'String','Ref','FontSize',10,'Linewidth',1.25,'TextColor','r')


circ_radius_Delta_L = 4.25/sqrt(3);
circ_radius_Delta_P = (4.25)/2;
Line_Circle_Delta= plot(((circ_radius_Delta_L-0)*sqrt(3)*cos(tc))+0,(circ_radius_Delta_L-0)*sqrt(3)*sin(tc)-25,'Color','k');
Phase_Circle_Delta= plot(((circ_radius_Delta_P-0)*cos(tc))+0,(circ_radius_Delta_P-0)*sin(tc)-25,'Color','k');

plot([2*circ_radius_Delta_P, -2*circ_radius_Delta_P],[-25 -25],'LineStyle',":",'LineWidth',1.5,'Color','k');
plot([0 0], [-25-(2*circ_radius_Delta_P), -25-(-2*circ_radius_Delta_P)],'LineStyle',":",'LineWidth',1.5,'Color','k');




annotation('textarrow',[0.57 0.54],[0.2175 0.2175],'String','Ref','FontSize',10,'Linewidth',1.25,'TextColor','r')

axis([-31 37 -30 7],'xy');
set(gca,'XColor', 'none','YColor','none')
%set(gca, 'Color', 'none')
%whitebg('w')
hf.Color = 'w';
set(gcf,'Color',[0.90 0.90 0.90])
hf.Color='w';


for j=1:1:2

    clearpoints(Vry_Star_Wave);
    clearpoints(Vyb_Star_Wave);
    clearpoints(Vbr_Star_Wave);

    clearpoints(Vry_Star_Wave_Shifted);
    clearpoints(Vyb_Star_Wave_Shifted);
    clearpoints(Vbr_Star_Wave_Shifted);

    clearpoints(Vry_Delta_Wave);
    clearpoints(Vyb_Delta_Wave);
    clearpoints(Vbr_Delta_Wave);

    clearpoints(Vry_Delta_Wave_Shifted);
    clearpoints(Vyb_Delta_Wave_Shifted);
    clearpoints(Vbr_Delta_Wave_Shifted);

    clearpoints(Vup_Wave);
    clearpoints(Vlow_Wave);
    clearpoints(VTotal_Wave);

    clearpoints(Iab_Star_Wave);
    clearpoints(Iab_Delta_Wave);
    clearpoints(Iab_Total_Wave);
    
for i=1:10:length(time)

    %%%%%%%%% Phase Voltage Star side
    clearpoints(Vrn_Phasor);
    clearpoints(Vrn_Phasor_cursor);

    addpoints(Vrn_Phasor,[0, circ_radius_star_P*cosd(angVr(i))+0],[0,circ_radius_star_P*sind(angVr(i))-0]);
    addpoints(Vrn_Phasor_cursor,circ_radius_star_P*cosd(angVr(i))+0,circ_radius_star_P*sind(angVr(i))-0);
    txt_vrn = text(1.05*(circ_radius_star_P-0)*cosd(angVr(i))+0.25, 1.05*(circ_radius_star_P-0)*sind(angVr(i))+0.25,'V_{RN}','FontSize',8); 


    clearpoints(Vyn_Phasor);
    clearpoints(Vyn_Phasor_cursor);

    addpoints(Vyn_Phasor,[0, circ_radius_star_P*cosd(angVy(i))+0],[0,circ_radius_star_P*sind(angVy(i))-0]);
    addpoints(Vyn_Phasor_cursor,circ_radius_star_P*cosd(angVy(i))+0,circ_radius_star_P*sind(angVy(i))-0);
    txt_vyn = text(1.05*(circ_radius_star_P-0)*cosd(angVy(i))+0.25, 1.05*(circ_radius_star_P-0)*sind(angVy(i))+0.25,'V_{YN}','FontSize',8); 

    clearpoints(Vbn_Phasor);
    clearpoints(Vbn_Phasor_cursor);

    addpoints(Vbn_Phasor,[0, circ_radius_star_P*cosd(angVb(i))+0],[0,circ_radius_star_P*sind(angVb(i))-0]);
    addpoints(Vbn_Phasor_cursor,circ_radius_star_P*cosd(angVb(i))+0,circ_radius_star_P*sind(angVb(i))-0);
    txt_vbn = text(1.05*(circ_radius_star_P-0)*cosd(angVb(i))+0.25, 1.05*(circ_radius_star_P-0)*sind(angVb(i))+0.25,'V_{BN}','FontSize',8); 
    
    %%%%%%%%%%% Line Voltages Star Side
    clearpoints(Vry_Phasor);
    clearpoints(Vry_Phasor_cursor);

    addpoints(Vry_Phasor,[0, circ_radius_star_L*sqrt(3)*cosd(angVry_Star(i))+0],[0,circ_radius_star_L*sqrt(3)*sind(angVry_Star(i))-0]);
    addpoints(Vry_Phasor_cursor,circ_radius_star_L*sqrt(3)*cosd(angVry_Star(i))+0,circ_radius_star_L*sqrt(3)*sind(angVry_Star(i))-0);
    txt_vry = text(1.05*(circ_radius_star_L-0)*sqrt(3)*cosd(angVry_Star(i))+0.25, 1.05*(circ_radius_star_L-0)*sqrt(3)*sind(angVry_Star(i))+0.25,'V_{RY}','FontSize',8); 

    clearpoints(Vyb_Phasor);
    clearpoints(Vyb_Phasor_cursor);

    addpoints(Vyb_Phasor,[0, circ_radius_star_L*sqrt(3)*cosd(angVyb_Star(i))+0],[0,circ_radius_star_L*sqrt(3)*sind(angVyb_Star(i))-0]);
    addpoints(Vyb_Phasor_cursor,circ_radius_star_L*sqrt(3)*cosd(angVyb_Star(i))+0,circ_radius_star_L*sqrt(3)*sind(angVyb_Star(i))-0);
    txt_vyb = text(1.05*(circ_radius_star_L-0)*sqrt(3)*cosd(angVyb_Star(i))+0.25, 1.05*(circ_radius_star_L-0)*sqrt(3)*sind(angVyb_Star(i))+0.25,'V_{YB}','FontSize',8); 


    clearpoints(Vbr_Phasor);
    clearpoints(Vbr_Phasor_cursor);

    addpoints(Vbr_Phasor,[0, circ_radius_star_L*sqrt(3)*cosd(angVbr_Star(i))+0],[0,circ_radius_star_L*sqrt(3)*sind(angVbr_Star(i))-0]);
    addpoints(Vbr_Phasor_cursor,circ_radius_star_L*sqrt(3)*cosd(angVbr_Star(i))+0,circ_radius_star_L*sqrt(3)*sind(angVbr_Star(i))-0);
    txt_vbr = text(1.05*(circ_radius_star_L-0)*sqrt(3)*cosd(angVbr_Star(i))+0.25, 1.05*(circ_radius_star_L-0)*sqrt(3)*sind(angVbr_Star(i))+0.25,'V_{BR}','FontSize',8); 

    %%%%%%%%%%%%%%%%%% Line Voltages Star Side 180 degrees shifted


    clearpoints(Vrys_Phasor);
    clearpoints(Vrys_Phasor_cursor);

    addpoints(Vrys_Phasor,[0, circ_radius_star_L*sqrt(3)*cosd(angVry_Star_Shifted(i))+0],[0,circ_radius_star_L*sqrt(3)*sind(angVry_Star_Shifted(i))-0]);
    addpoints(Vrys_Phasor_cursor,circ_radius_star_L*sqrt(3)*cosd(angVry_Star_Shifted(i))+0,circ_radius_star_L*sqrt(3)*sind(angVry_Star_Shifted(i))-0);
    txt_vrys = text(1.05*(circ_radius_star_L-0)*sqrt(3)*cosd(angVry_Star_Shifted(i))+0.25, 1.05*(circ_radius_star_L-0)*sqrt(3)*sind(angVry_Star_Shifted(i))+0.25,'V_{YR}','FontSize',8); 

    clearpoints(Vybs_Phasor);
    clearpoints(Vybs_Phasor_cursor);

    addpoints(Vybs_Phasor,[0, circ_radius_star_L*sqrt(3)*cosd(angVyb_Star_Shifted(i))+0],[0,circ_radius_star_L*sqrt(3)*sind(angVyb_Star_Shifted(i))-0]);
    addpoints(Vybs_Phasor_cursor,circ_radius_star_L*sqrt(3)*cosd(angVyb_Star_Shifted(i))+0,circ_radius_star_L*sqrt(3)*sind(angVyb_Star_Shifted(i))-0);
    txt_vybs = text(1.05*(circ_radius_star_L-0)*sqrt(3)*cosd(angVyb_Star_Shifted(i))+0.25, 1.05*(circ_radius_star_L-0)*sqrt(3)*sind(angVyb_Star_Shifted(i))+0.25,'V_{BY}','FontSize',8); 


    clearpoints(Vbrs_Phasor);
    clearpoints(Vbrs_Phasor_cursor);

    addpoints(Vbrs_Phasor,[0, circ_radius_star_L*sqrt(3)*cosd(angVbr_Star_Shifted(i))+0],[0,circ_radius_star_L*sqrt(3)*sind(angVbr_Star_Shifted(i))-0]);
    addpoints(Vbrs_Phasor_cursor,circ_radius_star_L*sqrt(3)*cosd(angVbr_Star_Shifted(i))+0,circ_radius_star_L*sqrt(3)*sind(angVbr_Star_Shifted(i))-0);
    txt_vbrs = text(1.05*(circ_radius_star_L-0)*sqrt(3)*cosd(angVbr_Star_Shifted(i))+0.25, 1.05*(circ_radius_star_L-0)*sqrt(3)*sind(angVbr_Star_Shifted(i))+0.25,'V_{RB}','FontSize',8); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Delta Line Voltage Phasors
    clearpoints(Vry_Delta_Phasor);
    clearpoints(Vry_Delta_Phasor_cursor);

    addpoints(Vry_Delta_Phasor,[0, circ_radius_Delta_L*sqrt(3)*cosd(angVry_Delta(i))+0],[-25,circ_radius_Delta_L*sqrt(3)*sind(angVry_Delta(i))-25]);
    addpoints(Vry_Delta_Phasor_cursor,circ_radius_Delta_L*sqrt(3)*cosd(angVry_Delta(i))+0,circ_radius_Delta_L*sqrt(3)*sind(angVry_Delta(i))-25);
    txt_vry_Delta = text(1.05*(circ_radius_Delta_L-0)*sqrt(3)*cosd(angVry_Delta(i))+0.25, 1.05*(circ_radius_Delta_L-0)*sqrt(3)*sind(angVry_Delta(i))+0.25-25,'V_{RY}','FontSize',8); 

    clearpoints(Vyb_Delta_Phasor);
    clearpoints(Vyb_Delta_Phasor_cursor);

    addpoints(Vyb_Delta_Phasor,[0, circ_radius_Delta_L*sqrt(3)*cosd(angVyb_Delta(i))+0],[-25,circ_radius_Delta_L*sqrt(3)*sind(angVyb_Delta(i))-25]);
    addpoints(Vyb_Delta_Phasor_cursor,circ_radius_Delta_L*sqrt(3)*cosd(angVyb_Delta(i))+0,circ_radius_Delta_L*sqrt(3)*sind(angVyb_Delta(i))-25);
    txt_vyb_Delta = text(1.05*(circ_radius_Delta_L-0)*sqrt(3)*cosd(angVyb_Delta(i))+0.25, 1.05*(circ_radius_Delta_L-0)*sqrt(3)*sind(angVyb_Delta(i))+0.25-25,'V_{YB}','FontSize',8); 


    clearpoints(Vbr_Delta_Phasor);
    clearpoints(Vbr_Delta_Phasor_cursor);

    addpoints(Vbr_Delta_Phasor,[0, circ_radius_Delta_L*sqrt(3)*cosd(angVbr_Delta(i))+0],[-25,circ_radius_Delta_L*sqrt(3)*sind(angVbr_Delta(i))-25]);
    addpoints(Vbr_Delta_Phasor_cursor,circ_radius_Delta_L*sqrt(3)*cosd(angVbr_Delta(i))+0,circ_radius_Delta_L*sqrt(3)*sind(angVbr_Delta(i))-25);
    txt_vbr_Delta = text(1.05*(circ_radius_Delta_L-0)*sqrt(3)*cosd(angVbr_Delta(i))+0.25, 1.05*(circ_radius_Delta_L-0)*sqrt(3)*sind(angVbr_Delta(i))+0.25-25,'V_{BR}','FontSize',8); 

    %%%%%%%%%%%%%%% Line Voltage Delta Shifted

    clearpoints(Vrys_Delta_Phasor);
    clearpoints(Vrys_Delta_Phasor_cursor);

    addpoints(Vrys_Delta_Phasor,[0, circ_radius_Delta_L*sqrt(3)*cosd(angVry_Delta_Shifted(i))+0],[-25,circ_radius_Delta_L*sqrt(3)*sind(angVry_Delta_Shifted(i))-25]);
    addpoints(Vrys_Delta_Phasor_cursor,circ_radius_Delta_L*sqrt(3)*cosd(angVry_Delta_Shifted(i))+0,circ_radius_Delta_L*sqrt(3)*sind(angVry_Delta_Shifted(i))-25);
    txt_vrys_Delta = text(1.05*(circ_radius_Delta_L-0)*sqrt(3)*cosd(angVry_Delta_Shifted(i))+0.25, 1.05*(circ_radius_Delta_L-0)*sqrt(3)*sind(angVry_Delta_Shifted(i))+0.25-25,'V_{YR}','FontSize',8); 

    clearpoints(Vybs_Delta_Phasor);
    clearpoints(Vybs_Delta_Phasor_cursor);

    addpoints(Vybs_Delta_Phasor,[0, circ_radius_Delta_L*sqrt(3)*cosd(angVyb_Delta_Shifted(i))+0],[-25,circ_radius_Delta_L*sqrt(3)*sind(angVyb_Delta_Shifted(i))-25]);
    addpoints(Vybs_Delta_Phasor_cursor,circ_radius_Delta_L*sqrt(3)*cosd(angVyb_Delta_Shifted(i))+0,circ_radius_Delta_L*sqrt(3)*sind(angVyb_Delta_Shifted(i))-25);
    txt_vybs_Delta = text(1.05*(circ_radius_Delta_L-0)*sqrt(3)*cosd(angVyb_Delta_Shifted(i))+0.25, 1.05*(circ_radius_Delta_L-0)*sqrt(3)*sind(angVyb_Delta_Shifted(i))+0.25-25,'V_{BY}','FontSize',8); 


    clearpoints(Vbrs_Delta_Phasor);
    clearpoints(Vbrs_Delta_Phasor_cursor);

    addpoints(Vbrs_Delta_Phasor,[0, circ_radius_Delta_L*sqrt(3)*cosd(angVbr_Delta_Shifted(i))+0],[-25,circ_radius_Delta_L*sqrt(3)*sind(angVbr_Delta_Shifted(i))-25]);
    addpoints(Vbrs_Delta_Phasor_cursor,circ_radius_Delta_L*sqrt(3)*cosd(angVbr_Delta_Shifted(i))+0,circ_radius_Delta_L*sqrt(3)*sind(angVbr_Delta_Shifted(i))-25);
    txt_vbrs_Delta = text(1.05*(circ_radius_Delta_L-0)*sqrt(3)*cosd(angVbr_Delta_Shifted(i))+0.25, 1.05*(circ_radius_Delta_L-0)*sqrt(3)*sind(angVbr_Delta_Shifted(i))+0.25-25,'V_{RB}','FontSize',8); 


    addpoints(Vry_Star_Wave,L(i),Vry_star(i));
    addpoints(Vyb_Star_Wave,L(i),Vyb_star(i));
    addpoints(Vbr_Star_Wave,L(i),Vbr_star(i));

    addpoints(Vry_Star_Wave_Shifted,L(i),-1*Vry_star(i));
    addpoints(Vyb_Star_Wave_Shifted,L(i),-1*Vyb_star(i));
    addpoints(Vbr_Star_Wave_Shifted,L(i),-1*Vbr_star(i));


    addpoints(Vry_Delta_Wave,L(i),Vry_delta(i));
    addpoints(Vyb_Delta_Wave,L(i),Vyb_delta(i));
    addpoints(Vbr_Delta_Wave,L(i),Vbr_delta(i));

    addpoints(Vry_Delta_Wave_Shifted,L(i),Vry_delta_shifted(i));
    addpoints(Vyb_Delta_Wave_Shifted,L(i),Vyb_delta_shifted(i));
    addpoints(Vbr_Delta_Wave_Shifted,L(i),Vbr_delta_shifted(i));

    
    addpoints(Vup_Wave,L1(i),V_up_scaled(i));
    addpoints(Vlow_Wave,L1(i),V_low_scaled(i));
    addpoints(VTotal_Wave,L(i),V_Tload_scaled(i));

    addpoints(Iab_Star_Wave,T(i),Iab_Star_Scaled(i));
    addpoints(Iab_Delta_Wave,T(i),Iab_Delta_Scaled(i));
    addpoints(Iab_Total_Wave,T(i),Iab_Total_Scaled(i))
    %%%%%%%%%%%%%%%%%%%%%% Color Changing code

    if Iab_Star(i)>=0.1 && Ibc_Star(i)<=-0.1
        set(D1,'Color',[0.9290 0.6940 0.1250]	);
        set(D2,'Color',[0.7 0.7 0.7]);
        set(D3,'Color',[0.7 0.7 0.7]);
        set(D4,'Color',[0.7 0.7 0.7]);
        set(D5,'Color',[0.7 0.7 0.7]);
        set(D6,'Color',[0.9290 0.6940 0.1250]	);

        set(UD1D3,'Color',[0.9290 0.6940 0.1250]	);
        set(UD3D5,'Color',[0.9290 0.6940 0.1250]	);
        set(UD4D6,'Color',[0.7 0.7 0.7]);
        set(UD6D2,'Color',[0.7 0.7 0.7]);

        set(ul1mp,'Color',[0.9290 0.6940 0.1250]	);
        set(ul2mp,'Color',[0.9290 0.6940 0.1250]	);
        set(ul3mp,'Color',[0.7 0.7 0.7]);

        set(tr1,'Color',[0.9290 0.6940 0.1250]	);
        set(tr2,'Color',[0.9290 0.6940 0.1250]	);
        set(tr3,'Color',[0.7 0.7 0.7]);

        set(l12,'Color',[0.9290 0.6940 0.1250]	);
        set(l23,'Color',[0.7 0.7 0.7]);
       
    end


    if Iab_Star(i)>=0.1 && Ica_Star(i)<=-0.10
        set(D1,'Color',[0.9290 0.6940 0.1250]	);
        set(D2,'Color',[0.9290 0.6940 0.1250]	);
        set(D3,'Color',[0.7 0.7 0.7]);
        set(D4,'Color',[0.7 0.7 0.7]);
        set(D5,'Color',[0.7 0.7 0.7]);
        set(D6,'Color',[0.7 0.7 0.7]);
            
        set(UD1D3,'Color',[0.9290 0.6940 0.1250]	);
        set(UD3D5,'Color',[0.9290 0.6940 0.1250]	);
        set(UD4D6,'Color',[0.7 0.7 0.7]);
        set(UD6D2,'Color',[0.9290 0.6940 0.1250]	);%
        set(ul1mp,'Color',[0.9290 0.6940 0.1250]	);
        set(ul2mp,'Color',[0.7 0.7 0.7]);
        set(ul3mp,'Color',[0.9290 0.6940 0.1250]	);%
        set(tr1,'Color',[0.9290 0.6940 0.1250]	);
        set(tr2,'Color',[0.7 0.7 0.7]);
        set(tr3,'Color',[0.9290 0.6940 0.1250]	);
        set(l12,'Color',[0.9290 0.6940 0.1250]	);
        set(l23,'Color',[0.9290 0.6940 0.1250]	);
        
    end


    if Ibc_Star(i)>=0.1 && Iab_Star(i)<=-0.1
        set(D1,'Color',[0.7 0.7 0.7]);
        set(D2,'Color',[0.7 0.7 0.7]);
        set(D3,'Color',[0.9290 0.6940 0.1250]	);
        set(D4,'Color',[0.9290 0.6940 0.1250]	);
        set(D5,'Color',[0.7 0.7 0.7]);
        set(D6,'Color',[0.7 0.7 0.7]);
            
        set(UD1D3,'Color',[0.7 0.7 0.7]);
        set(UD3D5,'Color',[0.9290 0.6940 0.1250]	);
        set(UD4D6,'Color',[0.9290 0.6940 0.1250]	);
        set(UD6D2,'Color',[0.7 0.7 0.7]);%

        set(ul1mp,'Color',[0.9290 0.6940 0.1250]	);
        set(ul2mp,'Color',[0.9290 0.6940 0.1250]	);
        set(ul3mp,'Color',[0.7 0.7 0.7]);%

        set(tr1,'Color',[0.9290 0.6940 0.1250]	);
        set(tr2,'Color',[0.9290 0.6940 0.1250]	);
        set(tr3,'Color',[0.7 0.7 0.7]);

        set(l12,'Color',[0.9290 0.6940 0.1250]	);
        set(l23,'Color',[0.7 0.7 0.7]);
    end
    if Ibc_Star(i)>=0.1 && Ica_Star(i)<=-0.1
        set(D1,'Color',[0.7 0.7 0.7]);
        set(D2,'Color',[0.9290 0.6940 0.1250]	);
        set(D3,'Color',[0.9290 0.6940 0.1250]	);
        set(D4,'Color',[0.7 0.7 0.7]);
        set(D5,'Color',[0.7 0.7 0.7]);
        set(D6,'Color',[0.7 0.7 0.7]);

        set(UD1D3,'Color',[0.7 0.7 0.7]);
        set(UD3D5,'Color',[0.9290 0.6940 0.1250]	);
        set(UD4D6,'Color',[0.7 0.7 0.7]);
        set(UD6D2,'Color',[0.9290 0.6940 0.1250]	);%

        set(ul1mp,'Color',[0.7 0.7 0.7]);
        set(ul2mp,'Color',[0.9290 0.6940 0.1250]	);
        set(ul3mp,'Color',[0.9290 0.6940 0.1250]	);%

        set(tr1,'Color',[0.7 0.7 0.7]);
        set(tr2,'Color',[0.9290 0.6940 0.1250]	);
        set(tr3,'Color',[0.9290 0.6940 0.1250]	);

        set(l12,'Color',[0.7 0.7 0.7]);
        set(l23,'Color',[0.9290 0.6940 0.1250]	);
    end
    if Ica_Star(i)>=0.1 && Iab_Star(i)<=-0.1
        set(D1,'Color',[0.7 0.7 0.7]);
        set(D2,'Color',[0.7 0.7 0.7]);
        set(D3,'Color',[0.7 0.7 0.7]);
        set(D4,'Color',[0.9290 0.6940 0.1250]	);
        set(D5,'Color',[0.9290 0.6940 0.1250]	);
        set(D6,'Color',[0.7 0.7 0.7]);

        set(UD1D3,'Color',[0.7 0.7 0.7]);
        set(UD3D5,'Color',[0.7 0.7 0.7]);
        set(UD4D6,'Color',[0.9290 0.6940 0.1250]	);
        set(UD6D2,'Color',[0.7 0.7 0.7]);%

        set(ul1mp,'Color',[0.9290 0.6940 0.1250]	);
        set(ul2mp,'Color',[0.7 0.7 0.7]);
        set(ul3mp,'Color',[0.9290 0.6940 0.1250]	);%

        set(tr1,'Color',[0.9290 0.6940 0.1250]	);
        set(tr2,'Color',[0.7 0.7 0.7]);
        set(tr3,'Color',[0.9290 0.6940 0.1250]	);

        set(l12,'Color',[0.9290 0.6940 0.1250]	);
        set(l23,'Color',[0.9290 0.6940 0.1250]	);
    end
    if Ica_Star(i)>=0.1 && Ibc_Star(i)<=-0.1
        set(D1,'Color',[0.7 0.7 0.7]);
        set(D2,'Color',[0.7 0.7 0.7]);
        set(D3,'Color',[0.7 0.7 0.7]);
        set(D4,'Color',[0.7 0.7 0.7]);
        set(D5,'Color',[0.9290 0.6940 0.1250]	);
        set(D6,'Color',[0.9290 0.6940 0.1250]	);

        set(UD1D3,'Color',[0.7 0.7 0.7]);
        set(UD3D5,'Color',[0.7 0.7 0.7]);
        set(UD4D6,'Color',[0.7 0.7 0.7]);
        set(UD6D2,'Color',[0.7 0.7 0.7]);%

        set(ul1mp,'Color',[0.7 0.7 0.7]);
        set(ul2mp,'Color',[0.9290 0.6940 0.1250]	);
        set(ul3mp,'Color',[0.9290 0.6940 0.1250]	);%

        set(tr1,'Color',[0.7 0.7 0.7]);
        set(tr2,'Color',[0.9290 0.6940 0.1250]	);
        set(tr3,'Color',[0.9290 0.6940 0.1250]	);

        set(l12,'Color',[0.7 0.7 0.7]);
        set(l23,'Color',[0.9290 0.6940 0.1250]	);
    end

 if Iab_Delta(i)>=0.1 && Ibc_Delta(i)<=-0.1
        set(D1L,'Color',[0.9290 0.6940 0.1250]	);
        set(D2L,'Color',[0.7 0.7 0.7]);
        set(D3L,'Color',[0.7 0.7 0.7]);
        set(D4L,'Color',[0.7 0.7 0.7]);
        set(D5L,'Color',[0.7 0.7 0.7]);
        set(D6L,'Color',[0.9290 0.6940 0.1250]	);
        
        set(BD1D3,'Color',[0.9290 0.6940 0.1250]	);
        set(BD3D5,'Color',[0.7 0.7 0.7]);
        set(BD4D6,'Color',[0.7 0.7 0.7]);
        set(BD6D2,'Color',[0.9290 0.6940 0.1250]	);%

        set(bl1mp,'Color',[0.9290 0.6940 0.1250]	);
        set(bl2mp,'Color',[0.9290 0.6940 0.1250]	);
        set(bl3mp,'Color',[0.7 0.7 0.7]);%

        set(trb1,'Color',[0.7 0.7 0.7]);
        set(trb2,'Color',[0.9290 0.6940 0.1250]	);
        set(trb3,'Color',[0.7 0.7 0.7]);

        set(lji2i1,'Color',[0.7 0.7 0.7]);
        set(lji3i2,'Color',[0.9290 0.6940 0.1250]	);
        set(lji1i3,'Color',[0.7 0.7 0.7]);
    end


    if Iab_Delta(i)>=0.1 && Ica_Delta(i)<=-0.10
         set(D1L,'Color',[0.9290 0.6940 0.1250]	);
        set(D2L,'Color',[0.9290 0.6940 0.1250]	);
        set(D3L,'Color',[0.7 0.7 0.7]);
        set(D4L,'Color',[0.7 0.7 0.7]);
        set(D5L,'Color',[0.7 0.7 0.7]);
        set(D6L,'Color',[0.7 0.7 0.7]);

        set(BD1D3,'Color',[0.9290 0.6940 0.1250]	);
        set(BD3D5,'Color',[0.7 0.7 0.7]);
        set(BD4D6,'Color',[0.7 0.7 0.7]);
        set(BD6D2,'Color',[0.7 0.7 0.7]);%

        set(bl1mp,'Color',[0.9290 0.6940 0.1250]	);
        set(bl2mp,'Color',[0.7 0.7 0.7]);
        set(bl3mp,'Color',[0.9290 0.6940 0.1250]	);%

        set(trb1,'Color',[0.7 0.7 0.7]);
        set(trb2,'Color',[0.7 0.7 0.7]);
        set(trb3,'Color',[0.9290 0.6940 0.1250]	);

        set(lji2i1,'Color',[0.7 0.7 0.7]);
        set(lji3i2,'Color',[0.7 0.7 0.7]);
        set(lji1i3,'Color',[0.9290 0.6940 0.1250]	);
    end
    if Ibc_Delta(i)>=0.1 && Iab_Delta(i)<=-0.1
        set(D1L,'Color',[0.7 0.7 0.7]);
        set(D2L,'Color',[0.7 0.7 0.7]);
        set(D3L,'Color',[0.9290 0.6940 0.1250]	);
        set(D4L,'Color',[0.9290 0.6940 0.1250]	);
        set(D5L,'Color',[0.7 0.7 0.7]);
        set(D6L,'Color',[0.7 0.7 0.7]);

        set(BD1D3,'Color',[0.7 0.7 0.7]);
        set(BD3D5,'Color',[0.7 0.7 0.7]);
        set(BD4D6,'Color',[0.9290 0.6940 0.1250]	);
        set(BD6D2,'Color',[0.9290 0.6940 0.1250]	);%

        set(bl1mp,'Color',[0.9290 0.6940 0.1250]	);
        set(bl2mp,'Color',[0.9290 0.6940 0.1250]	);
        set(bl3mp,'Color',[0.7 0.7 0.7]);%

        set(trb1,'Color',[0.7 0.7 0.7]);
        set(trb2,'Color',[0.9290 0.6940 0.1250]	);
        set(trb3,'Color',[0.7 0.7 0.7]);

        set(lji2i1,'Color',[0.7 0.7 0.7]);
        set(lji3i2,'Color',[0.9290 0.6940 0.1250]	);
        set(lji1i3,'Color',[0.7 0.7 0.7]);
    end
    if Ibc_Delta(i)>=0.1 && Ica_Delta(i)<=-0.1
        set(D1L,'Color',[0.7 0.7 0.7]);
        set(D2L,'Color',[0.9290 0.6940 0.1250]	);
        set(D3L,'Color',[0.9290 0.6940 0.1250]	);
        set(D4L,'Color',[0.7 0.7 0.7]);
        set(D5L,'Color',[0.7 0.7 0.7]);
        set(D6L,'Color',[0.7 0.7 0.7]);

        set(BD1D3,'Color',[0.7 0.7 0.7]);
        set(BD3D5,'Color',[0.7 0.7 0.7]);
        set(BD4D6,'Color',[0.7 0.7 0.7]);
        set(BD6D2,'Color',[0.7 0.7 0.7]);%

        set(bl1mp,'Color',[0.7 0.7 0.7]);
        set(bl2mp,'Color',[0.9290 0.6940 0.1250]	);
        set(bl3mp,'Color',[0.9290 0.6940 0.1250]	);%

        set(trb1,'Color',[0.9290 0.6940 0.1250]	);
        set(trb2,'Color',[0.7 0.7 0.7]);
        set(trb3,'Color',[0.7 0.7 0.7]);

        set(lji2i1,'Color',[0.9290 0.6940 0.1250]	);
        set(lji3i2,'Color',[0.7 0.7 0.7]);
        set(lji1i3,'Color',[0.7 0.7 0.7]);
    end
    if Ica_Delta(i)>=0.1 && Iab_Delta(i)<=-0.1
         set(D1L,'Color',[0.7 0.7 0.7]);
        set(D2L,'Color',[0.7 0.7 0.7]);
        set(D3L,'Color',[0.7 0.7 0.7]);
        set(D4L,'Color',[0.9290 0.6940 0.1250]	);
        set(D5L,'Color',[0.9290 0.6940 0.1250]	);
        set(D6L,'Color',[0.7 0.7 0.7]);

        set(BD1D3,'Color',[0.7 0.7 0.7]);
        set(BD3D5,'Color',[0.9290 0.6940 0.1250]	);
        set(BD4D6,'Color',[0.9290 0.6940 0.1250]	);
        set(BD6D2,'Color',[0.9290 0.6940 0.1250]	);%

        set(bl1mp,'Color',[0.9290 0.6940 0.1250]	);
        set(bl2mp,'Color',[0.7 0.7 0.7]);
        set(bl3mp,'Color',[0.9290 0.6940 0.1250]	);%

        set(trb1,'Color',[0.7 0.7 0.7]);
        set(trb2,'Color',[0.7 0.7 0.7]);
        set(trb3,'Color',[0.9290 0.6940 0.1250]	);

        set(lji2i1,'Color',[0.7 0.7 0.7]);
        set(lji3i2,'Color',[0.7 0.7 0.7]);
        set(lji1i3,'Color',[0.9290 0.6940 0.1250]	);
    end
    if Ica_Delta(i)>=0.1 && Ibc_Delta(i)<=-0.1
        set(D1L,'Color',[0.7 0.7 0.7]);
        set(D2L,'Color',[0.7 0.7 0.7]);
        set(D3L,'Color',[0.7 0.7 0.7]);
        set(D4L,'Color',[0.7 0.7 0.7]);
        set(D5L,'Color',[0.9290 0.6940 0.1250]	);
        set(D6L,'Color',[0.9290 0.6940 0.1250]	);

        set(BD1D3,'Color',[0.7 0.7 0.7]);
        set(BD3D5,'Color',[0.9290 0.6940 0.1250]	);
        set(BD4D6,'Color',[0.7 0.7 0.7]);
        set(BD6D2,'Color',[0.9290 0.6940 0.1250]	);%

        set(bl1mp,'Color',[0.7 0.7 0.7]);
        set(bl2mp,'Color',[0.9290 0.6940 0.1250]	);
        set(bl3mp,'Color',[0.9290 0.6940 0.1250]	);%

        set(trb1,'Color',[0.9290 0.6940 0.1250]	);
        set(trb2,'Color',[0.7 0.7 0.7]);
        set(trb3,'Color',[0.7 0.7 0.7]);

        set(lji2i1,'Color',[0.9290 0.6940 0.1250]	);
        set(lji3i2,'Color',[0.7 0.7 0.7]);
        set(lji1i3,'Color',[0.7 0.7 0.7]);
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    drawnow
    currFrame=getframe(hf);
    writeVideo(writeobj, currFrame)
    set(txt_vrn, 'String', ' ');
    set(txt_vyn, 'String', ' ');
    set(txt_vbn, 'String', ' ');
    set(txt_vry, 'String', ' ');
    set(txt_vyb, 'String', ' ');
    set(txt_vbr, 'String', ' ');
    set(txt_vrys, 'String', ' ');
    set(txt_vybs, 'String', ' ');
    set(txt_vbrs, 'String', ' ');
    set(txt_vry_Delta, 'String', ' ');
    set(txt_vyb_Delta, 'String', ' ');
    set(txt_vbr_Delta, 'String', ' ');
    set(txt_vrys_Delta, 'String', ' ');
    set(txt_vybs_Delta, 'String', ' ');
    set(txt_vbrs_Delta, 'String', ' ');

end
 end

% ah1 = axes('position',get(gca,'position'),'visible','off');
% lgd1 = legend(ah1,[Vry_Delta_Wave Vyb_Delta_Wave Vbr_Delta_Wave Vry_Delta_Wave_Shifted Vyb_Delta_Wave_Shifted Vbr_Delta_Wave_Shifted], 'V_{RY}', 'V_{YB}', 'V_{BR}','V_{YR}','V_{BY}','V_{RB}','Location','southeast','Orientation','horizontal','TextColor','blue','FontSize',10,'AutoUpdate','off');
% lgd1.NumColumns = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(txt_vrn, 'String', 'V_{Rn}');
    set(txt_vyn, 'String', 'V_{Yn}');
    set(txt_vbn, 'String', 'V_{Bn}');
    set(txt_vry, 'String', 'V_{RY}');
    set(txt_vyb, 'String', 'V_{YB}');
    set(txt_vbr, 'String', 'V_{BR}');
    set(txt_vrys, 'String', 'V_{YR}');
    set(txt_vybs, 'String', 'V_{BY}');
    set(txt_vbrs, 'String', 'V_{RB}');
    set(txt_vry_Delta, 'String', 'V_{RY}');
    set(txt_vyb_Delta, 'String', 'V_{YB}');
    set(txt_vbr_Delta, 'String', 'V_{BR}');
    set(txt_vrys_Delta, 'String', 'V_{YR}');
    set(txt_vybs_Delta, 'String', 'V_{BY}');
    set(txt_vbrs_Delta, 'String', 'V_{RB}');
 close(writeobj)
