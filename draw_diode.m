function [D1x1ref, D1x2ref, D1y1ref, D1y2ref, D1txref,D1tyref,D2x1ref, D2x2ref, D2y1ref, D2y2ref, D2txref,D2tyref]= draw_diode(xref,yref,mp)

D1rx1=[xref xref+1 xref+0.5 xref xref+0.5 xref+1 xref xref+1 xref+0.5 xref+0.5];
D1ry1=[yref yref yref yref-1 yref yref-1 yref-1 yref-1 yref yref+2];
D1rx2=[xref+0.5 xref+0.5];
D1ry2=[yref-mp yref-1];
D1txtxref=xref-1.225;
D1txtyref=yref-0.025;

%D2rx1=[xref xref+1 xref+0.5 xref xref+0.5 xref+1 xref xref+1 xref+0.5 xref+0.5];
%D2ry1=[yref-5 yref-5 yref-5 yref-1-5 yref-5 yref-1-5 yref-1-5 yref-1-5 yref-5 yref+2-5];
D2rx1=[xref xref+1 xref+0.5 xref xref+0.5 xref+1 xref xref+1 xref+0.5 xref+0.5];
D2ry1=[yref-5 yref-5 yref-5 yref-1-5 yref-5 yref-1-5 yref-1-5 yref-1-5 yref-5 yref-mp];
D2rx2=[xref+0.5 xref+0.5];
D2ry2=[yref-1-5-2 yref-1-5];
D2txtxref=xref-1.225;
D2txtyref=yref-0.025-5;


D1x1ref=D1rx1;
D1x2ref=D1rx2;
D1y1ref=D1ry1;
D1y2ref=D1ry2;
D1txref=D1rx1(:,1)-0.45;
D1tyref=D1ry1(:,1)+0.75;

D2x1ref=D2rx1;
D2x2ref=D2rx2;
D2y1ref=D2ry1;
D2y2ref=D2ry2;
% D2txref=D2txtxref;
% D2tyref=D2txtyref;
D2txref=D2rx1(:,1)-0.45;;
D2tyref=D2ry1(:,1)+0.5;

