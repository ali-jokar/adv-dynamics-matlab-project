clc, clear all, close all

global ma mb mc m1 m2 m3 L1 L2 L3 ra rb rc R g b k1 k2

ma = 2; mb = 3; mc = 1; m1 = 0.4; m2 = 0.3; m3 = 0.3; 
L1 = 0.7; L2 = 0.6; L3 = 0.6; ra = 0.4; rb = 0.5; rc = 0.3;
R = 3; g = 9.81; b = 0.8*L3; k1=2; k2=1;

options = odeset('maxstep',0.001);
z0 = [-pi/6;pi/3;pi/2;0.1;0.1;0.1;0.2;0.2;0.2;0;0;0];
[t,z] = ode45(@threeDisk_threePend,[0:0.01:10],z0,options);

figure
subplot(2,3,1);
plot(t,z(:,1))
xlabel("t(s)");
ylabel("q1(rad)");

subplot(2,3,2);
plot(t,z(:,2))
xlabel("t(s)");
ylabel("q2(rad)");

subplot(2,3,3);
plot(t,z(:,3))
xlabel("t(s)");
ylabel("q3(rad)");

subplot(2,3,4);
plot(t,z(:,4))
xlabel("t(s)");
ylabel("q4(rad)");

subplot(2,3,5);
plot(t,z(:,5))
xlabel("t(s)");
ylabel("q5(rad)");

subplot(2,3,6);
plot(t,z(:,6))
xlabel("t(s)");
ylabel("q6(rad)");

figure
subplot(2,3,1);
plot(t,z(:,7))
xlabel("t(s)");
ylabel("qdot_1(rad/s)");

subplot(2,3,2);
plot(t,z(:,8))
xlabel("t(s)");
ylabel("qdot_2(rad/s)");

subplot(2,3,3);
plot(t,z(:,9))
xlabel("t(s)");
ylabel("qdot_3(rad/s)");

subplot(2,3,4);
plot(t,z(:,10))
xlabel("t(s)");
ylabel("qdot_4(rad/s)");

subplot(2,3,5);
plot(t,z(:,11))
xlabel("t(s)");
ylabel("qdot_5(rad/s)");

subplot(2,3,6);
plot(t,z(:,12))
xlabel("t(s)");
ylabel("qdot_6(rad/s)");

E = t;
for i = 1:length(t)
    q1 = z(i,1); q2 = z(i,2); q3 = z(i,3); q4 = z(i,4); q5 = z(i,5); q6 = z(i,6); dq1 = z(i,7); dq2 = z(i,8); dq3 = z(i,9); dq4 = z(i,10); dq5 = z(i,11); dq6 = z(i,12);    
    E(i) =(m3*(dq2*cos(q2)*(R - rb) + L1*dq4*cos(q4) + L2*dq5*cos(q5) + (L3*dq6*cos(q6))/2)^2)/2 + (k1*((q1*cos(q1)*(R - ra) - q2*cos(q2)*(R - rb))^2 + (q1*sin(q1)*(R - ra) - q2*sin(q2)*(R - rb))^2))/2 + (k2*((q2*cos(q2)*(R - rb) - q3*cos(q3)*(R - rc))^2 + (q2*sin(q2)*(R - rb) - q3*sin(q3)*(R - rc))^2))/2 + (m2*(dq2*cos(q2)*(R - rb) + L1*dq4*cos(q4) + (L2*dq5*cos(q5))/2)^2)/2 + (m1*(dq2*cos(q2)*(R - rb) + (L1*dq4*cos(q4))/2)^2)/2 + (m3*(dq2*sin(q2)*(R - rb) + L1*dq4*sin(q4) + L2*dq5*sin(q5) + (L3*dq6*sin(q6))/2)^2)/2 + (m2*(dq2*sin(q2)*(R - rb) + L1*dq4*sin(q4) + (L2*dq5*sin(q5))/2)^2)/2 + (m1*(dq2*sin(q2)*(R - rb) + (L1*dq4*sin(q4))/2)^2)/2 + (L1^2*dq4^2*m1)/24 + (L2^2*dq5^2*m2)/24 + (L3^2*dq6^2*m3)/24 - g*m2*(L1*cos(q4) + (L2*cos(q5))/2 + cos(q2)*(R - rb)) + (3*dq1^2*ma*(R - ra)^2)/4 + (3*dq2^2*mb*(R - rb)^2)/4 + (3*dq3^2*mc*(R - rc)^2)/4 - g*m3*(L1*cos(q4) + L2*cos(q5) + (L3*cos(q6))/2 + cos(q2)*(R - rb)) - g*m1*((L1*cos(q4))/2 + cos(q2)*(R - rb)) - g*ma*cos(q1)*(R - ra) - g*mb*cos(q2)*(R - rb) - g*mc*cos(q3)*(R - rc);
 end
figure
plot(t,(E-E(1))/E(1)*100)
xlabel("t(s)");
ylabel("change of Energy");

q1 = z(:,1); q2 = z(:,2); q3 = z(:,3); q4 = z(:,4); q5 = z(:,5); q6 = z(:,6);
dq1 = z(:,7); dq2 = z(:,8); dq3 = z(:,9); dq4 = z(:,10); dq5 = z(:,11); dq6 = z(:,12);

figure
hold on

tt = 0:.1:2*pi+0.06;
xa=ra*cos(tt); ya=ra*sin(tt);
xb=rb*cos(tt); yb=rb*sin(tt);
xc=rc*cos(tt); yc=rc*sin(tt);

ttt = pi:0.01:2*pi;
xcy = R*cos(ttt); ycy = R*sin(ttt);

for i=1:length(t)
    Xa(i) = (R-ra)*sin(q1(i));
    Ya(i) = -(R-ra)*cos(q1(i));

    Xb(i) = (R-rb)*sin(q2(i));
    Yb(i) = -(R-rb)*cos(q2(i));

    Xc(i) = (R-rc)*sin(q3(i));
    Yc(i) = -(R-rc)*cos(q3(i));

    Xd(i) = Xb(i) + L1 *sin(q4(i));
    Yd(i) = Yb(i) - L1 *cos(q4(i));

    Xe(i) = Xd(i) + L2 * sin(q5(i));
    Ye(i) = Yd(i) - L2 * cos(q5(i));

    Xf(i) = Xe(i) + L2 * sin(q6(i));
    Yf(i) = Ye(i) - L2 * cos(q6(i));

    V1x=[Xb(i);Xd(i)];  V1y=[Yb(i);Yd(i)];
    V2x=[Xd(i);Xe(i)];  V2y=[Yd(i);Ye(i)];
    V3x=[Xe(i);Xf(i)];  V3y=[Ye(i);Yf(i)];

    V4x=[Xa(i);Xa(i)+ra*sin((R-ra)/ra*q1(i))];  V4y=[Ya(i);Ya(i)+ra*cos((R-ra)/ra*q1(i))];

    V5x=[Xb(i);Xb(i)+rb*sin((R-rb)/rb*q2(i))];  V5y=[Yb(i);Yb(i)+rb*cos((R-rb)/rb*q2(i))];

    V6x=[Xc(i);Xc(i)+rc*sin((R-rc)/rc*q3(i))];  V6y=[Yc(i);Yc(i)+rc*cos((R-rc)/rc*q3(i))];

    V7x=[Xa(i);Xb(i)];  V7y=[Ya(i);Yb(i)];
    V8x=[Xb(i);Xc(i)];  V8y=[Yb(i);Yc(i)];

    plot(xcy,ycy,'k-','linewidth',4)

    hold on

    plot(Xa(i)+xa,Ya(i)+ya,'r-','linewidth',4)
    plot(Xb(i)+xb,Yb(i)+yb,'r-','linewidth',4)
    plot(Xc(i)+xc,Yc(i)+yc,'r-','linewidth',4)

    plot(V1x,V1y,'b-','linewidth',4)

    plot(V2x,V2y,'g-','linewidth',4)

    plot(V3x,V3y,'m-','linewidth',4)

    plot(V4x,V4y,'r-','linewidth',4)

    plot(V5x,V5y,'r-','linewidth',4)

    plot(V6x,V6y,'r-','linewidth',4)

    plot(V7x,V7y,'c-','linewidth',4)

    plot(V8x,V8y,'c-','linewidth',4)

    plot(Xa(i),Ya(i),'yO','linewidth',6)

    plot(Xb(i),Yb(i),'yO','linewidth',6)

    plot(Xc(i),Yc(i),'yO','linewidth',6)

    plot(Xd(i),Yd(i),'yO','linewidth',6) 

    plot(Xe(i),Ye(i),'yO','linewidth',6) 
    
    str=['Time = ',num2str(t(i)),' s'];
    text(-R/5,-R/3,str,'fontsize',20,'fontweight','bold')
    axis equal 
    axis([-R*1.3 R*1.3 -R*1.3 0.3*R])
    hold off
    pause(0.01)
end

