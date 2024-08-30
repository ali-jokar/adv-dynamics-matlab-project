clc, clear all

syms t 
syms q1 q2 q3 q4 q5 q6
syms dq1 dq2 dq3 dq4 dq5 dq6
syms ma mb mc m1 m2 m3 L1 L2 L3 k1 k2 ra rb rc R g b alpha F

Ta = 1/2*3/2*ma*(R-ra)^2*dq1^2;
Tb = 1/2*3/2*mb*(R-rb)^2*dq2^2;
Tc = 1/2*3/2*mc*(R-rc)^2*dq3^2;

Vb = (R-rb)*dq2*[cos(q2);sin(q2);0];
Omega1 = dq4*[0;0;1];
rG1_b = L1/2*[sin(q4);-cos(q4);0];
VG1 = Vb + cross(Omega1,rG1_b);
IG1 = 1/12*m1*L1^2;
Tp1 = 1/2 *m1*VG1.'*VG1 + 1/2 *IG1*dq4^2;

Omega2 = dq5*[0;0;1];
rd_G1 = L1/2*[sin(q4);-cos(q4);0];
rG2_d = L2/2*[sin(q5);-cos(q5);0];
VG2 = VG1 + cross(Omega1,rd_G1) + cross(Omega2,rG2_d);
IG2 = 1/12*m2*L2^2;
Tp2 = 1/2 *m2*VG2.'*VG2 + 1/2 *IG2*dq5^2;

Omega3 = dq6*[0;0;1];
re_G2 = L2/2*[sin(q5);-cos(q5);0];
rG3_e = L3/2*[sin(q6);-cos(q6);0];
VG3 = VG2 + cross(Omega2,re_G2) + cross(Omega3,rG3_e);
IG3 = 1/12*m3*L3^2;
Tp3 = 1/2 *m3*VG3.'*VG3 + 1/2 *IG3*dq6^2;

% Kinetic Energy
T = Ta + Tb + Tc + Tp1 + Tp2 + Tp3;

% Potential Energy
ya = -(R-ra)*cos(q1);
yb = -(R-rb)*cos(q2);
yc = -(R-rc)*cos(q3);
yG1 = yb - L1/2*cos(q4);
yG2 = yG1 - L1/2*cos(q4) - L2/2*cos(q5);
yG3 = yG2 - L2/2*cos(q5) - L3/2*cos(q6);

V = ma*g*ya + mb*g*yb + mc*g*yc + m1*g*yG1 + m2*g*yG2 + m3*g*yG3 + 0.5*k1*(((R-ra)*q1*cos(q1)-(R-rb)*q2*cos(q2))^2 + ((R-ra)*q1*sin(q1)-(R-rb)*q2*sin(q2))^2 ) + 0.5*k2*(((R-rb)*q2*cos(q2)-(R-rc)*q3*cos(q3))^2 + ((R-rb)*q2*sin(q2)-(R-rc)*q3*sin(q3))^2 );

Lagrangian = T - V;


q = [q1;q2;q3;q4;q5;q6];
dq = [dq1;dq2;dq3;dq4;dq5;dq6];


% Generalized Force
Force = F*[cos(alpha);sin(alpha);0];
rEnd_G3 = (b-L3/2)*[sin(q6);-cos(q6);0];
VEnd = VG3 + cross(Omega3,rEnd_G3);
Power = Force.'*VEnd;
Q = simplify(jacobian(Power,dq))

dL_dq = jacobian(Lagrangian,q);
dL_ddq = jacobian(Lagrangian,dq);

% Mass Matrix
Mass = simplify(jacobian(dL_ddq,dq))
B = simplify(jacobian(dL_ddq,q)*dq + diff(dL_ddq.',t) - dL_dq.')

E = simplify(T + V)
