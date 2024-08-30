function dz = threeDisk_threePend(t,z)
dz = zeros(12,1);
dz(1:6) = z(7:12);

F = 0; alpha = pi*sin(0.2*t);

global ma mb mc m1 m2 m3 L1 L2 L3 ra rb rc R g b k1 k2
q1 = z(1); q2 = z(2); q3 = z(3); q4 = z(4); q5 = z(5); q6 = z(6); 
dq1 = z(7); dq2 = z(8); dq3 = z(9); dq4 = z(10); dq5 = z(11); dq6 = z(12); 

Q = [0, F*cos(alpha - q2)*(R - rb), 0, F*L1*cos(alpha - q4), F*L2*cos(alpha - q5), F*b*cos(alpha - q6)]';
 
Mass =[(3*ma*(R - ra)^2)/2,                                               0,                   0,                                               0,                                        0,                        0;
                  0,      ((R - rb)^2*(2*m1 + 2*m2 + 2*m3 + 3*mb))/2,                   0, (L1*cos(q2 - q4)*(R - rb)*(m1 + 2*m2 + 2*m3))/2, (L2*cos(q2 - q5)*(R - rb)*(m2 + 2*m3))/2, (L3*m3*cos(q2 - q6)*(R - rb))/2;
                  0,                                               0, (3*mc*(R - rc)^2)/2,                                               0,                                        0,                               0;
                  0, (L1*cos(q2 - q4)*(R - rb)*(m1 + 2*m2 + 2*m3))/2,                   0,                     (L1^2*(m1 + 3*m2 + 3*m3))/3,       (L1*L2*cos(q4 - q5)*(m2 + 2*m3))/2,       (L1*L3*m3*cos(q4 - q6))/2;
                  0,        (L2*cos(q2 - q5)*(R - rb)*(m2 + 2*m3))/2,                   0,              (L1*L2*cos(q4 - q5)*(m2 + 2*m3))/2,                     (L2^2*(m2 + 3*m3))/3,       (L2*L3*m3*cos(q5 - q6))/2;
                  0,                 (L3*m3*cos(q2 - q6)*(R - rb))/2,                   0,                       (L1*L3*m3*cos(q4 - q6))/2,                (L2*L3*m3*cos(q5 - q6))/2,                     (L3^2*m3)/3];
 

 
B =[
 
                                                                                                                                                                                                                                                                                                                                (k1*((q1*cos(q1)*(R - ra) - q2*cos(q2)*(R - rb))*(2*R - 2*ra)*(cos(q1) - q1*sin(q1)) + (q1*sin(q1)*(R - ra) - q2*sin(q2)*(R - rb))*(2*R - 2*ra)*(sin(q1) + q1*cos(q1))))/2 + g*ma*sin(q1)*(R - ra)
((R - rb)*(2*g*m1*sin(q2) + 2*g*m2*sin(q2) + 2*g*m3*sin(q2) + 2*g*mb*sin(q2) + 2*R*k1*q2 + 2*R*k2*q2 - 2*k1*q2*rb - 2*k2*q2*rb - 2*R*k1*q1*cos(q1 - q2) - 2*R*k2*q3*cos(q2 - q3) + 2*k1*q1*ra*cos(q1 - q2) + 2*k2*q3*rc*cos(q2 - q3) + L1*dq4^2*m1*sin(q2 - q4) + 2*L1*dq4^2*m2*sin(q2 - q4) + 2*L1*dq4^2*m3*sin(q2 - q4) + L2*dq5^2*m2*sin(q2 - q5) + 2*L2*dq5^2*m3*sin(q2 - q5) + L3*dq6^2*m3*sin(q2 - q6) - 2*R*k1*q1*q2*sin(q1 - q2) + 2*R*k2*q2*q3*sin(q2 - q3) + 2*k1*q1*q2*ra*sin(q1 - q2) - 2*k2*q2*q3*rc*sin(q2 - q3)))/2
                                                                                                                                                                                                                                                                                                                                g*mc*sin(q3)*(R - rc) - (k2*((q2*cos(q2)*(R - rb) - q3*cos(q3)*(R - rc))*(2*R - 2*rc)*(cos(q3) - q3*sin(q3)) + (q2*sin(q2)*(R - rb) - q3*sin(q3)*(R - rc))*(2*R - 2*rc)*(sin(q3) + q3*cos(q3))))/2
                                                                                                                                                                                                                 (L1*(g*m1*sin(q4) + 2*g*m2*sin(q4) + 2*g*m3*sin(q4) + dq2^2*m1*rb*sin(q2 - q4) + 2*dq2^2*m2*rb*sin(q2 - q4) + 2*dq2^2*m3*rb*sin(q2 - q4) + L2*dq5^2*m2*sin(q4 - q5) + 2*L2*dq5^2*m3*sin(q4 - q5) + L3*dq6^2*m3*sin(q4 - q6) - R*dq2^2*m1*sin(q2 - q4) - 2*R*dq2^2*m2*sin(q2 - q4) - 2*R*dq2^2*m3*sin(q2 - q4)))/2
                                                                                                                                                                                                                                                                                           (L2*(g*m2*sin(q5) + 2*g*m3*sin(q5) + dq2^2*m2*rb*sin(q2 - q5) + 2*dq2^2*m3*rb*sin(q2 - q5) - L1*dq4^2*m2*sin(q4 - q5) - 2*L1*dq4^2*m3*sin(q4 - q5) + L3*dq6^2*m3*sin(q5 - q6) - R*dq2^2*m2*sin(q2 - q5) - 2*R*dq2^2*m3*sin(q2 - q5)))/2
                                                                                                                                                                                                                                                                                                                                                                                                             -(L3*m3*(L1*dq4^2*sin(q4 - q6) - g*sin(q6) + L2*dq5^2*sin(q5 - q6) + R*dq2^2*sin(q2 - q6) - dq2^2*rb*sin(q2 - q6)))/2];
 

dz(7:12) = Mass\(Q-B);