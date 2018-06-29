function [c0, c1, c2, c3, c4] = ganiso (alpha_g, beta_g, gamma_g, we) 
%we=28.025*9.4*1e9;
gx=we*(2.00614/2);
gy=we*(2.00194/2);
gz=we*(2.00988/2);

ca = cosd(alpha_g);
cb = cosd(beta_g);
cg = cosd(gamma_g);
sa = sind(alpha_g);
sb = sind(beta_g);
sg = sind(gamma_g);

r11 = ca*cb*cg-sa*sg;
r12 = sa*cb*cg+ca*sg; 
r13 = -sb*cg;
r21 = -ca*cb*sg-sa*cg; 
r22 = -sa*cb*sg+ca*cg; 
r23 = sb*sg;
r31 = ca*sb; 
r32 = sa*sb;
r33 = cb;
c0 = 1/3*(gx+gy+gz);
c1 = 2*sqrt(2)/3*(gx*r11*r31+gy*r12*r32+gz*r13*r33);
c2 = 2*sqrt(2)/3*(gx*r21*r31+gy*r22*r32+gz*r23*r33);
c3 = 1/3*(gx*(r11^2-r21^2)+gy*(r12^2-r22^2)+gz*(r13^2-r23^2));
c4 = 2/3*(gx*r11*r21+gy*r22*r12+gz*r13*r23);
end