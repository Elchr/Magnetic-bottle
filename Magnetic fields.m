clc; clear; close all;

B0=10^(-6);
q=1.6*10^(-19);
m=1.67*10^(-27);
%m=9.1*10^(-31);
RE=6371000;
z0=3*10^7; %5*RE
E=100*10^6*1.6*10^(-19);
u=sqrt(2*E/m);
r=(2*m*u)/(B0*q);
rg=(m*u)/(B0*q);            
w=q*B0/m;
wmirr_th=(w*rg)/(2*z0)
T=abs(2*pi/w); 
dt=T/1000;
n=30;
total=n*T;
nstep=fix(total/dt);

xold=[r,0,0,0,0,m*u];
res=zeros(nstep,6);
x=zeros(nstep,1);
y=zeros(nstep,1);
time=zeros(nstep,1);


for i=1:1:nstep
   time(i)=i*dt;
   xnew(4)=xold(4)+((xold(5)^2-q^2*xold(1)^4*B0^2/4*exp(2*(xold(3)/z0)^2))/(m*xold(1)^3))*dt;
   xnew(5)=xold(5);
   xnew(6)=xold(6)+(q*B0*xold(3)*(xold(5)-q*xold(1)^2*B0*exp((xold(3)/z0)^2)/2)*exp((xold(3)/z0)^2)/(m*z0^2))*dt;
   xnew(1)=xold(1)+xnew(4)*dt/m;
   xnew(3)=xold(3)+xnew(6)*dt/m;
   xnew(2)=xold(2)+(xnew(5)-q*xnew(1)^2*B0/2*exp((xnew(3)/z0)^2))*dt/(m*xnew(1)^2);
   res(i,:)=[xnew(1),xnew(2),xnew(3),xnew(4),xnew(5),xnew(6)];
   
   x(i,1)=res(i,1)*cos(res(i,2));
   y(i,1)=res(i,1)*sin(res(i,2));
   
   xold=xnew;
end
rho_c=zeros(1000,n);
rho_m=zeros(n,1);
z_c=zeros(1000,n);
z_m=zeros(n,1);
phi_m=zeros(n,1);
x_m=zeros(n,1);
y_m=zeros(n,1);
rr=zeros(n,1);
aa=1;
bb=1000;
t_m=zeros(n,1);
for h=1:1:n
    
z_m(h)=sum(res(aa:bb,3))/1000;
phi_m(h)=sum(res(aa:bb,2))/1000;
rho_m(h)=sum(res(aa:bb,1))/1000;
x_m(h)=sum(x(aa:bb,1))/1000;
y_m(h)=sum(y(aa:bb,1))/1000;

rr(h)=sqrt(x_m(h)^2+y_m(h)^2);
%x_n(h)=rho_m(h)*cos(phi_m(h));
%y_n(h)=rho_m(h)*sin(phi_m(h));
t_m(h)=sum(time(aa:bb,1))/1000;

aa=aa+1000;
bb=bb+1000;
end
rho_m2=sum(res(1:n,1))/(n*1000);

figure(1)
plot(res(:,1),res(:,3),'.k')
xlabel('ñ (m)')
ylabel('z (m)')
title('ÄéÜãñáììá óôï åðßðåäï ñ-z ')
grid on
%hold on
%plot(rho_m,z_m,'.r')
%grid on
%hold off


figure(2)
plot3(res(:,1),res(:,2),res(:,3),'.k')
xlabel('ñ (m)')
ylabel('ö (deg)')
zlabel('z (m)')
title('Cylindrical coordinates - Position')
grid on

figure(3)
plot3(res(:,4),res(:,5),res(:,6),'.k')
xlabel('pñ (kg m/sec)')
ylabel('pö (kg m^2/sec)')
zlabel('pz (kg m/sec)')
title('Cylindrical coordinates - Momentum')
grid on

figure(4)
plot3(x(:,1),y(:,1),res(:,3),'.k')
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
title('Cartesian coordinates - Position')
grid on
hold on
plot3(x_m,y_m,z_m,'.r')
hold off

figure(5)
plot(time,res(:,3),'.k')
xlabel('t(sec)')
ylabel('z(m)')
title('Motion in magnetic bottle')
grid on

figure(6)
plot(abs(rho_m),z_m,'.k')
xlabel('<ñ>(m)')
ylabel('<z>(m)')
title('Average values per period - Cylindrical')
grid on

figure(7)
plot(rr,z_m,'.k')
xlabel('<ñ>(m)')
ylabel('<z>(m)')
title('Average values per period - Cylindrical')
grid on
figure(8)

plot(x_m,z_m,'.k')
xlabel('<x>(m)')
ylabel('<z>(m)')
title('Average values per period - Cartesian')
grid on

figure(9)
plot3(x_m,y_m,z_m,'.k')
grid on
xlabel('<x>(m)')
ylabel('<y>(m)')
zlabel('<z>(m)')
title('Average values per period - Cartesian')

logos=rg/(2*z0);
s=(4*5^6-3*5^5)^(-1/4)
test=4*5*RE*(1.3-0.56*s)/u

figure(10)
plot(t_m,z_m,'.k')
xlabel('t(sec)')
ylabel('<z>(m)')
title('Motion in magnetic bottle')
grid on
