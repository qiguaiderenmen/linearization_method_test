%% Defining Symbols and Variables
clc
clear all
syms x1 x2 x3 x4 u k1 k2 k3 k4 l lambda real
m= 0.1;
l=1;
g=10;
M = 1;
a= 1/(M+m);

%% Non-Linear State Space Model
% From Zak's Textbook refer section 1.7.4 for dynaical model
x1_dot = x2;
x2_dot = (-m*a*g*(sin(2*x3)/2)+a*l*sin(x3)*(x4^2)*m)/(1-m*a*(cos(x3)^2)) + ...
(a*u/(1-m*a*(cos(x3)^2)));
x3_dot = x4;
x4_dot = (g*sin(x3)-m*l*a*(x4^2)*(sin(2*x3)/2))/(l-m*l*a*(cos(x3)^2))-...
(-a*cos(x3)*u)/(l-m*l*a*(cos(x3)^2));
% Non linear State space model
X_dot = [x1_dot;x2_dot;x3_dot;x4_dot]

%% Linearization of Model at x = 0 and u = 0
AX = [x2; (-m*a*g*(sin(2*x3)/2)+a*l*sin(x3)*(x4^2)*m)/(1-m*a*(cos(x3)^2)); x4; (g*sin(x3)-m*l*a*(x4^2)*(sin(2*x3)/2))/(l-m*l*a*(cos(x3)^2))];
BU = [x2; (a*u/(1-m*a*(cos(x3)^2))); x4; (-a*cos(x3)*u)/(l-m*l*a*(cos(x3)^2))];
A = jacobian(AX,[x1;x2;x3;x4]);
B = jacobian(BU,[u]);
A = simplify(subs(A,[x1,x2,x3,x4],[0,0,0,0]))
B = simplify(subs(B,[x1,x2,x3,x4,u],[0,0,0,0,0]))

%% State Feedback Control Law
k = [k1 k2 k3 k4];
M = A-B*k;
M_lI = M - lambda*eye(4);
detMlI = det(M_lI);
eq1 = subs(detMlI, lambda, -1)==0;
eq2 = subs(detMlI, lambda, -2)==0;
eq3 = subs(detMlI, lambda, -1+1i)==0;
eq4 = subs(detMlI, lambda, -1-1i)==0;
eqn = [eq1, eq2, eq3, eq4];
[k1,k2,k3,k4] = solve(eqn, k);
K = [double(vpa(k1)) double(vpa(k2)) double(vpa(k3)) double(vpa(k4))]

% substituting u = -kx will solve for state feedback control law


%% Animation
tspan = linspace(0,10,250);
x0 = [0; 0; pi/4; 0];
[t,x] = ode45(@(t,x) invpen(t,x,double(A),double(B),K), tspan,x0);
figure
hold on
plot(t,x(:,1),'LineWidth',2);
plot(t,x(:,2),'LineWidth',2);
plot(t,x(:,3),'LineWidth',2);
plot(t,x(:,4),'LineWidth',2);
legend('x','x_dot','theta','theta_dot');
p.g = 10.0; % (m/s^2) gravity
p.l = 1; % (m) pendulum (pole) length
z = x';
pos = cartPolePosition(z,p);
x1 = pos(1,:);
y1 = pos(2,:);
x2 = pos(3,:);
y2 = pos(4,:);
padding = 0.2*p.l;
xLow = min(min(x1,x2)) - padding;
xUpp = max(max(x1,x2)) + padding;
yLow = min(min(y1,y2)) - padding;
yUpp = max(max(y1,y2)) + padding;
extents = [xLow,xUpp,yLow,yUpp];
% Create and clear a figure:
figure(2);
clf;
time = 0;
tic;
while time < t(end)
time = toc;
posDraw = interp1(t',pos',time')';
clf;
hold on;
drawCartPole(time,posDraw,extents);
drawnow;
end

drawCartPole(t(end),pos(:,end),extents);
drawCartPole(t(end),pos(:,end),extents);
drawCartPole(t(end),pos(:,end),extents);

function pos = cartPolePosition(z,p)
x = z(1,:); % Cart position (Not used in dynamics)
q = z(3,:);
l = p.l;
x1 = x;
y1 = zeros(size(x));
x2 = x1 + l*sin(q);
y2 = y1 + l*cos(q);
pos = [x1;y1;x2;y2];
end
function dx = invpen(t,x,A,B,K)
dx = (A-B*K)*x;
end
function drawCartPole(time,pos,extents)
x1 = pos(1);
y1 = pos(2);
x2 = pos(3);
y2 = pos(4);
title(sprintf('t = %2.2f%',time));
plot(extents(1:2),[0,0],'k-','LineWidth',2);
plot(x1, y1, 'ks','MarkerSize',30,'LineWidth',2);
plot([x1,x2], [y1, y2], 'k-','LineWidth',2);
plot(x2, y2, 'ko','MarkerSize',10,'LineWidth',2);
axis equal;
axis(extents);
axis off;
end
