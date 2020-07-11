function []= PD_control()
clc
clear all;
close all;

l1 = 0.3; l2 = 0.3; l3 = 0.3; m1 = 0.5; m2 = 0.5; m3 = 0.5; g = 9.8;

%% Inverse Kinematics

% we use geometric method to get IK. Equations below refers to Spong's
% textbook, p.88 & p.91.
x_init = [0.3, 0, 0.6];
xi = x_init(1);
yi = x_init(2);
zi = x_init(3);

q1 = atan(yi/xi);
D = (xi^2 + yi^2 + (zi-l1)^2 - l2^2 - l3^2) / (2*l2*l3);  % cos(q3)
q3 = atan2(sqrt(1-D^2), D);  % initial configure is elbow-down, according to DH table q3 should be positive
q2 = atan2(zi-l1, sqrt(xi^2 + yi^2)) - atan2(l3*sin(q3), l2+l3*cos(q3));

x0= [q1, q2, q3, 0, 0, 0]; % Initial Condition 

tf=8;

%% Solve the closed-loop system nonlinear differential equation (PlanarArmODE) via ode45
%%ode45 solves the differential equation and returns X with respect to T.
global torque 
global end_p
global err_p
torque=[];
end_p=[];
err_p=[];

[T,X] = ode45(@(t,x)planarArmODE(t,x),[0 tf],x0);


%% Animation

figure('Name', 'Animation of PD control with GP');

for i = 1:1:0.5*size(T,1)  % The full time interval will be too long, and most of the last parts are meaningless because the final angle chnages are limited
    
clf;
x_1 = 0; y_1 = 0; z_1 = l1;  % T01 last column
x_2 = l2*cos(X(i,1))*cos(X(i,2)); y_2 = l2*cos(X(i,2))*sin(X(i,1)); z_2 = l1+l2*sin(X(i,2));  % T02 last column
del4 = l3*cos(X(i,2)+X(i,3))+l2*cos(X(i,2));
x_3 = cos(X(i,1))*del4; y_3 = sin(X(i,1))*del4; z_3 = l1+l3*sin(X(i,2)+X(i,3))+l2*sin(X(i,2));  % T03 last column

gg(i) = x_3; hh(i) = y_3; ii(i) = z_3;  % store the values in arrays so that they will not be erased

plot3([0 x_1],[0 y_1],[0 z_1],'linewidth',4);

hold on;  % this 'hold on' command is really important. This was either mentioned by some website comments or John Dong in RBE 500. 
          % I don't know why it is essential, but without it the first link will not occur
          
plot3([x_1 x_2],[y_1 y_2],[z_1 z_2],'linewidth',4);
plot3([x_2 x_3],[y_2 y_3],[z_2 z_3],'linewidth',4);
axis([0 0.5 0 0.5  0 1]); 
pause(0.00001);  % difine the interval of the frames in the anime
end

plot3(gg, hh, ii,'color','y','linewidth',4)  % draw end-effector traj
hold on;
axis([0 0.5 0 0.5  0 1]);


%% Plot Data

figure('Name', 'plot of PD control with GP');

clf;
hold off
subplot(3,2,1);
plot(T, X(:,1),'r-');
hold on;
plot(T, X(:,2),'g--');
hold on
plot(T, X(:,3),'b.');
title('Joint angles vs time')

subplot(3,2,2);
plot(T, X(:,4),'r-');
hold on;
plot(T, X(:,5),'g--');
hold on
plot(T, X(:,6),'b.');
title('Velocities vs time')


subplot(3,2,3);
plot(T(1:size(T,1)), torque(1,1:size(T,1)),'r-' );  % the torque input converges much quicker than others
hold on
plot(T(1:size(T,1)), torque(2,1:size(T,1)),'g--');
hold on
plot(T(1:size(T,1)), torque(3,1:size(T,1)),'b.');
title('Control input')

subplot(3,2,4);
plot(T, end_p(1,1:size(T,1)),'r-' );
hold on
plot(T, end_p(2,1:size(T,1)),'g--');
hold on
plot(T, end_p(3,1:size(T,1)),'b.');
title('End-effector position vs time')

subplot(3,2,5);
plot(T, err_p(1,1:size(T,1)),'r-' );
hold on
plot(T, err_p(2,1:size(T,1)),'g--');
hold on
plot(T, err_p(3,1:size(T,1)),'b.');
title('Error vs time')

torque=[];

%% Definging Functions

    function dx = planarArmODE(t,x)
        x_final = [0, 0.3, 0.6];
        xf = x_final(1);
        yf = x_final(2);
        zf = x_final(3);

        qdes1 = atan(yf/xf);
        D_f = (xf^2 + yf^2 + (zf-l1)^2 - l2^2 - l3^2) / (2*l2*l3);  % cos(q3)
        qdes3 = atan2(-sqrt(1-D_f^2), D_f);  % initial configure is elbow-up, according to DH table q3 should be negative
        qdes2 = atan2(zf-l1, sqrt(xf^2 + yf^2)) - atan2(l3*sin(qdes3), l2+l3*cos(qdes3));

        theta_d=[qdes1;qdes2;qdes3]; % Desired Set-Point Position
        dtheta_d=[0;0;0]; % Desired velocity (Derivative of theta_d)
        ddtheta_d=[0;0;0];
        theta= x(1:3,1);
        dtheta= x(4:6,1);
        
        
        global Mmat Cmat
        Mmat = [(l2^2*(m2+m3+(m2+m3)*cos(2*x(2))) + l3^2*(m3+m3*cos(2*x(2)+2*x(3)))) /2 + l2*l3*m3*(cos(x(3))+cos(2*x(2)+x(3))), 0, 0; ...
                0, l2^2*(m2+m3)+l3^2*m3+2*l2*l3*m3*cos(x(3)), l3*m3*(l3+l2*cos(x(3))); ...
                0, l3*m3*(l3+l2*cos(x(3))), l3^2*m3];
        
        del1 = sin(2*x(2)+2*x(3));
        del2 = sin(2*x(2)+x(3));
        Cmat = [-x(4)*(x(5)*(l2^2*(m2+m3)*sin(2*x(2)) + l3^2*m3*del1 + 2*l2*l3*m3*del2) + x(6)*(l3^2*m3*del1 + l2*l3*m3*(del2+sin(x(3))))); ...
                x(4)^2/2*(l2^2*(m2+m3)*sin(2*x(2)) + l3^2*m3*del1 + 2*l2*l3*m3*del2) - x(6)*sin(x(3))*l2*l3*m3*(x(6) + 2*x(5)); ...
                l3*m3/2*(x(4)^2*(l3*del1+l2*sin(x(3))+l2*del2) + 2*x(5)^2*l2*sin(x(3)))];
        
        invM = inv(Mmat);
        invMC = invM*Cmat;
        
        tau = PDControl(theta_d,dtheta_d,ddtheta_d,theta,dtheta);
        del3 = l3*cos(x(2)+x(3))+l2*cos(x(2));
        endp = [cos(x(1))*del3; sin(x(1))*del3; l1+l3*sin(x(2)+x(3))+l2*sin(x(2))];
        err = x_final.' - endp;
            
        
        torque =[torque, tau];
        end_p =[end_p, endp];
        err_p =[err_p, err];
        dx=zeros(6,1);
        dx(1) = x(4); %dtheta1
        dx(2) = x(5); %dtheta2
        dx(3) = x(6); %dtheta3
        dx(4:6) = -invMC +invM*tau; % because ddot theta = -M^{-1}C + M^{-1} tau
    end


 function tau = PDControl(theta_d,dtheta_d,ddtheta_d,theta,dtheta)
        Kp=[20 0 0; 0 15 0; 0 0 30];
        Kv=[20 0 0; 0 15 0; 0 0 30];
        e=theta_d-theta; % position error
        de = dtheta_d - dtheta; % velocity error
        tau = Kp*e + Kv*de;
    end
    
disp('Finish.');

end
