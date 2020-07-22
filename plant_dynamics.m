    function dx = plant_dynamics(x, u)
    
        global l1 l2 l3 m1 m2 m3

        l1 = 0.3; l2 = 0.3; l3 = 0.3; m1 = 0.5; m2 = 0.5; m3 = 0.5; g = 9.8;
        
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
        invMC = invM * Cmat;
                
        
        dx=zeros(6,1);
        dx(1) = x(4); %dtheta1
        dx(2) = x(5); %dtheta2
        dx(3) = x(6); %dtheta3
        dx(4:6) = -invMC + invM * u; % because ddot theta = -M^{-1}C + M^{-1} tau
    end