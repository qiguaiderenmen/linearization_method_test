function [M, C, G] = MCGCalc()
    clc
    clear all
    
    pi = sym('pi', 'real');
    syms q1 q2 q3
    syms l1 l2 l3
    syms dq1 dq2 dq3 ddq1 ddq2 ddq3 real;
    syms m1 m2 m3 g real;
    
    a = [0 l2 l3];
    alpha = [pi/2 0 -pi/2];
    d = [l1 0 0];
    theta = [q1 q2 q3];
    
    T01 = dh2transMatrix(theta(1), d(1), alpha(1), a(1))
    T12 = dh2transMatrix(theta(2), d(2), alpha(2), a(2));
    T23 = dh2transMatrix(theta(3), d(3), alpha(3), a(3));
    
    T02 = simplify(T01 * T12)
    T03 = simplify(T02 * T23)
    
    double(pi);
    
    J_l1 = jacobian(T01(1:3,4),[q1]);
    J_l2 = jacobian(T02(1:3,4),[q1 q2]);
    J_l3 = jacobian(T03(1:3,4),[q1 q2 q3]);
    J_l1 = simplify(J_l1);
    J_l2 = simplify(J_l2);
    J_l3 = simplify(J_l3);
    
    v_m1 = J_l1 * dq1;
    v_m2 = J_l2 * [dq1 ; dq2];
    v_m3 = J_l3 * [dq1 ; dq2 ; dq3];
    v_m1 = simplify(v_m1);
    v_m2 = simplify(v_m2);
    v_m3 = simplify(v_m3);

    
    K1 = 0.5 * m1 * (v_m1.' * v_m1);
    K2 = 0.5 * m2 * (v_m2.' * v_m2);
    K3 = 0.5 * m3 * (v_m3.' * v_m3);
    K1 = simplify(K1);
    K2 = simplify(K2);
    K3 = simplify(K3);
    
    P1 = m1 * g * T01(3,4);
    P2 = m2 * g * T02(3,4);
    P3 = m3 * g * T03(3,4);
    P1 = simplify(P1);
    P2 = simplify(P2);
    P3 = simplify(P3);

    K = K1 + K2 + K3;
    P = P1 + P2 + P3;

    
    L = simplify( K - P );

    
    syms th1(t) th2(t) th3(t);
    
    A1 = diff(L,dq1);
    
    A1t = subs(A1,[q1 q2 q3 dq1 dq2 dq3], [th1 th2 th3 diff(th1(t),t) diff(th2(t),t) diff(th3(t), t)]);
    
    dA1t = diff(A1t, t);
    
    A1 = subs(dA1t,[th1 th2 th3 diff(th1(t),t) diff(th2(t),t) diff(th3(t),t) ...
        diff(th1(t),t,t) diff(th2(t),t,t) diff(th3(t), t,t)],[q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3]);
    
    B1 = diff(L,q1);
    
    A2 = diff(L,dq2);
    
    A2t = subs(A2,[q1 q2 q3 dq1 dq2 dq3], [th1 th2 th3 diff(th1(t),t) diff(th2(t),t) diff(th3(t), t)]);
    dA2t = diff(A2t, t);
    A2 = subs(dA2t,[th1 th2 th3 diff(th1(t),t) diff(th2(t),t) diff(th3(t),t) ...
        diff(th1(t),t,t) diff(th2(t),t,t) diff(th3(t), t,t)],[q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3]);
    B2 = diff(L,q2);
    A3 = diff(L,dq3);
    A3t = subs(A3,[q1 q2 q3 dq1 dq2 dq3], [th1 th2 th3 diff(th1(t),t) diff(th2(t),t) diff(th3(t), t)]);
    dA3t = diff(A3t, t);
    A3 = subs(dA3t,[th1 th2 th3 diff(th1(t),t) diff(th2(t),t) diff(th3(t),t) ...
        diff(th1(t),t,t) diff(th2(t),t,t) diff(th3(t), t,t)],[q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3]);
    B3 = diff(L,q3);
    
    Tau_l1 = A1 - B1;
    Tau_l2 = A2 - B2;
    Tau_l3 = A3 - B3;
    
    M11 = simplify(Tau_l1 - subs(Tau_l1,ddq1,0)) /ddq1;
    M12 = simplify(Tau_l1 - subs(Tau_l1,ddq2,0)) /ddq2;
    M13 = simplify(Tau_l1 - subs(Tau_l1,ddq3,0)) /ddq3;
    M21 = simplify(Tau_l2 - subs(Tau_l2,ddq1,0)) /ddq1;
    M22 = simplify(Tau_l2 - subs(Tau_l2,ddq2,0)) /ddq2;
    M23 = simplify(Tau_l2 - subs(Tau_l2,ddq3,0)) /ddq3;
    M31 = simplify(Tau_l3 - subs(Tau_l3,ddq1,0)) /ddq1;
    M32 = simplify(Tau_l3 - subs(Tau_l3,ddq2,0)) /ddq2;
    M33 = simplify(Tau_l3 - subs(Tau_l3,ddq3,0)) /ddq3;
    G11 = subs(Tau_l1, [dq1 dq2 dq3 ddq1 ddq2 ddq3], [0 0 0 0 0 0]);
    G21 = subs(Tau_l2, [dq1 dq2 dq3 ddq1 ddq2 ddq3], [0 0 0 0 0 0]);
    G31 = subs(Tau_l3, [dq1 dq2 dq3 ddq1 ddq2 ddq3], [0 0 0 0 0 0]);
    
    M = simplify([M11 M12 M13;
        M21 M22 M23;
        M31 M32 M33]);
    
    G = simplify([G11; G21; G31]);
    
    C11 = Tau_l1 - (M(1,:) * [ddq1 ddq2 ddq3].' + G11); 
    C21 = Tau_l2 - (M(2,:) * [ddq1 ddq2 ddq3].' + G21); 
    C31 = Tau_l3 - (M(3,:) * [ddq1 ddq2 ddq3].' + G31);
    
    C = simplify([C11; C21; C31]);
end