%% noligandequationsolver.m
% Owner: Team 6
% Date initiated: 11/01/2021
% Date last modified: 11/01/2021


%% Workspace initiation
clear, format short e, clf

%% Establishing constants, change based on new equations
kp          = 1;
kq          = 20;
kplus       = 81*1e9 * 1e-10;
kminus      = 1;
kbplus      = 1.2e10 * 1e-10;
kgplus      = 0.1;
kaplus      = 1;
kbminus     = 0.144;
kgminus     = 1e-4;
kaminus     = 1e18 * 1e-20;
mu          = 2;
muplus      = 1;
muminus     = muplus/mu;
thetamu     = 1;

Const = [kq, kminus, kp, kplus, muminus, thetamu, kaplus, muplus, kaminus, kbplus, kbminus, kgplus, kgminus];         
% C(1) = kq        C(2) = k_        C(3) = kp         C(4) = k+        
% C(5) = μ_        C(6) = θ{μ}        C(7) = ka+         C(8) = μ+    
% C(9) = ka-        C(10) = kb+        C(11) = kb-         C(12) = kg+
% C(13) = kg-    

R   = 265e-12 * 1e10;
Rz  = 0;
RG  = 0;
G   = 265e-12 * 1e10;
RzG = 0;
a   = 0;
by  = 0;
at  = 0;
tspan = logspace(-20, 20, 10000);        
yinit = [R, Rz, G, RG, RzG, a, at, by];

% x(1) = [R]        x(2) = [R*]        x(3) = [G]        x(4) = [RG]
% x(5) = [R*G]        x(6) = [α]        x(7) = [α†]        x(8) = [βγ]

%% Solving ODE system
DiffFileName = 'noligandequations';
DE = eval(sprintf('@(t,x,C) %s(t,x,C)', DiffFileName));
options = odeset('RelTol',1e-6,'AbsTol',1e-4);
[tout, yout] = ode15s(@(t,x) DE(t,x,Const), tspan, yinit, options);
%
semilogx(tout,yout(:,1), 'k-')
hold on;
semilogx(tout,yout(:,3), 'k--')
xlim([1e-3 1e4])
title('R and G vs. log(t)')
legend('R', 'G', 'Location', 'Best')
hold off;

figure
semilogx(tout,yout(:,2), 'k-')
hold on;
semilogx(tout,yout(:,5), 'k--')
xlim([1e-3 1e4])
legend('R*', 'R*G', 'Location', 'Best')
title('R* and R*G vs. log(t)')
hold off;

Rtot = yout(:,1)+yout(:,2)+yout(:,4)+yout(:,5);
Rztot = yout(:,2)+yout(:,5);
Gtot = yout(:,3)+yout(:,4)+yout(:,5)+yout(:,6)+yout(:,7);

figure
semilogx(tout, Rtot, 'k-')
hold on;
xlim([1e-3 1e4])
ylim([0 3.5])
xlabel('Time (s)')
ylabel('Concentration (1e-10 M)')
semilogx(tout, Rztot, 'k--')
semilogx(tout(1:20:end), Gtot(1:20:end), 'k.', 'MarkerSize', 3.5)
title('Total Concentrations vs. Time')
legend('Total CB1 Receptor', 'Total Activated Receptor', 'Total G Protein', 'Location', 'northwest')
hold off;

figure
semilogx(tout,yout(:,6), 'k-')
hold on;
semilogx(tout,yout(:,7), 'k--')
byp = yout(:,8)
semilogx(tout(1:20:end),byp(1:20:end), 'k.', 'MarkerSize', 3.5)
xlabel('Time (s)')
ylabel('Concentration (1e-10 M)')
xlim([1e-3 1e4])
ylim([0 2.2])
legend('α', 'α†','βγ', 'Location', 'northwest')
title('Concentration of α , α† , βγ vs. Time')
hold off;

% figure
% lRG = log(Rtot)-log(Gtot)
% plot(lRG, yout(:,7), 'k-')
% xlim([-4 10])
% hold on;
% title('log(Rtot/Gtot) vs. α†')
% hold off;

% close all

%EQUILIBRIUM CONCENTRATIONS (in 1e-10 M)
%alpha = 0.0899369
%alpha-cross = 1.1075
%beta-gamma = 1.19743
%Rtot = 2.65
%Gtot = 1.45257
%[RG] = 1.21251
%[R*G] = 0.111639
%[R*] = 0.0677131
%[G] = 0.128413
%[R] = 1.25813
