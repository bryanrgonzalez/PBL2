%% noligandequationsolver.m
% Owner: Team 6
% Date initiated: 11/01/2021
% Date last modified: 11/01/2021


%% Workspace initiation
clear, format short e, clf

%% Establishing constants, change based on new equations
kp      	= 1;
kq      	= 20;
kplus   	= 81*1e9 * 1e-10;
kminus  	= 1;
kbplus  	= 1.2e10 * 1e-10;
kgplus  	= 0.1;
kaplus  	= 1;
kbminus 	= 0.144;
kgminus 	= 1e-4;
kaminus 	= 1e18 * 1e-20;
mu      	= 2;
muplus  	= 1;
muminus 	= muplus/mu;
thetamu 	= 1;
kcon = 5.333333e-5;

Const = [kq, kminus, kp, kplus, muminus, thetamu, kaplus, muplus, kaminus, kbplus, kbminus, kgplus, kgminus,kcon];    	 
% C(1) = kq    	C(2) = k_    	C(3) = kp     	C(4) = k+   	 
% C(5) = μ_    	C(6) = θ{μ}    	C(7) = ka+     	C(8) = μ+    
% C(9) = ka-    	C(10) = kb+    	C(11) = kb-     	C(12) = kg+
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

% x(1) = [R]    	x(2) = [R*]    	x(3) = [G]    	x(4) = [RG]
% x(5) = [R*G]    	x(6) = [α]    	x(7) = [α†]    	x(8) = [βγ]

%% Solving ODE system
DiffFileName = 'internalizationequations';
DE = eval(sprintf('@(t,x,C) %s(t,x,C)', DiffFileName));
options = odeset('RelTol',1e-6,'AbsTol',1e-4);
[tout, yout] = ode15s(@(t,x) DE(t,x,Const), tspan, yinit, options);

Rtot = yout(:,1)+yout(:,2)+yout(:,4)+yout(:,5);
Gtot = yout(:,3)+yout(:,4)+yout(:,5)+yout(:,6)+yout(:,7);
ActivatedR = yout(:,2)+yout(:,5);
% %
semilogx(tout,yout(:,1), 'k-')
hold on;
semilogx(tout,yout(:,3), 'k--')
semilogx(tout,yout(:,4), 'k:')
xlim([1e-3 1e4])
title('R, G, RG vs. log(t)')
legend('R', 'G', 'RG', 'Location', 'Best')
grid('on')
xlabel('time (s)');
ylabel('Concentration (1e-10 M)');
hold off;

figure
semilogx(tout,yout(:,2), 'k-')
hold on;
semilogx(tout,yout(:,5), 'k--')
xlim([1e-3 1e4])
legend('R*', 'R*G', 'Location', 'Best')
title('R* and R*G vs. log(t)')
xlabel('time (s)');
grid('on')
ylabel('Concentration (1e-10 M)');
hold off;
%
% Rtot = yout(:,1)+yout(:,2)+yout(:,4)+yout(:,5);
% Gtot = yout(:,3)+yout(:,4)+yout(:,5)+yout(:,6)+yout(:,7);
%
figure
semilogx(tout, Rtot, 'k-')
hold on;
xlim([1e-3 1e4])
semilogx(tout, Gtot, 'k:')
title('Rtotal, Gtotal vs. log(t)')
legend('Rtotal', 'Gtotal', 'Location', 'Best')
grid('on')
xlabel('time (s)');
ylabel('Concentration (1e-10 M)');
ylim([0 5])
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

%yinit = [R, Rz, G, RG, RzG, a, at, by];

figure
semilogx(tout, Rtot, 'k-')
hold on;
semilogx(tout,ActivatedR, 'k--');
semilogx(tout(1:20:end), Gtot(1:20:end), 'k.', 'MarkerSize', 3.5)
gfree = yout(:,3);
semilogx(tout(1:20:end), gfree(1:20:end), 'ks', 'MarkerSize', 2)
legend('Total CB1 Receptor','Total Activated CB1 Receptor', 'Total G Protein', 'Total Free G Protein', 'Location', 'northwest')
title('Total Concentrations vs. Time')
grid('off')
xlim([1e-3 1e4])
ylim([0 3.5])
xlabel('Time (s)')
ylabel('Concentration (1e-10 M)')
hold off;

% figure
% lRG = log(Rtot)-log(Gtot)
% plot(lRG, yout(:,7), 'k-')
% xlim([-4 10])
% hold on;
% title('log(Rtot/Gtot) vs. α†')
% hold off;

% close all

%CONCENTRATIONS at t = 10000s (9963.22s) (in 1e-10 M)
%alpha = 0.100741
%alpha-cross = 0.799331
%beta-gamma = 0.900071
%Rtot = 1.55769
%Gtot = 1.14964
%[RG] = 0.868906
%[R*G] = 0.0801269
%[R*] = 0.0322045
%[G] = 0.200605
%[R] = 0.576453
