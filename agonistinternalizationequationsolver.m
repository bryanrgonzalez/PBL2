%Agonist-induced internalization solver
%% equationsolver.m
% Owner: Team 6
% Date initiated: 11/01/2021
% Date last modified: 11/01/2021


%% Workspace initiation
clear, format short e, clf
%kplus,kbplaus,kaminus (10^-20)
%% Establishing constants, change based on new equations
kp          = 1;
kq          = 20;
kf          = 35211267.61*1e-10;
kr          = 1;
kplus       = 81*1e9 * 1e-10;
kminus      = 1;
kbplus      = 1.2e10 * 1e-10;
kgplus      = 0.1;
kaplus      = 1;
kbminus     = 0.144;
kgminus     = 1e-4;
kaminus     = 1e18 * 1e-20;
zeta        = 1e3;
mu          = 2;
v           = 1;
zetaplus    = 1;
muplus      = 1;
vplus       = 1;
zetaminus   = 1e-3;
muminus     = muplus/mu;
vminus      = 1;
thetazeta   = 1;
thetamu     = 1;
thetav      = 1;
thetazetamu = 1;
thetazetav  = 1;
thetamuv    = 1;
thetaa      = 1;
kcon        = (5+1.0/3.0)*10^-5;
kint        = (4+5.0/6.0)*10^-4;
A           =3e-10*1e10; %M % 1e-15 mimics paper

Const = [kq, kminus, kp, kplus, muminus, thetamu, kaplus, muplus, kaminus, kbplus, kbminus,kr,kf,thetazeta,zetaminus,zetaplus,thetav,vminus,vplus,thetamuv,thetaa,thetazetav,thetazetamu, kgplus, kgminus,A, kcon, kint];         
% C(1) = kq        C(2) = k_        C(3) = kp         C(4) = k+        
% C(5) = μ_        C(6) = θ{μ}        C(7) = ka+         C(8) = μ+    
% C(9) = ka-        C(10) = kb+        C(11) = kb-         C(12) = kg+
% C(13) = kg-    
% kr=C(12);
% kf=C(13);
% thetazeta=C(14);
% zetaminus=C(15);
% zetaplus=C(16);
% thetav=C(17);
% vminus=C(18);
% vplus=C(19);
% thetamuv=C(20);
% thetaa=C(21);
% thetazetav=C(22);
% thetazetamu=C(23);
% kgplus = C(24);
% kgminus = C(25);
% kcon = C(26);
% kint = C(27);

tspan = logspace(-10, 10, 10000); 

%Updated Basal Levels with drug added at t=9963.22s
%Note: t=9963.22s resets to t=0
R = 0.369777;
G = 0.382346;
Rz = 0.0224164;
RzG = 0.0988111;
RG = 1.06666;
a = 0.116168;
at = 0.986013;
by = 1.10218;
AR=0;
ARG=0;
ARzG=0;
ARz=0;



yinit = [R, Rz, AR, ARz, G, RG, RzG, ARG, ARzG, a, at, by];

% R = x(1);
% Ra = x(2);
% AR = x(3);
% ARa = x(4);
% G = x(5);
% RG = x(6);
% RaG = x(7);
% ARG = x(8);
% ARaG = x(9);
% alpha = x(10);
% alphas = x(11);
% betagamma = x(12);
%% Solving ODE system
DiffFileName = 'agonistinternalizationequations';
DE = eval(sprintf('@(t,x,C) %s(t,x,C)', DiffFileName));
% options = odeset('RelTol',1e-6,'AbsTol',1e4);
[tout, yout] = ode15s(@(t,x) DE(t,x,Const), tspan, yinit);

% %Inactivated
% figure;
% semilogx(tout,yout(:,1), 'k-')
% hold on;
% semilogx(tout,yout(:,5), 'k--')
% semilogx(tout,yout(:,6), 'k:')
% xlim([1e-3 1e4])
% title('R, G, and RG vs. log(t)')
% legend('R', 'G', 'RG', 'Location', 'Best')
% grid('on')
% xlabel('time (s)');
% ylabel('Concentration (1e-10 M)');
% hold off;
% 
% %Activated
% figure;
% semilogx(tout,yout(:,2), 'k-')
% hold on;
% semilogx(tout,yout(:,7), 'k--')
% xlim([1e-6 1e4])
% title('R*, R*G, vs. log(t)')
% legend('R*', 'R*G', 'Location', 'Best')
% grid('on')
% xlabel('time (s)');
% ylabel('Concentration (1e-10 M)');
% hold off;


Rtot = yout(:,1)+yout(:,2)+yout(:,3)+yout(:,4)+yout(:,6)+yout(:,7)+yout(:,8)+yout(:,9);
Gtot = yout(:,5)+yout(:,6)+yout(:,7)+yout(:,8)+yout(:,9)+yout(:,11)+yout(:,10);
ActivatedR = yout(:,2)+yout(:,4);+yout(:,7)+yout(:,9);

yinit = [R, Rz, AR, ARz, G, RG, RzG, ARG, ARzG, a, at, by];

%Totals
figure
semilogx(tout, Rtot, 'k-')
hold on;
semilogx(tout, ActivatedR, 'k--');
semilogx(tout(1:20:end), Gtot(1:20:end), 'k.', 'MarkerSize', 3.5)
gfree = yout(:,5);
semilogx(tout(1:20:end), gfree(1:20:end), 'ks', 'MarkerSize', 2)
legend('Total CB1 Receptor','Total Activated CB1 Receptor', 'Total G Protein', 'Total Free G Protein', 'Location', 'northwest')
title('Total Concentrations vs. Time')
grid('off')
xlim([1e-3 1e4])
ylim([0 3.5])
xlabel('Time (s)')
ylabel('Concentration (1e-10 M)')
hold off;

%Alpha and Friends
figure
semilogx(tout,yout(:,10), 'k-')
hold on;
semilogx(tout,yout(:,11), 'k--')
byp = yout(:,12)
semilogx(tout(1:20:end),byp(1:20:end), 'k.', 'MarkerSize', 3.5)
xlabel('Time (s)')
ylabel('Concentration (1e-10 M)')
xlim([1e-3 1e4])
ylim([0 2.2])
legend('α', 'α†','βγ', 'Location', 'northwest')
title('Concentration of α , α† , βγ vs. Time')
hold off;

% %Ligand
% figure;
% semilogx(tout,yout(:,3), 'k-')
% hold on;
% semilogx(tout,yout(:,8), 'k--')
% title('AR, ARG, vs. log(t)')
% xlim([1e-3 1e4])
% legend('AR', 'ARG', 'Location', 'Best')
% grid('on')
% xlabel('time (s)');
% ylabel('Concentration (1e-10 M)');
% 
% figure;
% semilogx(tout,yout(:,4), 'k:')
% hold on;
% semilogx(tout, yout(:,9), 'k-.');
% title('AR*, AR*G vs. log(t)');
% xlim([1e-3 1e4])
% legend('AR*', 'AR*G', 'Location', 'Best');
% grid('on')
% xlabel('time (s)');
% ylabel('Concentration (1e-10 M)');



% figure
% lRG = log(Rtot)-log(Gtot)
% plot(lRG, yout(:,7), 'k-')
% xlim([-4 10])
% hold on;
% title('log(Rtot/Gtot) vs. α†')
% hold off;
