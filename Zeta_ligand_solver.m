%% equationsolver.m
% Owner: Team 6
% Date initiated: 11/01/2021
% Date last modified: 11/01/2021

%% Workspace initiation
clear, format short e, clf
%kplus,kbplaus,kaminus (10^-20)
%% Establishing constants, change based on new equations
kp      	= 1;
kq      	= 20;
kf      	= 35211267.61*1e-10;
kr      	= 1;
kplus   	= 81*1e9 * 1e-10;
kminus  	= 1;
kbplus  	= 1.2e10 * 1e-10;
kgplus  	= 0.1;
kaplus  	= 1;
kbminus 	= 0.144;
kgminus 	= 1e-4;
kaminus 	= 1e18 * 1e-20;
zeta    	= 1e3;
mu      	= 2;
v       	= 1;
zetaplus	= 1;
muplus  	= 1;
vplus   	= 1;
zetaminus   = 1e-3;
muminus 	= muplus/mu;
vminus  	= 1;
thetazeta   = 1;
thetamu 	= 1;
thetav  	= 1;
thetazetamu = 1;
thetazetav  = 1;
thetamuv	= 1;
thetaa  	= 1;
kcon    	= (5+1.0/3.0)*10^-5;
kint    	= (4+5.0/6.0)*10^-4;
A       	=3e-10*1e10; %M % 1e-15 mimics paper

Const = [kq, kminus, kp, kplus, muminus, thetamu, kaplus, muplus, kaminus, kbplus, kbminus,kr,kf,thetazeta,zetaminus,zetaplus,thetav,vminus,vplus,thetamuv,thetaa,thetazetav,thetazetamu, kgplus, kgminus,A, kcon, kint];    	 
% C(1) = kq    	C(2) = k_    	C(3) = kp     	C(4) = k+   	 
% C(5) = μ_    	C(6) = θ{μ}    	C(7) = ka+     	C(8) = μ+    
% C(9) = ka-    	C(10) = kb+    	C(11) = kb-     	C(12) = kg+
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
figure(1);
k = 1;
for i = 0.5:0.1:1.5;
[tout, yout] = ode15s(@(t,x) DE(t,x,[kq, kminus, kp, kplus, muminus, thetamu, kaplus, muplus, kaminus, kbplus, kbminus,kr,kf,thetazeta,1e-3,i,thetav,vminus,vplus,thetamuv,thetaa,thetazetav,thetazetamu, kgplus, kgminus,A,kcon,kint]), tspan, yinit);
semilogx(tout,yout(:,7)+yout(:,9),'color',rand(1,3));
hold on;
xlabel('Time (s)')
ylabel('Concentration (1e-10 M)')
xlim([1e-3 1e4])
Legend{k}=strcat('ζ= ', num2str(i*1000));
title('Parameter Sweep of ζ -', 'Total Concentration of Activated Receptor vs. Time with Internalization');
k = k+1;
end;
legend(Legend, 'Location', 'northwest');

figure(2);
r = 1;
for i = 0.5:0.1:1.5;
[tout, yout] = ode15s(@(t,x) DE(t,x,[kq, kminus, kp, kplus, muminus, thetamu, kaplus, muplus, kaminus, kbplus, kbminus,kr,kf,thetazeta,1e-3,i,thetav,vminus,vplus,thetamuv,thetaa,thetazetav,thetazetamu, kgplus, kgminus,A,kcon,kint]), tspan, yinit);
[x, y] = max(yout(:,11));
semilogx(tout,yout(:,11),'color',rand(1,3));
hold on;
xlabel('Time (s)')
ylabel('Concentration (1e-10 M)')
xlim([1e-3 1e4])
%zetarray(r) = i;
%maxZ(r) = x;
Legend2{r}=strcat('ζ= ', num2str(i*1000));
title('Parameter Sweep of ζ', 'Concentration of α† vs. Time with Internalization');
r = r+1;
end;
legend(Legend2, 'Location', 'northwest');
% figure(3)
% plot(zetarray,maxZ,'-x')

