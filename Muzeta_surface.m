%% equationsolver.m
% Owner: Team 6
% Date initiated: 11/01/2021
% Date last modified: 11/01/2021


%% Workspace initiation
clear, format short e, clf
%kplus,kbplaus,kaminus (10^-20)
%% Establishing constants, change based on new equations
kp      	= 1;
kq      	= 1000;
kf      	= 8.4e7;
kr      	= 0.37;
kplus   	= 3.6e7 * 1e-10;
kminus  	= 3e-3;
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
A           =1e-12*1e10; %M



Const = [kq, kminus, kp, kplus, muminus, thetamu, kaplus, muplus, kaminus, kbplus, kbminus,kr,kf,thetazeta,zetaminus,zetaplus,thetav,vminus,vplus,thetamuv,thetaa,thetazetav,thetazetamu, kgplus, kgminus,A];    	 
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

tspan = logspace(-10, 10, 10000); 

%THESE NEED TO BE BASAL LEVELS ADJUST WITH UPDATES
R = 1.97513e-10*1e10;
G = 1.52301e-10*1e10;
Rz = 0.00197945e-10*1e10;
RzG = 0.00432847e-10*1e10;
RG = 2.16856e-10*1e10;
a = 0.410409e-10*1e10;
at = 0.0436912e-10*1e10;
by = 0.454099e-10*1e10;
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
DiffFileName = 'ligandequations';
DE = eval(sprintf('@(t,x,C) %s(t,x,C)', DiffFileName));
% options = odeset('RelTol',1e-6,'AbsTol',1e4);

[i,j] = meshgrid(1.5:0.1:2.5);
j = j*1000-1000;
z = zeros(size(i));
c = 1;
d = 1;
for a = 1.5:0.1:2.5
c = 1;
    for b = 0.5:0.1:1.5
[tout, yout] = ode15s(@(t,x) DE(t,x,[kq, kminus, kp, kplus, 1, thetamu, kaplus, a, kaminus, kbplus, kbminus,kr,kf,thetazeta,1e-3,b,thetav,vminus,vplus,thetamuv,thetaa,thetazetav,thetazetamu, kgplus, kgminus,A]), tspan, yinit);
[y, x] = max(yout(:,11));
z(c,d) = y;
c = c+1;
    end
d = d+1;
end
surf(i,j,z)
xlabel("μ")
ylabel("ζ")
zlabel("Concentration of α† (M)")
title("Concentration of α† at its peak with varying μ and ζ")

