% ligandequations
function dxdt = ligandequations(t,x,C)
R = x(1);
Ra = x(2);
AR = x(3);
ARa = x(4);
G = x(5);
RG = x(6);
RaG = x(7);
ARG = x(8);
ARaG = x(9);
alpha = x(10);
alphas = x(11);
betagamma = x(12); 

kq=C(1);
kminus=C(2);
kp=C(3);
kplus=C(4);
muminus=C(5);
thetamu=C(6);
kaplus=C(7);
muplus=C(8);
kaminus=C(9);
kbplus=C(10);
kbminus=C(11);
kr=C(12);
kf=C(13);
thetazeta=C(14);
zetaminus=C(15);
zetaplus=C(16);
thetav=C(17);
vminus=C(18);
vplus=C(19);
thetamuv=C(20);
thetaa=C(21);
thetazetav=C(22);
thetazetamu=C(23);
kgplus = C(24);
kgminus = C(25);
A=C(26);

    dxdt = [
%dR/dt
kq.*Ra + kminus.*RG + kr.*AR - (kp + kplus.*G + kf.*A).*R;

%dRa/dt
kp.*R + (muminus.*thetamu.*kminus + kaplus).*RaG + thetazeta.*zetaminus.*kr.*ARa - (kq + thetamu.*muplus.*kplus.*G + thetazeta.*zetaplus.*kf.*A + kaminus.*alphas.*betagamma).*Ra;

%dAR/dt
kf.*A.*R + zetaminus.*kq.*ARa + thetav.*vminus.*kminus.*ARG - (kr + thetav.*vplus.*kplus.*G + zetaplus.*kp).*AR

%dARa/dt
zetaplus.*kp.*AR + thetazeta.*zetaplus.*kf.*A.*Ra + vminus.*(muminus.*thetamuv.*kminus + thetaa.*kaplus).*ARaG - (zetaminus.*(kq + thetazeta.*kr) + thetamuv.*muplus.*vplus.*kplus.*G + thetaa.*vplus.*kaminus.*alphas.*betagamma).*ARa;

%dG/dt
kminus.*(RG + thetamu.*muminus.*RaG + thetav.*vminus.*ARG + thetamuv.*muminus.*vminus.*ARaG) + kbplus.*alpha.*betagamma - (kplus.*(R + thetamu.*muplus.*Ra + vplus.*(thetav.*AR + thetamuv.*muplus.*ARa)) + kbminus).*G;

%dRG/dt
kplus.*R.*G + muminus.*kq.*RaG + vminus.*kr.*ARG - (kminus + vplus.*kf.*A + muplus.*kp).*RG;

%dRaG/dt
muplus.*kp.*RG + thetazetav.*zetaminus.*vminus.*kr.*ARaG + thetamu.*muplus.*kplus.*Ra.*G + kaminus.*alphas.*betagamma.*Ra - (muminus.*(kq + thetamu.*kminus) + thetazetav.*zetaplus.*vplus.*kf.*A + kaplus).*RaG;

%dARG/dt
thetav.*vplus.*kplus.*AR.*G + vplus.*kf.*A.*RG + thetazetamu.*zetaminus.*muminus.*kq.*ARaG - (vminus.*(kr + thetav.*kminus) + thetazetamu.*zetaplus.*muplus.*kp).*ARG;

%dARaG/dt
thetamuv.*muplus.*vplus.*kplus.*ARa.*G + thetazetamu.*zetaplus.*vplus.*kf.*A.*RaG + thetazetamu.*zetaplus.*muplus.*kp.*ARG + thetaa.*vplus.*kaminus.*alphas.*betagamma.*ARa - (thetamuv.*muminus.*vminus.*kminus + thetazetav.*zetaminus.*vminus.*kr + thetazetamu.*zetaminus.*muminus.*kq + thetaa.*vminus.*kaplus).*ARaG;

%dalpha/dt 
kgplus.*alphas + kbminus.*G - kbplus.*alpha.*betagamma - kgminus.*alpha;

%dalphas/dt
kaplus.*(thetaa.*vminus.*ARaG + RaG) + kgminus.*alpha - (kgplus + kaminus.*(Ra + thetaa.*vplus.*ARa).*betagamma).*alphas;

%dbetagamma/dt
kaplus.*(RaG + thetaa.*vminus.*ARaG) + kbminus.*G - kbplus.*alpha.*betagamma - kaminus.*(Ra + thetaa.*vplus.*ARa).*alphas.*betagamma;
];    
  
    
end


% x(1) = [R]        x(2) = [R*]        x(3) = [RG]        x(4) = [G]
% x(5) = [R*G]        x(6) = [α]        x(7) = [βγ]        x(8) = [α†] 
% x(9) = [A]        x(10) = [AR]        x(11) = [ARG]        x(12) = [AR*G]
% x(13) = [AR*]
% C(1) = kq        C(2) = k_        C(3) = kp         C(4) = k+         
% C(5) = μ_        C(6) = θ{μ}        C(7) = ka+          C(8) = μ+     
% C(9) = ka-         C(10) = kb+         C(11) = kb-     
% C(12) = kr        C(13) = kf        C(14) = θ{𝜁}        C(15) = 𝜁_
% C(16) = 𝜁 +         C(17) = θ{v}        C(18) = v_        C(19) = v+  
% C(20) = θ{μv}    C(21) =  θ{a}        C(22) = θ{𝜁v}        C(23) = θ{𝜁μ}
% C(24) = kg+         C(25) = kg- 
