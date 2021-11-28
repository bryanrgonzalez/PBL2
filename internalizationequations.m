function dxdt = internalizationequations(t,x,C)
R = x(1);
Ra = x(2);
G = x(3);
RG = x(4);
RaG = x(5);
alpha = x(6);
alphas = x(7);
betagamma = x(8);  

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
kgplus=C(12);
kgminus=C(13);
kcon = C(14);

    dxdt = [%dR/dt
            kq.*Ra+ kminus.*RG-(kp+kplus.*G).*R-kcon*R;
               
            %dR*/dt
            kp.*R+(muminus.*thetamu.*kminus+kaplus).*(RaG)-(kq+thetamu.*muplus.*kplus.*G+kaminus.*alphas.*betagamma).*Ra-kcon*Ra;
            
            %dG/dt
            kminus.*RG+thetamu.*muminus.*kminus.*RaG+kbplus.*alpha.*betagamma - (kplus.*(R+thetamu.*muplus.*Ra)+kbminus).*G + kcon*(RG+RaG);
            
            %dRG/dt
            kplus.*R.*G+muminus.*kq.*RaG-(kminus+muplus.*kp).*RG-kcon*RG;
            
            %dR*G/dt
            muplus.*kp.*RG+thetamu.*muplus.*kplus.*Ra.*G+kaminus.*alphas.*betagamma.*Ra-(muminus.*(kq+thetamu.*kminus)+kaplus).*RaG-kcon*RaG;
            
            %d(alpha)/dt
            kgplus.*alphas + kbminus.*G - kbplus.*alpha.*betagamma - kgminus.*alpha;
            
            %d(alpha-prime)/dt
            kaplus.*RaG + kgminus.*alpha - (kgplus+kaminus.*Ra.*betagamma).*alphas;
            
            %d(beta-gamma)/dt
            kaplus.*RaG + kbminus.*G - kbplus.*alpha.*betagamma - kaminus.*Ra.*alphas.*betagamma];
           
end

% x(1) = [R]        x(2) = [R*]        x(3) = [G]        x(4) = [RG]
% x(5) = [R*G]        x(6) = [α]        x(7) = [α†]        x(8) = [βγ]
%
% C(1) = kq        C(2) = k_        C(3) = kp         C(4) = k+        
% C(5) = μ_        C(6) = θ{μ}        C(7) = ka+         C(8) = μ+    
% C(9) = ka-    C(10) = kb+    C(11) = kb-         C(12) = kg+
% C(13) = kg-   C(14) = kcon
