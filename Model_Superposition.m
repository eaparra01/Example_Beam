%
%   Dynamic Response to Harmonic Transverse Excitation of Cantilever
%   Euler-Bernoulli Beam Carrying a Point Mass Method.
%   This method uses moda shape and modal coordinates functions which are
%   calculated by seperation of varianbles method and Laplace
%   Transformation with Clamped-Free Configuration
%   Calculations for this program by Edgar Andres Parra Ricaurte
%

% clear
syms x k L xm M E Inertia omega omega_xm rho Area u U real


L = 0.732;
E = 70e+9;
height =  0.015;
width = 0.00302;
Area = width*height;
Inertia = (height*(width^3))/12;
rho = 2700 ;
mass = rho*Area*L;
xm = L ;                  % Mass position in tne beam
M = 0.58;                 % Mass


A = cos(k*L)+cosh(k*L);
B = sin(k*L)+sinh(k*L);
C = sin(k*L-k*xm)+sinh(k*L-k*xm);


D = sinh(k*L)-sin(k*L);
F = cos(k*L)+cosh(k*L);
G = cos(k*xm-k*L)+cosh(k*xm-k*L);

R = -(B*G-C*F)/(B*D-A*F);
S = (A*G-C*D)/(B*D-A*F);

fx = M*(k^4)*((E*Inertia)/(rho*Area))*(R*(cosh(k*xm) - cos(k*xm)) + ...
     S*(sinh(k*xm) - sin(k*xm))) - 2*E*Inertia*(k^3);

v=1;
values = 16;
ks1=zeros(values,1);
for n = 1:values
    result = vpasolve(fx == 0, k ,n);
    if isempty( result) ~= true
        ks1(v,1) = abs(result);
        v = v+1;
    end
end

ks1 = sort(double(ks1));
v=1;
for n = 1:length(ks1)
    if(n==1)
        ks(v,1) = ks1(n);
        v = v+1;
    elseif(v>1 && ks1(n-1)~=ks1(n))
        ks(v,1) = ks1(n);
        v = v+1;
    end
end

omega = sqrt((ks.^(4))*E*Inertia/(rho*Area))
hz = omega/(2*pi)
