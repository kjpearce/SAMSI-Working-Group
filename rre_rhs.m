function [T,y] = rre_rhs(params,t, s0, e0)
%
% ODE15s solution the Reaction Rate Equation for 
% the cooperative enzyme kinetics problem
%
% programmer: Alen Alexanderian 
%
% input: 
%    params --- vector of rate parameters 
%
%         params(1) = k1
%         params(2) = k1_minus
%         params(3) = k2
%         params(4) = k3
%         params(5) = k3_minus
%         params(6) = k4
%
%    s0 --- initial substrate concentration
%    e0 --- initial enzyme concentration
%  
%    output:
%         t --- vector of time steps
%         y --- state vector y(t)
%    note: the state vector is y = [p e s c1 c2]'
%          p      --- product
%          e      --- enzyme 
%          s      --- substrate
%          c1, c2 --- complex species


% initial state
yzero = zeros(5, 1);
yzero(1) = 0;   
yzero(2) = e0;   
yzero(3) = s0;
yzero(4) = 0;
yzero(5) = 0; 

options = odeset('AbsTol',1e-10, 'RelTol', 1e-10, 'NonNegative',1);

frhs = @(t, y)(rre_kinetics_rhs(t, y, params));
 
[T,y] = ode15s(frhs, t, yzero, options);

%-------------- sub-function ----------
function yprime = rre_kinetics_rhs(t, y, params)
yprime = zeros(5, 1);

% rate constants
k1       = params(1)^2;
k1_minus = params(2)^2;
k2       = params(3)^2;
k3       = params(4)^2;
k3_minus = params(5)^2;
k4       = params(6)^2;

% system species
p  = y(1);
e  = y(2);
s  = y(3);
c1 = y(4);
c2 = y(5);

yprime(1) = k2 * c1 + k4 * c2; 

yprime(2) = (k1_minus + k2) * c1 - k1 * s * e;

yprime(3) = -k1 * s * e + k1_minus * c1 - k3 * s * c1 ...
            + k3_minus * c2; 

yprime(4) = k1 * s * e - (k1_minus + k2) * c1 - k3 * s * c1 ... 
            + (k4 + k3_minus) * c2; 

yprime(5) = k3 * s * c1 - (k4 + k3_minus) * c2; 

end
end