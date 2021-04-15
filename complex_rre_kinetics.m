
%
%          SIR_rhs_complex
%
  function dy = complex_rre_kinetics(t,y,params)


% rate constants
k1       = params(1);
k1_minus = params(2);
k2       = params(3);
k3       = params(4);
k3_minus = params(5);
k4       = params(6);

% system species
p  = y(1);
e  = y(2);
s  = y(3);
c1 = y(4);
c2 = y(5);

%          SIR_rhs_complex

%
  dy = [ k2 * c1 + k4 * c2;
        (k1_minus + k2) * c1 - k1 * s * e;
        -k1 * s * e + k1_minus * c1 - k3 * s * c1 + k3_minus * c2; 
        k1 * s * e - (k1_minus + k2) * c1 - k3 * s * c1 + (k4 + k3_minus) * c2;
        k3 * s * c1 - (k4 + k3_minus) * c2];


