
%
%          SIR_rhs
%
  function dy = rre_senseq(t,y,params)
dy = zeros(35,1);

k1       = params(1);
km1      = params(2);
k2       = params(3);
k3       = params(4);
km3      = params(5);
k4       = params(6);

P = y(1);
E = y(2);
S = y(3);
C1 = y(4);
C2 = y(5);

P_k1 = y(6);
E_k1 = y(7);
S_k1 = y(8);
C1_k1 = y(9);
C2_k1 = y(10);

P_km1 = y(11);
E_km1 = y(12);
S_km1 = y(13);
C1_km1 = y(14);
C2_km1 = y(15);

P_k2 = y(16);
E_k2 = y(17);
S_k2 = y(18);
C1_k2 = y(19);
C2_k2 = y(20);

P_k3 = y(21);
E_k3 = y(22);
S_k3 = y(23);
C1_k3 = y(24);
C2_k3 = y(25);

P_km3 = y(26);
E_km3 = y(27);
S_km3 = y(28);
C1_km3 = y(29);
C2_km3 = y(30);

P_k4 = y(31);
E_k4 = y(32);
S_k4 = y(33);
C1_k4 = y(34);
C2_k4 = y(35);


Svec = [P_k1; E_k1; S_k1; C1_k1; C2_k1; P_km1; E_km1; S_km1; C1_km1; C2_km1; P_k2; E_k2; S_k2; C1_k2; C2_k2; P_k3; E_k3; S_k3; C1_k3; C2_k3; P_km3; E_km3; S_km3; C1_km3; C2_km3; P_k4; E_k4; S_k4; C1_k4; C2_k4] ; 


J = [ 0 0 0 k2 k4 ; 0 -k1*S -k1*E (km1+k2) 0 ; 0 -k1*S (-k1*E-k3*C1) (km1-k3*S) km3 ; 0 k1*S (k1*E-k3*C1) (-(km1+k2)-k3*S) (k4+km3) ; 0 0 k3*C1 k3*S -(k4+km3) ] ;
Der = blkdiag(J,J,J,J,J,J); 

Grad_k1 = [ 0; S*E; -S*E; S*E; 0 ];
Grad_km1 = [ 0; C1; C1; -C1; 0 ];
Grad_k2 = [ C1; C1; 0; -C1; 0 ];
Grad_k3 = [ 0; 0; -S*C1; -S*C1; S*C1 ];
Grad_km3 = [ 0; 0; C2; C2; -C2];
Grad_k4 = [ C2; 0; 0; C2; -C2];

Grad = [Grad_k1; Grad_km1; Grad_k2; Grad_k3; Grad_km3; Grad_k4 ];

Sen = Der*Svec + Grad;


%%% Construct right-hand side
dy = [ k2 * C1 + k4 * C2; %dP
        (km1 + k2) * C1 - k1 * S * E; %dE
        -k1 * S * E + km1 * C1 - k3 * S * C1 + km3 * C2; %dS
        k1 * S * E - (km1 + k2) * C1 - k3 * S * C1 + (k4 + km3) * C2; %dC1
        k3 * S * C1 - (k4 + km3) * C2; %dC2
        Sen];





