
clear
close all

%%%% Cooperative Enzyme Kinetics Example (2 Intermed Complexes) : GSA Working Group
%%%% Kate Pearce



%%%% nonzero initial species conc: units in molar concentration

s0 = 5e-7;
e0 = 2e-7;


%%%% "optimal" parameter values for toy problem

k1 = 3e5; 
k1_minus = 1e-3;  
k2 = 0.1; 
k3 = 9e5;
k3_minus = 1e-2;
k4 = 0.45;


%%%% initial values for ode solver
%%%% state vector: y = [p; e; s; c1; c2]
Y0 = [0; e0; s0; 0; 0];

params = [k1 k1_minus k2 k3 k3_minus k4]';
%%%% model solution: p = y(:,1)
[t, y] = rre_kinetics(params, s0, e0);
p = y(:,1);


%%%% Compute the sensitivities (d product/d par) using complex-step derivative approximations

h = 1e-3;
tfinal = 100;
tspan = 0 : 0.01 : tfinal; 
odeoptions = odeset('AbsTol',1e-10, 'RelTol', 1e-10, 'NonNegative',1);

k1_complex = complex(k1,h);
params = [k1_complex k1_minus k2 k3 k3_minus k4]';
[t,Y] = ode15s(@complex_rre_kinetics,tspan,Y0,odeoptions,params);
p_k1_comp = imag(Y(:,1))/h;

k1_minus_complex = complex(k1_minus,h);
params = [k1 k1_minus_complex k2 k3 k3_minus k4]';
[t,Y] = ode15s(@complex_rre_kinetics,tspan,Y0,odeoptions,params);
p_km1_comp = imag(Y(:,1))/h;

k2_complex = complex(k2,h);
params = [k1 k1_minus k2_complex k3 k3_minus k4]';
[t,Y] = ode15s(@complex_rre_kinetics,tspan,Y0,odeoptions,params);
p_k2_comp = imag(Y(:,1))/h;

k3_complex = complex(k3,h);
params = [k1 k1_minus k2 k3_complex k3_minus k4]';
[t,Y] = ode15s(@complex_rre_kinetics,tspan,Y0,odeoptions,params);
p_k3_comp = imag(Y(:,1))/h;

k3_minus_complex = complex(k3_minus,h);
params = [k1 k1_minus k2 k3 k3_minus_complex k4]';
[t,Y] = ode15s(@complex_rre_kinetics,tspan,Y0,odeoptions,params);
p_km3_comp = imag(Y(:,1))/h;

k4_complex = complex(k4,h);
params = [k1 k1_minus k2 k3 k3_minus k4_complex]';
[t,Y] = ode15s(@complex_rre_kinetics,tspan,Y0,odeoptions,params);
p_k4_comp = imag(Y(:,1))/h;


%%%%%% Compute the sensitivities (d product/d par) using finite difference approximations

k1_pert = k1+h;
params = [k1_pert k1_minus k2 k3 k3_minus k4]';
[t,Y] = ode15s(@complex_rre_kinetics,tspan,Y0,odeoptions,params); %%% complex step RHS is exactly same as fin diff RHS; NOT true for sens. eqs RHS 
p_k1_fin = Y(:,1);
p_k1_fin = (p_k1_fin - p)/h;

k1_minus_pert = k1_minus+h;
params = [k1 k1_minus_pert k2 k3 k3_minus k4]';
[t,Y] = ode15s(@complex_rre_kinetics,tspan,Y0,odeoptions,params);  
p_km1_fin = Y(:,1);
p_km1_fin = (p_km1_fin - p)/h;

k2_pert = k2+h;
params = [k1 k1_minus k2_pert k3 k3_minus k4]';
[t,Y] = ode15s(@complex_rre_kinetics,tspan,Y0,odeoptions,params);  
p_k2_fin = Y(:,1);
p_k2_fin = (p_k2_fin - p)/h;

k3_pert = k3+h;
params = [k1 k1_minus k2 k3_pert k3_minus k4]';
[t,Y] = ode15s(@complex_rre_kinetics,tspan,Y0,odeoptions,params);  
p_k3_fin = Y(:,1);
p_k3_fin = (p_k3_fin - p)/h;

k3_minus_pert = k3_minus+h;
params = [k1 k1_minus k2 k3 k3_minus_pert k4]';
[t,Y] = ode15s(@complex_rre_kinetics,tspan,Y0,odeoptions,params); 
p_km3_fin = Y(:,1);
p_km3_fin = (p_km3_fin - p)/h;

k4_pert = k4+h;
params = [k1 k1_minus k2 k3 k3_minus k4_pert]';
[t,Y] = ode15s(@complex_rre_kinetics,tspan,Y0,odeoptions,params);
p_k4_fin = Y(:,1);
p_k4_fin = (p_k4_fin - p)/h;




%%%% Compute the sensitivities (d product/d par) using sensitivity
%%%% equations 
clear Y0

Y0 = zeros(35,1);
Y0(2,1) = e0; 
Y0(3,1) = s0; 

k1 = 3e5; 
k1_minus = 1e-3;  
k2 = 0.1; 
k3 = 9e5;
k3_minus = 1e-2;
k4 = 0.45;

params = [k1 k1_minus k2 k3 k3_minus k4]';

tfinal = 100;
tspan = [0 : 0.01 : tfinal]; 
odeoptions = odeset('AbsTol',1e-10, 'RelTol', 1e-10, 'NonNegative',1);

[t, Y] = ode15s(@rre_senseq,tspan,Y0,odeoptions,params);

%%% extract states and sensitivities
P = Y(:,1);
E = Y(:,2);
S = Y(:,3);
C1 = Y(:,4);
C2 = Y(:,5);

P_k1 = Y(:,6);
E_k1 = Y(:,7);
S_k1 = Y(:,8);
C1_k1 = Y(:,9);
C2_k1 = Y(:,10);

P_km1 = Y(:,11);
E_km1 = Y(:,12);
S_km1 = Y(:,13);
C1_km1 = Y(:,14);
C2_km1 = Y(:,15);

P_k2 = Y(:,16);
E_k2 = Y(:,17);
S_k2 = Y(:,18);
C1_k2 = Y(:,19);
C2_k2 = Y(:,20);

P_k3 = Y(:,21);
E_k3 = Y(:,22);
S_k3 = Y(:,23);
C1_k3 = Y(:,24);
C2_k3 = Y(:,25);

P_km3 = Y(:,26);
E_km3 = Y(:,27);
S_km3 = Y(:,28);
C1_km3 = Y(:,29);
C2_km3 = Y(:,30);

P_k4 = Y(:,31);
E_k4 = Y(:,32);
S_k4 = Y(:,33);
C1_k4 = Y(:,34);
C2_k4 = Y(:,35);

figsens_k1 = figure();
plot(t, p_k1_comp, '-r', t, p_k1_fin, '--k', t, P_k1, '-.b', 'linewidth', 3)
set(gca,'Fontsize',[22]);
xlabel('Time')
ylabel('P_{k_1}')
legend('Complex-Step','Finite Diff','Sens Eq','Location','SouthEast')

figsens_km1 = figure();
plot(t, p_km1_comp, '-r', t, p_km1_fin, '--k', t, P_km1, '-.b', 'linewidth', 3)
set(gca,'Fontsize',[22]);
xlabel('Time')
ylabel('P_{k_{-1}}')
legend('Complex-Step','Finite Diff','Sens Eq','Location','SouthEast')

figsens_k2 = figure();
plot(t, p_k2_comp, '-r', t, p_k2_fin, '-k', t, P_k2, '-.b', 'linewidth', 3)
set(gca,'Fontsize',[22]);
xlabel('Time')
ylabel('P_{k_2}')
legend('Complex-Step','Finite Diff','Sens Eq','Location','SouthEast')

figsens_k3 = figure();
plot(t, p_k3_comp, '-r', t, p_k3_fin, '-k', t, P_k3, '-.b', 'linewidth', 3)
set(gca,'Fontsize',[22]);
xlabel('Time')
ylabel('P_{k_3}')
legend('Complex-Step','Finite Diff','Sens Eq','Location','SouthEast')

figsens_km3 = figure();
plot(t, p_km3_comp, '-r', t, p_km3_fin, '-k', t, P_km3, '-.b', 'linewidth', 3)
set(gca,'Fontsize',[22]);
xlabel('Time')
ylabel('P_{{k}_{-3}}')
legend('Complex-Step','Finite Diff','Sens Eq','Location','SouthEast')

figsens_k4 = figure();
plot(t, p_k4_comp, '-r', t, p_k4_fin, '-k', t, P_k4, '-.b', 'linewidth', 3)
set(gca,'Fontsize',[22]);
xlabel('Time')
ylabel('P_{k_4}')
legend('Complex-Step','Finite Diff','Sens Eq','Location','SouthEast')




%%% relative error forward diff/complex-step 
% ForCompErr_k1 = max(abs(p_k1_comp-p_k1_fin))/max(abs(p_k1_comp));
% ForCompErr_km1 = max(abs(p_km1_comp-p_km1_fin))/max(abs(p_km1_comp));
% ForCompErr_k2 = max(abs(p_k2_comp-p_k2_fin))/max(abs(p_k2_comp));
% ForCompErr_k3 = max(abs(p_k3_comp-p_k3_fin))/max(abs(p_k3_comp));
% ForCompErr_km3 = max(abs(p_km3_comp-p_km3_fin))/max(abs(p_km3_comp));
% ForCompErr_k4 = max(abs(p_k4_comp-p_k4_fin))/max(abs(p_k4_comp));
k1_ErrForComp = max(abs(abs(p_k1_comp)-abs(p_k1_fin)))/max(abs(p_k1_comp));
k1_minus_ErrForComp = max(abs(abs(p_km1_comp)-abs(p_km1_fin)))/max(abs(p_km1_comp));
k2_ErrForComp = max(abs(abs(p_k2_comp)-abs(p_k2_fin)))/max(abs(p_k2_comp));
k3_ErrForComp = max(abs(abs(p_k3_comp)-abs(p_k3_fin)))/max(abs(p_k3_comp));
k3_minus_ErrForComp = max(abs(abs(p_km3_comp)-abs(p_km3_fin)))/max(abs(p_km3_comp));
k4_ErrForComp = max(abs(abs(p_k4_comp)-abs(p_k4_fin)))/max(abs(p_k4_comp));

% SensCompErr_k1 = max(abs(p_k1_comp - P_k1))/max(abs(p_k1_comp));
% SensCompErr_km1 = max(abs(p_km1_comp - P_km1))/max(abs(p_km1_comp));
% SensCompErr_k2 = max(abs(p_k2_comp - P_k2))/max(abs(p_k2_comp));
% SensCompErr_k3 = max(abs(p_k3_comp - P_k3))/max(abs(p_k3_comp));
% SensCompErr_km3 = max(abs(p_km3_comp - P_km3))/max(abs(p_km3_comp));
% SensCompErr_k4 = max(abs(p_k4_comp - P_k4))/max(abs(p_k4_comp));
k1_ErrSensComp = max(abs(abs(p_k1_comp) - abs(P_k1)))/max(abs(p_k1_comp));
k1_minus_ErrSensComp = max(abs(abs(p_km1_comp) - abs(P_km1)))/max(abs(p_km1_comp));
k2_ErrSensComp = max(abs(abs(p_k2_comp) - abs(P_k2)))/max(abs(p_k2_comp));
k3_ErrSensComp = max(abs(abs(p_k3_comp) - abs(P_k3)))/max(abs(p_k3_comp));
k3_minus_ErrSensComp = max(abs(abs(p_km3_comp) - abs(P_km3)))/max(abs(p_km3_comp));
k4_ErrSensComp = max(abs(abs(p_k4_comp) - abs(P_k4)))/max(abs(p_k4_comp));


k1_ErrSensFor = max(abs(p_k1_fin - P_k1))/max(abs(P_k1));
k1_minus_ErrSensFor = max(abs(p_km1_fin - P_km1))/max(abs(P_km1));
k2_ErrSensFor = max(abs(p_k2_fin - P_k2))/max(abs(P_k2));
k3_ErrSensFor = max(abs(p_k3_fin - P_k3))/max(abs(P_k3));
k3_minus_ErrSensFor = max(abs(p_km3_fin - P_km3))/max(abs(P_km3));
k4_ErrSensFor = max(abs(p_k4_fin - P_k4))/max(abs(P_k4));

clearvars k1_pert k2_pert k3_pert k4_pert k1_minus_pert k1_complex k1_minus_complex k2_complex k3_complex k3_minus_complex k3_minus_pert k4_complex 