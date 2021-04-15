%%%%% Comparison of Finite Differences, Complex-Step Approximations, and
%%%%% Sensitivity Equations for GSA Working Group 'Toy Problem'

%%%%% Kate Pearce

%%%%% Cooperative Enzyme Model with num intermed complex = 2



clear
close all

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
Y0 = [0; e0; s0; 0; 0]; %%% for model soln, finite diff, comp step

Y0_sens = zeros(35,1); %%% for sens eq initial conditions
Y0_sens(2,1) = e0; 
Y0_sens(3,1) = s0; 


params = [k1 k1_minus k2 k3 k3_minus k4]';
%%%% model solution: p = y(:,1)
[t, y] = rre_kinetics(params, s0, e0);
p = y(:,1);


tfinal = 100;
tspan = 0 : 0.01 : tfinal; 
odeoptions = odeset('AbsTol',1e-10, 'RelTol', 1e-10, 'NonNegative',1);

ordmag = 17;
hvec = zeros(ordmag,1);

FinDiff = struct;
CompStep = struct;
SensEqs = struct;

for i=1:ordmag
    
    h = 10^(-(i-1)); %%% order of magnitude 
    
    hvec(i) = h;
    %%%% Compute the sensitivities (d product/d par) using complex-step derivative approximations
    
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
    
    CompStep(i).k1 = p_k1_comp;
    CompStep(i).km1 = p_km1_comp;
    CompStep(i).k2 = p_k2_comp;
    CompStep(i).k3 = p_k3_comp;
    CompStep(i).km3 = p_km3_comp;
    CompStep(i).k4 = p_k4_comp;

    

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

    FinDiff(i).k1 = p_k1_fin;
    FinDiff(i).km1 = p_km1_fin;
    FinDiff(i).k2 = p_k2_fin;
    FinDiff(i).k3 = p_k3_fin;
    FinDiff(i).km3 = p_km3_fin;
    FinDiff(i).k4 = p_k4_fin;
    


    %%%% Compute the sensitivities (d product/d par) using sensitivity
    %%%% equations 

    [t, Y] = ode15s(@rre_senseq,tspan,Y0_sens,odeoptions,params);

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
    
    SensEqs(i).k1 = P_k1;
    SensEqs(i).km1 = P_km1;
    SensEqs(i).k2 = P_k2;
    SensEqs(i).k3 = P_k3;
    SensEqs(i).km3 = P_km3;
    SensEqs(i).k4 = P_k4;
    
    
    
    %%%%% Plot all sensitivities together
    
    figsens_k1 = figure();
    plot(t, p_k1_comp, '-r', t, p_k1_fin, '-k', t, P_k1, '-b', 'linewidth', 3)
    set(gca,'Fontsize',[22]);
    xlabel('Time')
    ylabel('P_{k_1}')
    legend('Complex-Step','Finite Diff','Sens Eq','Location','SouthEast')

    figsens_km1 = figure();
    plot(t, p_km1_comp, '-r', t, p_km1_fin, '-k', t, P_km1, '-b', 'linewidth', 3)
    set(gca,'Fontsize',[22]);
    xlabel('Time')
    ylabel('P_{k_{-1}}')
    legend('Complex-Step','Finite Diff','Sens Eq','Location','SouthEast')

    figsens_k2 = figure();
    plot(t, p_k2_comp, '-r', t, p_k2_fin, '-k', t, P_k2, '-b', 'linewidth', 3)
    set(gca,'Fontsize',[22]);
    xlabel('Time')
    ylabel('P_{k_2}')
    legend('Complex-Step','Finite Diff','Sens Eq','Location','SouthEast')

    figsens_k3 = figure();
    plot(t, p_k3_comp, '-r', t, p_k3_fin, '-k', t, P_k3, '-b', 'linewidth', 3)
    set(gca,'Fontsize',[22]);
    xlabel('Time')
    ylabel('P_{k_3}')
    legend('Complex-Step','Finite Diff','Sens Eq','Location','SouthEast')

    figsens_km3 = figure();
    plot(t, p_km3_comp, '-r', t, p_km3_fin, '-k', t, P_km3, '-b', 'linewidth', 3)
    set(gca,'Fontsize',[22]);
    xlabel('Time')
    ylabel('P_{{k}_{-3}}')
    legend('Complex-Step','Finite Diff','Sens Eq','Location','SouthEast')

    figsens_k4 = figure();
    plot(t, p_k4_comp, '-r', t, p_k4_fin, '-k', t, P_k4, '-b', 'linewidth', 3)
    set(gca,'Fontsize',[22]);
    xlabel('Time')
    ylabel('P_{k_4}')
    legend('Complex-Step','Finite Diff','Sens Eq','Location','SouthEast')
    
    

    %%%%% Compute relative errors between forward diff and complex
    %%%%% step
    ForCompErr_k1 = max(abs(p_k1_comp-p_k1_fin))/max(abs(p_k1_comp));
    ForCompErr_km1 = max(abs(p_km1_comp-p_km1_fin))/max(abs(p_km1_comp));
    ForCompErr_k2 = max(abs(p_k2_comp-p_k2_fin))/max(abs(p_k2_comp));
    ForCompErr_k3 = max(abs(p_k3_comp-p_k3_fin))/max(abs(p_k3_comp));
    ForCompErr_km3 = max(abs(p_km3_comp-p_km3_fin))/max(abs(p_km3_comp));
    ForCompErr_k4 = max(abs(p_k4_comp-p_k4_fin))/max(abs(p_k4_comp));
    
    %%%%% Compute relative errors between sens eqs and complex
    %%%%% step
    SensCompErr_k1 = max(abs(p_k1_comp - P_k1))/max(abs(p_k1_comp));
    SensCompErr_km1 = max(abs(p_km1_comp - P_km1))/max(abs(p_km1_comp));
    SensCompErr_k2 = max(abs(p_k2_comp - P_k2))/max(abs(p_k2_comp));
    SensCompErr_k3 = max(abs(p_k3_comp - P_k3))/max(abs(p_k3_comp));
    SensCompErr_km3 = max(abs(p_km3_comp - P_km3))/max(abs(p_km3_comp));
    SensCompErr_k4 = max(abs(p_k4_comp - P_k4))/max(abs(p_k4_comp));
    
    %%%%% Compute relative errors between sens eqs and forward diff
    SensForErr_k1 = max(abs(p_k1_fin - P_k1))/max(abs(P_k1));
    SensForErr_km1 = max(abs(p_km1_fin - P_km1))/max(abs(P_km1));
    SensForErr_k2 = max(abs(p_k2_fin - P_k2))/max(abs(P_k2));
    SensForErr_k3 = max(abs(p_k3_fin - P_k3))/max(abs(P_k3));
    SensForErr_km3 = max(abs(p_km3_fin - P_km3))/max(abs(P_km3));
    SensForErr_k4 = max(abs(p_k4_fin - P_k4))/max(abs(P_k4));
    
    
    
    %%% store errors
    CompStep(i).k1_senserr = SensCompErr_k1;
    CompStep(i).km1_senserr = SensCompErr_km1;
    CompStep(i).k2_senserr = SensCompErr_k2;
    CompStep(i).k3_senserr = SensCompErr_k3;
    CompStep(i).km3_senserr = SensCompErr_km3;
    CompStep(i).k4_senserr = SensCompErr_k4;
    
    
    CompStep(i).k1_forerr = ForCompErr_k1;
    CompStep(i).km1_forerr = ForCompErr_km1;
    CompStep(i).k2_forerr = ForCompErr_k2;
    CompStep(i).k3_forerr = ForCompErr_k3;
    CompStep(i).km3_forerr = ForCompErr_km3;
    CompStep(i).k4_forerr = ForCompErr_k4;
    
    SensEqs(i).k1_forerr = SensForErr_k1;
    SensEqs(i).km1_forerr = SensForErr_km1;
    SensEqs(i).k2_forerr = SensForErr_k2;
    SensEqs(i).k3_forerr = SensForErr_k3;
    SensEqs(i).km3_forerr = SensForErr_km3;
    SensEqs(i).k4_forerr = SensForErr_k4;
    
end

