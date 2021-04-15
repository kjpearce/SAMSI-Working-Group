%%%%% Comparison of Finite Differences, Complex-Step Approximations, and
%%%%% Sensitivity Equations for GSA Working Group 

%%%%% Author: Kate Pearce

%%%%% Cooperative Enzyme Model with num intermed complex = 2:
%%%%% Reaction System:

%%%%%       kp1    kp2
%%%%% S + E <-> C1 -> P + E
%%%%%       km1

%%%%%         kp3    kp4
%%%%% S + C1 <-> C2 -> P + E
%%%%%         km3

clear
close all

%%% step size(s) for comparisons: can be a vector
step = 1e-10;


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

paramnames = {'k^{+}_1','k^{-}_1','k^{+}_2','k^{+}_3','k^{-}_3','k^{+}_4'};

%%%% IC for ode solver(s)
%%%% state vector: y = [p; e; s; c1; c2]

%%% for model soln, finite diff, comp step
Y0 = [0; e0; s0; 0; 0]; 

%%% for sens eq 
Y0_sens = zeros(35,1); 
Y0_sens(2,1) = e0; 
Y0_sens(3,1) = s0; 

%%% nominal parameter vector
params = [k1; k1_minus; k2; k3; k3_minus; k4];

%%%% model solution at nominal vals: p = y(:,1) from Alen's code
[~, y] = rre_kinetics(params, s0, e0);
p = y(:,1);


tfinal = 100;
tspan = 0 : 0.01 : tfinal; 
odeoptions = odeset('AbsTol',1e-10, 'RelTol', 1e-10);


%%% how many step sizes
numsteps = length(step);

for stepind = 1:numsteps

    %%% step size
    h = step(stepind);
    
    numpar = length(params);
    numpts = length(tspan);

    %%%% Compute the sensitivities (d product/d par) using complex-step derivative approximations
    Chi_CS = zeros(numpts, numpar);

    for parind = 1:numpar
        %%% temp parameter vector to perturb
        partemp = params;

        %%% parameter to perturb
        par_real = partemp(parind); 
        %%% complex perturbation
        par_pert = complex(par_real,h); 

        %%% full parameter vector with complex pert
        partemp(parind) = par_pert; 

        [~,Y] = ode15s(@complex_rre_kinetics,tspan,Y0,odeoptions,partemp);

        %%% Y(:,1) model soln for reaction product 
        sens = imag(Y(:,1))/h;

        %%% store parameter sensitivity in corresp column
        Chi_CS(:,parind) = sens;
    end





    %%%%%% Compute the sensitivities (d product/d par) using finite difference approximations
    
    Chi_FD = zeros(numpts, numpar);
    
    for parind = 1:numpar
        %%% temp parameter vector to perturb
        partemp = params;

        %%% parameter to perturb
        par_real = partemp(parind); 
        %%% finite diff perturbation
        par_pert = par_real + h; 

        %%% full parameter vector w pert
        partemp(parind) = par_pert; 

        %%% RHS the same for FD as CS
        [~,Y] = ode15s(@complex_rre_kinetics,tspan,Y0,odeoptions,partemp);

        %%% Y(:,1) model soln for reaction product 
        sens = (Y(:,1)-p)./h;

        %%% store parameter sensitivity in corresp column
        Chi_FD(:,parind) = sens;
    end
    
    
    %%%% Compute the sensitivities (d product/d par) using sensitivity
    %%%% equations 
    [~, Y] = ode15s(@rre_senseq,tspan,Y0_sens,odeoptions,params);

    %%% extract sens for product
    P_k1 = Y(:,6);
    P_km1 = Y(:,11);
    P_k2 = Y(:,16);
    P_k3 = Y(:,21);
    P_km3 = Y(:,26);
    P_k4 = Y(:,31);


    Chi_SE = [P_k1 P_km1 P_k2 P_k3 P_km3 P_k4];

%     %%%%% Plot all sensitivities together
%     for parind = 1:numpar
%         figure()
%         plot(tspan, Chi_CS(:,parind), '-k',tspan, Chi_FD(:,parind), '-.r', tspan, Chi_SE(:,parind), ':b','LineWidth',3)
%         set(gca, 'FontSize', 22)
%         xlabel('Time')
%         ylabel(sprintf('Sensitivity wrt %s',paramnames{parind}))
%         legend('Complex-Step','Forward Diff','Sensivity Eq','Location','SouthEast')
%     end

    %%%%% Plot only complex-step and sensitivity equations together
    for parind = 1:numpar
        figure()
        plot(tspan, Chi_CS(:,parind), '-k', tspan, Chi_SE(:,parind), ':b','LineWidth',3)
        set(gca, 'FontSize', 22)
        xlabel('Time')
        ylabel(sprintf('Sensitivity wrt %s',paramnames{parind}))
        legend('Complex-Step','Sensivity Eq','Location','SouthEast')
    end


%     %%%%% Compute relative errors between forward diff and complex
%     %%%%% step
%     ForCompErr_k1 = max(abs(p_k1_comp-p_k1_fin))/max(abs(p_k1_comp));
%     ForCompErr_km1 = max(abs(p_km1_comp-p_km1_fin))/max(abs(p_km1_comp));
%     ForCompErr_k2 = max(abs(p_k2_comp-p_k2_fin))/max(abs(p_k2_comp));
%     ForCompErr_k3 = max(abs(p_k3_comp-p_k3_fin))/max(abs(p_k3_comp));
%     ForCompErr_km3 = max(abs(p_km3_comp-p_km3_fin))/max(abs(p_km3_comp));
%     ForCompErr_k4 = max(abs(p_k4_comp-p_k4_fin))/max(abs(p_k4_comp));
% 
%     %%%%% Compute relative errors between sens eqs and complex
%     %%%%% step
%     SensCompErr_k1 = max(abs(p_k1_comp - P_k1))/max(abs(p_k1_comp));
%     SensCompErr_km1 = max(abs(p_km1_comp - P_km1))/max(abs(p_km1_comp));
%     SensCompErr_k2 = max(abs(p_k2_comp - P_k2))/max(abs(p_k2_comp));
%     SensCompErr_k3 = max(abs(p_k3_comp - P_k3))/max(abs(p_k3_comp));
%     SensCompErr_km3 = max(abs(p_km3_comp - P_km3))/max(abs(p_km3_comp));
%     SensCompErr_k4 = max(abs(p_k4_comp - P_k4))/max(abs(p_k4_comp));
% 
%     %%%%% Compute relative errors between sens eqs and forward diff
%     SensForErr_k1 = max(abs(p_k1_fin - P_k1))/max(abs(P_k1));
%     SensForErr_km1 = max(abs(p_km1_fin - P_km1))/max(abs(P_km1));
%     SensForErr_k2 = max(abs(p_k2_fin - P_k2))/max(abs(P_k2));
%     SensForErr_k3 = max(abs(p_k3_fin - P_k3))/max(abs(P_k3));
%     SensForErr_km3 = max(abs(p_km3_fin - P_km3))/max(abs(P_km3));
%     SensForErr_k4 = max(abs(p_k4_fin - P_k4))/max(abs(P_k4));

end

