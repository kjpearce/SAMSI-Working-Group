% 
% clear
% close all
% 
% 
% % species units in molar concentration
% 
% s0 = 5e-7;
% e0 = 2e-7;
% 
% k1 = 3e5; 
% k1_minus = 1e-3;  
% k2 = 0.1; 
% k3 = 9e5;
% k3_minus = 1e-2;
% k4 = 0.45;
% 
% 
% 
% params = [k1 k1_minus k2 k3 k3_minus k4]';
% [t y] = rre_kinetics(sqrt(params), s0, e0);
% 
% % %    note: the state vector is y = [p e s c1 c2]'
% 
% p = y(:, 1); 
% e = y(:, 2); 
% s = y(:, 3);
% c1 = y(:, 4);
% c2 = y(:, 5);
% 
% figure()
% plot(t, s, t, p, t, e, 'linewidth', 2);
% legend('substrate', 'product', 'enzyme');
% xlabel('time');
% set(gca, 'fontsize', 20);
% 
% figure()
% plot(t, c1, t, c2, 'linewidth', 2);
% legend('complex 1', 'complex 2');
% xlabel('time');
% set(gca, 'fontsize', 20);
% 
% 
% 
% 
% %%%%% generate data 
% down = 100; %%% downsample model response by 100
% timevec = t(1:down:end);
% numpts = length(timevec);
% pvec = zeros(1,numpts);
% 
% %%%% create fixed but known variance 
% varn = 1e-16;
% 
% for pt = 1:numpts
%     pvec(pt) = p(1+down*(pt-1)); %%% could add noise here if we wanted
% end
% 
% %%%% error
% noise = sqrt(varn)*randn(size(pvec));
% 
% %%%% DATAVEC IS OUR VECTOR WITH ADDED NOISE
% datavec = pvec + noise;
% datavec(1)=0;
% 
% q0 = [1e5, 1, 1, 1e5, 1, 1];
% qinit = sqrt(q0);
% 
% [qopt, fmin] = fminroutine(datavec, timevec, s0, e0, qinit);

[T,OptCurves] = rre_rhs(qopt,t,s0,e0);
soln_p = OptCurves(:,1);

s0 = 5e-7;
e0 = 2e-7;

h = 1e-16;

k1 = qopt(1);
km1 = qopt(2);

k1_complex = complex(k1,h);
params = [gamma r delta_complex];
[t,Y] = ode45(@SIR_rhs_complex,t_vals,Y0,ode_options,params);
S_delta = imag(Y(:,1))/h;
I_delta = imag(Y(:,2))/h;

gamma_complex = complex(gamma,h);
params = [gamma_complex r delta];
[t,Y] = ode45(@SIR_rhs_complex,t_vals,Y0,ode_options,params);
S_gamma = imag(Y(:,1))/h;
I_gamma = imag(Y(:,2))/h;

r_complex = complex(r,h);
params = [gamma r_complex delta];
[t,Y] = ode45(@SIR_rhs_complex,t_vals,Y0,ode_options,params);
S_r = imag(Y(:,1))/h;
I_r = imag(Y(:,2))/h;

%
% Compute the relative errors.
% 

S_delta_err = max(abs(S_delta_sen-S_delta))/max(abs(S_delta_sen));
I_delta_err = max(abs(I_delta_sen-I_delta))/max(abs(I_delta_sen));
S_gamma_err = max(abs(S_gamma_sen-S_gamma))/max(abs(S_gamma_sen));
I_gamma_err = max(abs(I_gamma_sen-I_gamma))/max(abs(I_gamma_sen));
S_r_err = max(abs(S_r_sen-S_r))/max(abs(S_r_sen));
I_r_err = max(abs(I_r_sen-I_r))/max(abs(I_r_sen));








Nparams = 6;
SolnMatrix = zeros(numpts,Nparams);
eps = 0.005;

h = 1e-4;
for par = 1:Nparams
    qoptcomplex = qopt;
    qoptcomplex(par) = qoptcomplex(par) + 1i*h;
    NewSolCurvesComplex = rre_kinetics(qoptcomplex,s0,e0);
    SolnMatrixComplex(:,par)= NewSolCurvesComplex(:,1);
end

% 
% for par=1:Nparams
%     qoptpert = qopt;
%     qoptpert(par) = qoptpert(par)*(1+eps);
%     NewSolCurves = rre_kinetics(qoptpert,s0,e0);
%     SolnMatrix(:,par) = NewSolCurves(:,1);
% end
% 
% 
% %%%%% complex step 
% 
% 
% 
% 
% 
%     
% %%%%%% build finite diff approx: scaled length(timevec) x 6 sensitivity matrix dpdq, evaluated
% %%%%%% at timevec points
% dpdq = zeros(length(timevec),length(params));
% for i = 1:length(timevec)
%     dpdq(i,1)= (p_k1(timevec(i)) - p(timevec(i)))/(eps*k1);
%     dpdq(i,2)= (p_km1(timevec(i)) - p(timevec(i)))/(eps*k1_minus);
%     dpdq(i,3)= (p_k2(timevec(i)) - p(timevec(i)))/(eps*k2);
%     dpdq(i,4)= (p_k3(timevec(i)) - p(timevec(i)))/(eps*k3);
%     dpdq(i,5)= (p_km3(timevec(i)) - p(timevec(i)))/(eps*k3_minus);
%     dpdq(i,6)= (p_k4(timevec(i)) - p(timevec(i)))/(eps*k4);
% 
%     dpdq(i,1)=dpdq(i,1)*0.2*k1;
%     dpdq(i,2)=dpdq(i,2)*0.2*k1_minus;
%     dpdq(i,3)=dpdq(i,3)*0.2*k2;
%     dpdq(i,4)=dpdq(i,4)*0.2*k3; 
%     dpdq(i,5)=dpdq(i,5)*0.2*k3_minus;
%     dpdq(i,6)=dpdq(i,6)*0.2*k4;
% 
% end
% 
% %%%%%% construct Fisher information matrix
% FIM = dpdq'*dpdq;
% 
% %%%%%% eigenvec/vals of FIM
% [V, lambda] = eig(FIM);
% 
% logEvals = zeros(1,length(params));
% for i=1:length(params)
%     logEvals(i)=log10(abs(lambda(i,i)));
% end
% 
% %%%%%%% order the eigenvalues smallest to largest
% [log_sort,ind] = sort(logEvals);
% V_sort = V(:,ind);
% 
% 
% 
% 
% figure()
% bar(log_sort)
% ylabel('|Eigenvalue| (log scale)')
% xticks(1:1:6)
% xticklabels({'\lambda_1','\lambda_2','\lambda_3','\lambda_4','\lambda_5','\lambda_6'})
% title('Eigenvalues of FIM')
% ax = gca;
% ax.FontSize = 20;
% 
% figure()
% bar(V_sort(:,1),'grouped')
% title('First eigenvector')
% xticks(1:1:6)
% xticklabels({'k_1','k_{-1}','k_2','k_3','k_{-3}','k_4'})
% ax = gca;
% ax.FontSize = 20;
% 
% figure()
% bar(V_sort(:,2),'grouped')
% title('Second eigenvector')
% xticks(1:1:6)
% xticklabels({'k_1','k_{-1}','k_2','k_3','k_{-3}','k_4'})
% ax = gca;
% ax.FontSize = 20;
% 
% 
