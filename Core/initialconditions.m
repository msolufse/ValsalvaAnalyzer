function Init = initialconditions(pars,data)
% function initialconditions.m
% Input: model parameters and data structure
% Output: initial conditions for differential equations

% Model parameters
alpha  = pars(1); 
beta   = pars(2);
A      = pars(3);
    
K_P    = pars(4);
K_b    = pars(5);
K_pb   = pars(6);
K_pr   = pars(7);
K_s    = pars(8);
    
q_w    = pars(15);
q_pb   = pars(16);
q_pr   = pars(17);
q_s    = pars(18);
    
xi_w   = pars(19);
xi_pb  = pars(20);
xi_pr  = pars(21);
xi_s   = pars(22);
    
HI     = pars(23);
    
% Calculate initial conditions ensuring that the model starts in steady
% state
Pbar      = alpha * data.SPbar + (1 - alpha) * data.PPbar; 
dPdt_bar  = alpha * data.dSBPdt_bar + (1 - alpha) * data.dPPdt_bar; 
dPadt_bar = dPdt_bar - data.dPthdt_bar; 
    
Pc_m_0 = Pbar + K_P * dPdt_bar; 
Pa_m_0 = (Pbar - data.Pthbar) + K_P * dPadt_bar; % Should there be a K_P
    
y_c_0 = -q_w * (Pc_m_0 - xi_w); 
y_a_0 = -q_w * (Pa_m_0 - xi_w); 
    
e_wc_0 = 1 - sqrt((1 + exp(y_c_0))/(A + exp(y_c_0))); 
e_wa_0 = 1 - sqrt((1 + exp(y_a_0))/(A + exp(y_a_0))); 
    
e_bc_0 = K_b * e_wc_0; 
e_ba_0 = K_b * e_wa_0;  
    
n_0   = beta*(e_wc_0 - e_bc_0) + (1 - beta)*(e_wa_0 - e_ba_0);
    
G_pb_0 = 1/(1 + exp(-q_pb*(n_0 - xi_pb)));
G_pr_0 = 1/(1 + exp(q_pr*(data.Pthbar - xi_pr)));
G_s_0  = 1/(1 + exp(q_s*(n_0 - xi_s))); 
    
T_pb_0 = K_pb*G_pb_0; 
T_pr_0 = K_pr*G_pr_0;
T_s_0  = K_s*G_s_0; 
    
% Output vector with initial conditions
Init = [Pc_m_0; Pa_m_0; e_bc_0; e_ba_0; T_pb_0; T_pr_0; T_s_0; HI];

end % function initialconditions.m % 