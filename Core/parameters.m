function [pars,data] = parameters(data)
% function parameters.m
% Input: data structure
% Outout: list of model parameters and augmented data structure
% Description: sets nominal model parameters

% Initial mean values
SPbar  = data.SPbar;
HminR  = data.HminR; 
HmaxR  = data.HmaxR; 
Hbar   = data.Hbar; 

% Model parameters
D  = 0.75; 
B  = 0.4; 
A  = 5;         

% Gains
K_P  = 1; 
K_b  = .1;
K_pb = 1;
K_pr = 1; 
K_s  = 1; 

% Time scales
tau_P  = 1;  
tau_b  = 1;              
tau_pb = 0.5;          
tau_pr = 4;
tau_s  = 10;          
tau_H  = .1;

% Sigmoid shifts
q_w  = .04;         
q_pb = 10;          
q_pr = 1; 
q_s  = 10;   

% Patient specific parameters
xi_w  = D * SPbar + (1 - D) * data.PPbar;           
xi_pr = data.Pthbar;  

% Intrinsic HR
HI     = Hbar;
HIa    = data.Hbara; 
iH_Ia  = data.i_HR4e_min; 
iH_Iaa = data.i_HR4e_Max; 
iH_I   = data.i_HR4e_min;  
tchar  = data.Tdata(floor(iH_I));

% Maximal HR
HM  = 1.5*max(data.HRdata);    
H_s = (1/K_s)*(HM/HI - 1); 

% At end of expiration and inspiration
Gpr_ss = 1/(1 + exp(q_pr*(data.Pthbar - xi_pr)));

Tpr_ss = K_pr*Gpr_ss; 
H_pr = (HmaxR - HminR)/HI/Tpr_ss ;

% Calculate sigmoid shifts
Pc_ss  = xi_w; 
Pa_ss  = xi_w - data.Pthbar; 

ewc_ss = 1 - sqrt((1 + exp(-q_w*(Pc_ss - xi_w)))/(A + exp(-q_w*(Pc_ss - xi_w)))); 
ewa_ss = 1 - sqrt((1 + exp(-q_w*(Pa_ss - xi_w)))/(A + exp(-q_w*(Pa_ss - xi_w)))); 

ebc_ss = K_b*ewc_ss; 
eba_ss = K_b*ewa_ss;  

n_ss   = B*(ewc_ss - ebc_ss) + (1 - B)*(ewa_ss - eba_ss);

Tpb_ss = .8;
Ts_ss  = .2; 

% Steady-state sigmoid shifts 
xi_pb = n_ss + log(K_pb/Tpb_ss - 1)/q_pb;  
xi_s  = n_ss - log(K_s/Ts_ss - 1)/q_s;   

H_pb = (H_pr*Tpr_ss + H_s*Ts_ss)/Tpb_ss;

% Parameter vector (output)

pars = [D; B;                                   % Convex combination parameters 1-2
    A;                                          % Wall strain parameter 3 
    K_P; K_b; K_pb; K_pr; K_s;                  % Gains 4-8
    tau_P; tau_b; tau_pb; tau_pr; tau_s; tau_H; % Time Constants 9-14
    q_w; q_pb; q_pr; q_s;                       % Sigmoid Steepnesses 15-18
    xi_w; xi_pb; xi_pr; xi_s;                   % Sigmoid Shifts 19-22
    HI; H_pb; H_pr; H_s;                        % Heart Rate Parameters 23-26
    HIa; tchar;];                               % Mean HR after VM and index for time at which it changes 27
               
pars_names = {'$\alpha$', '$\beta$', '$A$', ...
    '$K_P$','$K_b$','$K_{pb}$','$K_{pr}$','$K_s$', ...
    '$\tau_P$','$\tau_b$','$\tau_{pb}$','$\tau_{pr}$','$\tau_s$','$\tau_H$',...
    '$q_w$','$q_{pb}$','$q_{pr}$','$q_{s}$', ...
    '$\xi_w$','$\xi_{pb}$','$\xi_{pr}$','$\xi_{s}$', ...
    '$H_I$','$H_{pb}$','$H_{pr}$','$H_{s}$', ...
    '$H_{Ia}$','$tchar$'};

% Parameter bounds

% Vary nominal parameters by +/- 50%
lb  = pars/10; 
ub  = pars*10;

% B - Convex combination
lb([1,2])  = .01;                       
ub([1,2])  = 1;
lb([11])   = 0.0001;

% Log scaled outputs
pars = log(pars);
lb   = log(lb);
ub   = log(ub); 
data.lb = lb;
data.ub = ub;
data.pars_names = pars_names; 

end  % function parameters.m %