clearvars
clear -all
clc

% Schmalz, Joseph



%% ---------------------------------
% ---------- Run Parameters --------
duration =  150*1000; %ms
step = 0.5; %ms
N_E = 500; 

% applied external current to initiate the seizure
Iapp = 200;

% ---------------------------------
% ---------------------------------

%%
% ---------------------------------------
% ------ FTSTS pulse parameters --------
pulse_width = 1; %ms
pulse_interval = 10; %ms
pulse_amplitude = 2*1000; %1500 2000 3000 %ms
pulse_phase_dif = 0*3.14;% radians
aE = -1; % polarity
aI = 1;  % polarity
% ---------------------------------------
% ---------------------------------------


%%

% network parameters 
C = 100; % 
g_L = 4; % 
E_L = -57; %   
E_E = 0;  % 
E_K = -90;  % 
f0 = 0.002;  
beta = 1.5;   % 
tau_ref = 5;  % 
tau_E = 15;  % 
tau_I = 15;  % 
tau_phi = 100;  % 
phi_0 = -55;  % 
delta_phi = 2.5;  % 
tau_Cl = 5*1000;   % 
V_d = 0.2357;   %   
Cl_ineq = 6; % 
Cl_out = 110;  % 
tau_K = 5*1000;  % 
delta_K = 40;  % 
sigma_E = 0.02; % 
sigma_I = 0.03;  % 
tau_STDP = 15; % 
F = 96500; %  
fmax = 1/(5);  %  

eta = 1*10^-3;% 
a_LTP = 1;
a_LTD = -1*a_LTP;




x0(1) = C;
x0(2) = g_L;
x0(3) = E_L ;
x0(4) = E_E;
x0(5) = E_K ;
x0(6) = f0 ;
x0(7) = beta ;
x0(8) = tau_ref ;
x0(9) = tau_E;
x0(10) = tau_I;
x0(11) = tau_phi;
x0(12) = phi_0 ;
x0(13) = delta_phi ;
x0(14) = tau_Cl ;
x0(15) = V_d ;
x0(16) = Cl_ineq ;
x0(17) = Cl_out ;
x0(18) = tau_K ;
x0(19) = delta_K ;
x0(20) = sigma_E ;
x0(21) = sigma_I ;
x0(22) = tau_STDP ;
x0(23) = F ;
x0(24) = fmax ;
x0(25) = eta ;
x0(26) = a_LTP ;
x0(27) = a_LTD ;



% noise parameters
sigma_current = 20;    % unit: pA
tau_x_current = 200; % spatial constant unit: neuron index
tau_t_current = 15; % time unit: ms

PN(1) = sigma_current;
PN(2) = tau_x_current;
PN(3) = tau_t_current;

% ICs
V0 = -67;
IC =  zeros(N_E,18);
IC(:,1) = V0;
IC(:,2) = Cl_ineq;
IC(:,3) = 0;
IC(:,4) = -55;
IC(:,5) = 0;
IC(:,6) = 0;
IC(:,7) = 0;
IC(:,8) = 0;
IC(:,9) = 0;
IC(:,10) = -67;
IC(:,11) = Cl_ineq;
IC(:,12) = 0;
IC(:,13) = -55;
IC(:,14)  = 0;
IC(:,15)  = 0;
IC(:,16)  = 0;
IC(:,17)  = 0;
IC(:,18)  = 0;



% weight matrix
[W_EE, W_EI, W_IE, W_II] = make_weights(sigma_E,sigma_I,N_E);
W0 = ones(N_E,N_E);

% FTSTS parameter to pass to model
FTSTS(1) = pulse_width;
FTSTS(2) = pulse_interval;
FTSTS(3) = pulse_amplitude;
FTSTS(4) = pulse_phase_dif;
FTSTS(5) = aE;
FTSTS(6) = aI;



%%


% run
tic
[R_t,R_spE,R_spI,R_V,R_Vi,R_s_E,R_s_I,R_Cl_in,R_W,W,R_W2,Wee] = model_seizure(duration,step,N_E,x0,W_EE, W_EI, W_IE, W_II,PN,IC,Iapp,W0,FTSTS);
toc


% save('./Standard-FTSTS-Results/R_1500_spE.mat','R_spE')
% save('./Standard-FTSTS-Results/R_1500_spI.mat','R_spI')
% save('./Standard-FTSTS-Results/W_EE_t_1500.mat','R_W')
% save('./Standard-FTSTS-Results/W_EI_t_1500.mat','R_W2')
% save('./Standard-FTSTS-Results/W_EE_1500.mat','Wee')
% save('./Standard-FTSTS-Results/W_EI_1500.mat','W')

%% plots


%%
% ---------------------------
% raster plot

figure(1)
plot(R_t/1000,R_spE','k.','MarkerSize',2)
ylim([0.9,N_E+0.1])
ylabel('Neuron Index of E Population')
xlabel('Time (ms)')
xlim([0,150])
%



%%
% ---------------------------
% plot the average weight over time

ts = [20:20:duration]./1000;

figure(2)
subplot(2,1,1)
plot(ts,R_W)
title('E-to-I')
xlim([0,duration/1000])
xlabel('Time (sec)')
ylabel('Neuron Index')

subplot(2,1,2)
plot(ts,R_W2)
title('E-to-E')
xlim([0,duration/1000])
xlabel('Time (sec)')
ylabel('Neuron Index')


% ---------------------------
% heat map of the final weights 


figure(4)
colormap('jet');   % set colormap
imagesc(W_EE.*Wee-W_EE);       
colorbar;
title('Percent Change in Weight (\DeltaW_{ee})')
ylabel('Excitatory Neuron Index')
xlabel('Excitatory Neuron Index')


figure(5)
colormap('jet');   % set colormap
imagesc(W_EI.*W-W_EI);        
colorbar;
title('Change in Conductance')
xlabel('Inhibitory Neuron Index')
ylabel('Excitatory Neuron Index')



