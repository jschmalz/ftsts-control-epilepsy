function [R_t,R_spE,R_spI,R_V,R_Vi,R_s_E,R_s_I,R_Cl_in,R_W,W,R_W2,Wee] = model_seizure(duration,step,N_E,x0,W_EE, W_EI, W_IE, W_II,PN,IC,Iapp,W0,FTSTS)

% record
R_t = zeros(1,floor(duration/step));
R_spE = zeros(N_E,floor(duration/step));
R_spI = zeros(N_E,floor(duration/step));
R_V = zeros(N_E,floor(duration/step));
R_Vi = zeros(N_E,floor(duration/step));
R_s_E = zeros(N_E,floor(duration/step));
R_s_I = zeros(N_E,floor(duration/step));
R_Cl_in = zeros(N_E,floor(duration/step));

IC_t = 0;

len_segments = 20; 
segment_duration = floor(duration/len_segments); 

R_W = zeros(1,len_segments);
R_W2 = zeros(1,len_segments);

% plastic weight

Wp = ones(N_E,N_E,2);
Apre = zeros(N_E,N_E,2);
Apost = zeros(N_E,N_E,2);
% 
Wp2 = ones(N_E,N_E,2);
Apre2 = zeros(N_E,N_E,2);
Apost2 = zeros(N_E,N_E,2);

% save(,'W')
% save('weightee_500_100_long300.mat','Wee')
% load('weight_500_110_decW.mat')
% load('weightee_500_110_decW.mat')

% load('weight_500_100_long500.mat')
% load('weightee_500_100_long500.mat')

Wp(:,:,1) = W0;
Wp(:,:,2) = W0;
Wp2(:,:,1) = W0;
Wp2(:,:,2) = W0;

%% Pink noise generation
sigma_current = PN(1);
tau_x_current = PN(2);
tau_t_current = PN(3);
dtI = step;
nI = [N_E,1];
Iz = zeros(nI);

NI = floor(duration/step);
Iext = pink_noise(sigma_current,tau_x_current,tau_t_current,dtI,Iz,nI);

% FTSTS parameters
pulse_width = FTSTS(1); 
pulse_interval = FTSTS(2);
pulse_amplitude = FTSTS(3);
pulse_phase_dif = FTSTS(4);% degrees
aE = FTSTS(5);
aI = FTSTS(6);


% % remove neurons
% V_zeros = ones(N_E,1);
% I = 


%calc pulse
[U,Ue] = FTSTS_pulse2(duration,step,pulse_width,pulse_interval,pulse_amplitude,pulse_phase_dif,aE,aI);

figure(21)
plot([step:step:duration],Ue,'k',[step:step:duration],U,'r')
legend('E','I')
xlim([100,150])

% W_EE_nodiag = W_EE - diag(diag(W_EE));

% disp('Percent Complete (%)')
% run
for i = 1:segment_duration
    
    N_stim = len_segments/step; % number of steps in segment duration
    a = (i-1)*N_stim +1;
    b = i*N_stim;
    
    [t,V,s_E,s_I,sp_E,Vi,si_E,si_I,sp_I,Cl_in,IC,IC_t,Wp,Apre,Apost,Wp2,Apre2,Apost2,Iext] = ode_seizure(len_segments,step,...
        N_E,x0,...
        W_EE, W_EI, W_IE, W_II,...g
        PN,IC_t,IC,Iapp,...
        Wp,Apre,Apost,...
        Wp2,Apre2,Apost2,...
        Iext,Ue(1,a:b),U(1,a:b));
    
    Iext = pink_noise(sigma_current,tau_x_current,tau_t_current,dtI,Iext,nI);
    
    
    
    R_t(1,a:b) = t(1,1:end-1);
    R_V(:,a:b) = V(:,1:end-1);
    R_Vi(:,a:b) = Vi(:,1:end-1);
    R_spE(:,a:b) = sp_E(:,1:end-1);
    R_spI(:,a:b) = sp_I(:,1:end-1);
    R_s_E(:,a:b) = s_E(:,1:end-1);
    R_s_I(:,a:b) = s_I(:,1:end-1);
    R_Cl_in(:,a:b) = Cl_in(:,1:end-1);
    R_W(:,i) = mean(sum(W_EI.*squeeze(Wp(:,:,end)),1));
    R_W2(:,i) = mean(sum(W_EE.*squeeze(Wp2(:,:,end)),1));
    
    if mod(floor(i*len_segments),floor(0.25*duration)) == 0
        disp(100*i/segment_duration)
    end
end



W = squeeze(Wp(:,:,end));
Wee = squeeze(Wp2(:,:,end));

end