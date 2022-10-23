function [t,V,s_E,s_I,sp_E,Vi,si_E,si_I,sp_I,Cl_in,FVs,FVt,Wp,Apre,Apost,Wp2,Apre2,Apost2,Iext] = ode_seizure(duration,step,N_E,x0,W_EE, W_EI, W_IE, W_II,PN,IC_t,IC,Iapp,Wp,Apre,Apost,Wp2,Apre2,Apost2,Iext,Ue,U)
% parameters
C = x0(1); 
g_L = x0(2); 
E_L = x0(3); 
E_E = x0(4);  
E_K = x0(5);  
f0 = x0(6);  
beta = x0(7);   
tau_ref = x0(8);  
tau_E = x0(9);  
tau_I = x0(10);  
tau_phi = x0(11);  
phi_0 = x0(12);
delta_phi = x0(13);  
tau_Cl = x0(14);  
V_d = x0(15);   
Cl_ineq = x0(16); 
Cl_out = x0(17);  
tau_K = x0(18);  
delta_K = x0(19);  
sigma_E = x0(20);
sigma_I = x0(21); 
tau_STDP = x0(22);
F = x0(23); 
fmax = x0(24); 
eta = x0(25) ;
a_LTP = x0(26);
a_LTD = x0(27);

% initalize vectors
V = zeros(N_E,floor(duration/step));
Cl_in = zeros(N_E,floor(duration/step));
g_K = zeros(N_E,floor(duration/step));
phi = zeros(N_E,floor(duration/step));
I_Cl = zeros(N_E,floor(duration/step));
Tref = zeros(N_E,floor(duration/step));
s_E = zeros(N_E,floor(duration/step));
s_I = zeros(N_E,floor(duration/step));
t = zeros(1,floor(duration/step));
sp_E = zeros(N_E,floor(duration/step));


Vi = zeros(N_E,floor(duration/step));
Cl_ini = zeros(N_E,floor(duration/step));
gi_K = zeros(N_E,floor(duration/step));
phii = zeros(N_E,floor(duration/step));
I_Cli = zeros(N_E,floor(duration/step));
Trefi = zeros(N_E,floor(duration/step));
si_E = zeros(N_E,floor(duration/step));
si_I = zeros(N_E,floor(duration/step));
sp_I = zeros(N_E,floor(duration/step));

Wsave = zeros(N_E,N_E,floor(duration/step));

t(1,1) = IC_t;
V(:,1) = IC(:,1);
Cl_in(:,1) = IC(:,2);
g_K(:,1) = IC(:,3);
phi(:,1) = IC(:,4);
I_Cl(:,1) = IC(:,5);
Tref(:,1) = IC(:,6);
s_E(:,1) = IC(:,7);
s_I(:,1) = IC(:,8);
sp_E(:,1) = IC(:,9);
Vi(:,1) = IC(:,10);
Cl_ini(:,1) = IC(:,11);
gi_K(:,1) = IC(:,12);
phii(:,1) = IC(:,13);
I_Cli(:,1) = IC(:,14);
Trefi(:,1) = IC(:,15);
si_E(:,1) = IC(:,16);
si_I(:,1) = IC(:,17);
sp_I(:,1) = IC(:,18);

%% Pink noise generation
sigma_current = PN(1);
tau_x_current = PN(2);
tau_t_current = PN(3);
dtI = step;
nI = [N_E,1];
Iz = zeros(nI);

% W_EE_nodiag = W_EE - diag(diag(W_EE));
%%
for i = 1:floor(duration/step)
    
    
    Iext = pink_noise(sigma_current,tau_x_current,tau_t_current,dtI,Iext,nI);
    
    Iapp_zeros = ones(N_E,1);
%     Na = floor(N_E/2 - 0.025*N_E);
%     Nb = floor(N_E/2 + 0.025*N_E);
%     Iapp_zeros(floor(350):floor(150),1) = 1;
    
    % E inputs
    I_s = Iext;
    I_FTSTS = Ue(1,i)*Iapp_zeros*(t(i)>10000)*(t(i)<=15000);
    Na = floor(N_E/2 - 0.025*N_E);
    Nb = floor(N_E/2 + 0.025*N_E);
    I_d  = zeros(N_E,1);
    I_d(Na:Nb,1) = Iapp*(t(i)>1100)*(t(i)<=4100) ;% + I_FTSTS(floor(Na-50):floor(Nb+50),1);
    I = I_s + I_d + I_FTSTS;
    g_E = s_E(:,i);
    g_I = s_I(:,i);
    
    % I inputs
    Ii_s = Iext;
    Ii_FTSTS = U(1,i)*Iapp_zeros*(t(i)>10000)*(t(i)<=15000);
    Ii_d  = zeros(N_E,1);
    Ii_d(Na:Nb,1) = Iapp*(t(i)>1100)*(t(i)<=4100) ;%  + Ii_FTSTS(floor(Na-50):floor(Nb+50),1);
    Ii = Ii_s + Ii_d + Ii_FTSTS;
    gi_E = si_E(:,i);
    gi_I = si_I(:,i);
    
    % E pre-calculations
    f(:,1) = (f0.*exp(((V(:,i))-phi(:,i))./beta) ).*(Tref(:,i)<0);
    E_Cl = -26.7*log(Cl_out./Cl_in(:,i));
    I_Cl(:,i) = g_I.*(V(:,i)-E_Cl);
    
    % I pre-calculations
    fi(:,1) = (f0.*exp(((Vi(:,i))-phii(:,i))./beta) ).*(Trefi(:,i)<0);
    Ei_Cl = -26.7*log(Cl_out./Cl_ini(:,i));
    I_Cli(:,i) = gi_I.*(Vi(:,i)-Ei_Cl);
    
    % E  spikes?
    spike = (f(:,1).*step) > rand(N_E,1);
    
    % I  spikes?
    spikei = (fi(:,1).*step) > rand(N_E,1);
    
    % ODEs    
    dVdt = (1/C).*( g_L.*(E_L-V(:,i)) + (g_E./fmax).*(E_E-V(:,i)) + (g_I./fmax).*(E_Cl - V(:,i)) + (g_K(:,i)./fmax).*(E_K-V(:,i)) + I ); %fmax
    dphidt = ( phi_0 - phi(:,i)  )./tau_phi;
    dCl_indt = I_Cl(:,i)./(V_d.*F) + (Cl_ineq - Cl_in(:,i))./tau_Cl; % neg sign here
    dg_Kdt = (1./tau_K).*( -g_K(:,i)   );
    
    dVidt = (1/C).*( g_L.*(E_L-Vi(:,i)) + (gi_E./fmax).*(E_E-Vi(:,i)) + (gi_I./fmax).*(E_Cl - Vi(:,i)) + (gi_K(:,i)./fmax).*(E_K-Vi(:,i)) + Ii );
    dphiidt = ( phi_0 - phii(:,i)  )./tau_phi;
    dCl_inidt = I_Cli(:,i)./(V_d.*F) + (Cl_ineq - Cl_ini(:,i))./tau_Cl; % neg sign here
    dgi_Kdt = (1./tau_K).*( -gi_K(:,i)   );
    
% --- update synaptic variables ---
    % dxdt update
    dApostdt = -Apost(:,:,end)/tau_STDP;
    dApredt = -Apre(:,:,end)/tau_STDP;
     
    % x update (note only save last 60 seconds)
    % shift everything
    Apost(:,:,1:end-1) = Apost(:,:,2:end);
    Apre(:,:,1:end-1) = Apre(:,:,2:end);
    Wp(:,:,1:end-1) = max(min(Wp(:,:,2:end),2),0);
    Wp(:,:,end-1) =  Wp(:,:,end-1) - diag(diag( Wp(:,:,end-1)))  + diag(ones(N_E,1));
    
%     max(min(Wee,2),0)
    
    % plasticity update - "on_pre"
%     Apre(:,:,end) = Apre(:,:,end-1) + 0.005*spike.*(W_EI~=0) + step*dApredt;
%     Wp(:,:,end) = Wp(:,:,end-1) + eta*a_LTD*Apost(:,:,end-1).*spike.*(W_EI~=0);
    Apre(:,:,end) = Apre(:,:,end-1) + spike + step*dApredt;
    Wp(:,:,end) = Wp(:,:,end-1) + eta*a_LTD*Apost(:,:,end-1).*spike;
     
     
    % plasticity update - "on_post"
%     Apost(:,:,end) = Apost(:,:,end-1) + 0.005*spikei'.*(W_EI~=0)  + step*dApostdt;
%     Wp(:,:,end) = Wp(:,:,end) + eta*a_LTP*Apre(:,:,end-1).*spikei'.*(W_EI~=0);
    Apost(:,:,end) = Apost(:,:,end-1) + spikei'  + step*dApostdt;
    Wp(:,:,end) = Wp(:,:,end) + eta*a_LTP*Apre(:,:,end-1).*spikei';
%     Wprec(i+1) = mean(mean(Wp(:,:,end),1));
    
% --- update synaptic variables ---
    % dxdt update
    dApost2dt = -Apost2(:,:,end)/tau_STDP;
    dApre2dt = -Apre2(:,:,end)/tau_STDP;
     
    % x update (note only save last 60 seconds)
    % shift everything
    Apost2(:,:,1:end-1) = Apost2(:,:,2:end);
    Apre2(:,:,1:end-1) = Apre2(:,:,2:end);
    Wp2(:,:,1:end-1) = max(min(Wp2(:,:,2:end),2),0);
    Wp2(:,:,end-1) =  Wp2(:,:,end-1) - diag(diag( Wp2(:,:,end-1))) + diag(ones(N_E,1));
    
    % plasticity update - "on_pre"
%     Apre2(:,:,end) = Apre2(:,:,end-1) + 0.005*spike.*(W_EE~=0) + step*dApre2dt;
%     Wp2(:,:,end) = Wp2(:,:,end-1) + eta*a_LTD*Apost2(:,:,end-1).*spike.*(W_EE~=0);
    Apre2(:,:,end) = Apre2(:,:,end-1) + spike + step*dApre2dt;
    Wp2(:,:,end) = Wp2(:,:,end-1) + eta*a_LTD*Apost2(:,:,end-1).*spike;
     
   %%%%% opposite plasiticyt rule %%%%
     
    % plasticity update - "on_post"
%     Apost2(:,:,end) = Apost2(:,:,end-1) + 0.005*spikei'.*(W_EE~=0)  + step*dApost2dt;
%     Wp2(:,:,end) = Wp2(:,:,end) + eta*a_LTP*Apre2(:,:,end-1).*spikei'.*(W_EE~=0);
    Apost2(:,:,end) = Apost2(:,:,end-1) + spike'  + step*dApost2dt;
    Wp2(:,:,end) = Wp2(:,:,end) + eta*a_LTP*Apre2(:,:,end-1).*spike';
%     Wprec2(i+1) = mean(mean(Wp2(:,:,end),1));
% %     
%  
    
    % update
    V(:,i+1) = V(:,i) + step.*dVdt;
    phi(:,i+1) = phi(:,i) + ( spike ).*delta_phi  + step.*dphidt;
    Cl_in(:,i+1) = Cl_in(:,i) + step.*dCl_indt;
    g_K(:,i+1) = g_K(:,i) +( spike ).*delta_K./tau_K  + step.*dg_Kdt;
    Tref(:,i+1) = Tref(:,i) - step;
%     TAP(:,i+1) = TAP(:,i) - step;
    
    Vi(:,i+1) = Vi(:,i) + step.*dVidt;
    phii(:,i+1) = phii(:,i) + (spikei).*delta_phi  + step.*dphiidt;
    Cl_ini(:,i+1) = Cl_ini(:,i) + step.*dCl_inidt;
    gi_K(:,i+1) = gi_K(:,i) +( spikei ).*delta_K./tau_K  + step.*dgi_Kdt;
    Trefi(:,i+1) = Trefi(:,i) - step;
   
    % synaptic input
    ds_Edt = (-s_E(:,i)  )./tau_E ; 
    ds_Idt = (-s_I(:,i)   )./tau_I ;

    dsi_Edt = (-si_E(:,i)  )./tau_E ; 
    dsi_Idt = (-si_I(:,i)   )./tau_I ;
    
    
    s_E(:,i+1) = s_E(:,i) +  (spike'*(W_EE.*squeeze(Wp2(:,:,end-1))))'./tau_E  + step.*ds_Edt; %
    s_I(:,i+1) = s_I(:,i) +  (spikei'*(W_IE))'./tau_E  + step.*ds_Idt; % 
    
    si_E(:,i+1) = si_E(:,i) +  (spike'*(W_EI.*squeeze(Wp(:,:,end-1))))'./tau_E  + step.*dsi_Edt; %
%     si_E(:,i+1) = si_E(:,i) +  (spike'*(W_EI))'./tau_E  + step.*dsi_Edt; %
    si_I(:,i+1) = si_I(:,i) +  (spikei'*(W_II))'./tau_E  + step.*dsi_Idt; % 
    
    % E spikes
    
    
   
    
    % look for spikes
    sp_E(:,i+1) = [1:1:N_E]'.*spike;%( V(:,i) >= phi(:,i) ).*(Tref(:,i)<=0); 
    sp_I(:,i+1) = [1:1:N_E]'.*spikei;
    
    
    Tref(:,i+1) = Tref(:,i+1).*( spike ~= 1 ) + tau_ref.*spike;%( V(:,i) >= phi(:,i) ).*(Tref(:,i)<=0);
%     TAP(:,i+1) = TAP(:,i+1).*( spike ~= 1 ) + tau_AP.*spike;
    
    Trefi(:,i+1) = Trefi(:,i+1).*( spikei ~= 1 ) + tau_ref.*spikei;
    
 
    V(:,i) = V(:,i).*( spike ~= 1 ) + 0.5*(40+V(:,max(i-1,1))).*( spike );
    V(:,i+1) = V(:,i+1).*( spike ~=1 ) + (V(:,max(i-1,1))-20).*( spike );
    
    Vi(:,i) = Vi(:,i).*( spikei ~= 1 ) + 0.5*(40+Vi(:,max(i-1,1))).*( spikei );
    Vi(:,i+1) = Vi(:,i+1).*( spikei ~=1 ) + (Vi(:,max(i-1,1))-20).*( spikei );
    
    
    t(1,i+1) = t(1,i) + step;
    
%     % save weight to look at
%     Wsave(:,:,i) = Wp(:,:,end-1);
end

% calculate finals values (FVs)
FVs = zeros(N_E,18);
FVt = t(end);

FVs(:,1) = V(:,end);
FVs(:,2) = Cl_in(:,end);
FVs(:,3) = g_K(:,end);
FVs(:,4) = phi(:,end);
FVs(:,5) = I_Cl(:,end);
FVs(:,6) = Tref(:,end);
FVs(:,7) = s_E(:,end);
FVs(:,8) = s_I(:,end);
FVs(:,9) = sp_E(:,end);
FVs(:,10) = Vi(:,end);
FVs(:,11) = Cl_ini(:,end);
FVs(:,12) = gi_K(:,end);
FVs(:,13) = phii(:,end);
FVs(:,14) = I_Cli(:,end);
FVs(:,15) = Trefi(:,end);
FVs(:,16) = si_E(:,end);
FVs(:,17) = si_I(:,end);
FVs(:,18) = sp_I(:,end);




end