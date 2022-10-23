function [Ui,Ue] =FTSTS_pulse2(duration,step,pulse_width,pulse_interval,pulse_amplitude,pulse_phase_dif,aE,aI)

U = zeros(1,floor(duration/step));
Tpulse = (pulse_width*2+pulse_interval);
Npulse = floor(Tpulse/step);
deltaT = (2*pulse_width)*pulse_phase_dif/(2*3.14);  % completely out-of-phase reference: phi_I - phi_E = pi
I = [1:1:floor(duration/step)-Npulse]-1;

for i = 1:Npulse
   
   U(1,I+i) =  -pulse_amplitude*(mod(I+i,Npulse)>=1).*(mod(I+i,Npulse)<1+floor(pulse_width/step)) ...
          + pulse_amplitude*(mod(I+i,Npulse)>=1+floor(pulse_width/step)).*(mod(I+i,Npulse)<1+2*floor(pulse_width/step));
    
end


% shift phase
Nphase = floor(abs(deltaT)/step);
Ui = zeros(size(U));
if deltaT>0
    Ui(1,Nphase+1:end) = aI*U(1,1:end-Nphase);
    Ui(1,1:Nphase) = aI*U(1,end-Nphase+1:end);
elseif deltaT<0
    Ui(1,1:end-Nphase) = aI*U(1,Nphase+1:end);
    Ui(1,end-Nphase+1:end) = aI*U(1,1:Nphase);
else
    Ui = aI*U;
end

% E FTSTS
Ue = aE*U;


end