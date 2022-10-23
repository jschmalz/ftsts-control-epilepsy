clearvars
clc
clear -all

duration = 1000;
step = 0.5;

% FTSTS pulse parameters
pulse_width = 1; 
pulse_interval = 10;
pulse_amplitude = 1000;
pulse_phase_dif = -12*3.14;% degrees


aE = 1; 
aI = -1; 


FTSTS(1) = pulse_width;
FTSTS(2) = pulse_interval;
FTSTS(3) = pulse_amplitude;
FTSTS(4) = pulse_phase_dif;
FTSTS(5) = aE;
FTSTS(6) = aI;


%calc pulse
[Ui,Ue] = FTSTS_pulse2(duration,step,pulse_width,pulse_interval,pulse_amplitude,pulse_phase_dif,aE,aI);



figure(21)
plot([step:step:duration],Ue,'k',[step:step:duration],Ui,'r')
legend('E','I')
xlim([100,150])



