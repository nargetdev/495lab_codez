% ME495 Lab 1
% Last Updated: Winter 2011
% Example program for Frequency response of flexible shaft
% Transfer function and frequency response using system parameters
K_m 			= Km_calculated;		% N-m/A or V/(rad/sec)
R_m 			= Rm;		% Ohms
L_m 			= 0.002;	% Henry
J_m 			= 3.8e-5;	% Kg-m^2
D_flywheel		= 0.137;	% m
L_flywheel		= 0.0127;	% m
rho_flywheel	= 7755;		% kg/m^3
D_shaft			= 3.18e-3;	% m
L_shaft			= 0.305;	% m
G_shaft			= 7.31e10;	% N/m^2

% Computed parameters
J_f = rho_flywheel*pi*D_flywheel^4*L_flywheel/32;
K  = pi*G_shaft*D_shaft^4/(32*L_shaft);
J_c = J_m + J_f;

% Key in the damping terms you obtained in the following two lines
% B_s 		= 0.000264622410542826;		% Type in your value here
% B_b         = 8.47E-4;		% Type in your value here
% B_m         = 0.0132;        % Type in your value here
% B_tr_given     	= B_s + B_b;
% B_tl_given        = B_s + B_m;
B_tr_given=B_tr;
B_tl_given=B_tl;

num = [K_m*B_s,  K_m*K];
den = [R_m*J_c*J_f,   R_m*J_c*B_tr+J_f*(K_m^2+B_tl_given*R_m), ...
       	R_m*(K*(J_c+J_f)+B_b*B_s)+(K_m^2+R_m*B_m)*B_tr_given, K*(K_m^2+(B_m+B_b)*R_m)];

% import the data file obtained from Labview
freq_test=2*pi*freqsweep_freq(:,1);  
   % data inside [..] should be the test frequency in Hertz

mag_test = freqsweep_bodemag(:,1); 
%mag_test = mag_test*0.5083;
%mag_test = mag_test/3.6;
   % data inside [..] should be output Bode magnitude in dB
		   
phase_test = freqsweep_bodephase(:,1);
   % data inside [..] should be output phase lag in degree

[mag_model, phase_model] = bode(num, den, freq_test);
figure(1)
plot(freq_test/(2*pi), (mag_test-20), 'ro', freq_test/(2*pi), 20*log10(mag_model),'b-')
xlabel('w (Hz)')
ylabel('| theta_2 dot (rad/sec) / V_i_n (volts) | in dB')
title('Line: model    circles: test results')
figure(2)
plot(freq_test/(2*pi),phase_test, 'ro', freq_test/(2*pi), phase_model,'b-')
xlabel('w (Hz)')
ylabel('Phase of (theta_2 dot - V_i_n) in deg')
title('Line: model    circles: test results')
