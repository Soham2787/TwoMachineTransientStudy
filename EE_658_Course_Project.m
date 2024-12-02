% NAME: Soham Deepak Koyande
% Roll No.: 21D170040
% Date: 29th Nov 2024
% Course Code: EE 658
% Course Name: Power System Dynamics and Control
% Course Instructor: Prof. Anil Kulkarni
% COURSE PROJECT: Simulating a system of 2 identical synchronous machine
%                 connected by a transmission line

%% Parameters to be changed

clear
clc
close all

% initial system parameters
R_L_1_0 = 10;               % initial load resistance on bus 1
R_L_1   = R_L_1_0;

R_L_2_0 = 0.5263;           % initial load resistance on bus 2
R_L_2   = R_L_2_0;

r_line_0 = 0.1;             % initial line resitance      
r_line = r_line_0;

x_line_0 = 0.3;             % initial line impedance
x_line = x_line_0;

z_complex_line_0 = r_line + i*x_line;


% simulation step size and time duration
h = 0.0001;                 % numerical integration time step
Ts = 1;                     % checking steady state
Tend1 = 10;                 % load change/ fault at this time
Tend2 = 30;                 % final time


% updating the system parameters to perform tests (you can update these)
R_L_1_new = 10;             % new resistance at bus 1 updated at Tend1
R_L_2_new = 0.5263;         % new resistance at bus 2 updated at Tend1
r_line_new = 0.1;           % new line resistance
x_line_new = 0.3;           % new line impedance
governer_1_enabled = 1;     % governer of machine 1 enabled(1)/ disabled(0)
governer_2_enabled = 1;     % governer of machine 2 enabled(1)/ disabled(0)
Tfault_remove_period = 0.08;% time duration of fault in system


% after fault parameters of system (same as initial, by default)
R_L_1_new2 = R_L_1_0;       % load resistance at bus 1 post fault clearance
R_L_2_new2 = R_L_2_0;       % load resistance at bus 2 post fault clearance
r_line_new2 = r_line_0;     % line resistance post fault clearance
x_line_new2 = x_line_0;     % line impedance post fault clearance


fault = 0;
if R_L_1_new == 0 | R_L_2_new == 0 | sqrt((r_line_new^2) + (x_line_new^2)) == 0 | R_L_1_new > 100 | R_L_2_new > 100 | r_line ~= 0.1 | x_line ~= 0.3
    fault = 1;
end

Tfault_remove = Tend1 + Tfault_remove_period;

%% Machine Parameters (identical for both machines):

Ra      = 0.0;              % Armature resistance (in p.u.)
H       = 3;                % Inertia Constant (in MJ/MVA)

xd      = 2;                % d-axis synchronous reactance (in p.u.)
xd_d    = 0.32;             % d-axis transient reactance (in p.u.)
xd_dd   = 0.2;              % d-axis sub-transient reactance (in p.u.)
Tdo_d   = 5;                % d-axis transient open circuit time constant (in seconds)
Tdo_dd  = 0.05;             % d-axis sub-transient open circuit time constant (in seconds)

xq      = 1.9;              % q-axis synchronous reactance (in p.u.)
xq_d    = 0.75;             % q-axis transient reactance (in p.u.)
xq_dd   = 0.2;              % q-axis sub-transient reactance (in p.u.)
Tqo_d   = 1;                % q-axis transient open circuit time constant (in seconds)
Tqo_dd  = 0.05;             % q-axis sub-transient open circuit time constant (in seconds)


% Assume Tdcdd = Tddd

% Finding short circuit time constants
% To determine Tdd and Tddd from Td0d and Td0dd
a=(1- xd./xd_d + xd./xd_dd);
b=-(Tdo_d + Tdo_dd);
c=(xd_dd./xd_d).*Tdo_d.*Tdo_dd;
Td_dd1= (-b + sqrt(b.*b - 4*a.*c))./(2*a);
Td_dd2= (-b - sqrt(b.*b - 4*a.*c))./(2*a);
Td_dd= min(Td_dd1,Td_dd2);
Td_d = Tdo_d.*Tdo_dd.*(xd_dd./xd)./Td_dd;
% To determine Tqd and Tqdd from Tq0d and Tq0dd
a=(1- xq./xq_d + xq./xq_dd);
b=-(Tqo_d + Tqo_dd);
c=(xq_dd./xq_d).*Tqo_d.*Tqo_dd;
Tq_dd1= (-b + sqrt(b.*b - 4*a.*c))./(2*a);
Tq_dd2= (-b - sqrt(b.*b - 4*a.*c))./(2*a);
Tq_dd= min(Tq_dd1,Tq_dd2);
Tq_d = Tqo_d.*Tqo_dd.*(xq_dd./xq)./Tq_dd;


% Exciter Parameters:
ka      = 200;          % Exciter transfer function gain
Ta      = 0.02;         % Exciter time parameter (in seconds)
Efd_max = 6;            % Max. limit of field voltage (in p.u.)
Efd_min = -6;           % Min. limit of field voltage (in p.u.)

% Governor Parameters:
Tg1     = 2;            % in seconds
Tg2     = 6;            % in seconds
Kg      = 20;           % unitless
Pm_max  = 1.2; %1.105;  % in p.u.
Pm_min  = 0.6;          % in p.u.

W_b = 2*pi*50;          % base speed 
W_0 = W_b;              % system speed
W_ref = W_b;            % reference speed


%% initial condition calculations

P_g_2_0 = 1;            % machine 2 initial active power injection (given)          

V_mag_1_0 = 1;          % bus 1 initial voltage magnitude (given)
V_phase_1_0 = 0;        % bus 1 initial voltage phase (given)
V_1 = 1;

V_mag_2_0 = 1;          % bus 2 initial voltage magnitude (given)


% Solving the load-flow problem to find bus 2 voltage phasor and line
% current

% Admittance matrix
Y_1 = 1/R_L_1_0;
Y_2 = 1/R_L_2_0;
Y_1_2 = 1/z_complex_line_0;

Y_BUS = [Y_1 + Y_1_2 , -Y_1_2;
         -Y_1_2     , Y_2 + Y_1_2];

% guess
V_2 = V_mag_2_0;
tol   = 1e-6;
iter_max = 1000;
converged = false;

% Gauss Siedel iteration
for iter = 1:iter_max
    V_2_old = V_mag_2_0;

    Q_g_2_0 = -imag(conj(V_2)*((Y_BUS(2,1)*V_1) + (Y_BUS(2,2)*V_2)));

    S_2 = P_g_2_0 + i*Q_g_2_0;

    V_2 = ((conj(S_2))/(conj(V_2)*Y_BUS(2,2))) - ((1/Y_BUS(2,2))*(Y_BUS(2,1)*conj(V_1)));

    V_2 = V_mag_2_0*(V_2/abs(V_2));

    if abs(V_2 - V_2_old) < tol
        converged = true;
        break;
    end
end

V_complex_1_0 = V_mag_1_0*(cos(V_phase_1_0) + i*sin(V_phase_1_0));
V_complex_2_0 = V_2;
V_phase_2_0 = angle(V_2);

% initial line current
I_complex_line_0 = (V_complex_1_0 - V_complex_2_0)/(z_complex_line_0);

I_mag_line_0 = abs(I_complex_line_0);
I_phase_line_0 = angle(I_complex_line_0);



%% Machine initial conditions

% current through resistances
I_L_1_0 = V_complex_1_0/R_L_1_0;            % current through bus 1 resistance
I_L_2_0 = V_complex_2_0/R_L_2_0;            % current through bus 2 resistance

% current injection by machine 1
I_complex_1_0 = I_complex_line_0 + I_L_1_0; 
I_mag_1_0 = abs(I_complex_1_0);
I_phase_1_0 = angle(I_complex_1_0);

% current injection by machine 2
I_complex_2_0 = I_L_2_0 - I_complex_line_0;
I_mag_2_0 = abs(I_complex_2_0);
I_phase_2_0 = angle(I_complex_2_0);

% power angle of machines 1 and 2 respectively
delta_1_0 = angle(V_complex_1_0 + i*(xq)*(I_complex_1_0));
delta_2_0 = angle(V_complex_2_0 + i*(xq)*(I_complex_2_0));

% dq0 frame voltage of machines 1 and 2 respectively
V_d_1_0 = V_mag_1_0*(sin(V_phase_1_0 - delta_1_0));
V_d_2_0 = V_mag_2_0*(sin(V_phase_2_0 - delta_2_0));

V_q_1_0 = V_mag_1_0*(cos(V_phase_1_0 - delta_1_0));
V_q_2_0 = V_mag_2_0*(cos(V_phase_2_0 - delta_2_0));

% dq0 frame currents of machines 1 and 2 respectively
I_d_1_0 = I_mag_1_0*(sin(I_phase_1_0 - delta_1_0));
I_d_2_0 = I_mag_2_0*(sin(I_phase_2_0 - delta_2_0));

I_q_1_0 = I_mag_1_0*(cos(I_phase_1_0 - delta_1_0));
I_q_2_0 = I_mag_2_0*(cos(I_phase_2_0 - delta_2_0));

% field voltage of machine 1 and 2 respectively
E_fd_1_0 = abs(V_complex_1_0 + i*(xq)*(I_complex_1_0)) - ((xd - xq)*I_d_1_0);
E_fd_2_0 = abs(V_complex_2_0 + i*(xq)*(I_complex_2_0)) - ((xd - xq)*I_d_2_0);

% dq0 frame (d-axis and q-axis) flux magnitudes
psi_d_1_0 = V_q_1_0;
psi_d_2_0 = V_q_2_0;

psi_q_1_0 = -V_d_1_0;
psi_q_2_0 = -V_d_2_0;

psi_H_1_0 = psi_d_1_0;
psi_H_2_0 = psi_d_2_0;

psi_F_1_0 = psi_d_1_0 + ((xd_d/(xd - xd_d))*E_fd_1_0);
psi_F_2_0 = psi_d_2_0 + ((xd_d/(xd - xd_d))*E_fd_2_0);

psi_K_1_0 = psi_q_1_0;
psi_K_2_0 = psi_q_2_0;

psi_G_1_0 = psi_q_1_0;
psi_G_2_0 = psi_q_2_0;

% machine rotor speeds 
omega_1_0 = W_b;
omega_2_0 = W_b;

% machine mechanical power inputs
Pm_1_0 = (psi_d_1_0*I_q_1_0) - (psi_q_1_0*I_d_1_0);
Pm_2_0 = (psi_d_2_0*I_q_2_0) - (psi_q_2_0*I_d_2_0);

% Exciter states
Xe_1_0 = E_fd_1_0;
Xe_2_0 = E_fd_2_0;

% reference voltages
V_ref_1_0 = V_mag_1_0 + (E_fd_1_0/ka);
V_ref_2_0 = V_mag_2_0 + (E_fd_2_0/ka);

% speed governor initial state values
Xg_1_0 = 0;
Xg_2_0 = 0;

% defining vectors for each machine's initial condition (to be updated
% during the simulation)

InitialCondition_1 = [V_d_1_0, V_q_1_0, I_d_1_0, I_q_1_0, psi_d_1_0, psi_q_1_0, E_fd_1_0, psi_H_1_0, psi_F_1_0, psi_K_1_0, psi_G_1_0, omega_1_0, Pm_1_0, Xe_1_0, V_ref_1_0, delta_1_0, Xg_1_0]';
InitialCondition_2 = [V_d_2_0, V_q_2_0, I_d_2_0, I_q_2_0, psi_d_2_0, psi_q_2_0, E_fd_2_0, psi_H_2_0, psi_F_2_0, psi_K_2_0, psi_G_2_0, omega_2_0, Pm_2_0, Xe_2_0, V_ref_2_0, delta_2_0, Xg_2_0]';
InitialCondition   = [InitialCondition_1, InitialCondition_2];

% Displaying the initial conditions
% Define row names as the variables
row_names = { ...
    'V_d', 'V_q', 'I_d', 'I_q', ...
    'psi_d', 'psi_q', 'E_fd', 'psi_H', ...
    'psi_F', 'psi_K', 'psi_G', 'omega', ...
    'Pm', 'Xe', 'V_ref', 'delta', 'Xg' ...
};

% Define column names
col_names = {'Index', 'InitialCondition_1', 'InitialCondition_2'};

% Create the table
InitialCondition = table( ...
    [1:17]', InitialCondition_1, InitialCondition_2, ...
    'RowNames', row_names, ...
    'VariableNames', col_names ...
);

% Display the table
disp(InitialCondition);


%%         [Vd0; Vq0; id0; iq0; psid0; psiq0; Efd0; psiH0; psiF0; psiK0; psiG0; omega0; Pm0; Xe0; Vref0; delta0; Xg0]
%%         [1  ; 2  ; 3  ; 4  ; 5    ; 6    ; 7   ; 8    ; 9    ; 10   ; 11   ; 12    ; 13 ; 14 ; 15   ; 16    ; 17  ]

% starting the simulation code

X_old_1 = InitialCondition_1;
X_old_2 = InitialCondition_2;

for k = 1:round(Ts/h)

    t = (k)*h;

    V_d_1(k)   = X_old_1(1);
    V_d_2(k)   = X_old_2(1);

    V_q_1(k)   = X_old_1(2);
    V_q_2(k)   = X_old_2(2);

    V_rms_1(k) = sqrt((V_d_1(k))^2 + (V_q_1(k))^2);
    V_rms_2(k) = sqrt((V_d_2(k))^2 + (V_q_2(k))^2);


    I_d_1(k)   = X_old_1(3);
    I_d_2(k)   = X_old_2(3);

    I_q_1(k)   = X_old_1(4);
    I_q_2(k)   = X_old_2(4);

    psi_d_1(k) = X_old_1(5);
    psi_d_2(k) = X_old_2(5);

    psi_q_1(k) = X_old_1(6);
    psi_q_2(k) = X_old_2(6);

    E_fd_1(k)  = X_old_1(7);
    E_fd_2(k)  = X_old_2(7);

    psi_H_1(k) = X_old_1(8);
    psi_H_2(k) = X_old_2(8);

    psi_F_1(k) = X_old_1(9);
    psi_F_2(k) = X_old_2(9);

    psi_K_1(k) = X_old_1(10);
    psi_K_2(k) = X_old_2(10);

    psi_G_2(k) = X_old_1(11);
    psi_G_2(k) = X_old_2(11);

    omega_1(k) = X_old_1(12);
    omega_2(k) = X_old_2(12);

    Pm_1(k) = X_old_1(13);
    Pm_2(k) = X_old_2(13);

    Xe_1(k) = X_old_1(14);
    Xe_2(k) = X_old_2(14);

    V_ref_1(k) = X_old_1(15);
    V_ref_2(k) = X_old_2(15);

    delta_1(k) = X_old_1(16);
    delta_2(k) = X_old_2(16);

    Xg_1(k) = X_old_1(17);
    Xg_2(k) = X_old_2(17);

    Ck = sqrt(2/3)*[cos(W_0*t)           , sin(W_0*t)           , sqrt(1/2);
                    cos(W_0*t - (2*pi/3)), sin(W_0*t - (2*pi/3)), sqrt(1/2);
                    cos(W_0*t + (2*pi/3)), sin(W_0*t + (2*pi/3)), sqrt(1/2)];

    V_dqo_1 = [V_d_1(k), V_q_1(k), 0];
    V_dqo_2 = [V_d_2(k), V_q_2(k), 0];

    I_dqo_1 = [I_d_1(k), I_q_1(k), 0];
    I_dqo_2 = [I_d_2(k), I_q_2(k), 0];

    psi_dqo_1 = [psi_d_1(k), psi_q_1(k), 0];
    psi_dqo_2 = [psi_d_2(k), psi_q_2(k), 0];

    V_abc_1 = Ck*V_dqo_1';
    V_abc_2 = Ck*V_dqo_2'; 

    I_abc_1 = Ck*I_dqo_1';
    I_abc_2 = Ck*I_dqo_2'; 

    psi_abc_1 = Ck*psi_dqo_1';
    psi_abc_2 = Ck*psi_dqo_2'; 

    V_a_1(k) = V_abc_1(1);
    V_b_1(k) = V_abc_1(2);
    V_c_1(k) = V_abc_1(3);

    I_a_1(k) = I_abc_1(1);
    I_b_1(k) = I_abc_1(2);
    I_c_1(k) = I_abc_1(3);

    psi_a_1(k) = psi_abc_1(1);
    psi_b_1(k) = psi_abc_1(2);
    psi_c_1(k) = psi_abc_1(3);

    V_a_2(k) = V_abc_2(1);
    V_b_2(k) = V_abc_2(2);
    V_c_2(k) = V_abc_2(3);

    I_a_2(k) = I_abc_2(1);
    I_b_2(k) = I_abc_2(2);
    I_c_2(k) = I_abc_2(3);

    psi_a_2(k) = psi_abc_2(1);
    psi_b_2(k) = psi_abc_2(2);
    psi_c_2(k) = psi_abc_2(3);

    P_line(k) = Pm_1(k) - (((V_rms_1(k))^2)/R_L_1);
    
    Pe_1(k) = (psi_d_1(k)*I_q_1(k)) - (psi_q_1(k)*I_d_1(k));
    Pe_2(k) = (psi_d_2(k)*I_q_2(k)) - (psi_q_2(k)*I_d_2(k));

    V_rms_1_3phase(k) = V_rms_1(k);
    V_rms_2_3phase(k) = V_rms_2(k);

    theta_1(k) = delta_1(k) + asin(V_d_1(k)/V_rms_1_3phase(k));
    theta_2(k) = delta_2(k) + asin(V_d_2(k)/V_rms_2_3phase(k));


        
end
%%


for k=round((Ts/h))+1:round(Tend1/h)
    t = k*h;
    R_L_1 = R_L_1_0;
    R_L_2 = R_L_2_0;


    % differential equations

    % equation (1)
    psi_H_1(k) = X_old_1(8) + h*(1/Td_dd)*(-X_old_1(8) + X_old_1(5));
    X_old_1(8) = psi_H_1(k);

    psi_H_2(k) = X_old_2(8) + h*(1/Td_dd)*(-X_old_2(8) + X_old_2(5));
    X_old_2(8) = psi_H_2(k);

    % equation (2)
    psi_F_1(k) = X_old_1(9) + h*(1/Td_d)*(-X_old_1(9) + X_old_1(5) + ((xd_d/(xd - xd_d))*X_old_1(7)));
    X_old_1(9) = psi_F_1(k);

    psi_F_2(k) = X_old_2(9) + h*(1/Td_d)*(-X_old_2(9) + X_old_2(5) + ((xd_d/(xd - xd_d))*X_old_2(7)));
    X_old_2(9) = psi_F_2(k);

    % equation (3)
    psi_K_1(k) = X_old_1(10) + h*(1/Tq_dd)*(-X_old_1(10) + X_old_1(6));
    X_old_1(10) = psi_K_1(k);

    psi_K_2(k) = X_old_2(10) + h*(1/Tq_dd)*(-X_old_2(10) + X_old_2(6));
    X_old_2(10) = psi_K_2(k);

    % equation (4)
    psi_G_1(k) = X_old_1(11) + h*(1/Tq_d)*(-X_old_1(11) + X_old_1(6));
    X_old_1(11) = psi_G_1(k);

    psi_G_2(k) = X_old_2(11) + h*(1/Tq_d)*(-X_old_2(11) + X_old_2(6));
    X_old_2(11) = psi_G_2(k);



     % Mechanical Equations

    delta_1(k) = X_old_1(16) + h*(X_old_1(12) - W_0);    % equation (9)
    X_old_1(16) = delta_1(k);

    omega_1(k) = X_old_1(12) + h*(W_b/(2*H))*(X_old_1(13) - ((X_old_1(5)*X_old_1(4)) - (X_old_1(6)*X_old_1(3))));  % equation (10)
    X_old_1(12) = omega_1(k);

    delta_2(k) = X_old_2(16) + h*(X_old_2(12) - W_0);    % equation (9)
    X_old_2(16) = delta_2(k);
    
    omega_2(k) = X_old_2(12) + h*(W_b/(2*H))*(X_old_2(13) - ((X_old_2(5)*X_old_2(4)) - (X_old_2(6)*X_old_2(3))));  % equation (10)
    X_old_2(12) = omega_2(k);

    


    % Exciter

    % Machine 1
    Xe_1(k)  = X_old_1(14) + h*(1/Ta)*(-X_old_1(14) + ka*(X_old_1(15) - V_rms_1(k-1)));

    if Xe_1(k) < Efd_min
        E_fd_1(k) = Efd_min;

    elseif Xe_1(k) > Efd_max
        E_fd_1(k) = Efd_max;

    elseif Xe_1(k) <= Efd_max && Xe_1(k) >= Efd_min
        E_fd_1(k) = Xe_1(k);

    end
    X_old_1(14) = Xe_1(k);
    X_old_1(7) = E_fd_1(k);

    % Machine 2
    Xe_2(k)  = X_old_2(14) + h*(1/Ta)*(-X_old_2(14) + ka*(X_old_2(15) - V_rms_2(k-1)));

    if Xe_2(k) < Efd_min
        E_fd_2(k) = Efd_min;

    elseif Xe_2(k) > Efd_max
        E_fd_2(k) = Efd_max;

    elseif Xe_2(k) <= Efd_max && Xe_2(k) >= Efd_min
        E_fd_2(k) = Xe_2(k);

    end
    X_old_2(14) = Xe_2(k);
    X_old_2(7) = E_fd_2(k);

    



    % Speed Governer
    if governer_1_enabled == 1
        Xg_1(k) = X_old_1(17) + h*(1/Tg2)*(-X_old_1(17) + ((W_ref - X_old_1(12))/W_b));
        Pm_1(k) = InitialCondition_1(13) + Kg*((1 - (Tg1/Tg2))*Xg_1(k) + ((Tg1/Tg2)*((W_ref - X_old_1(12))/W_b)));

        if Pm_1(k) < Pm_min
            Pm_1(k) = Pm_min;

        elseif Pm_1(k) > Pm_max
            Pm_1(k) = Pm_max;

        elseif Pm_1(k) <= Pm_max && Pm_1(k) >= Pm_min
            Pm_1(k) = Pm_1(k);

        end

    elseif governer_1_enabled == 0
        Pm_1(k) = InitialCondition_1(13);
        Xg_1(k) = 0;

    end
    X_old_1(17) = Xg_1(k);
    X_old_1(13) = Pm_1(k);

    if governer_2_enabled ==1
        Xg_2(k) = X_old_2(17) + h*(1/Tg2)*(-X_old_2(17) + ((W_ref - X_old_2(12))/W_b));
        Pm_2(k) = InitialCondition_2(13) + Kg*((1 - (Tg1/Tg2))*Xg_2(k) + ((Tg1/Tg2)*((W_ref - X_old_2(12))/W_b)));

        if Pm_2(k) < Pm_min
            Pm_2(k) = Pm_min;

        elseif Pm_2(k) > Pm_max
            Pm_2(k) = Pm_max;

        elseif Pm_2(k) <= Pm_max && Pm_2(k) >= Pm_min
            Pm_2(k) = Pm_2(k);

        end

    elseif governer_2_enabled == 0
        Pm_2(k) = InitialCondition_2(13);
        Xg_2(k) = 0;

    end
    X_old_2(17) = Xg_2(k);
    X_old_2(13) = Pm_2(k);


    % algebraic equations

    C_KG_delta_1 = (((xq_d - xq_dd)/xq_d)*X_old_1(10) + ((xq - xq_d)/xq)*(xq_dd/xq_d)*X_old_1(11))*cos(X_old_1(16));
    S_HF_delta_1 = (((xd_d - xd_dd)/xd_d)*X_old_1(8) + ((xd - xd_d)/xd)*(xd_dd/xd_d)*X_old_1(9))*sin(X_old_1(16));

    S_KG_delta_1 = (((xq_d - xq_dd)/xq_d)*X_old_1(10) + ((xq - xq_d)/xq)*(xq_dd/xq_d)*X_old_1(11))*sin(X_old_1(16));
    C_HF_delta_1 = (((xd_d - xd_dd)/xd_d)*X_old_1(8) + ((xd - xd_d)/xd)*(xd_dd/xd_d)*X_old_1(9))*cos(X_old_1(16));

    C_KG_delta_2 = (((xq_d - xq_dd)/xq_d)*X_old_2(10) + ((xq - xq_d)/xq)*(xq_dd/xq_d)*X_old_2(11))*cos(X_old_2(16));
    S_HF_delta_2 = (((xd_d - xd_dd)/xd_d)*X_old_2(8) + ((xd - xd_d)/xd)*(xd_dd/xd_d)*X_old_2(9))*sin(X_old_2(16));

    S_KG_delta_2 = (((xq_d - xq_dd)/xq_d)*X_old_2(10) + ((xq - xq_d)/xq)*(xq_dd/xq_d)*X_old_2(11))*sin(X_old_2(16));
    C_HF_delta_2 = (((xd_d - xd_dd)/xd_d)*X_old_2(8) + ((xd - xd_d)/xd)*(xd_dd/xd_d)*X_old_2(9))*cos(X_old_2(16));



    b = [C_KG_delta_1 - S_HF_delta_1;
        S_KG_delta_1 + C_HF_delta_1;
        C_KG_delta_2 - S_HF_delta_2;
        S_KG_delta_2 + C_HF_delta_2;
        0                      ;
        0                      ;
        0                      ;
        0                      ;
        0                      ;
        0                      ;
        0                      ;
        0                      ;
        0                      ;
        0                      ;];

    A = [0     , 1     , 0     , 0     , 0      , -xq_dd    , 0         , 0         , 0             , 0             , 0             , 0             , 0                 , 0                 ;
        1      , 0     , 0     , 0     , -xd_dd , 0         , 0         , 0         , 0             , 0             , 0             , 0             , 0                 , 0                 ;
        0      , 0     , 0     , 1     , 0      , 0         , 0         , -xq_dd    , 0             , 0             , 0             , 0             , 0                 , 0                 ;
        0      , 0     , 1     , 0     , 0      , 0         , -xd_dd    , 0         , 0             , 0             , 0             , 0             , 0                 , 0                 ;
        0      , 1     , 0     , 0     , Ra     , 0         , 0         , 0         , 1             , 0             , 0             , 0             , 0                 , 0                 ;
        -1     , 0     , 0     , 0     , 0      , Ra        , 0         , 0         , 0             , 1             , 0             , 0             , 0                 , 0                 ;
        0      , 0     , 0     , 1     , 0      , 0         , Ra        , 0         , 0             , 0             , 1             , 0             , 0                 , 0                 ;
        0      , 0     , -1    , 0     , 0      , 0         , 0         , Ra        , 0             , 0             , 0             , 1             , 0                 , 0                 ;
        0      , 0     , 0     , 0     , 0      , -R_L_1    , 0         , 0         , 0             , 1             , 0             , 0             , 0                 , R_L_1             ;
        0      , 0     , 0     , 0     , -R_L_1 , 0         , 0         , 0         , 1             , 0             , 0             , 0             , R_L_1             , 0                 ;
        0      , 0     , 0     , 0     , 0      , 0         , 0         , -R_L_2    , 0             , 0             , 0             , 1             , 0                 , -R_L_2            ;
        0      , 0     , 0     , 0     , 0      , 0         , -R_L_2    , 0         , 0             , 0             , 1             , 0             , -R_L_2            , 0                 ;
        0      , 0     , 0     , 0     , 0      , 0         , 0         , 0         , W_b/x_line    , 0             , -W_b/x_line   , 0             , -r_line*W_b/x_line, -W_b              ;
        0      , 0     , 0     , 0     , 0      , 0         , 0         , 0         , 0             , W_b/x_line    , 0             , -W_b/x_line   , W_b               , -r_line*W_b/x_line    ];

    Y = A\b;

    % bringing it back to machine reference

    psi_d_1(k)   = [1 0]*[cos(X_old_1(16)), -sin(X_old_1(16));
        sin(X_old_1(16)),  cos(X_old_1(16))]*Y(1:2);

    X_old_1(5) = psi_d_1(k);

    psi_q_1(k)  = [0 1]*[cos(X_old_1(16)), -sin(X_old_1(16));
        sin(X_old_1(16)),  cos(X_old_1(16))]*Y(1:2);

    X_old_1(6) = psi_q_1(k);


    psi_d_2(k)   = [1 0]*[cos(X_old_2(16)), -sin(X_old_2(16));
        sin(X_old_2(16)),  cos(X_old_2(16))]*Y(3:4);

    X_old_2(5) = psi_d_2(k);

    psi_q_2(k)   = [0 1]*[cos(X_old_2(16)), -sin(X_old_2(16));
        sin(X_old_2(16)),  cos(X_old_2(16))]*Y(3:4);

    X_old_2(6) = psi_q_2(k);



    I_d_1(k)     = [1 0]*[cos(X_old_1(16)), -sin(X_old_1(16));
        sin(X_old_1(16)),  cos(X_old_1(16))]*Y(5:6);

    X_old_1(3) = I_d_1(k);

    I_q_1(k)     = [0 1]*[cos(X_old_1(16)), -sin(X_old_1(16));
        sin(X_old_1(16)),  cos(X_old_1(16))]*Y(5:6);

    X_old_1(4) = I_q_1(k);


    I_d_2(k)     = [1 0]*[cos(X_old_2(16)), -sin(X_old_2(16));
        sin(X_old_2(16)),  cos(X_old_2(16))]*Y(7:8);

    X_old_2(3) = I_d_2(k);

    I_q_2(k)     = [0 1]*[cos(X_old_2(16)), -sin(X_old_2(16));
        sin(X_old_2(16)),  cos(X_old_2(16))]*Y(7:8);

    X_old_2(4) = I_q_2(k);



    V_d_1(k)     = [1 0]*[cos(X_old_1(16)), -sin(X_old_1(16));
        sin(X_old_1(16)),  cos(X_old_1(16))]*Y(9:10);

    X_old_1(1) = V_d_1(k);

    V_q_1(k)     = [0 1]*[cos(X_old_1(16)), -sin(X_old_1(16));
        sin(X_old_1(16)),  cos(X_old_1(16))]*Y(9:10);

    X_old_1(2) = V_q_1(k);

    V_rms_1(k) = sqrt((V_d_1(k))^2 + (V_q_1(k))^2);

    V_d_2(k)     = [1 0]*[cos(X_old_2(16)), -sin(X_old_2(16));
        sin(X_old_2(16)),  cos(X_old_2(16))]*Y(11:12);

    X_old_2(1) = V_d_2(k);

    V_q_2(k)     = [0 1]*[cos(X_old_2(16)), -sin(X_old_2(16));
        sin(X_old_2(16)),  cos(X_old_2(16))]*Y(11:12);

    X_old_2(2) = V_q_2(k);

    V_rms_2(k) = sqrt((V_d_2(k))^2 + (V_q_2(k))^2);

   

    

    % converting to abc frame

    Ck = sqrt(2/3)*[cos(W_0*t)           , sin(W_0*t)           ,   sqrt(1/2);
        cos(W_0*t - (2*pi/3)), sin(W_0*t - (2*pi/3)),   sqrt(1/2);
        cos(W_0*t + (2*pi/3)), sin(W_0*t + (2*pi/3)),   sqrt(1/2)];

    V_dqo_1 = [V_d_1(k), V_q_1(k), 0];
    V_dqo_2 = [V_d_2(k), V_q_2(k), 0];

    I_dqo_1 = [I_d_1(k), I_q_1(k), 0];
    I_dqo_2 = [I_d_2(k), I_q_2(k), 0];

    psi_dqo_1 = [psi_d_1(k), psi_q_1(k), 0];
    psi_dqo_2 = [psi_d_2(k), psi_q_2(k), 0];

    V_abc_1 = Ck*V_dqo_1';
    V_abc_2 = Ck*V_dqo_2';

    I_abc_1 = Ck*I_dqo_1';
    I_abc_2 = Ck*I_dqo_2';

    psi_abc_1 = Ck*psi_dqo_1';
    psi_abc_2 = Ck*psi_dqo_2';

    V_a_1(k) = V_abc_1(1);
    V_b_1(k) = V_abc_1(2);
    V_c_1(k) = V_abc_1(3);

    I_a_1(k) = I_abc_1(1);
    I_b_1(k) = I_abc_1(2);
    I_c_1(k) = I_abc_1(3);

    psi_a_1(k) = psi_abc_1(1);
    psi_b_1(k) = psi_abc_1(2);
    psi_c_1(k) = psi_abc_1(3);

    V_a_2(k) = V_abc_2(1);
    V_b_2(k) = V_abc_2(2);
    V_c_2(k) = V_abc_2(3);

    I_a_2(k) = I_abc_2(1);
    I_b_2(k) = I_abc_2(2);
    I_c_2(k) = I_abc_2(3);

    psi_a_2(k) = psi_abc_2(1);
    psi_b_2(k) = psi_abc_2(2);
    psi_c_2(k) = psi_abc_2(3);

    P_line(k) = Pm_1(k) - ((V_rms_1(k)^2)/R_L_1);
    
    Pe_1(k) = (psi_d_1(k)*I_q_1(k)) - (psi_q_1(k)*I_d_1(k));
    Pe_2(k) = (psi_d_2(k)*I_q_2(k)) - (psi_q_2(k)*I_d_2(k));

    V_rms_1_3phase(k) = V_rms_1(k);
    V_rms_2_3phase(k) = V_rms_2(k);

    theta_1(k) = delta_1(k) + asin(V_d_1(k)/V_rms_1_3phase(k));
    theta_2(k) = delta_2(k) + asin(V_d_2(k)/V_rms_2_3phase(k));
end

if fault == 0
        
    main_title = 'Load'
    for k=round((Tend1/h))+1:round(Tend2/h)
        t = k*h;
        R_L_1 = R_L_1_new;
        R_L_2 = R_L_2_new;
        r_line = r_line_new;
        x_line = x_line_new;


        % differential equations

        % equation (1)
        psi_H_1(k) = X_old_1(8) + h*(1/Td_dd)*(-X_old_1(8) + X_old_1(5));
        X_old_1(8) = psi_H_1(k);

        psi_H_2(k) = X_old_2(8) + h*(1/Td_dd)*(-X_old_2(8) + X_old_2(5));
        X_old_2(8) = psi_H_2(k);

        % equation (2)
        psi_F_1(k) = X_old_1(9) + h*(1/Td_d)*(-X_old_1(9) + X_old_1(5) + ((xd_d/(xd - xd_d))*X_old_1(7)));
        X_old_1(9) = psi_F_1(k);

        psi_F_2(k) = X_old_2(9) + h*(1/Td_d)*(-X_old_2(9) + X_old_2(5) + ((xd_d/(xd - xd_d))*X_old_2(7)));
        X_old_2(9) = psi_F_2(k);

        % equation (3)
        psi_K_1(k) = X_old_1(10) + h*(1/Tq_dd)*(-X_old_1(10) + X_old_1(6));
        X_old_1(10) = psi_K_1(k);

        psi_K_2(k) = X_old_2(10) + h*(1/Tq_dd)*(-X_old_2(10) + X_old_2(6));
        X_old_2(10) = psi_K_2(k);

        % equation (4)
        psi_G_1(k) = X_old_1(11) + h*(1/Tq_d)*(-X_old_1(11) + X_old_1(6));
        X_old_1(11) = psi_G_1(k);

        psi_G_2(k) = X_old_2(11) + h*(1/Tq_d)*(-X_old_2(11) + X_old_2(6));
        X_old_2(11) = psi_G_2(k);



        % Mechanical Equations

        delta_1(k) = X_old_1(16) + h*(X_old_1(12) - W_0);    % equation (9)
        X_old_1(16) = delta_1(k);

        omega_1(k) = X_old_1(12) + h*(W_b/(2*H))*(X_old_1(13) - ((X_old_1(5)*X_old_1(4)) - (X_old_1(6)*X_old_1(3))));  % equation (10)
        X_old_1(12) = omega_1(k);

        delta_2(k) = X_old_2(16) + h*(X_old_2(12) - W_0);    % equation (9)
        X_old_2(16) = delta_2(k);

        omega_2(k) = X_old_2(12) + h*(W_b/(2*H))*(X_old_2(13) - ((X_old_2(5)*X_old_2(4)) - (X_old_2(6)*X_old_2(3))));  % equation (10)
        X_old_2(12) = omega_2(k);




        % Exciter

        % Machine 1
        Xe_1(k)  = X_old_1(14) + h*(1/Ta)*(-X_old_1(14) + ka*(X_old_1(15) - V_rms_1(k-1)));

        if Xe_1(k) < Efd_min
            E_fd_1(k) = Efd_min;

        elseif Xe_1(k) > Efd_max
            E_fd_1(k) = Efd_max;

        elseif Xe_1(k) <= Efd_max && Xe_1(k) >= Efd_min
            E_fd_1(k) = Xe_1(k);

        end
        X_old_1(14) = Xe_1(k);
        X_old_1(7) = E_fd_1(k);

        % Machine 2
        Xe_2(k)  = X_old_2(14) + h*(1/Ta)*(-X_old_2(14) + ka*(X_old_2(15) - V_rms_2(k-1)));

        if Xe_2(k) < Efd_min
            E_fd_2(k) = Efd_min;

        elseif Xe_2(k) > Efd_max
            E_fd_2(k) = Efd_max;

        elseif Xe_2(k) <= Efd_max && Xe_2(k) >= Efd_min
            E_fd_2(k) = Xe_2(k);

        end
        X_old_2(14) = Xe_2(k);
        X_old_2(7) = E_fd_2(k);





        % Speed Governer
        if governer_1_enabled == 1
            Xg_1(k) = X_old_1(17) + h*(1/Tg2)*(-X_old_1(17) + ((W_ref - X_old_1(12))/W_b));
            Pm_1(k) = InitialCondition_1(13) + Kg*((1 - (Tg1/Tg2))*Xg_1(k) + ((Tg1/Tg2)*((W_ref - X_old_1(12))/W_b)));

            if Pm_1(k) < Pm_min
                Pm_1(k) = Pm_min;

            elseif Pm_1(k) > Pm_max
                Pm_1(k) = Pm_max;

            elseif Pm_1(k) <= Pm_max && Pm_1(k) >= Pm_min
                Pm_1(k) = Pm_1(k);

            end

        elseif governer_1_enabled == 0
            Pm_1(k) = InitialCondition_1(13);
            Xg_1(k) = 0;

        end
        X_old_1(17) = Xg_1(k);
        X_old_1(13) = Pm_1(k);

        if governer_2_enabled ==1
            Xg_2(k) = X_old_2(17) + h*(1/Tg2)*(-X_old_2(17) + ((W_ref - X_old_2(12))/W_b));
            Pm_2(k) = InitialCondition_2(13) + Kg*((1 - (Tg1/Tg2))*Xg_2(k) + ((Tg1/Tg2)*((W_ref - X_old_2(12))/W_b)));

            if Pm_2(k) < Pm_min
                Pm_2(k) = Pm_min;

            elseif Pm_2(k) > Pm_max
                Pm_2(k) = Pm_max;

            elseif Pm_2(k) <= Pm_max && Pm_2(k) >= Pm_min
                Pm_2(k) = Pm_2(k);

            end

        elseif governer_2_enabled == 0
            Pm_2(k) = InitialCondition_2(13);
            Xg_2(k) = 0;

        end
        X_old_2(17) = Xg_2(k);
        X_old_2(13) = Pm_2(k);


        % algebraic equations

        C_KG_delta_1 = (((xq_d - xq_dd)/xq_d)*X_old_1(10) + ((xq - xq_d)/xq)*(xq_dd/xq_d)*X_old_1(11))*cos(X_old_1(16));
        S_HF_delta_1 = (((xd_d - xd_dd)/xd_d)*X_old_1(8) + ((xd - xd_d)/xd)*(xd_dd/xd_d)*X_old_1(9))*sin(X_old_1(16));

        S_KG_delta_1 = (((xq_d - xq_dd)/xq_d)*X_old_1(10) + ((xq - xq_d)/xq)*(xq_dd/xq_d)*X_old_1(11))*sin(X_old_1(16));
        C_HF_delta_1 = (((xd_d - xd_dd)/xd_d)*X_old_1(8) + ((xd - xd_d)/xd)*(xd_dd/xd_d)*X_old_1(9))*cos(X_old_1(16));

        C_KG_delta_2 = (((xq_d - xq_dd)/xq_d)*X_old_2(10) + ((xq - xq_d)/xq)*(xq_dd/xq_d)*X_old_2(11))*cos(X_old_2(16));
        S_HF_delta_2 = (((xd_d - xd_dd)/xd_d)*X_old_2(8) + ((xd - xd_d)/xd)*(xd_dd/xd_d)*X_old_2(9))*sin(X_old_2(16));

        S_KG_delta_2 = (((xq_d - xq_dd)/xq_d)*X_old_2(10) + ((xq - xq_d)/xq)*(xq_dd/xq_d)*X_old_2(11))*sin(X_old_2(16));
        C_HF_delta_2 = (((xd_d - xd_dd)/xd_d)*X_old_2(8) + ((xd - xd_d)/xd)*(xd_dd/xd_d)*X_old_2(9))*cos(X_old_2(16));



        b = [C_KG_delta_1 - S_HF_delta_1;
            S_KG_delta_1 + C_HF_delta_1;
            C_KG_delta_2 - S_HF_delta_2;
            S_KG_delta_2 + C_HF_delta_2;
            0                      ;
            0                      ;
            0                      ;
            0                      ;
            0                      ;
            0                      ;
            0                      ;
            0                      ;
            0                      ;
            0                      ;];

        A = [0     , 1     , 0     , 0     , 0      , -xq_dd    , 0         , 0         , 0             , 0             , 0             , 0             , 0                 , 0                 ;
            1      , 0     , 0     , 0     , -xd_dd , 0         , 0         , 0         , 0             , 0             , 0             , 0             , 0                 , 0                 ;
            0      , 0     , 0     , 1     , 0      , 0         , 0         , -xq_dd    , 0             , 0             , 0             , 0             , 0                 , 0                 ;
            0      , 0     , 1     , 0     , 0      , 0         , -xd_dd    , 0         , 0             , 0             , 0             , 0             , 0                 , 0                 ;
            0      , 1     , 0     , 0     , Ra     , 0         , 0         , 0         , 1             , 0             , 0             , 0             , 0                 , 0                 ;
            -1     , 0     , 0     , 0     , 0      , Ra        , 0         , 0         , 0             , 1             , 0             , 0             , 0                 , 0                 ;
            0      , 0     , 0     , 1     , 0      , 0         , Ra        , 0         , 0             , 0             , 1             , 0             , 0                 , 0                 ;
            0      , 0     , -1    , 0     , 0      , 0         , 0         , Ra        , 0             , 0             , 0             , 1             , 0                 , 0                 ;
            0      , 0     , 0     , 0     , 0      , -R_L_1    , 0         , 0         , 0             , 1             , 0             , 0             , 0                 , R_L_1             ;
            0      , 0     , 0     , 0     , -R_L_1 , 0         , 0         , 0         , 1             , 0             , 0             , 0             , R_L_1             , 0                 ;
            0      , 0     , 0     , 0     , 0      , 0         , 0         , -R_L_2    , 0             , 0             , 0             , 1             , 0                 , -R_L_2            ;
            0      , 0     , 0     , 0     , 0      , 0         , -R_L_2    , 0         , 0             , 0             , 1             , 0             , -R_L_2            , 0                 ;
            0      , 0     , 0     , 0     , 0      , 0         , 0         , 0         , W_b/x_line    , 0             , -W_b/x_line   , 0             , -r_line*W_b/x_line, -W_b              ;
            0      , 0     , 0     , 0     , 0      , 0         , 0         , 0         , 0             , W_b/x_line    , 0             , -W_b/x_line   , W_b               , -r_line*W_b/x_line    ];

        Y = A\b;

        % bringing it back to machine reference

        psi_d_1(k)   = [1 0]*[cos(X_old_1(16)), -sin(X_old_1(16));
            sin(X_old_1(16)),  cos(X_old_1(16))]*Y(1:2);

        X_old_1(5) = psi_d_1(k);

        psi_q_1(k)  = [0 1]*[cos(X_old_1(16)), -sin(X_old_1(16));
            sin(X_old_1(16)),  cos(X_old_1(16))]*Y(1:2);

        X_old_1(6) = psi_q_1(k);


        psi_d_2(k)   = [1 0]*[cos(X_old_2(16)), -sin(X_old_2(16));
            sin(X_old_2(16)),  cos(X_old_2(16))]*Y(3:4);

        X_old_2(5) = psi_d_2(k);

        psi_q_2(k)   = [0 1]*[cos(X_old_2(16)), -sin(X_old_2(16));
            sin(X_old_2(16)),  cos(X_old_2(16))]*Y(3:4);

        X_old_2(6) = psi_q_2(k);



        I_d_1(k)     = [1 0]*[cos(X_old_1(16)), -sin(X_old_1(16));
            sin(X_old_1(16)),  cos(X_old_1(16))]*Y(5:6);

        X_old_1(3) = I_d_1(k);

        I_q_1(k)     = [0 1]*[cos(X_old_1(16)), -sin(X_old_1(16));
            sin(X_old_1(16)),  cos(X_old_1(16))]*Y(5:6);

        X_old_1(4) = I_q_1(k);


        I_d_2(k)     = [1 0]*[cos(X_old_2(16)), -sin(X_old_2(16));
            sin(X_old_2(16)),  cos(X_old_2(16))]*Y(7:8);

        X_old_2(3) = I_d_2(k);

        I_q_2(k)     = [0 1]*[cos(X_old_2(16)), -sin(X_old_2(16));
            sin(X_old_2(16)),  cos(X_old_2(16))]*Y(7:8);

        X_old_2(4) = I_q_2(k);



        V_d_1(k)     = [1 0]*[cos(X_old_1(16)), -sin(X_old_1(16));
            sin(X_old_1(16)),  cos(X_old_1(16))]*Y(9:10);

        X_old_1(1) = V_d_1(k);

        V_q_1(k)     = [0 1]*[cos(X_old_1(16)), -sin(X_old_1(16));
            sin(X_old_1(16)),  cos(X_old_1(16))]*Y(9:10);

        X_old_1(2) = V_q_1(k);

        V_rms_1(k) = sqrt((V_d_1(k))^2 + (V_q_1(k))^2);

        V_d_2(k)     = [1 0]*[cos(X_old_2(16)), -sin(X_old_2(16));
            sin(X_old_2(16)),  cos(X_old_2(16))]*Y(11:12);

        X_old_2(1) = V_d_2(k);

        V_q_2(k)     = [0 1]*[cos(X_old_2(16)), -sin(X_old_2(16));
            sin(X_old_2(16)),  cos(X_old_2(16))]*Y(11:12);

        X_old_2(2) = V_q_2(k);

        V_rms_2(k) = sqrt((V_d_2(k))^2 + (V_q_2(k))^2);





        % converting to abc frame

        Ck = sqrt(2/3)*[cos(W_0*t)           , sin(W_0*t)           ,   sqrt(1/2);
            cos(W_0*t - (2*pi/3)), sin(W_0*t - (2*pi/3)),   sqrt(1/2);
            cos(W_0*t + (2*pi/3)), sin(W_0*t + (2*pi/3)),   sqrt(1/2)];

        V_dqo_1 = [V_d_1(k), V_q_1(k), 0];
        V_dqo_2 = [V_d_2(k), V_q_2(k), 0];

        I_dqo_1 = [I_d_1(k), I_q_1(k), 0];
        I_dqo_2 = [I_d_2(k), I_q_2(k), 0];

        psi_dqo_1 = [psi_d_1(k), psi_q_1(k), 0];
        psi_dqo_2 = [psi_d_2(k), psi_q_2(k), 0];

        V_abc_1 = Ck*V_dqo_1';
        V_abc_2 = Ck*V_dqo_2';

        I_abc_1 = Ck*I_dqo_1';
        I_abc_2 = Ck*I_dqo_2';

        psi_abc_1 = Ck*psi_dqo_1';
        psi_abc_2 = Ck*psi_dqo_2';

        V_a_1(k) = V_abc_1(1);
        V_b_1(k) = V_abc_1(2);
        V_c_1(k) = V_abc_1(3);

        I_a_1(k) = I_abc_1(1);
        I_b_1(k) = I_abc_1(2);
        I_c_1(k) = I_abc_1(3);

        psi_a_1(k) = psi_abc_1(1);
        psi_b_1(k) = psi_abc_1(2);
        psi_c_1(k) = psi_abc_1(3);

        V_a_2(k) = V_abc_2(1);
        V_b_2(k) = V_abc_2(2);
        V_c_2(k) = V_abc_2(3);

        I_a_2(k) = I_abc_2(1);
        I_b_2(k) = I_abc_2(2);
        I_c_2(k) = I_abc_2(3);

        psi_a_2(k) = psi_abc_2(1);
        psi_b_2(k) = psi_abc_2(2);
        psi_c_2(k) = psi_abc_2(3);

        P_line(k) = Pm_1(k) - ((V_rms_1(k)^2)/R_L_1);

        Pe_1(k) = (psi_d_1(k)*I_q_1(k)) - (psi_q_1(k)*I_d_1(k));
        Pe_2(k) = (psi_d_2(k)*I_q_2(k)) - (psi_q_2(k)*I_d_2(k));

        V_rms_1_3phase(k) = V_rms_1(k);
        V_rms_2_3phase(k) = V_rms_2(k);

        theta_1(k) = delta_1(k) + asin(V_d_1(k)/V_rms_1_3phase(k));
        theta_2(k) = delta_2(k) + asin(V_d_2(k)/V_rms_2_3phase(k));
    end
end

if fault == 1
    for k=round((Tend1/h))+1:round(Tfault_remove/h)
        t = k*h;
        R_L_1 = R_L_1_new;
        R_L_2 = R_L_2_new;
        r_line = r_line_new;
        x_line = x_line_new;


        % differential equations

        % equation (1)
        psi_H_1(k) = X_old_1(8) + h*(1/Td_dd)*(-X_old_1(8) + X_old_1(5));
        X_old_1(8) = psi_H_1(k);

        psi_H_2(k) = X_old_2(8) + h*(1/Td_dd)*(-X_old_2(8) + X_old_2(5));
        X_old_2(8) = psi_H_2(k);

        % equation (2)
        psi_F_1(k) = X_old_1(9) + h*(1/Td_d)*(-X_old_1(9) + X_old_1(5) + ((xd_d/(xd - xd_d))*X_old_1(7)));
        X_old_1(9) = psi_F_1(k);

        psi_F_2(k) = X_old_2(9) + h*(1/Td_d)*(-X_old_2(9) + X_old_2(5) + ((xd_d/(xd - xd_d))*X_old_2(7)));
        X_old_2(9) = psi_F_2(k);

        % equation (3)
        psi_K_1(k) = X_old_1(10) + h*(1/Tq_dd)*(-X_old_1(10) + X_old_1(6));
        X_old_1(10) = psi_K_1(k);

        psi_K_2(k) = X_old_2(10) + h*(1/Tq_dd)*(-X_old_2(10) + X_old_2(6));
        X_old_2(10) = psi_K_2(k);

        % equation (4)
        psi_G_1(k) = X_old_1(11) + h*(1/Tq_d)*(-X_old_1(11) + X_old_1(6));
        X_old_1(11) = psi_G_1(k);

        psi_G_2(k) = X_old_2(11) + h*(1/Tq_d)*(-X_old_2(11) + X_old_2(6));
        X_old_2(11) = psi_G_2(k);



        % Mechanical Equations

        delta_1(k) = X_old_1(16) + h*(X_old_1(12) - W_0);    % equation (9)
        X_old_1(16) = delta_1(k);

        omega_1(k) = X_old_1(12) + h*(W_b/(2*H))*(X_old_1(13) - ((X_old_1(5)*X_old_1(4)) - (X_old_1(6)*X_old_1(3))));  % equation (10)
        X_old_1(12) = omega_1(k);

        delta_2(k) = X_old_2(16) + h*(X_old_2(12) - W_0);    % equation (9)
        X_old_2(16) = delta_2(k);

        omega_2(k) = X_old_2(12) + h*(W_b/(2*H))*(X_old_2(13) - ((X_old_2(5)*X_old_2(4)) - (X_old_2(6)*X_old_2(3))));  % equation (10)
        X_old_2(12) = omega_2(k);




        % Exciter

        % Machine 1
        Xe_1(k)  = X_old_1(14) + h*(1/Ta)*(-X_old_1(14) + ka*(X_old_1(15) - V_rms_1(k-1)));

        if Xe_1(k) < Efd_min
            E_fd_1(k) = Efd_min;

        elseif Xe_1(k) > Efd_max
            E_fd_1(k) = Efd_max;

        elseif Xe_1(k) <= Efd_max && Xe_1(k) >= Efd_min
            E_fd_1(k) = Xe_1(k);

        end
        X_old_1(14) = Xe_1(k);
        X_old_1(7) = E_fd_1(k);

        % Machine 2
        Xe_2(k)  = X_old_2(14) + h*(1/Ta)*(-X_old_2(14) + ka*(X_old_2(15) - V_rms_2(k-1)));

        if Xe_2(k) < Efd_min
            E_fd_2(k) = Efd_min;

        elseif Xe_2(k) > Efd_max
            E_fd_2(k) = Efd_max;

        elseif Xe_2(k) <= Efd_max && Xe_2(k) >= Efd_min
            E_fd_2(k) = Xe_2(k);

        end
        X_old_2(14) = Xe_2(k);
        X_old_2(7) = E_fd_2(k);





        % Speed Governer
        if governer_1_enabled == 1
            Xg_1(k) = X_old_1(17) + h*(1/Tg2)*(-X_old_1(17) + ((W_ref - X_old_1(12))/W_b));
            Pm_1(k) = InitialCondition_1(13) + Kg*((1 - (Tg1/Tg2))*Xg_1(k) + ((Tg1/Tg2)*((W_ref - X_old_1(12))/W_b)));

            if Pm_1(k) < Pm_min
                Pm_1(k) = Pm_min;

            elseif Pm_1(k) > Pm_max
                Pm_1(k) = Pm_max;

            elseif Pm_1(k) <= Pm_max && Pm_1(k) >= Pm_min
                Pm_1(k) = Pm_1(k);

            end

        elseif governer_1_enabled == 0
            Pm_1(k) = InitialCondition_1(13);
            Xg_1(k) = 0;

        end
        X_old_1(17) = Xg_1(k);
        X_old_1(13) = Pm_1(k);

        if governer_2_enabled ==1
            Xg_2(k) = X_old_2(17) + h*(1/Tg2)*(-X_old_2(17) + ((W_ref - X_old_2(12))/W_b));
            Pm_2(k) = InitialCondition_2(13) + Kg*((1 - (Tg1/Tg2))*Xg_2(k) + ((Tg1/Tg2)*((W_ref - X_old_2(12))/W_b)));

            if Pm_2(k) < Pm_min
                Pm_2(k) = Pm_min;

            elseif Pm_2(k) > Pm_max
                Pm_2(k) = Pm_max;

            elseif Pm_2(k) <= Pm_max && Pm_2(k) >= Pm_min
                Pm_2(k) = Pm_2(k);

            end

        elseif governer_2_enabled == 0
            Pm_2(k) = InitialCondition_2(13);
            Xg_2(k) = 0;

        end
        X_old_2(17) = Xg_2(k);
        X_old_2(13) = Pm_2(k);


        % algebraic equations

        C_KG_delta_1 = (((xq_d - xq_dd)/xq_d)*X_old_1(10) + ((xq - xq_d)/xq)*(xq_dd/xq_d)*X_old_1(11))*cos(X_old_1(16));
        S_HF_delta_1 = (((xd_d - xd_dd)/xd_d)*X_old_1(8) + ((xd - xd_d)/xd)*(xd_dd/xd_d)*X_old_1(9))*sin(X_old_1(16));

        S_KG_delta_1 = (((xq_d - xq_dd)/xq_d)*X_old_1(10) + ((xq - xq_d)/xq)*(xq_dd/xq_d)*X_old_1(11))*sin(X_old_1(16));
        C_HF_delta_1 = (((xd_d - xd_dd)/xd_d)*X_old_1(8) + ((xd - xd_d)/xd)*(xd_dd/xd_d)*X_old_1(9))*cos(X_old_1(16));

        C_KG_delta_2 = (((xq_d - xq_dd)/xq_d)*X_old_2(10) + ((xq - xq_d)/xq)*(xq_dd/xq_d)*X_old_2(11))*cos(X_old_2(16));
        S_HF_delta_2 = (((xd_d - xd_dd)/xd_d)*X_old_2(8) + ((xd - xd_d)/xd)*(xd_dd/xd_d)*X_old_2(9))*sin(X_old_2(16));

        S_KG_delta_2 = (((xq_d - xq_dd)/xq_d)*X_old_2(10) + ((xq - xq_d)/xq)*(xq_dd/xq_d)*X_old_2(11))*sin(X_old_2(16));
        C_HF_delta_2 = (((xd_d - xd_dd)/xd_d)*X_old_2(8) + ((xd - xd_d)/xd)*(xd_dd/xd_d)*X_old_2(9))*cos(X_old_2(16));



        b = [C_KG_delta_1 - S_HF_delta_1;
            S_KG_delta_1 + C_HF_delta_1;
            C_KG_delta_2 - S_HF_delta_2;
            S_KG_delta_2 + C_HF_delta_2;
            0                      ;
            0                      ;
            0                      ;
            0                      ;
            0                      ;
            0                      ;
            0                      ;
            0                      ;
            0                      ;
            0                      ;];

        A = [0     , 1     , 0     , 0     , 0      , -xq_dd    , 0         , 0         , 0             , 0             , 0             , 0             , 0                 , 0                 ;
            1      , 0     , 0     , 0     , -xd_dd , 0         , 0         , 0         , 0             , 0             , 0             , 0             , 0                 , 0                 ;
            0      , 0     , 0     , 1     , 0      , 0         , 0         , -xq_dd    , 0             , 0             , 0             , 0             , 0                 , 0                 ;
            0      , 0     , 1     , 0     , 0      , 0         , -xd_dd    , 0         , 0             , 0             , 0             , 0             , 0                 , 0                 ;
            0      , 1     , 0     , 0     , Ra     , 0         , 0         , 0         , 1             , 0             , 0             , 0             , 0                 , 0                 ;
            -1     , 0     , 0     , 0     , 0      , Ra        , 0         , 0         , 0             , 1             , 0             , 0             , 0                 , 0                 ;
            0      , 0     , 0     , 1     , 0      , 0         , Ra        , 0         , 0             , 0             , 1             , 0             , 0                 , 0                 ;
            0      , 0     , -1    , 0     , 0      , 0         , 0         , Ra        , 0             , 0             , 0             , 1             , 0                 , 0                 ;
            0      , 0     , 0     , 0     , 0      , -R_L_1    , 0         , 0         , 0             , 1             , 0             , 0             , 0                 , R_L_1             ;
            0      , 0     , 0     , 0     , -R_L_1 , 0         , 0         , 0         , 1             , 0             , 0             , 0             , R_L_1             , 0                 ;
            0      , 0     , 0     , 0     , 0      , 0         , 0         , -R_L_2    , 0             , 0             , 0             , 1             , 0                 , -R_L_2            ;
            0      , 0     , 0     , 0     , 0      , 0         , -R_L_2    , 0         , 0             , 0             , 1             , 0             , -R_L_2            , 0                 ;
            0      , 0     , 0     , 0     , 0      , 0         , 0         , 0         , W_b/x_line    , 0             , -W_b/x_line   , 0             , -r_line*W_b/x_line, -W_b              ;
            0      , 0     , 0     , 0     , 0      , 0         , 0         , 0         , 0             , W_b/x_line    , 0             , -W_b/x_line   , W_b               , -r_line*W_b/x_line    ];

        Y = A\b;

        % bringing it back to machine reference

        psi_d_1(k)   = [1 0]*[cos(X_old_1(16)), -sin(X_old_1(16));
            sin(X_old_1(16)),  cos(X_old_1(16))]*Y(1:2);

        X_old_1(5) = psi_d_1(k);

        psi_q_1(k)  = [0 1]*[cos(X_old_1(16)), -sin(X_old_1(16));
            sin(X_old_1(16)),  cos(X_old_1(16))]*Y(1:2);

        X_old_1(6) = psi_q_1(k);


        psi_d_2(k)   = [1 0]*[cos(X_old_2(16)), -sin(X_old_2(16));
            sin(X_old_2(16)),  cos(X_old_2(16))]*Y(3:4);

        X_old_2(5) = psi_d_2(k);

        psi_q_2(k)   = [0 1]*[cos(X_old_2(16)), -sin(X_old_2(16));
            sin(X_old_2(16)),  cos(X_old_2(16))]*Y(3:4);

        X_old_2(6) = psi_q_2(k);



        I_d_1(k)     = [1 0]*[cos(X_old_1(16)), -sin(X_old_1(16));
            sin(X_old_1(16)),  cos(X_old_1(16))]*Y(5:6);

        X_old_1(3) = I_d_1(k);

        I_q_1(k)     = [0 1]*[cos(X_old_1(16)), -sin(X_old_1(16));
            sin(X_old_1(16)),  cos(X_old_1(16))]*Y(5:6);

        X_old_1(4) = I_q_1(k);


        I_d_2(k)     = [1 0]*[cos(X_old_2(16)), -sin(X_old_2(16));
            sin(X_old_2(16)),  cos(X_old_2(16))]*Y(7:8);

        X_old_2(3) = I_d_2(k);

        I_q_2(k)     = [0 1]*[cos(X_old_2(16)), -sin(X_old_2(16));
            sin(X_old_2(16)),  cos(X_old_2(16))]*Y(7:8);

        X_old_2(4) = I_q_2(k);



        V_d_1(k)     = [1 0]*[cos(X_old_1(16)), -sin(X_old_1(16));
            sin(X_old_1(16)),  cos(X_old_1(16))]*Y(9:10);

        X_old_1(1) = V_d_1(k);

        V_q_1(k)     = [0 1]*[cos(X_old_1(16)), -sin(X_old_1(16));
            sin(X_old_1(16)),  cos(X_old_1(16))]*Y(9:10);

        X_old_1(2) = V_q_1(k);

        V_rms_1(k) = sqrt((V_d_1(k))^2 + (V_q_1(k))^2);

        V_d_2(k)     = [1 0]*[cos(X_old_2(16)), -sin(X_old_2(16));
            sin(X_old_2(16)),  cos(X_old_2(16))]*Y(11:12);

        X_old_2(1) = V_d_2(k);

        V_q_2(k)     = [0 1]*[cos(X_old_2(16)), -sin(X_old_2(16));
            sin(X_old_2(16)),  cos(X_old_2(16))]*Y(11:12);

        X_old_2(2) = V_q_2(k);

        V_rms_2(k) = sqrt((V_d_2(k))^2 + (V_q_2(k))^2);





        % converting to abc frame

        Ck = sqrt(2/3)*[cos(W_0*t)           , sin(W_0*t)           ,   sqrt(1/2);
            cos(W_0*t - (2*pi/3)), sin(W_0*t - (2*pi/3)),   sqrt(1/2);
            cos(W_0*t + (2*pi/3)), sin(W_0*t + (2*pi/3)),   sqrt(1/2)];

        V_dqo_1 = [V_d_1(k), V_q_1(k), 0];
        V_dqo_2 = [V_d_2(k), V_q_2(k), 0];

        I_dqo_1 = [I_d_1(k), I_q_1(k), 0];
        I_dqo_2 = [I_d_2(k), I_q_2(k), 0];

        psi_dqo_1 = [psi_d_1(k), psi_q_1(k), 0];
        psi_dqo_2 = [psi_d_2(k), psi_q_2(k), 0];

        V_abc_1 = Ck*V_dqo_1';
        V_abc_2 = Ck*V_dqo_2';

        I_abc_1 = Ck*I_dqo_1';
        I_abc_2 = Ck*I_dqo_2';

        psi_abc_1 = Ck*psi_dqo_1';
        psi_abc_2 = Ck*psi_dqo_2';

        V_a_1(k) = V_abc_1(1);
        V_b_1(k) = V_abc_1(2);
        V_c_1(k) = V_abc_1(3);

        I_a_1(k) = I_abc_1(1);
        I_b_1(k) = I_abc_1(2);
        I_c_1(k) = I_abc_1(3);

        psi_a_1(k) = psi_abc_1(1);
        psi_b_1(k) = psi_abc_1(2);
        psi_c_1(k) = psi_abc_1(3);

        V_a_2(k) = V_abc_2(1);
        V_b_2(k) = V_abc_2(2);
        V_c_2(k) = V_abc_2(3);

        I_a_2(k) = I_abc_2(1);
        I_b_2(k) = I_abc_2(2);
        I_c_2(k) = I_abc_2(3);

        psi_a_2(k) = psi_abc_2(1);
        psi_b_2(k) = psi_abc_2(2);
        psi_c_2(k) = psi_abc_2(3);

        P_line(k) = Pm_1(k) - ((V_rms_1(k)^2)/R_L_1);

        Pe_1(k) = (psi_d_1(k)*I_q_1(k)) - (psi_q_1(k)*I_d_1(k));
        Pe_2(k) = (psi_d_2(k)*I_q_2(k)) - (psi_q_2(k)*I_d_2(k));

        V_rms_1_3phase(k) = V_rms_1(k);
        V_rms_2_3phase(k) = V_rms_2(k);

        theta_1(k) = delta_1(k) + asin(V_d_1(k)/V_rms_1_3phase(k));
        theta_2(k) = delta_2(k) + asin(V_d_2(k)/V_rms_2_3phase(k));
    end

    for k=round((Tfault_remove/h))+1:round(Tend2/h)
        t = k*h;
        R_L_1 = R_L_1_new2;
        R_L_2 = R_L_2_new2;
        r_line = r_line_new2;
        x_line = x_line_new2;


        % differential equations

        % equation (1)
        psi_H_1(k) = X_old_1(8) + h*(1/Td_dd)*(-X_old_1(8) + X_old_1(5));
        X_old_1(8) = psi_H_1(k);

        psi_H_2(k) = X_old_2(8) + h*(1/Td_dd)*(-X_old_2(8) + X_old_2(5));
        X_old_2(8) = psi_H_2(k);

        % equation (2)
        psi_F_1(k) = X_old_1(9) + h*(1/Td_d)*(-X_old_1(9) + X_old_1(5) + ((xd_d/(xd - xd_d))*X_old_1(7)));
        X_old_1(9) = psi_F_1(k);

        psi_F_2(k) = X_old_2(9) + h*(1/Td_d)*(-X_old_2(9) + X_old_2(5) + ((xd_d/(xd - xd_d))*X_old_2(7)));
        X_old_2(9) = psi_F_2(k);

        % equation (3)
        psi_K_1(k) = X_old_1(10) + h*(1/Tq_dd)*(-X_old_1(10) + X_old_1(6));
        X_old_1(10) = psi_K_1(k);

        psi_K_2(k) = X_old_2(10) + h*(1/Tq_dd)*(-X_old_2(10) + X_old_2(6));
        X_old_2(10) = psi_K_2(k);

        % equation (4)
        psi_G_1(k) = X_old_1(11) + h*(1/Tq_d)*(-X_old_1(11) + X_old_1(6));
        X_old_1(11) = psi_G_1(k);

        psi_G_2(k) = X_old_2(11) + h*(1/Tq_d)*(-X_old_2(11) + X_old_2(6));
        X_old_2(11) = psi_G_2(k);



        % Mechanical Equations

        delta_1(k) = X_old_1(16) + h*(X_old_1(12) - W_0);    % equation (9)
        X_old_1(16) = delta_1(k);

        omega_1(k) = X_old_1(12) + h*(W_b/(2*H))*(X_old_1(13) - ((X_old_1(5)*X_old_1(4)) - (X_old_1(6)*X_old_1(3))));  % equation (10)
        X_old_1(12) = omega_1(k);

        delta_2(k) = X_old_2(16) + h*(X_old_2(12) - W_0);    % equation (9)
        X_old_2(16) = delta_2(k);

        omega_2(k) = X_old_2(12) + h*(W_b/(2*H))*(X_old_2(13) - ((X_old_2(5)*X_old_2(4)) - (X_old_2(6)*X_old_2(3))));  % equation (10)
        X_old_2(12) = omega_2(k);




        % Exciter

        % Machine 1
        Xe_1(k)  = X_old_1(14) + h*(1/Ta)*(-X_old_1(14) + ka*(X_old_1(15) - V_rms_1(k-1)));

        if Xe_1(k) < Efd_min
            E_fd_1(k) = Efd_min;

        elseif Xe_1(k) > Efd_max
            E_fd_1(k) = Efd_max;

        elseif Xe_1(k) <= Efd_max && Xe_1(k) >= Efd_min
            E_fd_1(k) = Xe_1(k);

        end
        X_old_1(14) = Xe_1(k);
        X_old_1(7) = E_fd_1(k);

        % Machine 2
        Xe_2(k)  = X_old_2(14) + h*(1/Ta)*(-X_old_2(14) + ka*(X_old_2(15) - V_rms_2(k-1)));

        if Xe_2(k) < Efd_min
            E_fd_2(k) = Efd_min;

        elseif Xe_2(k) > Efd_max
            E_fd_2(k) = Efd_max;

        elseif Xe_2(k) <= Efd_max && Xe_2(k) >= Efd_min
            E_fd_2(k) = Xe_2(k);

        end
        X_old_2(14) = Xe_2(k);
        X_old_2(7) = E_fd_2(k);





        % Speed Governer
        if governer_1_enabled == 1
            Xg_1(k) = X_old_1(17) + h*(1/Tg2)*(-X_old_1(17) + ((W_ref - X_old_1(12))/W_b));
            Pm_1(k) = InitialCondition_1(13) + Kg*((1 - (Tg1/Tg2))*Xg_1(k) + ((Tg1/Tg2)*((W_ref - X_old_1(12))/W_b)));

            if Pm_1(k) < Pm_min
                Pm_1(k) = Pm_min;

            elseif Pm_1(k) > Pm_max
                Pm_1(k) = Pm_max;

            elseif Pm_1(k) <= Pm_max && Pm_1(k) >= Pm_min
                Pm_1(k) = Pm_1(k);

            end

        elseif governer_1_enabled == 0
            Pm_1(k) = InitialCondition_1(13);
            Xg_1(k) = 0;

        end
        X_old_1(17) = Xg_1(k);
        X_old_1(13) = Pm_1(k);

        if governer_2_enabled ==1
            Xg_2(k) = X_old_2(17) + h*(1/Tg2)*(-X_old_2(17) + ((W_ref - X_old_2(12))/W_b));
            Pm_2(k) = InitialCondition_2(13) + Kg*((1 - (Tg1/Tg2))*Xg_2(k) + ((Tg1/Tg2)*((W_ref - X_old_2(12))/W_b)));

            if Pm_2(k) < Pm_min
                Pm_2(k) = Pm_min;

            elseif Pm_2(k) > Pm_max
                Pm_2(k) = Pm_max;

            elseif Pm_2(k) <= Pm_max && Pm_2(k) >= Pm_min
                Pm_2(k) = Pm_2(k);

            end

        elseif governer_2_enabled == 0
            Pm_2(k) = InitialCondition_2(13);
            Xg_2(k) = 0;

        end
        X_old_2(17) = Xg_2(k);
        X_old_2(13) = Pm_2(k);


        % algebraic equations

        C_KG_delta_1 = (((xq_d - xq_dd)/xq_d)*X_old_1(10) + ((xq - xq_d)/xq)*(xq_dd/xq_d)*X_old_1(11))*cos(X_old_1(16));
        S_HF_delta_1 = (((xd_d - xd_dd)/xd_d)*X_old_1(8) + ((xd - xd_d)/xd)*(xd_dd/xd_d)*X_old_1(9))*sin(X_old_1(16));

        S_KG_delta_1 = (((xq_d - xq_dd)/xq_d)*X_old_1(10) + ((xq - xq_d)/xq)*(xq_dd/xq_d)*X_old_1(11))*sin(X_old_1(16));
        C_HF_delta_1 = (((xd_d - xd_dd)/xd_d)*X_old_1(8) + ((xd - xd_d)/xd)*(xd_dd/xd_d)*X_old_1(9))*cos(X_old_1(16));

        C_KG_delta_2 = (((xq_d - xq_dd)/xq_d)*X_old_2(10) + ((xq - xq_d)/xq)*(xq_dd/xq_d)*X_old_2(11))*cos(X_old_2(16));
        S_HF_delta_2 = (((xd_d - xd_dd)/xd_d)*X_old_2(8) + ((xd - xd_d)/xd)*(xd_dd/xd_d)*X_old_2(9))*sin(X_old_2(16));

        S_KG_delta_2 = (((xq_d - xq_dd)/xq_d)*X_old_2(10) + ((xq - xq_d)/xq)*(xq_dd/xq_d)*X_old_2(11))*sin(X_old_2(16));
        C_HF_delta_2 = (((xd_d - xd_dd)/xd_d)*X_old_2(8) + ((xd - xd_d)/xd)*(xd_dd/xd_d)*X_old_2(9))*cos(X_old_2(16));



        b = [C_KG_delta_1 - S_HF_delta_1;
            S_KG_delta_1 + C_HF_delta_1;
            C_KG_delta_2 - S_HF_delta_2;
            S_KG_delta_2 + C_HF_delta_2;
            0                      ;
            0                      ;
            0                      ;
            0                      ;
            0                      ;
            0                      ;
            0                      ;
            0                      ;
            0                      ;
            0                      ;];

        A = [0     , 1     , 0     , 0     , 0      , -xq_dd    , 0         , 0         , 0             , 0             , 0             , 0             , 0                 , 0                 ;
            1      , 0     , 0     , 0     , -xd_dd , 0         , 0         , 0         , 0             , 0             , 0             , 0             , 0                 , 0                 ;
            0      , 0     , 0     , 1     , 0      , 0         , 0         , -xq_dd    , 0             , 0             , 0             , 0             , 0                 , 0                 ;
            0      , 0     , 1     , 0     , 0      , 0         , -xd_dd    , 0         , 0             , 0             , 0             , 0             , 0                 , 0                 ;
            0      , 1     , 0     , 0     , Ra     , 0         , 0         , 0         , 1             , 0             , 0             , 0             , 0                 , 0                 ;
            -1     , 0     , 0     , 0     , 0      , Ra        , 0         , 0         , 0             , 1             , 0             , 0             , 0                 , 0                 ;
            0      , 0     , 0     , 1     , 0      , 0         , Ra        , 0         , 0             , 0             , 1             , 0             , 0                 , 0                 ;
            0      , 0     , -1    , 0     , 0      , 0         , 0         , Ra        , 0             , 0             , 0             , 1             , 0                 , 0                 ;
            0      , 0     , 0     , 0     , 0      , -R_L_1    , 0         , 0         , 0             , 1             , 0             , 0             , 0                 , R_L_1             ;
            0      , 0     , 0     , 0     , -R_L_1 , 0         , 0         , 0         , 1             , 0             , 0             , 0             , R_L_1             , 0                 ;
            0      , 0     , 0     , 0     , 0      , 0         , 0         , -R_L_2    , 0             , 0             , 0             , 1             , 0                 , -R_L_2            ;
            0      , 0     , 0     , 0     , 0      , 0         , -R_L_2    , 0         , 0             , 0             , 1             , 0             , -R_L_2            , 0                 ;
            0      , 0     , 0     , 0     , 0      , 0         , 0         , 0         , W_b/x_line    , 0             , -W_b/x_line   , 0             , -r_line*W_b/x_line, -W_b              ;
            0      , 0     , 0     , 0     , 0      , 0         , 0         , 0         , 0             , W_b/x_line    , 0             , -W_b/x_line   , W_b               , -r_line*W_b/x_line    ];

        Y = A\b;

        % bringing it back to machine reference

        psi_d_1(k)   = [1 0]*[cos(X_old_1(16)), -sin(X_old_1(16));
            sin(X_old_1(16)),  cos(X_old_1(16))]*Y(1:2);

        X_old_1(5) = psi_d_1(k);

        psi_q_1(k)  = [0 1]*[cos(X_old_1(16)), -sin(X_old_1(16));
            sin(X_old_1(16)),  cos(X_old_1(16))]*Y(1:2);

        X_old_1(6) = psi_q_1(k);


        psi_d_2(k)   = [1 0]*[cos(X_old_2(16)), -sin(X_old_2(16));
            sin(X_old_2(16)),  cos(X_old_2(16))]*Y(3:4);

        X_old_2(5) = psi_d_2(k);

        psi_q_2(k)   = [0 1]*[cos(X_old_2(16)), -sin(X_old_2(16));
            sin(X_old_2(16)),  cos(X_old_2(16))]*Y(3:4);

        X_old_2(6) = psi_q_2(k);



        I_d_1(k)     = [1 0]*[cos(X_old_1(16)), -sin(X_old_1(16));
            sin(X_old_1(16)),  cos(X_old_1(16))]*Y(5:6);

        X_old_1(3) = I_d_1(k);

        I_q_1(k)     = [0 1]*[cos(X_old_1(16)), -sin(X_old_1(16));
            sin(X_old_1(16)),  cos(X_old_1(16))]*Y(5:6);

        X_old_1(4) = I_q_1(k);


        I_d_2(k)     = [1 0]*[cos(X_old_2(16)), -sin(X_old_2(16));
            sin(X_old_2(16)),  cos(X_old_2(16))]*Y(7:8);

        X_old_2(3) = I_d_2(k);

        I_q_2(k)     = [0 1]*[cos(X_old_2(16)), -sin(X_old_2(16));
            sin(X_old_2(16)),  cos(X_old_2(16))]*Y(7:8);

        X_old_2(4) = I_q_2(k);



        V_d_1(k)     = [1 0]*[cos(X_old_1(16)), -sin(X_old_1(16));
            sin(X_old_1(16)),  cos(X_old_1(16))]*Y(9:10);

        X_old_1(1) = V_d_1(k);

        V_q_1(k)     = [0 1]*[cos(X_old_1(16)), -sin(X_old_1(16));
            sin(X_old_1(16)),  cos(X_old_1(16))]*Y(9:10);

        X_old_1(2) = V_q_1(k);

        V_rms_1(k) = sqrt((V_d_1(k))^2 + (V_q_1(k))^2);

        V_d_2(k)     = [1 0]*[cos(X_old_2(16)), -sin(X_old_2(16));
            sin(X_old_2(16)),  cos(X_old_2(16))]*Y(11:12);

        X_old_2(1) = V_d_2(k);

        V_q_2(k)     = [0 1]*[cos(X_old_2(16)), -sin(X_old_2(16));
            sin(X_old_2(16)),  cos(X_old_2(16))]*Y(11:12);

        X_old_2(2) = V_q_2(k);

        V_rms_2(k) = sqrt((V_d_2(k))^2 + (V_q_2(k))^2);





        % converting to abc frame

        Ck = sqrt(2/3)*[cos(W_0*t)           , sin(W_0*t)           ,   sqrt(1/2);
            cos(W_0*t - (2*pi/3)), sin(W_0*t - (2*pi/3)),   sqrt(1/2);
            cos(W_0*t + (2*pi/3)), sin(W_0*t + (2*pi/3)),   sqrt(1/2)];

        V_dqo_1 = [V_d_1(k), V_q_1(k), 0];
        V_dqo_2 = [V_d_2(k), V_q_2(k), 0];

        I_dqo_1 = [I_d_1(k), I_q_1(k), 0];
        I_dqo_2 = [I_d_2(k), I_q_2(k), 0];

        psi_dqo_1 = [psi_d_1(k), psi_q_1(k), 0];
        psi_dqo_2 = [psi_d_2(k), psi_q_2(k), 0];

        V_abc_1 = Ck*V_dqo_1';
        V_abc_2 = Ck*V_dqo_2';

        I_abc_1 = Ck*I_dqo_1';
        I_abc_2 = Ck*I_dqo_2';

        psi_abc_1 = Ck*psi_dqo_1';
        psi_abc_2 = Ck*psi_dqo_2';

        V_a_1(k) = V_abc_1(1);
        V_b_1(k) = V_abc_1(2);
        V_c_1(k) = V_abc_1(3);

        I_a_1(k) = I_abc_1(1);
        I_b_1(k) = I_abc_1(2);
        I_c_1(k) = I_abc_1(3);

        psi_a_1(k) = psi_abc_1(1);
        psi_b_1(k) = psi_abc_1(2);
        psi_c_1(k) = psi_abc_1(3);

        V_a_2(k) = V_abc_2(1);
        V_b_2(k) = V_abc_2(2);
        V_c_2(k) = V_abc_2(3);

        I_a_2(k) = I_abc_2(1);
        I_b_2(k) = I_abc_2(2);
        I_c_2(k) = I_abc_2(3);

        psi_a_2(k) = psi_abc_2(1);
        psi_b_2(k) = psi_abc_2(2);
        psi_c_2(k) = psi_abc_2(3);

        P_line(k) = Pm_1(k) - ((V_rms_1(k)^2)/R_L_1);

        Pe_1(k) = (psi_d_1(k)*I_q_1(k)) - (psi_q_1(k)*I_d_1(k));
        Pe_2(k) = (psi_d_2(k)*I_q_2(k)) - (psi_q_2(k)*I_d_2(k));

        V_rms_1_3phase(k) = V_rms_1(k);
        V_rms_2_3phase(k) = V_rms_2(k);

        theta_1(k) = delta_1(k) + asin(V_d_1(k)/V_rms_1_3phase(k));
        theta_2(k) = delta_2(k) + asin(V_d_2(k)/V_rms_2_3phase(k));
    end
end 





%-----------------------------------------------------------------------------------------------------







%%

% ________________________________________________
% PLOTS
% ________________________________________________

% psi_H_1                       psi_H_2
% psi_F_1                       psi_F_2
% psi_K_1                       psi_K_2
% psi_G_1                       psi_G_2
% PSI_a1, PSI_b1, PSI_c1        PSI_a2, PSI_b2, PSI_c2            
% PSId1, PSIq1                  PSId2, PSIq2
% I_a1, I_b1, I_c1              I_a2, I_b2, I_c2
% Id1, Iq1                      Id2, Iq2
% V_a1, V_b1, V_c1              V_a2, V_b2, V_c2
% Vd1, Vq1                      Vd2, Vq2
% delta_1                       delta_2
% omega_1                       omega_2
% P_M1, XG_1                    P_M2, XG_2
% Efd1, XE_1, V_rms1            Efd2, XE_2, V_rms2


% ________________________________________________
% PLOTS
% ________________________________________________

% Define time vector
t = 0: h: Tend2-h;

% Create main figure and panel for buttons
f = figure;
p = uipanel("Title", "Plots", "Position", [.01, .05, .08, .9]);

% Define push buttons for each plot type
b1 = uicontrol("Parent", p, "Style", "pushbutton", "String", "Ia", "Position", [15, 640, 120, 30], "Callback", @(src, event) subplot2(t, I_a_1, I_a_2, 'a - phase current', 'I_a (machine 1)', 'I_a (machine 2)'));
b2 = uicontrol("Parent", p, "Style", "pushbutton", "String", "Idq", "Position", [15, 600, 120, 30], "Callback", @(src, event) subplot4(t, I_d_1, I_q_1, I_d_2, I_q_2, 'Currents (d-q)', 'I_d (machine 1)', 'I_q (machine 1)', 'I_d (machine 2)', 'I_q (machine 2)'));
b3 = uicontrol("Parent", p, "Style", "pushbutton", "String", "Va", "Position", [15, 560, 120, 30], "Callback", @(src, event) subplot2(t, V_a_1, V_a_2, 'a-phase voltage', 'V_a (machine 1)', 'V_a (machine 2)'));
b4 = uicontrol("Parent", p, "Style", "pushbutton", "String", "Vdq", "Position", [15, 520, 120, 30], "Callback", @(src, event) subplot4(t, V_d_1, V_q_1, V_d_2, V_q_2, 'Voltage (d-q)', 'V_d (machine 1)', 'V_q (machine 1)', 'V_d (machine 2)', 'V_q (machine 2)'));
b5 = uicontrol("Parent", p, "Style", "pushbutton", "String", "delta", "Position", [15, 480, 120, 30], "Callback", @(src, event) subplot2(t, delta_1, delta_2, 'Delta (in rad)', 'Delta (machine 1)', 'Delta (machine 2)'));
b6 = uicontrol("Parent", p, "Style", "pushbutton", "String", "W (all)", "Position", [15, 440, 120, 30], "Callback", @(src, event) plot3(t, omega_1./(2*pi), omega_2./(2*pi), ((omega_1 + omega_2)/2)./(2*pi), "Angular Speeds \omega (in rad/s)", "\omega (in rad/s)", "\omega 1", "\omega 2", "\omega COI"));
b7 = uicontrol("Parent", p, "Style", "pushbutton", "String", "Power (all)", "Position", [15, 400, 120, 30], "Callback", @(src, event) plot3(t, Pm_1, Pm_2, P_line, "Power Transients", "Power (pu)", "Pm 1", "Pm 2", "P line"));
b8 = uicontrol("Parent", p, "Style", "pushbutton", "String", "Efd", "Position", [15, 360, 120, 30], "Callback", @(src, event) subplot4(t, E_fd_1, Xe_1, E_fd_2, Xe_2, 'Field Voltage (pu)', 'Efd (machine 1)', 'State (Xe) (machine 1)', 'Efd (machine 2)', 'State (Xe) (machine 2)'));
b9 = uicontrol("Parent", p, "Style", "pushbutton", "String", "mech_diff", "Position", [15, 320, 120, 30], "Callback", @(src, event) subplot2(t, (delta_1 - delta_2)*(180/pi), omega_1 - omega_2 , 'Relative Mechanical Quantities','Rotor Angle Difference ''\delta 1 - \delta 2 (in degrees)', 'Diff in Rotor Angular Speeds ''\omega 1 - \omega 2 (in rad/s)'));
b10= uicontrol("Parent", p, "Style", "pushbutton", "String", "Power(E & M)", "Position", [15, 280, 120, 30], "Callback", @(src, event) subplot4(t, Pm_1, Pe_1, Pm_2, Pe_2, 'Mech & Elec Power (pu)', 'Mech Power Machine 1', 'Elec Power Machine 1', 'Mech Power Machine 2', 'Elec Power Machine 2'));
b11= uicontrol("Parent", p, "Style", "pushbutton", "String", "Bus (V & theta)", "Position", [15, 240, 120, 30], "Callback", @(src, event) subplot4(t, V_rms_1_3phase, V_rms_2_3phase, theta_1.*(180/pi), theta_2.*(180/pi), 'Bus Voltage and Angle', 'Voltage Mag Bus 1', 'Voltage Mag Bus 2', 'Voltage Angle Bus 1', 'Voltage Angle Bus 2'));
b12= uicontrol("Parent", p, "Style", "pushbutton", "String", "bus_diff", "Position", [15, 200, 120, 30], "Callback", @(src, event) subplot2(t, (V_rms_1_3phase - V_rms_2_3phase), (theta_2 - theta_1)*(180/pi) , 'Bus Voltage and Angle Diff','Bus Voltage Difference (V1 - V2) (pu)', 'Bus Angle Diff ''\theta 2 - \theta 1 (in degrees)'));
b13 = uicontrol("Parent", p, "Style", "pushbutton", "String", "Power Elec(all)", "Position", [15, 160, 120, 30], "Callback", @(src, event) plot3(t, Pe_1, Pe_2, P_line, "Power Transients", "Elec Power (pu)", "Pe 1", "Pe 2", "P line"));
% Function for one plot
function [] = plot1(time, data1, titleText, yLabel)
    figure;
    
    plot(time, data1, 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel(yLabel);
    title(titleText);
    grid on
end

% Function for two vertically stacked plots
function [] = subplot2(time, data1, data2, titleText, yLabel1, yLabel2)
    figure;
    tiledlayout(2, 1);

    % First plot
    nexttile;
    plot(time, data1, 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel(yLabel1);
    title(titleText);
    grid on;

    % Second plot
    nexttile;
    plot(time, data2, 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel(yLabel2);
    grid on;
end

% Function for three data series in same plot
function [] = plot3(time, data1, data2, data3, titleText, yLabel, yLabel1, yLabel2, yLabel3)
    figure;
    plot(time, data1, 'LineWidth', 2, 'Color', 'b', 'LineStyle', '-');
    hold on
    plot(time, data2, 'LineWidth', 2, 'Color', 'g', 'LineStyle', '-');
    plot(time, data3, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '-');
    xlabel('Time (s)');
    ylabel(yLabel);
    legend(yLabel1, yLabel2, yLabel3);
    title(titleText);
    grid on
    legend
end

% Function for four data series in two stacked plots
function [] = subplot4(time, data1, data2, data3, data4, titleText, yLabel1, yLabel2, yLabel3, yLabel4)
    figure;
    tiledlayout(2, 1);

    % Top plot with two lines
    nexttile;
    plot(time, data1, 'LineWidth', 2); hold on;
    plot(time, data2, 'LineWidth', 1.5); hold off;
    xlabel('Time (s)');
    ylabel([yLabel1, ', ', yLabel2]);
    legend(yLabel1, yLabel2);
    title(titleText);
    grid on;

    % Bottom plot with two lines
    nexttile;
    plot(time, data3, 'LineWidth', 2); hold on;
    plot(time, data4, 'LineWidth', 1.5); hold off;
    xlabel('Time (s)');
    ylabel([yLabel3, ', ', yLabel4]);
    legend(yLabel3, yLabel4);
    grid on;
end