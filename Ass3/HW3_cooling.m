%% Create a file which reads in CEA data and performs the tasks of HW3
%
% hola Adam
%
% 31/10/17  Adam, Max, Ivan, Alberto



clc
clear all
close all
format Long

% Known constants, given for problem
g0=9.80665;
OF = 0.98; % the O/F mixture ratio
Tw = 700; % Kelvin
T_cool = 298; %Kelvin
density_cool = 1010; %kg/m3
Cp_cool = 3000; % J/K/kg
mu_cool = 9.7e-4; %Pa/s
K_cool = 0.32; % W/K/m

%% READ AND PREPROCESS CEA OUTPUT FILE

% open the required data files and scan for required info into data arrays
data = fopen('HW3_single.plt');
C = textscan(data, '%f%f%f%f%f%f%f%f%f%f%f%f%f','commentStyle', '#');
fclose(data);

% Convert to matrices, easier to extract data
col_len = length(C{1});
clean_data = cell2mat(C);
n_cases = size(clean_data, 1)/3;
F = [10e3 1000e3];
n_thrust = length(F);

[Pc, Pc_Pt, Pc_Pe, M_t, M_e, Cf, Cf_t, Viscosity, Pr_n, Pr_fr, Cp, gamma, Isp, Isp_vac, Tc, Ae_At] = deal(zeros(n_cases,1));

% Extract only the data we need, and have the correct SI units
for n = 1:n_cases
    
    Pc(n) = clean_data(-2+3*n, 1)*10^5;
    Pc_Pt(n) = clean_data(-1+3*n, 2);
    Pc_Pe(n) = clean_data(3*n, 2);
    M_t(n) = clean_data(-1+3*n, 3);
    M_e(n) = clean_data(3*n, 3);
    Ae_At(n) = clean_data(3*n, 4);
    Cf(n) = clean_data(3*n, 5);     % Converted to Vacuum value below
    Viscosity(n) = clean_data(-2+3*n, 6)*10^-4; % convert from mPoise to Pa
    Pr(n) = clean_data(-2+3*n,7);
    Pr_fr(n) = clean_data(-2+3*n, 8);
    Cp(n) = clean_data(-2+3*n, 9)*10^3;  % Convert from KJ/Kg to KG             THIS WAS OUR PROBLEM
    gamma(n) = clean_data(-2+3*n, 10);
    Isp(n) = clean_data(3*n, 11)/g0;        % convert into Seconds
    Isp_vac(n) = clean_data(3*n, 12)/g0;
    Tc(n) = clean_data(-2+3*n, 13);     %Kelvin
    
end

M = M_t;
% Pr = Pr_fr;
h = zeros(n_cases, n_thrust);
q = zeros(n_cases, n_thrust);
At = zeros(n_cases, n_thrust);

%% CALCULATE Q FOR A SPECIFIC RANGE OF Pc AND Thrust
% Looking at the HOT GAS side first
for k = 1:n_thrust
    for i = 1:n_cases
        Taw(i) = Tc(i);
        sigma = 1/((0.5*Tw/Tc(i)*(1+(gamma(i)-1)/2*M(i)^2)+0.5)^0.68*(1+(gamma(i)-1)/2*M(i)^2)^0.12);
        Cf(i) = Cf(i) + 1/Pc_Pe(i)*Ae_At(i);          % Correction for an adapted nozzle to vacuum
        At(i,k) = F(1,k)/(Cf(i)*Pc(i));  %Area throat
        Dt(i,k) = 2*sqrt(At(i)/pi); %Diameter throat
        Dc(i,k) = 3*Dt(k);      % Assume chamber diameter is 3 times that of throat
        Ac(i,k) = (pi/4)*Dc(k)^2;  % Area of Chamber
        Lc(i,k) = 4*Dc(k);      % Assume cham length is 4 times the cham diameter
        Cstar(i) = (Isp(i)*g0)/Cf(i);
        Rc(i,k) = Dt(i,k);
        %h(i, k) = 0.026 * (Viscosity(i)^0.2 * Cp(i) / (Dt^0.2 * Pr(i)^0.6)) * (Pc(i)/Cstar)^0.8 * (At/A)^0.9 * (Dt/Rc)^0.1 * sigma;
        h(i, k) = 0.026 * (Viscosity(i)^0.2 * Cp(i) / (Dt(k)^0.2 * Pr(i)^0.6)) * (Pc(i)/Cstar)^0.8 * (At(k)/At(k))^0.9 * (Dt(k)/Rc(k))^0.1;
        q(i, k) = h(i, k) * (Taw(i) - Tw);
        
        h_hg(i,k) = 0.026*((Viscosity(i)^0.2*Cp(i))/(Dt(i)^0.2*Pr(i)^0.6))*(Pc(i)/Cstar(i))^0.8 *(Dt(i)/Rc(i))^0.1 * (At(i,k)/At(i,k))^0.9;
        q_hg(i,k) = h_hg(i,k)*(Taw(i) - Tw);
        
        M_dot_tot(i,k) = F(1,k)/(Isp(i));   % The total mass flow rate
        M_dot_cool(i,k) = M_dot_tot(i,k)/(1+OF);             % The propellant/coolant flow rate found with O/f ratio
        
    end
 end


%% Looking at the COOLANT side now
% need to input the N and coolant pipe diameter from another script



Pr_cool = (mu_cool*Cp_cool)/K_cool;
D_cool = 0.002;
h_c = 0.029*((Cp_cool*mu_cool^0.2)/Pr_cool^(2/3)) *(M_dot_cool/D_cool^0.2)*(T_cool/Tw)^0.55;  %% WONKEY AS FREAK


return






%% PLOT THE RESULT

C = [0.9047 0.1918 0.1988; 0.2941 0.5447 0.7494; 0.3718 0.7176 0.3612;...
    1.0000 0.5482 0.1000; 0.8650 0.8110 0.4330; 0.6859 0.4035 0.2412]; %Colors used in plots


% title('Heat flux on the throttle ( q_{t} ) vs. Chamber Pressure ( P_{C} ) for different Thrusts')
% ylabel('q_{t} [MW/m^{2}]')
% xlabel('P_{C} [bar]')
% legend('Thrust = 10^{-1}', 'Thrust = 10^{0}','Thrust = 10^{1}', 'Thrust = 10^{2}', 'Thrust = 10^{3}');
% xlim([5 250])
% ylim([1 1000])
% grid on


