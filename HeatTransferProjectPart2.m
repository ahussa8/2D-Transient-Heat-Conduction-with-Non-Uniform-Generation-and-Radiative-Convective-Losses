clc; clear; close all

L = 1.0;
dx = 1/(10-1); dy = dx;
k   = 2200;
rho = 3515;
cp  = 520;
a   = k/(rho*cp);
G_inc = 800;
G_ref = 120;
sigma = 5.670374419e-8;

T_L_fun      = @(t) 300 + 10*sin(2*pi*t/5);
T_T_fun      = @(t) 300 + 20*cos(2*pi*t/8);
T_R_fun      = @(t) 300;
q_bottom_fun = @(t) 200*(1 + sin(2*pi*t/6));

t_eval = 5;
T_L      = T_L_fun(t_eval);
T_T      = T_T_fun(t_eval);
T_R      = T_R_fun(t_eval);
q_bottom = q_bottom_fun(t_eval);
T_B      = T_L + dy*(q_bottom/k);
T_inf    = T_T;
u        = 5;
alpha_s  = (G_inc - G_ref)/G_inc;

Ts_range = T_T:20:1000;
n        = numel(Ts_range);

epsilon_vals = zeros(1,n);
h_vals       = zeros(1,n);
q_conv       = zeros(1,n);
q_rad        = zeros(1,n);
J_rad        = zeros(1,n);
q_net        = zeros(1,n);

fprintf('  Ts (°C)    ε     α_s    J_rad (W/m²)\n');
fprintf('--------------------------------------\n');
for i = 1:n
    Ts = Ts_range(i);
    eps_i = emissivity_vs_T(Ts);
    epsilon_vals(i) = eps_i;

    T_f   = 0.5*(Ts + T_inf);
    mu    = dynamic_viscosity(T_f);
    rho_a = air_density(T_f);
    cp_a  = specific_heat_capacity(T_f);
    k_air = thermal_conductivity(T_f);
    Pr    = mu*cp_a / k_air;

    Re_L = rho_a * u * L / mu;
    if Re_L < 5e5
        Nu = 0.664 * sqrt(Re_L) * Pr^(1/3);
    else
        Nu = 0.037 * Re_L^0.8 * Pr^(1/3);
    end

    h_vals(i) = Nu * k_air / L;
    q_rad(i)  = eps_i * sigma * (Ts^4 - T_inf^4);
    q_conv(i)= h_vals(i) * (Ts - T_inf);
    J_rad(i)  = eps_i * sigma * Ts^4 + (1 - eps_i) * G_inc;
    q_net(i)  = alpha_s * G_inc - q_rad(i) - q_conv(i);

    fprintf('%8.1f   %6.3f  %6.3f   %8.2f\n', Ts-273.15, eps_i, alpha_s, J_rad(i));
end
fprintf('\n');

T_C = Ts_range - 273.15;
Results = table( ...
    Ts_range', T_C', epsilon_vals', repmat(alpha_s,n,1), h_vals', ...
    q_conv', q_rad', J_rad', q_net', ...
    'VariableNames',{'T_s_K','T_s_C','epsilon','alpha_s','h','q_conv','q_rad','J_rad','q_net'});
disp(Results);

figure;
plot(Ts_range-273.15, q_net, 'LineWidth',2);
xlabel('Surface Temperature T_s (°C)');
ylabel('Net Heat Flux q_{net} (W/m²)');
title('Net Heat Flux vs Surface Temperature');
grid on;

figure;
plot(Ts_range-273.15, q_conv, 'r-',  'LineWidth',1.5); hold on;
plot(Ts_range-273.15, q_rad,  'b--', 'LineWidth',1.5);
xlabel('Surface Temperature T_s (°C)');
ylabel('Heat Flux (W/m²)');
title('Convective vs Radiative Losses');
legend('q_{conv}','q_{rad}','Location','northwest');
grid on; hold off;

h_pl = 6.626e-34; c_pl = 3e8; kB = 1.381e-23;
lam_um = logspace(-1,2,2000); lam_m = lam_um*1e-6;
T_bb   = [1500 2000 2500 3000 3500 4000 4500 5000];

figure; hold on;
for Tval = T_bb
    B = (2*h_pl*c_pl^2) ./ (lam_m.^5) ./ (exp(h_pl*c_pl./(lam_m*kB*Tval)) - 1);
    plot(lam_um, B, 'DisplayName', sprintf('T = %d K', Tval));
end
set(gca,'XScale','log');
xlabel('Wavelength (\mum)');
ylabel('Spectral Radiance (W/m²/sr/\mum)');
title('Blackbody Spectral Radiance');
legend('Location','northeast');
grid on; hold off;

eps_plot = [0.5 0.8];
T_gray   = 3000;

figure; hold on;
for epsv = eps_plot
    B_g = epsv * (2*h_pl*c_pl^2) ./ (lam_m.^5) ./ (exp(h_pl*c_pl./(lam_m*kB*T_gray)) - 1);
    plot(lam_um, B_g, 'DisplayName', sprintf('\\epsilon = %.1f', epsv));
end
xlabel('Wavelength (\mum)');
ylabel('Spectral Radiance (W/m²/sr/\mum)');
title('Graybody Spectral Radiance at 3000 K');
legend('Location','northeast');
grid on; hold off;

function eps = emissivity_vs_T(T)
    eps = 0.10 + 1e-4*(T - 300);
    eps = max(min(eps,0.95),0.05);
end

function mu = dynamic_viscosity(T)
    C1 = 1.458e-6; S = 110.4;
    mu = C1 * T^(3/2) / (T + S);
end

function rho = air_density(T)
    P = 101325; R = 287.05;
    rho = P / (R * T);
end

function cp = specific_heat_capacity(T)
    cp = 1005 + 0.1*(T - 300);
end

function k = thermal_conductivity(T)
    k = 0.0262 + 0.0001*(T - 300);
end
