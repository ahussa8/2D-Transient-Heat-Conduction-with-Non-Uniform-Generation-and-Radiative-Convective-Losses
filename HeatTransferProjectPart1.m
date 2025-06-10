close all
clear all
clc

% grid and time
nx = 10;  ny = nx;
nt = 1400;
x  = linspace(0,1,nx);
y  = linspace(0,1,ny);
dx = x(2)-x(1);
dy = dx;
dt = 1e-3;

% solver params
tolerance = 1e-1;
Gauss_Seidel_iteration = 1;

% material (diamond)
k   = 2200;
rho = 3515;
cp  = 520;
a   = k/(rho*cp);
alpha = 1.1*(dt/dx^2);

% initial field
T = 300*ones(nx,ny);
T_old    = T;
T_initial = T;

% define your time‐dependent boundary functions here:
T_L_fun      = @(t) 300 + 20*sin(2*pi*t/2);    % left edge oscillates ±20 K over 2 s 
T_T_fun      = @(t) 300 + 10*cos(2*pi*t/1.5);  % top edge oscillates ±10 K over 1.5 s
T_R_fun      = @(t) 300;                       % right edge held constant
q_bottom_fun = @(t) 500*(1+sin(2*pi*t/3));     % bottom flux oscillates between 0 and 1000 W/m²

% define source blocks (unchanged)
Q_val1 = 15000;
Q_val2 = 10000;
src1 = false(nx,ny);
src2 = false(nx,ny);
src1(1:5,1:5) = true;
src2(6:8,6:8) = true;
Qmask = src1*Q_val1 + src2*Q_val2;

figure(1)
for k = 1:nt
    t = k*dt;
    % update time‐dependent boundaries BEFORE each time step
    T_L      = T_L_fun(t);
    T_T      = T_T_fun(t);
    T_R      = T_R_fun(t);
    q_bottom = q_bottom_fun(t);
    T_B      = T_R + dy*(q_bottom/k);  % Neumann at bottom

    % reapply BCs to both T and T_old (so solver sees them)
    T(2:end-1,1)   = T_L;      % left
    T_old(2:end-1,1) = T_L;
    T(1,2:end-1)   = T_T;      % top
    T_old(1,2:end-1) = T_T;
    T(2:end-1,end) = T_R;      % right
    T_old(2:end-1,end) = T_R;
    T(end,2:end-1) = T_B;      % bottom
    T_old(end,2:end-1) = T_B;

    % corners (average or inherited)
    T(1,1)     = (T_L+T_T)/2; 
    T_old(1,1) = T(1,1);
    T(1,end)   = T_T;  T_old(1,end)   = T_T;
    T(end,1)   = T_L;  T_old(end,1)   = T_L;
    T(end,end) = T(end,end-1); 
    T_old(end,end) = T(end,end);

    % Gauss–Seidel loop
    error = inf;
    while error > tolerance
        i = 2:nx-1; j = 2:ny-1;
        T(i,j) = ( T_initial(i,j).*(1-4*alpha) + ...
                   alpha*( T(i-1,j) + T_old(i+1,j) + ...
                           T(i,j+1) + T_old(i,j-1) ) ...
                 ) + Qmask(i,j)*dt;
        error = max(abs(T_old(:)-T(:)));
        T_old = T;
        Gauss_Seidel_iteration = Gauss_Seidel_iteration + 1;
    end

    T_initial = T;

    % plot
    contourf(x,y,T,20,'LineColor','none')
    colorbar, colormap(jet)
    set(gca,'ydir','reverse')
    xlabel('X (m)'), ylabel('Y (m)')
    title(sprintf('t=%.3f s — GS iters: %d', t, Gauss_Seidel_iteration))
    drawnow
end
