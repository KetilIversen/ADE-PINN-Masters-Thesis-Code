%Clear everything
clc
clear all 
close all

% Initialize model domains
Lx = 2;     % Default 2      (Length of x domain)
Lt = 1;     % Default 1      (Length of t domain)
Nx = 256;   % Default 256    (Gridpoints in x)
Nt = 2000;  % Default 10000  (Computational timesteps)
Ntp = 100;  % Default 100    (Saved timesteps)

% Initialize model step sizes and grid
dx = Lx/(Nx-1);
dt = Lt/(Nt-1);
x = linspace(-Lx/2,Lx/2,Nx);
t = linspace(0,Lt,Nt);
timesteps = linspace(0,Lt,Ntp);

%% Initialize PDE parameters
V_x = 1;       % Default 1       (Velocity x)
D = 0.25/pi;   % Default 0.25/pi (Diffusion rate)

%% Source function
f(Nx,1)  =  0;
lambda_1 = -0.8; % Default -0.8 (Source location x)
lambda_2 =  5.0; % Default  5.0 (Source stength)
S =  100;        % Default  100 (Source spread)
for i=2:Nx-1
    f(i,1) = lambda_2*exp(-S*(((x(i) - lambda_1)).^2));
end

%% Initialize solution
u = zeros(Nx,Ntp);
u_new = sparse(u(:,1));
u_old = sparse(u(:,1));

%% Iterate w/ RK4
tic
E = SparseE(Nx);
e = ones(Nx,1);
snapstep = ceil(Nt/Ntp);
snap = 1;
for n=1:Nt-1
    k1 = C_Derivative(u_old,         dx,D,V_x,f,Nx,E,e);
    k2 = C_Derivative(u_old+k1*dt/2, dx,D,V_x,f,Nx,E,e);
    k3 = C_Derivative(u_old+k2*dt/2, dx,D,V_x,f,Nx,E,e);
    k4 = C_Derivative(u_old+k3*dt,   dx,D,V_x,f,Nx,E,e);
    u_new = u_old + (1/6)*dt*(k1+2*k2+2*k3+k4);
    u_old = u_new;
    
    if mod(n,snapstep)==0
        u(:,snap+1) = u_new;
        snap = snap + 1;
    end
    if mod(n,10)==0
        disp(['Finished with iteration: ',num2str(n)])
    end
end
toc
%% Plot solution
DiscretePlotter(x,u,Ntp,Lt);

%% Functions
% RK4
function k = C_Derivative(u_vec,dx,D,V_x,f_vec,Nx,E,e)  
    
    alpha_x = D/(dx^2) - V_x/(2*dx);
    beta    = - (2*D/(dx^2));
    gamma_x = D/(dx^2) + V_x/(2*dx);
    
    row = [gamma_x*e, beta*e, alpha_x*e];
    d = [-1, 0, 1];
    A = spdiags(row,d,Nx,Nx);
    A = A.*E;
    
    k = A*u_vec + f_vec;
    
end
% Sparse E
function E = SparseE(Nx)
    E = zeros(Nx,Nx);
    for i = 2:Nx-1
       E(i,i-1) = 1;
       E(i,i)   = 1;
       E(i,i+1) = 1;
    end
    E = sparse(E);
end
% Plotter
function DiscretePlotter(x,u,Ntp,Lt)
    
    figure('DefaultAxesFontSize',18)
    x_width=1200 ;y_width=600;
    set(gcf, 'Position', [0 0 x_width y_width]);
    set(gcf, 'defaultAxesTickLabelInterpreter','latex')
    set(gcf, 'defaulttextinterpreter','latex')
    set(gcf, 'defaultLegendInterpreter','latex')
    
    t0 = 0; t1 = 1*Lt/5; t2 = 2*Lt/5; t3 = 3*Lt/5; t4 = 4*Lt/5; t5 = 5*Lt/5;
    N0 = ceil(t0*Ntp + 1); N1 = ceil(t1*Ntp/Lt); N2 = ceil(t2*Ntp/Lt);
    N3 = ceil(t3*Ntp/Lt); N4 = ceil(t4*Ntp/Lt); N5 = floor(t5*Ntp/Lt);

    subplot(2,3,1)
    plot(x,u(:,N0),'LineWidth',3)
    legend('$u(x,t)$')
    xlabel('$x$')
    ylabel('$u$')
    ylim([-0.5 1])
    title('$t$ = ' + string(t0))
    grid on

    subplot(2,3,2)
    plot(x,u(:,N1),'LineWidth',3)
    xlabel('$x$')
    ylabel('$u$')
    ylim([-0.5 1])
    title('$t$ = ' + string(t1))
    grid on

    subplot(2,3,3)
    plot(x,u(:,N2),'LineWidth',3)
    xlabel('$x$')
    ylabel('$u$')
    ylim([-0.5 1])
    title('$t$ = ' + string(t2))
    grid on

    subplot(2,3,4)
    plot(x,u(:,N3),'LineWidth',3)
    xlabel('$x$')
    ylabel('$u$')
    ylim([-0.5 1])
    title('$t$ = ' + string(t3))
    grid on

    subplot(2,3,5)
    plot(x,u(:,N4),'LineWidth',3)
    xlabel('$x$')
    ylabel('$u$')
    ylim([-0.5 1])
    title('$t$ = ' + string(t4))
    grid on

    subplot(2,3,6)
    plot(x,u(:,N5),'LineWidth',3)
    xlabel('$x$')
    ylabel('$u$')
    ylim([-0.5 1])
    title('$t$ = ' + string(t5))
    grid on
end
