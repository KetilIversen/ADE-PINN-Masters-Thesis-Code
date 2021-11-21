close all
clear all

% Initialize model domains
Lx = 4;     % Default 4     (Length of x domain)
Ly = 4;     % Default 4     (Length of y domain)
Lt = 2;     % Default 2     (Length of t domain)
Nx = 100;   % Default 100   (Gridpoints in x)
Ny = 100;   % Default 100   (Gridpoints in y)
Nt = 20000; % Default 20000 (Computational timesteps)
Ntp = 2000; % Default 2000  (Saved timesteps)

% Initialize model step sizes and grid
dx = Lx/(Nx-1);
dy = Ly/(Ny-1);
dt = Lt/(Nt-1);
x = linspace(-Lx/2,Lx/2,Nx);
y = linspace(-Ly/2,Ly/2,Ny);
t = linspace(0,Lt,Nt);
timesteps = linspace(0,Lt,Ntp);
[X,Y] = meshgrid(x,y);

%% Initialize PDE parameters
D = 0.2;                            % Default 0.2 (Diffusion rate)
V_x = @(t) 2.0;                     % Default 2.0 (Velocity x)           
V_y = @(t) 3.0*cos(3.0*pi*t) + 1.0; % Default 3.0*cos(3.0*pi*t) + 1.0 (Velocity y)

% Initialize solution
u(Nx,Ny,Ntp) = 0;
u_new = u(:,:,1);
u_old = u(:,:,1);

%% Source function
f(Nx,Ny) = 0;
lambda_1 = -1.4;    % Default -1.4 (Source location x)
lambda_2 = -1.0;    % Default -1.0 (Source location y)
lambda_3 =  1.0;    % Default  1.0 (Source strength)
S = 25.0;            % Default  4.0  (Source spread)
for i=2:Nx-1
    for j=2:Ny-1
        f(i,j) = lambda_3*exp(-S*(((x(i) - lambda_1)).^2+(y(j) - lambda_2).^2));
    end
end

%% Flatten solution
u_vec = reshape(u,Nx*Ny,Ntp);
u_new = sparse(reshape(u_new,Nx*Ny,1));
u_old = sparse(reshape(u_old,Nx*Ny,1));
f_vec = sparse(reshape(f,Nx*Ny,1));

%% Iterate w/ RK4
tic
E = SparseE(Nx,Ny);
e = ones(Nx*Ny,1);
snapstep = ceil(Nt/Ntp);
snap = 1;
for n=1:Nt-1    
    k1 = C_Derivative(t(n)       ,u_old          ,dx,dy,D,V_x,V_y,f_vec,Nx,Ny,E,e);
    k2 = C_Derivative(t(n)+ dt/2 ,u_old+k1*dt/2  ,dx,dy,D,V_x,V_y,f_vec,Nx,Ny,E,e);
    k3 = C_Derivative(t(n)+ dt/2 ,u_old+k2*dt/2  ,dx,dy,D,V_x,V_y,f_vec,Nx,Ny,E,e);
    k4 = C_Derivative(t(n)+ dt   ,u_old+k3*dt    ,dx,dy,D,V_x,V_y,f_vec,Nx,Ny,E,e);
    u_new = u_old + (1/6)*dt*(k1+2*k2+2*k3+k4);
    u_old = u_new;
    
    if mod(n,snapstep)==0
        u_vec(:,snap+1) = u_new;
        snap = snap + 1;
    end
    if mod(n,10)==0
        disp(['Finished with iteration: ',num2str(n)])
    end
end

%% Unflatten solution
u = reshape(u_vec,Nx,Ny,Ntp);
toc
%% Plot solution
DiscretePlotter(X,Y,u,0,0.25,timesteps,Ntp,Lt)

%% Functions
% RK4
function k = C_Derivative(t,u_vec,dx,dy,D,V_x,V_y,f_vec,Nx,Ny,E,e)  
    
    alpha_y = D/(dy^2) - V_y(t)/(2*dy);
    alpha_x = D/(dx^2) - V_x(t)/(2*dx);
    beta  = - (2*D/(dx^2) + 2*D/(dy^2));
    gamma_x = D/(dx^2) + V_x(t)/(2*dx);
    gamma_y = D/(dy^2) + V_y(t)/(2*dy);
    
    row = [gamma_y*e, gamma_x*e, beta*e, alpha_x*e, alpha_y*e];
    d = [-Nx, -1, 0, 1, Nx];
    A = spdiags(row,d,Nx*Ny,Nx*Ny);
    A = A.*E;
    
    k = A*u_vec + f_vec;
    
end
% Sparse E
function E = SparseE(Nx,Ny)
    E = zeros(Nx*Ny,Nx*Ny);
    for i=(Nx+2):(Nx*Ny-(Nx+1))
        if mod(i,Nx)~= 0 && mod(i,Nx) ~= 1
            E(i,i) = 1;
            E(i,i-1) = 1;
            E(i,i+1) = 1;
            E(i,i-Nx) = 1;
            E(i,i+Nx) = 1;
        end
    end
    E = sparse(E);
end
% Plotter
function DiscretePlotter(X,Y,u,zlimL,zlimU,timesteps,Ntp,Lt)
    % Discrete time plots
    figure('DefaultAxesFontSize',18)
    x_width=1200 ;y_width=600;
    set(gcf, 'Position', [0 0 x_width y_width]);
    set(gcf, 'defaultAxesTickLabelInterpreter','latex')
    set(gcf, 'defaulttextinterpreter','latex')
    set(gcf, 'defaultLegendInterpreter','latex')
    
    t0 = 0; t1 = Lt/5; t2 = 2*Lt/5; t3 = 3*Lt/5; t4 = 4*Lt/5; t5 = 5*Lt/5;
    N0 = ceil((t0*Ntp + 1)/Lt); N1 = ceil((t1*Ntp)/Lt); N2 = ceil((t2*Ntp)/Lt);
    N3 = ceil((t3*Ntp)/Lt); N4 = ceil((t4*Ntp)/Lt); N5 = floor((t5*Ntp)/Lt);

    subplot(2,3,1)
    contourf(X,Y,u(:,:,N0)', 12)
    colormap('turbo')
    colorbar
    hold on
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$u$')
    caxis([zlimL zlimU])
    title(['$t$ = ' num2str(timesteps(N0))])
    grid on

    subplot(2,3,2)
    contourf(X,Y,u(:,:,N1)', 12)
    colormap('turbo')
    colorbar
    hold on
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$u$')
    caxis([zlimL zlimU])
    title(['$t$ = ' num2str(timesteps(N1))])
    grid on

    subplot(2,3,3)
    contourf(X,Y,u(:,:,N2)', 12)
    colormap('turbo')
    colorbar
    hold on
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$u$')
    caxis([zlimL zlimU])
    title(['$t$ = ' num2str(timesteps(N2))])
    grid on

    subplot(2,3,4)
    contourf(X,Y,u(:,:,N3)', 12)
    colormap('turbo')
    colorbar
    hold on
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$u$')
    caxis([zlimL zlimU])
    title(['$t$ = ' num2str(timesteps(N3))])
    grid on

    subplot(2,3,5)
    contourf(X,Y,u(:,:,N4)', 12)
    colormap('turbo')
    colorbar
    hold on
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$u$')
    caxis([zlimL zlimU])
    title(['$t$ = ' num2str(timesteps(N4))])
    grid on

    subplot(2,3,6)
    contourf(X,Y,u(:,:,N5)', 12)
    colormap('turbo')
    colorbar
    hold on
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$u$')
    caxis([zlimL zlimU])
    title(['$t$ = ' num2str(timesteps(N5))])
    grid on
end


