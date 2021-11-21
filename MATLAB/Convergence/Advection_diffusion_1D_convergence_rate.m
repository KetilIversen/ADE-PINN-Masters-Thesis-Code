%Clear everything
clc
%clear all 
close all

Nt = [1000,1000,1000,1000,2000,8000,32000];
Nx = [16,32,64,128,256,512,1024];
Nt_max = 128000;
Nx_max = 2048;
dx = 2./(Nx-1);

u_exact = ADE_Solver(Nx_max,Nt_max);

for i = 1:length(Nx)

    xq = linspace(-1,1,Nx(i));
    x  = linspace(-1,1,Nx_max);
    t  = linspace(0,1,100);

    [Tq,Xq] = meshgrid(t,xq);
    [T,X] = meshgrid(t,x);

    u = ADE_Solver(Nx(i),Nt(i));

    u_interp = interp2(Tq,Xq,u,T,X);
    
    Err(i) = L2NormErrorRel(u_interp,u_exact);
    
%     if i > 1
%         
%         Err(i-1) = L2NormErrorRel(V_old,V);
%         V_old = V;
%        
%     else
%         
%         V_old = V;
%         
%     end
 
%     figure
%     contourf(T,X,V);
%     figure
%     contourf(Tq,Xq,u);
    
end

%% Plot solution
%DiscretePlotter(xq,u,100,1);

% Error plot over step size
figure('DefaultAxesFontSize',18)
x_width=1180 ;y_width=550;
set(gcf, 'Position', [0 0 x_width y_width]);
set(gcf, 'defaultAxesTickLabelInterpreter','latex')
set(gcf, 'defaulttextinterpreter','latex')
set(gcf, 'defaultLegendInterpreter','latex')
title('Order of convergence')

o1 = (Err(1)/(dx(1)))*dx;
o2 = (Err(1)/(dx(1)^2))*dx.^2;
o3 = (Err(1)/(dx(1)^3))*dx.^3;
loglog(dx,Err,'*-','LineWidth',3)
hold on
loglog(dx,o1,'k--','LineWidth',2)
hold on
loglog(dx,o2,'k-.','LineWidth',2)
hold on
loglog(dx,o3,'k:','LineWidth',2)
hold off
legend('Relative $L^2$-error','$\mathcal{O}(dx)$','$\mathcal{O}(dx^2)$','$\mathcal{O}(dx^3)$','Location','northwest');
xlabel('$dx$')
ylabel('$L^2$-error')
grid on


function Err = L2NormErrorRel(u,u_star)
    u_vec = reshape(u,[],1);
    u_star_vec = reshape(u_star,[],1);
    
    Err = norm(u_vec-u_star_vec,2)/(norm(u_star_vec,2));
    
end


function u = ADE_Solver(Nx,Nt)

    % Initialize model domains
    Lx = 2;     % Default 2      (Length of x domain)
    Lt = 1;     % Default 1      (Length of t domain)
    %Nx = 256;   % Default 256    (Gridpoints in x)
    %Nt = 8000;  % Default 10000  (Computational timesteps)
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

end

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
