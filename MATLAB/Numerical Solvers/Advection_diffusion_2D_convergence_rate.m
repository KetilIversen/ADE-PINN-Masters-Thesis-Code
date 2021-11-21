%close all
%clear all

Nt = [200,200,400,2000,2000,2000,2000];
N = [10,20,40,80,100,120,150];
Nt_max = 4000;
N_max = 200;
dx = 4./(N-1);

u_exact = ADE_Solver(N_max,Nt_max);
%u_exact = u_exact(:,:,end);

for i = 1:length(N)

    xq = linspace(-2,2,N(i));
    yq = linspace(-2,2,N(i));
    tq = linspace(0,2,Nt(i));
    x  = linspace(-2,2,N_max);
    y  = linspace(-2,2,N_max);
    t  = linspace(0,2,Nt_max);

    [Xq,Yq,Tq] = meshgrid(xq,yq,tq);
    [X,Y,T] = meshgrid(x,y,t);

    u = ADE_Solver(N(i),Nt(i));

    u_interp = interp3(Xq,Yq,Tq,u,X,Y,T);
    
    Err(i) = L2NormErrorRel(u_interp,u_exact);
    
end

%%
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


function u = ADE_Solver(N,Nt)

    % Initialize model domains
    Lx = 4;     % Default 4     (Length of x domain)
    Ly = 4;     % Default 4     (Length of y domain)
    Lt = 2;     % Default 2     (Length of t domain)
    Nx = N;   % Default 100   (Gridpoints in x)
    Ny = N;   % Default 100   (Gridpoints in y)
    %Nt = 20000; % Default 20000 (Computational timesteps)
    Ntp = Nt; % Default 2000  (Saved timesteps)

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
    S = 4.0;            % Default 4.0  (Source spread)
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

end

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

function Err = L2NormErrorRel(u,u_star)
    u_vec = reshape(u,[],1);
    u_star_vec = reshape(u_star,[],1);
    
    Err = norm(u_vec-u_star_vec,2)/(norm(u_star_vec,2));
    
end

