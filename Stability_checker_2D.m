close all

Nx = 100;
Ny = 100;
Nt = 2000;
Lx = 4;
Ly = 4;
Lt = 2;
D = 0.2;
V_x = 2.0;
V_y = 3.0+1.0;
dx = Lx/(Nx-1);
dy = Ly/(Ny-1);
dt = Lt/(Nt-1);

alpha_y = D/(dy^2) - V_y/(2*dy);
alpha_x = D/(dx^2) - V_x/(2*dx);
beta  = - (2*D/(dx^2) + 2*D/(dy^2));
gamma_x = D/(dx^2) + V_x/(2*dx);
gamma_y = D/(dy^2) + V_y/(2*dy);


A = zeros(Nx*Ny,Nx*Ny);
for i=(Nx+2):(Nx*Ny-(Nx+1))
    if mod(i,Nx)~= 0 && mod(i,Nx) ~= 1
        A(i,i) = beta;
        A(i,i-1) = gamma_x;
        A(i,i+1) = alpha_x;
        A(i,i-Nx) = gamma_y;
        A(i,i+Nx) = alpha_y;
    end
end
tic
lambda = eig(A);
toc
%%
%%
StabilityPlotter(lambda,dt)

% figure
% spy(A);


%% Functions
function RK4Plotter()
    % #### Comparing domain of stability of Runge Kutta method of order 2,3,4,6
    % Specify x range and number of points
    x0 = -3;
    x1 = 3;
    Nx = 301;
    % Specify y range and number of points
    y0 = -3;
    y1 = 3;
    Ny = 301;
    % Construct mesh
    xv = linspace(x0,x1,Nx);
    yv = linspace(y0,y1,Ny);
    [x,y] = meshgrid(xv,yv);
    % Calculate z
    z = x + 1i*y;
    % 2nd order Runge-Kutta growth factor
    g1 = 1 + z + 0.5*z.^2;
    % 3rd order Runge-Kutta growth factor
    g2 = 1 + z + 1/2*z.^2 + 1/6*z.^3;
    % 4th order Runge-Kutta growth factor
    g3 = 1 + z + 1/2*z.^2 + 1/6*z.^3 + 1/24*z.^4;
    % 6th order Runge-Kutta growth factor
    %g4 = 1 + z + 1/2*z.^2 + 1/6*z.^3 + 1/24*z.^4+1/120*z.^5+1/720*z.^6;
    % Calculate magnitude of g
    gmag1 = abs(g1);
    gmag2 = abs(g2);
    gmag3 = abs(g3);
    %gmag4 = abs(g4);
    % Plot contours of gmag
    contour(x,y,gmag1,[1 1],'-','LineWidth',2,'Color',[0.8500, 0.3250, 0.0980]);
    hold on;
    contour(x,y,gmag2,[1 1],'-','LineWidth',2,'Color',[0.9290, 0.6940, 0.1250]);
    contour(x,y,gmag3,[1 1],'-','LineWidth',2,'Color',[0.4940, 0.1840, 0.5560]);
    %contour(x,y,gmag4,[1 1],'b-','LineWidth',2);
    %axis([x0,x1,y0,y1]);
end

function CirclePlotter(x,y,r)
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    plot(xunit, yunit,'LineWidth',2,'Color',[0, 0.4470, 0.7410]);
end

function StabilityPlotter(lambda,dt)

    figure('DefaultAxesFontSize',18)
    x_width=1000 ;y_width=750;
    set(gcf, 'Position', [0 0 x_width y_width]);
    set(gcf, 'defaultAxesTickLabelInterpreter','latex')
    set(gcf, 'defaulttextinterpreter','latex')
    set(gcf, 'defaultLegendInterpreter','latex')
    title('Scaled eigenvalues of 1D ADE finite difference scheme')
    
    subplot(1,2,1)
    hold on
    %scatter(real(lambda*dt), imag(lambda*dt), 160, '.','MarkerEdgeColor',[0.4660 0.6740 0.1880])
    scatter(real(lambda*dt), imag(lambda*dt), 160, 'k.')
    CirclePlotter(-1,0,1);
    RK4Plotter();
    axis([-4.0,2.0,-3,3])
    axis('square');
    xlabel('Real($\lambda\Delta t$)');
    ylabel('Imag($\lambda\Delta t$)');
    grid on;
    legend('$\lambda\Delta t$','F-Euler','RK2','RK3','RK4')
    hold off
    
    subplot(1,2,2)
    hold on
    %scatter(real(lambda*dt), imag(lambda*dt), 160, '.','MarkerEdgeColor',[0.4660 0.6740 0.1880])
    scatter(real(lambda*dt), imag(lambda*dt), 160, 'k.')
    CirclePlotter(-1,0,1);
    RK4Plotter();
    axis([-0.1,0.1,-0.1,0.1])
    axis('square');
    xlabel('Real($\lambda\Delta t$)');
    ylabel('Imag($\lambda\Delta t$)');
    grid on;
    legend('$\lambda\Delta t$','F-Euler','RK2','RK3','RK4')
    hold off
    
end