%clear all
close all
clc

savegif = 0;
filename_cont = 'contour.gif';
filename_surf = 'surface.gif';
frameskips = 20;

Lx = 4;
Ly = 4;
Lt = 2;
Nx = 100;
Ny = 100;
Nt = 2000;
Ntp = 2000; % Number of saved time intervals
%Nm = 45;   % Number of measurement stations
%Nu = 52000; % Number of dense measurements
timesteps = linspace(0,Lt,Ntp);

Epochs_Adam = 10000;

x = linspace(-Lx/2,Lx/2,Nx);
y = linspace(-Ly/2,Ly/2,Ny);

A = 1.0;
Spread = 4.0;
Sx = -1.4; % Default -1.4
Sy = -1.0;  % Default 0.6
for i=1:Nx
    for j=1:Ny
        f(i,j) = A*exp(-Spread*(((x(i) - Sx)).^2+(y(j) - Sy).^2)); 
    end
end


x = linspace(-Lx/2,Lx/2,Nx);
y = linspace(-Ly/2,Ly/2,Ny);
t = linspace(0,Lt,Nt);
[X,Y] = meshgrid(x,y);

u_predT = permute(u_pred,[2 1 3]);
u = u_predT;
u_trueT = permute(u_true,[2 1 3]);
u_star = u_trueT;

L2ErrAbs = L2NormErrorAbs(u,u_star);
L2ErrRel = L2NormErrorRel(u,u_star);

%% Error
%ErrorPlotter(t,err_norm);
%% Training history
History(loss_history,Epochs_Adam)
%% Plotting
%AnimatedColorPlotter(X,Y,u,timesteps,Ntp,savegif,filename_cont,frameskips);
%% Contour animation
%AnimatedSurfPlotter(X,Y,u,timesteps,Ntp,Lx,Ly,savegif,filename_surf,frameskips)
%% Discrete plots
DiscretePlotter(X,Y,u,0,0.25,timesteps,Ntp,Lt);
DiscretePlotter(X,Y,u_star,0,0.25,timesteps,Ntp,Lt);
DiscretePlotter(X,Y,abs(u-u_star),0,0.025,timesteps,Ntp,Lt);
%% Tracer plot
%TracerPlot(l1_history,l2_history,Sx,Sy,X,Y,u,f,0,0.2)
%% Training data plotter
%TrainingDataPlotter(X_u_train,Nm,Ntp)
TrainingDataPlotterDense(X_u_train,X_f)
%% Measurment stations
%MeasurementStations(X_u_train,u_train,u_train_pred,Ntp,Nm)
%MeasurementStationsNoPred(X_u_train,u_train,Ntp,Nm)
%% Plotting functions
function AnimatedColorPlotter(X,Y,u,timesteps,Ntp,savegif,filename_cont,frameskips)
    figure('DefaultAxesFontSize',18)
    for i=1:frameskips:Ntp
        %figure('DefaultAxesFontSize',18)
        contourf(X,Y,u(:,:,i)', 12)
        colormap('turbo')
        colorbar
        xlabel('X')
        ylabel('Y')
        title(['Time = ' num2str(timesteps(i)) 's'])
        pause(0.02);
        
        if savegif==1
             frame = getframe(1);
             im = frame2im(frame);
             [imind, cm] = rgb2ind(im,256);
             if i==1
                 imwrite(imind,cm,filename_cont,'gif','Loopcount',inf,'DelayTime', 0.01);
             else
                 imwrite(imind,cm,filename_cont,'gif','WriteMode','append','DelayTime',0.01);
             end
        end
    end
end
function AnimatedSurfPlotter(X,Y,u,timesteps,Ntp,Lx,Ly,savegif,filename_surf,frameskips)
    figure('DefaultAxesFontSize',18)
    for i=1:frameskips:Ntp
        %figure('DefaultAxesFontSize',18)
        surf(X,Y,u(:,:,i)')
        shading flat
        colorbar
        axis([-Lx/2 Lx/2 -Ly/2 Ly/2 0 0.5]);
        title(['Time = ' num2str(timesteps(i)) 's'])
        pause(0.02);
        
        if savegif==1
            frame = getframe(2);
            im = frame2im(frame);
            [imind, cm] = rgb2ind(im,256);
            if i==1
             imwrite(imind,cm,filename_surf,'gif','Loopcount',inf,'DelayTime', 0.01);
            else
             imwrite(imind,cm,filename_surf,'gif','WriteMode','append','DelayTime',0.01);
            end
        end
    end
end
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
function ErrorPlotter(t,err_norm)
    figure('DefaultAxesFontSize',18)
    x_width=800 ;y_width=400;
    set(gcf, 'Position', [0 0 x_width y_width]);
    set(gcf, 'defaultAxesTickLabelInterpreter','latex')
    set(gcf, 'defaulttextinterpreter','latex')
    set(gcf, 'defaultLegendInterpreter','latex')
    
    %semilogy(t, err_norm, 'LineWidth', 2)
    plot(t, err_norm, 'LineWidth', 2)
    title('Error')
    xlabel('$t$')
    ylabel('$\frac{||u(\cdot,\cdot,t)-u^*(\cdot,\cdot,t)||_2}{||u^*(\cdot,\cdot,t)||_2}$')
    grid on
   
end
function History(loss_history, Epochs_Adam)
    figure('DefaultAxesFontSize',18)
    x_width=1200 ;y_width=450;
    set(gcf, 'Position', [0 0 x_width y_width]);
    set(gcf, 'defaultAxesTickLabelInterpreter','latex')
    set(gcf, 'defaulttextinterpreter','latex')
    set(gcf, 'defaultLegendInterpreter','latex')
    
    hold on
    epochs = 1:1:length(loss_history);
    [min_loss, epoch] = min(loss_history);
    semilogy(epochs(1:Epochs_Adam),loss_history(1:Epochs_Adam), 'LineWidth', 2)
    semilogy(epochs(Epochs_Adam:end),loss_history(Epochs_Adam:end), 'LineWidth', 2)
    semilogy(epoch,min_loss,'ks','MarkerSize',16,'LineWidth',2)
    xlabel('Epoch')
    ylabel('Loss')
    xlim([0,length(loss_history)*1.1])
    title('Loss history')
    legend('ADAM', 'LBFG-S')
    grid on
    set(gca,'YScale','log')
    hold off

end
function TrainingDataPlotter(X_u_train,Nm,Ntp)
    data = Nm*Ntp;
    figure('DefaultAxesFontSize',18)
    scatter3(X_u_train(1:data,1), X_u_train(1:data,2), X_u_train(1:data,3),'x')
    hold on
    scatter3(X_u_train(data+1:end,1), X_u_train(data+1:end,2), X_u_train(data+1:end,3))
    xlabel('x')
    ylabel('y')
    zlabel('t')
    %hold off
end
function TrainingDataPlotterDense(X_u_train,X_f)
    figure('DefaultAxesFontSize',18)
    x_width=800 ;y_width=400;
    set(gcf, 'Position', [0 0 x_width y_width]);
    set(gcf, 'defaultAxesTickLabelInterpreter','latex')
    set(gcf, 'defaulttextinterpreter','latex')
    set(gcf, 'defaultLegendInterpreter','latex')
    
    hold on
    scatter3(X_f(:,1), X_f(:,2), X_f(:,3),'x','MarkerFaceAlpha',0.6)
    scatter3(X_u_train(:,1), X_u_train(:,2), X_u_train(:,3),'o','filled')
    legend('Collocation points', 'Data','Interpreter','latex')
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$t$')
    title('Training data')
    hold off
    view([-45,45])
end
function MeasurementStations(X_u_train,u_train,u_train_pred,Ntp,Nm)
    w = 4;
    h = 4;
    if mod(Nm,w*h) ~= 0
        k = floor(double(Nm)/(w*h)) + 1;
    else
        k = Nm/(w*h);
    end
    for j=(1:k)
        figure('DefaultAxesFontSize',12)
        if j == k
            r = mod(Nm,w*h);
        else
            r = w*h;
        end
        for i=1:r
            n = i + (j-1)*w*h;
            startpoint = (n-1)*Ntp+1;
            range = startpoint + int64((1:Ntp-1));
            subplot(w,h,double(i))
            plot(X_u_train(range,3),u_train(range))
            hold on
            plot(X_u_train(range,3),u_train_pred(range),'--','LineWidth',2)
            hold off
            legend({'Train u', 'Pred u'},'Location','northwest')
            xlabel('t')
            ylabel('u')
            title(['Measurement station # ' num2str(n)])
            grid on
        end
    end
end
function MeasurementStationsNoPred(X_u_train,u_train,Ntp,Nm)
    w = 4;
    h = 4;
    if mod(Nm,w*h) ~= 0
        k = floor(double(Nm)/(w*h)) + 1;
    else
        k = Nm/(w*h);
    end
    for j=(1:k)
        figure('DefaultAxesFontSize',12)
        if j == k
            r = mod(Nm,w*h);
        else
            r = w*h;
        end
        for i=1:r
            n = i + (j-1)*w*h;
            startpoint = (n-1)*Ntp+1;
            range = startpoint + int64((1:Ntp-1));
            subplot(w,h,double(i))
            plot(X_u_train(range,3),u_train(range))
            xlabel('t')
            ylabel('u')
            title(['Measurement station # ' num2str(n)])
            grid on
        end
    end
end


function Err = L2NormErrorAbs(u,u_star)
    u_vec = reshape(u,[],1);
    u_star_vec = reshape(u_star,[],1);
    
    %N = length(u_vec);
    N = 1;
    
    Err = norm(u_vec-u_star_vec,2)/N;
    
end

function Err = L2NormErrorRel(u,u_star)
    u_vec = reshape(u,[],1);
    u_star_vec = reshape(u_star,[],1);
    
    %N = length(u_vec);
    N = 1;
    
    Err = norm(u_vec-u_star_vec,2)/(N*norm(u_star_vec,2));
    
end
