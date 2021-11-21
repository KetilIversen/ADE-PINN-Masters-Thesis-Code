%clear all
close all
clc

savegif = 0;
filename_cont = 'contour.gif';
filename_surf = 'surface.gif';
frameskips = 20;

Lx = x(end)-x(1);
Lt = t(end)-t(1);
Nx = length(x);
Nt = length(t);
%Ntp = 2000; % Number of saved time intervals
%Nm = 45;   % Number of measurement stations
%Nu = 5200; % Number of dense measurements

Epochs_Adam = 15000;

Sx = -0.8;
A = 5.0;

timesteps = linspace(0,Lt,Nt);
%x = linspace(-Lx/2,Lx/2,Nx);

%x = linspace(-Lx/2,Lx/2,Nx);
%y = linspace(-Ly/2,Ly/2,Ny);
%t = linspace(0,Lt,Nt);
%[X,Y] = meshgrid(x,y);

u_predT = permute(u_pred,[2 1]);
u = u_predT;
u_trueT = permute(u_true,[2 1]);
u_star = u_trueT;

%err_norm = L2_norm_error(u,u_star,Nt);
L2ErrRel = L2NormErrorRel(u,u_star);

%% Error
%ErrorPlotter(t,err_norm);
%% Training history
History(loss_history,l1_history,l2_history,Sx,A,Epochs_Adam)
%% Plotting
ColorPlotter(X,T,u,X_u_train);
%% Discrete plots
DiscretePlotter(x,t,u,u_star,Nt,Lt,-0.1,1.0);
%% Tracer plot
%TracerPlot(l1_history,l2_history,Sx,Sy,X,Y,u,f,0,0.2)
%% Training data plotter
%TrainingDataPlotter(X_u_train,Nm,Ntp)
%TrainingDataPlotterDense(X_u_train,Nu)
%% Measurment stations
%MeasurementStations(X_u_train,u_train,u_train_pred,Ntp,Nm)
%MeasurementStationsNoPred(X_u_train,u_train,Ntp,Nm)
%% Plotting functions
function ColorPlotter(X,T,u,X_u_train)
    figure('DefaultAxesFontSize',18)
    x_width=800 ;y_width=400;
    set(gcf, 'Position', [0 0 x_width y_width]);
    set(gcf, 'defaultAxesTickLabelInterpreter','latex')
    set(gcf, 'defaulttextinterpreter','latex')
    set(gcf, 'defaultLegendInterpreter','latex')
    hold on
    %contourf(T,X,u(:,:), 12)
    %scatter(X_f(:,2),X_f(:,1),'x','MarkerFaceAlpha',0.6)
    scatter(X_u_train(:,2),X_u_train(:,1),'o','filled')
    legend('Data','Interpreter','latex')
    %colormap('turbo')
    %colorbar
    xlabel('$t$')
    ylabel('$x$')
    title('Training data')
    hold off
end
function DiscretePlotter(x,t,u,u_star,Ntp,Lt,ulimL,ulimU)
    % Discrete time plots
    figure('DefaultAxesFontSize',18)
    x_width=800 ;y_width=400;
    set(gcf, 'Position', [0 0 x_width y_width]);
    set(gcf, 'defaultAxesTickLabelInterpreter','latex')
    set(gcf, 'defaulttextinterpreter','latex')
    set(gcf, 'defaultLegendInterpreter','latex')
    
    t0 = 0; t1 = Lt/5; t2 = 2*Lt/5; t3 = 3*Lt/5; t4 = 4*Lt/5; t5 = 5*Lt/5;
    N0 = ceil((t0*Ntp + 1)/Lt); N1 = ceil((t1*Ntp)/Lt); N2 = ceil((t2*Ntp)/Lt);
    N3 = ceil((t3*Ntp)/Lt); N4 = ceil((t4*Ntp)/Lt); N5 = floor((t5*Ntp)/Lt);

    subplot(2,3,1)
    hold on
    plot(x,u(N0,:),'LineWidth',4,'Color',[0 0.4470 0.7410]);
    plot(x,u_star(N0,:),'--','LineWidth',4,'Color',[0.8500 0.3250 0.0980]);
    xlabel('$x$')
    ylabel('$u$')
    ylim([ulimL,ulimU])
    legend('Prediction','Reference')
    title(['$t$ = ' num2str(t(N0))])
    grid on
    hold off
    
    subplot(2,3,2)
    hold on
    plot(x,u(N1,:),'LineWidth',4,'Color',[0 0.4470 0.7410]);
    plot(x,u_star(N1,:),'--','LineWidth',4,'Color',[0.8500 0.3250 0.0980]);
    hold off
    xlabel('$x$')
    ylabel('$u$')
    ylim([ulimL,ulimU])
    title(['$t$ = ' num2str(t(N1))])
    grid on

    subplot(2,3,3)
    hold on
    plot(x,u(N2,:),'LineWidth',4,'Color',[0 0.4470 0.7410]);
    plot(x,u_star(N2,:),'--','LineWidth',4,'Color',[0.8500 0.3250 0.0980]);
    hold off
    xlabel('$x$')
    ylabel('$u$')
    ylim([ulimL,ulimU])
    title(['$t$ = ' num2str(t(N2))])
    grid on

    subplot(2,3,4)
    hold on
    plot(x,u(N3,:),'LineWidth',4,'Color',[0 0.4470 0.7410]);
    plot(x,u_star(N3,:),'--','LineWidth',4,'Color',[0.8500 0.3250 0.0980]);
    hold off
    xlabel('$x$')
    ylabel('$u$')
    ylim([ulimL,ulimU])
    title(['$t$ = ' num2str(t(N3))])
    grid on

    subplot(2,3,5)
    hold on
    plot(x,u(N4,:),'LineWidth',4,'Color',[0 0.4470 0.7410]);
    plot(x,u_star(N4,:),'--','LineWidth',4,'Color',[0.8500 0.3250 0.0980]);
    hold off
    xlabel('$x$')
    ylabel('$u$')
    ylim([ulimL,ulimU])
    title(['$t$ = ' num2str(t(N4))])
    grid on

    subplot(2,3,6)
    hold on
    plot(x,u(N5,:),'LineWidth',4,'Color',[0 0.4470 0.7410]);
    plot(x,u_star(N5,:),'--','LineWidth',4,'Color',[0.8500 0.3250 0.0980]);
    hold off
    xlabel('$x$')
    ylabel('$u$')
    ylim([ulimL,ulimU])
    title(['$t$ = ' num2str(t(N5))])
    grid on
    
end
function ErrorPlotter(t,err_norm)
    figure('DefaultAxesFontSize',18)
    x_width=800 ;y_width=400;
    set(gcf, 'Position', [0 0 x_width y_width]);
    set(gcf, 'defaultAxesTickLabelInterpreter','latex')
    set(gcf, 'defaulttextinterpreter','latex')
    set(gcf, 'defaultLegendInterpreter','latex')
    
    semilogy(t, err_norm, 'LineWidth', 2)
    %plot(t, err_norm, 'LineWidth', 2)
    title('Error')
    xlabel('$t$')
    ylabel('$\frac{||u(\cdot,t)-u^*(\cdot,t)||_2}{||u^*(\cdot,t)||_2}$')
    grid on

end
function History(loss_history,l1_history,l2_history,Sx,A,Epochs_Adam)
    figure('DefaultAxesFontSize',18)
    x_width=1200 ;y_width=900;
    set(gcf, 'Position', [0 0 x_width y_width]);
    set(gcf, 'defaultAxesTickLabelInterpreter','latex')
    set(gcf, 'defaulttextinterpreter','latex')
    set(gcf, 'defaultLegendInterpreter','latex')
    
    subplot(2,1,1)
    hold on
    plot(l1_history, 'LineWidth', 2)
    plot(l2_history, 'LineWidth', 2)
    yline(Sx,'--', 'LineWidth', 2)
    yline(A ,'--', 'LineWidth', 2)
    hold off
    xlabel('Epoch')
    ylabel('$\lambda_{i}$')
    title('Predicted values for $\lambda_{1}$ \& $\lambda_{2}$')
    legend('$\lambda_{1}$','$\lambda_{2}$','Location','northwest')
    grid on
    
    subplot(2,1,2)
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
function TracerPlot(l1_history,l2_history,Sx,Sy,X,Y,u,f,zlimL,zlimU)
    figure('DefaultAxesFontSize',18)
    fontsize = 14;
    hold on
    colormap('turbo')
    colorbar
    contourf(X,Y,u(:,:,end)', 12)
    %contour(X,Y,f', 12, 'k--')
    plot(l1_history,l2_history,'k--','LineWidth',2)
    plot(Sx,Sy,'ks','MarkerSize',18,'LineWidth',2)
    plot(l1_history(1),l2_history(1),'ko','MarkerSize',18,'LineWidth',2)
    plot(l1_history(end),l2_history(end),'kx','MarkerSize',18,'LineWidth',2)
    xlabel('x')
    ylabel('y')
    zlabel('u')
    caxis([zlimL zlimU])
    title(['Predicted source path'])
    grid on
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
function TrainingDataPlotterDense(X_u_train,Nu)
    data = Nu;
    figure('DefaultAxesFontSize',18)
    scatter3(X_u_train(1:data,1), X_u_train(1:data,2), X_u_train(1:data,3),'x')
    hold on
    scatter3(X_u_train(data+1:end,1), X_u_train(data+1:end,2), X_u_train(data+1:end,3))
    xlabel('x')
    ylabel('y')
    zlabel('t')
    %hold off
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
%% Error calc
% function err_norm = L2_norm_error(u,u_star,Nt)
%     err_norm(Nt) = 0;
%     for i=1:Nt
%         err_norm(i) = norm(u(i,:)-u_star(i,:),2)/norm(u_star(i,:),2); 
%     end
%     
%     total_rel_err = norm(u-u_star,2)/norm(u_star,2)
% 
% end

function Err = L2NormErrorRel(u,u_star)
    u_vec = reshape(u,[],1);
    u_star_vec = reshape(u_star,[],1);
    
    %N = length(u_vec);
    N = 1;
    
    Err = norm(u_vec-u_star_vec,2)/(N*norm(u_star_vec,2));
    
end