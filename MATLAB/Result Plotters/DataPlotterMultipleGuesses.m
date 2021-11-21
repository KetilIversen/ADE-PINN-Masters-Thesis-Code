close all
%clear all

Lx = 4;     % Default 4
Ly = 4;     % Default 4
Nx = 100;   % Default 100
Ny = 100;   % Default 100
x = linspace(-Lx/2,Lx/2,Nx);
y = linspace(-Ly/2,Ly/2,Ny);
[X,Y] = meshgrid(x,y);
% Continous source
Sx = -1.4;    % Default -0.8
Sy = -1.0;     % Default 0.0
A = 1.0;     % Default 1.0
Spread = 4.0; % Default 4.0
for i=1:Nx
    for j=1:Ny
        f(i,j) = A*exp(-Spread*(((x(i) - Sx)).^2+(y(j) - Sy).^2));
    end
end

zlimL = 0.0;
zlimU = 0.25; % Default 0.25

maxiter = 28000;

u_trueT = permute(u_true,[2 1 3]);
u = u_trueT;

non_conv = 3;

[paths, points] =  PathsandPoints();

[l1_history_ext, l2_history_ext, l3_history_ext, loss_history_ext] = Historyloader(maxiter,paths,points);

for i=1:size(l1_history_ext,1)
   %c(i,:) = rand(1,3); %Random colour
   c(i,:) = [0.9,0.9,0.9];
end
for i=0:non_conv-1
   c(end-i,:) = [1 0 0]; 
end

l1_pred = l1_history_ext(1:end-non_conv,:);
l2_pred = l2_history_ext(1:end-non_conv,:);
l3_pred = l3_history_ext(1:end-non_conv,:);

l1_mean = mean(l1_pred);
l2_mean = mean(l2_pred);
l3_mean = mean(l3_pred);
l1_std = std(l1_pred);
l2_std = std(l2_pred);
l3_std = std(l3_pred);

disp('l1 mean' + string(l1_mean(end)));
disp('l1 std'  + string(l1_std(end)));
disp('l2 mean' + string(l2_mean(end)));
disp('l2 std'  + string(l2_std(end)));
disp('l3 mean' + string(l3_mean(end)));
disp('l3 std'  + string(l3_std(end)));

%LambdaHistory(l1_history_ext,l2_history_ext,l3_history_ext,loss_history_ext,Sx,Sy,A,c,non_conv)
TracerPlot(l1_history_ext,l2_history_ext,Sx,Sy,X,Y,u,f,zlimL,zlimU,c,l1_mean,l2_mean,l1_std,l2_std)
MeanStdPlotter(l1_mean,l2_mean,l3_mean,l1_std,l2_std,l3_std,Sx,Sy,A)
function MeanStdPlotter(l1_mean,l2_mean,l3_mean,l1_std,l2_std,l3_std,Sx,Sy,A)
    figure('DefaultAxesFontSize',18)
    x_width=2000 ;y_width=600;
    set(gcf, 'Position', [0 0 x_width y_width]);
    set(gcf, 'defaultAxesTickLabelInterpreter','latex')
    set(gcf, 'defaulttextinterpreter','latex')
    set(gcf, 'defaultLegendInterpreter','latex')
    
    epoch = linspace(1,length(l1_mean),length(l1_mean));
    step = 150;
    hold on
    errorbar(epoch(1:step:end), l1_mean(1:step:end), l1_std(1:step:end), 'LineWidth', 2)
    errorbar(epoch(1:step:end), l2_mean(1:step:end), l2_std(1:step:end), 'LineWidth', 2)
    errorbar(epoch(1:step:end), l3_mean(1:step:end), l3_std(1:step:end), 'LineWidth', 2)
    yline(Sy,'--', 'LineWidth', 2)
    yline(Sx,'--', 'LineWidth', 2)
    yline(A ,'--', 'LineWidth', 2)
    hold off
    xlim([1, size(l1_mean,2)+10])
    xlabel('Epoch')
    ylabel('$\lambda_{i}$')
    title('Predicted values for $\lambda_{1}$, $\lambda_{2}$ \& $\lambda_{3}$')
    legend('$\lambda_{1}$','$\lambda_{2}$','$\lambda_{3}$','Location','northwest')
    grid on
    
end
function LambdaHistory(l1_history_ext,l2_history_ext,l3_history_ext,loss_history_ext,Sx,Sy,A,c,non_conv)
    figure('DefaultAxesFontSize',18)
    x_width=1200 ;y_width=900;
    set(gcf, 'Position', [0 0 x_width y_width]);
    set(gcf, 'defaultAxesTickLabelInterpreter','latex')
    set(gcf, 'defaulttextinterpreter','latex')
    set(gcf, 'defaultLegendInterpreter','latex')
    
    subplot(2,1,1)
    hold on
    for i=1:size(l1_history_ext,1)-non_conv
        plot(l1_history_ext(i,:), 'LineWidth', 2, 'Color', c(i,:))
        plot(l2_history_ext(i,:), 'LineWidth', 2, 'Color', c(i,:))
        plot(l3_history_ext(i,:), 'LineWidth', 2, 'Color', c(i,:))
    end
    yline(Sy,'--', 'LineWidth', 2)
    yline(Sx,'--', 'LineWidth', 2)
    yline(A ,'--', 'LineWidth', 2)
    hold off
    xlabel('Epoch')
    ylabel('$\lambda_{i}$')
    title('Predicted values for $\lambda_{1}$, $\lambda_{2}$ \& $\lambda_{3}$')

    
    grid on
    
    subplot(2,1,2)
    hold on
    for i=1:size(loss_history_ext,1)-non_conv
        semilogy(loss_history_ext(i,:), 'LineWidth', 2, 'Color', c(i,:))
    end
    xlabel('Epoch')
    ylabel('Loss')
    xlim([0,length(loss_history_ext(1,:))*1.1])
    title('Total objective function loss')
    grid on
    set(gca,'YScale','log')
    hold off

end
function TracerPlot(l1_history_ext,l2_history_ext,Sx,Sy,X,Y,u,f,zlimL,zlimU,c,l1_mean,l2_mean,l1_std,l2_std)
    figure('DefaultAxesFontSize',18)
    x_width=2000 ;y_width=600;
    set(gcf, 'Position', [0 0 x_width y_width]);
    set(gcf, 'defaultAxesTickLabelInterpreter','latex')
    set(gcf, 'defaulttextinterpreter','latex')
    set(gcf, 'defaultLegendInterpreter','latex')
    
    hold on
    %colormap('gray')
    colormap('turbo')
    colorbar
    contourf(X,Y,u(:,:,end)', 12)
    %contourf(X,Y,f',12)
    for i=1:size(l1_history_ext,1)
        plot(l1_history_ext(i,:),l2_history_ext(i,:),'--','LineWidth',2, 'Color', c(i,:))
        plot(l1_history_ext(i,1),l2_history_ext(i,1),'o','MarkerSize',18,'LineWidth',2, 'Color', c(i,:))
        plot(l1_history_ext(i,end),l2_history_ext(i,end),'x','MarkerSize',18,'LineWidth',2, 'Color', c(i,:))
    end
    plot(Sx,Sy,'ks','MarkerSize',18,'LineWidth',2)
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$u$')
    caxis([zlimL zlimU])
    title(['Predicted Source Paths'])
    grid on
    hold off
%     
%     subplot(1,2,2)
%     hold on
%     %colormap('gray')
%     %colorbar
%     %contour(X,Y,u(:,:,end)', 12)
%     for i=1:size(l1_history_ext,1)
%         plot(l1_history_ext(i,end),l2_history_ext(i,end),'x','MarkerSize',18,'LineWidth',2, 'Color', c(i,:))
%     end
%     plot(Sx,Sy,'ks','MarkerSize',18,'LineWidth',2)
%     plot(l1_mean(end),l2_mean(end),'k*','MarkerSize',18,'LineWidth',2)
%     
%     th = 0:pi/50:2*pi;
%     xunit = l1_std(end) * cos(th) + l1_mean(end);
%     yunit = l2_std(end) * sin(th) + l2_mean(end);
%     plot(xunit, yunit,'k--','LineWidth',2);
% 
%     xlabel('$x$')
%     ylabel('$y$')
%     zlabel('$u$')
%     xlim([-1.45,-1.35])
%     ylim([-1.1,-0.9])
%     caxis([zlimL zlimU])
%     title(['Final predictions'])
%     grid on
    
end

function [l1_history_ext, l2_history_ext, l3_history_ext, loss_history_ext] = Historyloader(maxiter,paths,points)
    
    for i=1:length(paths)
        path = paths(i);
        point = points(:,i)';

        [l1_row, l2_row, l3_row, loss_row] = LoaderFunc(path,point,maxiter);

        l1_history_ext(i,:) = l1_row;
        l2_history_ext(i,:) = l2_row;
        l3_history_ext(i,:) = l3_row;
        loss_history_ext(i,:) = loss_row;
    end

end

function [l1_history_ext_row, l2_history_ext_row, l3_history_ext_row, loss_history_ext_row] = LoaderFunc(path,point,maxiter)
    
    l1_loader = load(path, 'l1_history');
    l2_loader = load(path, 'l2_history');
    l3_loader = load(path, 'l3_history');
    loss_loader = load(path, 'loss_history');
    
    l1_history = l1_loader.l1_history;
    l1_history = [l1_history, l1_history(end)*ones(1,maxiter-length(l1_history))];
    
    l2_history = l2_loader.l2_history;
    l2_history = [l2_history, l2_history(end)*ones(1,maxiter-length(l2_history))];
    
    l3_history = l3_loader.l3_history;
    l3_history = [l3_history, l3_history(end)*ones(1,maxiter-length(l3_history))];
    
    loss_history = loss_loader.loss_history;
    loss_history = [loss_history, loss_history(end)*ones(1,maxiter-length(loss_history))];
    
    l1_history_ext_row = [point(1),  l1_history];
    l2_history_ext_row = [point(2),  l2_history];
    l3_history_ext_row = [point(3),  l3_history];
    loss_history_ext_row = loss_history;

end

function [paths, points] = PathsandPoints()

    paths(1) = "C:\Users\Ketil\Documents\MATLAB\Master\Results in thesis final\Multiple guesses\PINNresultsInverse(0.0,0.0).mat";
    paths(2) = "C:\Users\Ketil\Documents\MATLAB\Master\Results in thesis final\Multiple guesses\PINNresultsInverse(0.0,-1.0).mat";
    paths(3) = "C:\Users\Ketil\Documents\MATLAB\Master\Results in thesis final\Multiple guesses\PINNresultsInverse(0.0,-2.0).mat";
    paths(4) = "C:\Users\Ketil\Documents\MATLAB\Master\Results in thesis final\Multiple guesses\PINNresultsInverse(0.5,-0.5).mat";
    paths(5) = "C:\Users\Ketil\Documents\MATLAB\Master\Results in thesis final\Multiple guesses\PINNresultsInverse(1.0,-1.0).mat";
    paths(6) = "C:\Users\Ketil\Documents\MATLAB\Master\Results in thesis final\Multiple guesses\PINNresultsInverse(-0.5,0.5).mat";
    paths(7) = "C:\Users\Ketil\Documents\MATLAB\Master\Results in thesis final\Multiple guesses\PINNresultsInverse(-0.5,-0.5).mat";
    paths(8) = "C:\Users\Ketil\Documents\MATLAB\Master\Results in thesis final\Multiple guesses\PINNresultsInverse(-1.0,0.0).mat";
    paths(9) = "C:\Users\Ketil\Documents\MATLAB\Master\Results in thesis final\Multiple guesses\PINNresultsInverse(-1.0,1.0).mat";
    paths(10) = "C:\Users\Ketil\Documents\MATLAB\Master\Results in thesis final\Multiple guesses\PINNresultsInverse(-2.0,0.0).mat";
    
    %Non-convergent
    paths(11) = "C:\Users\Ketil\Documents\MATLAB\Master\Results in thesis final\Multiple guesses\PINNresultsInverse(0.5,0.5).mat";
    paths(12) = "C:\Users\Ketil\Documents\MATLAB\Master\Results in thesis final\Multiple guesses\PINNresultsInverse(1.0,1.0).mat";
    paths(13) = "C:\Users\Ketil\Documents\MATLAB\Master\Results in thesis final\Multiple guesses\PINNresultsInverse(2.0,0.0).mat";
    
    points(:,1) = [0.0,0.0,0.5];
    points(:,2) = [0.0,-1.0,0.5];
    points(:,3) = [0.0,-2.0,0.5];
    points(:,4) = [0.5,-0.5,0.5];
    points(:,5) = [1.0,-1.0,0.5];
    points(:,6) = [-0.5,0.5,0.5];
    points(:,7) = [-0.5,-0.5,0.5];
    points(:,8) = [-1.0,0.0,0.5];
    points(:,9) = [-1.0,1.0,0.5];
    points(:,10) = [-2.0,0.0,0.5];
    
    points(:,11) = [0.5,0.5,0.5];
    points(:,12) = [1.0,1.0,0.5];
    points(:,13) = [2.0,0.0,0.5];

end
