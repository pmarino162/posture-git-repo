clear; clc; clf; close all

%% Setup saveFig   
    saveFig = true;
    saveDir = 'C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Posture\Conferences\NCM 2022\Figures\Neural Space';
    set(0, 'DefaultFigureRenderer', 'painters');
    
%% Setup colormap    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat');
    rcmap = customRainbow;
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    pcmap = orli;

%% Define Trajectory
    t = 0:0.01:1;
    numTimesteps = length(t);
    traj = 0.3*sin(t*2*pi);
    traj2 = .1*sin(3*t);
    figure
    plot(t,traj)
    traj = [t',traj',traj2'];
    trans = [-.5 0 0];
    thetaZ = 110;
    thetaY = 20;
    thetaX = -20;
    Rz = [cosd(thetaZ) -sind(thetaZ) 0; sind(thetaZ) cosd(thetaZ) 0; 0 0 1];
    Ry = [cosd(thetaY) 0 sind(thetaY); 0 1 0; -sind(thetaY) 0 cosd(thetaY)];
    Rx = [1 0 0; 0 cosd(thetaX) -sind(thetaX); 0 sind(thetaX) cosd(thetaX)];
    R = Rz*Ry*Rx;
    traj = (traj+trans)*R;
    figure
    plot3(traj(:,1),traj(:,2),traj(:,3))
    xlabel('x'); ylabel('y'); zlabel('z')
    grid on
    axis square
    
%% Fit PC Plane
    [coeff,score,latent,tsquared,explained,mu] = pca(traj);
    planeColor = [.8,.8,.8];
    d = 0.25;

%% Plot/Animate
    trajColor = 4;
    PC1color = 7;
    PC2color = 8;
    PC3color = 5;
    fps = 50;

    fs = 14;
    posture = 1;
    axLim = .75;
    xLim = [-axLim axLim]; yLim = [-axLim axLim]; zLim = [-axLim axLim];
    ticks = [-.5 0 .5];
    viewAngle = [-31.0260,9.8615];
    
%% Empty space 
    f = figure
    plot3(0,0,0,'-','Color',rcmap(trajColor,:),'LineWidth',2);
    xlabel('Neuron 1 FR'); ylabel('Neuron 2 FR'); zlabel('Neuron 3 FR')
    grid on
    ax = gca;
    ax.XLim = xLim; ax.YLim = yLim; ax.ZLim = zLim;
    view(viewAngle);
    xticks(ticks); yticks(ticks); zticks(ticks);
    xticklabels({}); yticklabels({}); zticklabels({});
    set(gca,'fontname','arial'); set(gca,'fontsize',fs)
    if saveFig
        saveas(gcf,fullfile(saveDir,'EmptySpace.svg'));
        saveas(gcf,fullfile(saveDir,'EmptySpace.fig'));
    end
    
%% Initial point
    f = figure
    plot3(traj(1,1),traj(1,2),traj(1,3),'.','Color',rcmap(trajColor,:),'MarkerSize',25);
    xlabel('Neuron 1 FR'); ylabel('Neuron 2 FR'); zlabel('Neuron 3 FR')
    grid on
    view(viewAngle);
    ax = gca;
    ax.XLim = xLim; ax.YLim = yLim; ax.ZLim = zLim;
    xticks(ticks); yticks(ticks); zticks(ticks);
    xticklabels({}); yticklabels({}); zticklabels({});
    set(gca,'fontname','arial'); set(gca,'fontsize',fs)
    if saveFig
        saveas(gcf,fullfile(saveDir,'InitialPoint.svg'));
        saveas(gcf,fullfile(saveDir,'InitialPoint.fig'));
    end
    
%% Trajectory animation
    fileName = fullfile(saveDir,'trajectory.mp4');
    v = VideoWriter(fileName); v.FrameRate = fps; open(v);
    f = figure; set(gcf,'color','white'); hold on
    %Plot first point
    plot3(traj(1,1),traj(1,2),traj(1,3),'.','Color',rcmap(trajColor,:),'MarkerSize',25);
    view(viewAngle)
    xlabel('Neuron 1 FR'); ylabel('Neuron 2 FR'); zlabel('Neuron 3 FR')
    grid on
    ax = gca;
    ax.XLim = xLim; ax.YLim = yLim; ax.ZLim = zLim;
    xticks(ticks); yticks(ticks); zticks(ticks);
    xticklabels({}); yticklabels({}); zticklabels({});
    set(gca,'fontname','arial'); set(gca,'fontsize',fs)
    %Animate trajectory    
    frameList = [1:length(t)+1];    
    for frame = frameList
        if frame ~= 1
            plot3(traj(1:frame-1,1),traj(1:frame-1,2),traj(1:frame-1,3),'Color',rcmap(trajColor,:),'LineWidth',2);
            plot3(traj(1,1),traj(1,2),traj(1,3),'.','MarkerSize',25,'Color',rcmap(trajColor,:),'MarkerFaceColor',rcmap(trajColor,:));
            if frame == frameList(end)
                plot3(traj(end,1),traj(end,2),traj(end,3),'>','MarkerSize',9,'Color',rcmap(trajColor,:),'MarkerFaceColor',rcmap(trajColor,:));
            end
        end
        drawnow 
        M(frame) = getframe(gcf); writeVideo(v,M(frame));
    end
    close(v);
        

%% Add in plane
    f = figure
    plot3(traj(:,1),traj(:,2),traj(:,3),'-','Color',rcmap(trajColor,:),'LineWidth',2);
    hold on;
    plot3(traj(1,1),traj(1,2),traj(1,3),'.','Color',rcmap(trajColor,:),'MarkerSize',25);
    plot3(traj(end,1),traj(end,2),traj(end,3),'>','MarkerFaceColor',rcmap(trajColor,:),'MarkerEdgeColor',rcmap(trajColor,:),'MarkerSize',9);
    
    n = cross(coeff(:,1),coeff(:,2));
    r0 = mu;
    
    orgX = .25; orgY = .6;
    orgZ = (dot(n,r0)-n(1)*orgX-n(2)*orgY)/n(3);
    quiver3(orgX,orgY,orgZ,-coeff(1,1),-coeff(2,1),-coeff(3,1),'LineWidth',2,'Color',rcmap(PC1color,:))
    quiver3(orgX,orgY,orgZ,-coeff(1,2),-coeff(2,2),-coeff(3,2),'LineWidth',2,'Color',rcmap(PC2color,:))
    quiver3(orgX,orgY,orgZ,coeff(1,3),coeff(2,3),coeff(3,3),'LineWidth',2,'Color',rcmap(PC3color,:))
    
    [surfX,surfY] = meshgrid(xLim(1):d:xLim(2),yLim(1):d:yLim(2));
    for row = 1:size(surfX,1)
        for col = 1:size(surfX,2)
            zX = surfX(row,col);
            zY = surfY(row,col);
            surfZ(row,col) = (dot(n,r0)-n(1)*zX-n(2)*zY)/n(3);
        end
    end
    s = surf(surfX,surfY,surfZ,'EdgeColor','none','FaceColor',planeColor,'FaceAlpha',0.5);
    
    
    xlabel('Neuron 1 FR'); ylabel('Neuron 2 FR'); zlabel('Neuron 3 FR')
    grid on
    ax = gca;
    view(viewAngle)
    ax.XLim = xLim; ax.YLim = yLim; ax.ZLim = zLim;
    xticks(ticks); yticks(ticks); zticks(ticks);
    xticklabels({}); yticklabels({}); zticklabels({});
    set(gca,'fontname','arial'); set(gca,'fontsize',fs)
    if saveFig
        saveas(gcf,fullfile(saveDir,'NeuralSpaceWPlane.svg'));
        saveas(gcf,fullfile(saveDir,'NeuralSpaceWPlane.fig'));
    end
    
%% Plot in PC12 plane; shade
    f = figure
    pcProj = (traj-mu)*-coeff;
    plot3(pcProj(:,1),pcProj(:,2),pcProj(:,3),'-','Color',rcmap(trajColor,:),'LineWidth',2);
    hold on;
    plot3(pcProj(1,1),pcProj(1,2),pcProj(1,3),'.','Color',rcmap(trajColor,:),'MarkerSize',25);
    plot3(pcProj(end,1),pcProj(end,2),pcProj(end,3),'^','MarkerFaceColor',rcmap(trajColor,:),'MarkerEdgeColor',rcmap(trajColor,:),'MarkerSize',9);
    ax = gca;
    xlim = ax.XLim; ylim = ax.YLim;
    [surfX,surfY] = meshgrid([xlim(1),xlim(2)],[ylim(1),ylim(2)]);
    surfZ = zeros(size(surfX))-0.1;
    surf(surfX,surfY,surfZ,'EdgeColor','none','FaceColor',planeColor,'FaceAlpha',0.5);
     plot3(pcProj(:,1),pcProj(:,2),pcProj(:,3),'-','Color',rcmap(trajColor,:),'LineWidth',2);
    hold on;
    plot3(pcProj(1,1),pcProj(1,2),pcProj(1,3),'.','Color',rcmap(trajColor,:),'MarkerSize',25);
    plot3(pcProj(end,1),pcProj(end,2),pcProj(end,3),'^','MarkerFaceColor',rcmap(trajColor,:),'MarkerEdgeColor',rcmap(trajColor,:),'MarkerSize',9);
    
    orgX = xlim(1); orgY = ylim(1); orgZ = 0;
    quiver3(orgX,orgY,orgZ,.75,0,0,'LineWidth',2,'Color',rcmap(PC1color,:))
    quiver3(orgX,orgY,orgZ,0,.75,0,'LineWidth',2,'Color',rcmap(PC2color,:))
    
    xticklabels({}); yticklabels({}); zticklabels({});
    
    xlabel(['Dim 1'])
    ylabel(['Dim 2'])
    set(gca,'fontname','arial')
    set(gca,'fontsize',fs)
    view([0 90])
    ax = gca;
    ax.XLim = xlim; ax.YLim = ylim;
    
    if saveFig
        saveas(gcf,fullfile(saveDir,'PCProj.svg'));
        saveas(gcf,fullfile(saveDir,'PCProj.fig'));
    end

%% Plot in PC23 plane; shade
    f = figure
    pcProj = (traj-mu)*-coeff;
    plot3(pcProj(:,1),pcProj(:,2),pcProj(:,3),'-','Color',rcmap(trajColor,:),'LineWidth',2);
    hold on;
    plot3(pcProj(1,1),pcProj(1,2),pcProj(1,3),'.','Color',rcmap(trajColor,:),'MarkerSize',25);
    plot3(pcProj(end,1),pcProj(end,2),pcProj(end,3),'^','MarkerFaceColor',rcmap(trajColor,:),'MarkerEdgeColor',rcmap(trajColor,:),'MarkerSize',9);
    ax = gca;
    xlim = ax.XLim; ylim = ax.YLim;

    
    orgX = xlim(1); orgY = ylim(1); orgZ = 0;
    quiver3(orgX,orgY,orgZ,0,.75,0,'LineWidth',2,'Color',rcmap(PC2color,:))
    quiver3(orgX,orgY,orgZ,0,0,.75,'LineWidth',2,'Color',rcmap(PC3color,:))
    
    xticklabels({}); yticklabels({}); zticklabels({});
    
    ylabel('Dim 2')
    zlabel('Dim 3')
    set(gca,'fontname','arial')
    set(gca,'fontsize',fs)
    view([90 0])
    ax = gca;
    ax.XLim = xlim; ax.YLim = ylim;
    
    if saveFig
        saveas(gcf,fullfile(saveDir,'PCProj2.svg'));
        saveas(gcf,fullfile(saveDir,'PCProj2.fig'));
    end    