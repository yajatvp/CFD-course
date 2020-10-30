clear 
close;

%% Initialization
dt=0.2; %time step in sec
tmax=3000; %simulation time limit
nx=32; ny=32; %number of grid points
dx=1/nx;dy=1/ny;
al=9.7*10^(-5); %aluminum thermal diffusivity
r=2*al*dt/(dx*dx);
T(1:1:nx+1,1:1:ny+1,1:1:tmax/dt+1)=200; %Temperature spatio-temporal array initialization

%% CFD
sx=al*dt/(dx)^2;sy=al*dt/(dy)^2;
res(1)=0;
for nt=2:1:tmax/dt+1
    % boundary conditions
    T(:,1,nt)=250+(190-250).*(0:1:nx).*dx;
    T(:,ny+1,nt)=300+(330-300).*(0:1:nx).*dx;
    T(1,:,nt)=250+(300-250).*(0:1:ny).*dy;
    T(nx+1,:,nt)=190+(330-190).*(0:1:ny).*dy;
    
    % finite diff scheme
    for y=2:1:ny
        for x=2:1:nx
            T(x,y,nt)=T(x,y,nt-1)+sx*(T(x,y+1,nt-1)+T(x,y-1,nt-1)-2*T(x,y,nt-1))...
                +sy*(T(x+1,y,nt-1)+T(x-1,y,nt-1)-2*T(x,y,nt-1));
        end
    end
%     res(nt)=mean(mean(T(:,:,nt)-T(:,:,nt-1))); %residual avg.
%     if res(nt)<10^(-4)
%         tim=(nt-1)*dt;
%         break
%     end
end

%% POST PROCESSING
nt=1+100/dt;
figure;set(gcf, 'Position', get(0, 'Screensize'));
pcolor((0:1:nx).*dx,(0:1:ny).*dy,T(:,:,nt)');shading flat;axis equal;caxis([190, 330]);
colormap(gca,'hot');c = colorbar;c.Label.String = 'Temperature [K]';
xlim([0,1]);ylim([0,1]);
ylabel('y [m]');xlabel('x [m]');title(['t=',num2str((nt-1)*dt),' s']);

nt=1+500/dt;
figure;set(gcf, 'Position', get(0, 'Screensize'));
pcolor((0:1:nx).*dx,(0:1:ny).*dy,T(:,:,nt)');shading flat;axis equal;caxis([190, 330]);
colormap(gca,'hot');c = colorbar;c.Label.String = 'Temperature [K]';
xlim([0,1]);ylim([0,1]);
ylabel('y [m]');xlabel('x [m]');title(['t=',num2str((nt-1)*dt),' s']);

nt=1+3000/dt;
figure;set(gcf, 'Position', get(0, 'Screensize'));
pcolor((0:1:nx).*dx,(0:1:ny).*dy,T(:,:,nt)');shading flat;axis equal;caxis([190, 330]);
colormap(gca,'hot');c = colorbar;c.Label.String = 'Temperature [K]';
xlim([0,1]);ylim([0,1]);
ylabel('y [m]');xlabel('x [m]');title(['t=',num2str((nt-1)*dt),' s']);
%%
for nt=1:10/dt:tmax/dt+1
    figure;set(gcf, 'Position', get(0, 'Screensize'));
    pcolor((0:1:nx).*dx,(0:1:ny).*dy,T(:,:,nt)');shading flat;axis equal;caxis([190, 330]);
    colormap(gca,'hot');c = colorbar;c.Label.String = 'Temperature [K]';
    xlim([0,1]);ylim([0,1]);
    ylabel('y [m]');xlabel('x [m]');title(['t=',num2str((nt-1)*dt),' s']);
    w = waitforbuttonpress;
    close
end

%%
figure;set(gcf, 'Position', get(0, 'Screensize'));
contour((0:1:nx).*dx,(0:1:ny).*dy,T(:,:,end)',15);shading flat;axis equal;caxis([190, 320]);
colormap(gca,'copper');c = colorbar;c.Label.String = 'Temperature [K]';
xlim([0,1]);ylim([0,1]);
ylabel('y [m]');xlabel('x [m]');title(['t=',num2str((nt-1)*dt),' s']);
%%
figure;set(gcf, 'Position', get(0, 'Screensize'));
plot([100,256,400,2500,74*74,88*88,10000,12100,14400],...
    [267.1391,267.1456,267.1473,267.1495,267.1497,267.1497,267.1499,267.1498,267.1498],'-o');grid on;
ylabel('Temperature at midpoint - (0.5, 0.5) m  [K]');xlabel('Total grid points - n_xn_y');
