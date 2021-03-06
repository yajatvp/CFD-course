%% CODE #2: CN scheme with factorization �
clear
close;
%% Initialization
dt=0.2; %time step in sec
tmax=3000; %simulation time limit
%NX=6; NY=6; %number of grid points including boundary points
%nx=NX-1; ny=NY-1;
nx=32;ny=32;
dx=1/(nx-1);dy=1/(ny-1);
al=9.7*10^(-5); %aluminum thermal diffusivity
r=al*dt/(dx*dx);
T(1:1:nx,1:1:ny,1:1:tmax/dt+1)=200; %Temperature spatio-temporal array initialization
%% CFD - WITH FACTORIZATION
%%% THIS CODE IS DESIGNED ONLY FOR EQUAL NUMBER OF X & Y NODES - BECAUSE OF
%%% THE SYMMETRY OF THE PROBLEM. THIS IS TO CLARIFY THAT THE STUDENT IS ABLE TO SOLVE IT FOR UNEQUAL NUMBER AS WELL;
%%% THIS PROBLEM DID NOT REQUIRE SO IT IS NOT ATTEMPTED.
A1=diag((1+r).*ones((nx-2),1));
b=(-r/2).*ones(-1+(nx-2),1);
A2=diag(b,1);
A3=diag(b,-1);
A=A1+A2+A3;
for nt=2:1:tmax/dt+1
% boundary conditions
T(:,1,nt)=250+(190-250).*(0:1:nx-1).*dx;
T(:,ny,nt)=300+(330-300).*(0:1:nx-1).*dx;
T(1,:,nt)=250+(300-250).*(0:1:ny-1).*dy;
T(nx,:,nt)=190+(330-190).*(0:1:ny-1).*dy;
end
for nt=2:1:tmax/dt+1
% crank-nicholson - To compute T(x,y,nt)
nt
clear del1 del2
for y=2:1:ny-1
for x=2:1:nx-1
B(x-1)=r*((T(x+1,y,nt-1)+T(x-1,y,nt-1)-2*T(x,y,nt-1))+(T(x,y+1,nt-1)+T(x,y-1,nt-1)-2*T(x,y,nt-1)));
end
del1(:,y-1)=A\B';
end
for x=2:1:nx-1
del2(x-1,:)=A\del1(x-1,:)';
end
% Calculating T at nt timestep
T(2:nx-1,2:ny-1,nt)=T(2:nx-1,2:ny-1,nt-1)+del2(:,:);
end
%% POST PROCESSING
nt=1+100/dt;
figure;set(gcf, 'Position', get(0, 'Screensize'));
contourf((0:1:nx-1).*dx,(0:1:ny-1).*dy,T(:,:,nt)');shading flat;axis equal;caxis([190, 330]);
colormap(gca,'hot');c = colorbar;c.Label.String = 'Temperature [K]';
xlim([0,1]);ylim([0,1]);
ylabel('y [m]');xlabel('x [m]');title(['t=',num2str((nt-1)*dt),' s']);
nt=1+500/dt;
figure;set(gcf, 'Position', get(0, 'Screensize'));
contourf((0:1:nx-1).*dx,(0:1:ny-1).*dy,T(:,:,nt)');shading flat;axis equal;caxis([190, 330]);
colormap(gca,'hot');c = colorbar;c.Label.String = 'Temperature [K]';
xlim([0,1]);ylim([0,1]);
ylabel('y [m]');xlabel('x [m]');title(['t=',num2str((nt-1)*dt),' s']);
nt=1+3000/dt;
figure;set(gcf, 'Position', get(0, 'Screensize'));
contourf((0:1:nx-1).*dx,(0:1:ny-1).*dy,T(:,:,nt)');shading flat;axis equal;caxis([190, 330]);
colormap(gca,'hot');c = colorbar;c.Label.String = 'Temperature [K]';
xlim([0,1]);ylim([0,1]);
ylabel('y [m]');xlabel('x [m]');title(['t=',num2str((nt-1)*dt),' s']);