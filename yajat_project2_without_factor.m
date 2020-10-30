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
%% CFD - WITHOUT FACTORIZATION
%%% THIS CODE IS DESIGNED ONLY FOR EQUAL NUMBER OF X & Y NODES - BECAUSE OF
%%% THE SYMMETRY OF THE PROBLEM. THIS IS TO CLARIFY THAT THE STUDENT IS ABLE TO SOLVE IT FOR UNEQUAL NUMBER AS WELL;
%%% THIS PROBLEM DID NOT REQUIRE SO IT IS NOT ATTEMPTED.
A1=diag((1+2*r).*ones((nx-2)^2,1));
b=(-r/2).*ones(-1+(nx-2)^2,1);
for i=1:size(b)
if mod(i,nx-2)==0
b(i)=0;
end
end
A2=diag(b,1);
A3=diag(b,-1);
A4=diag((-r/2).*ones(-nx+2+(nx-2)^2,1),nx-2);
A5=diag((-r/2).*ones(-nx+2+(nx-2)^2,1),-nx+2);
A=A1+A2+A3+A4+A5;
for nt=2:1:tmax/dt+1
% boundary conditions
T(:,1,nt)=250+(190-250).*(0:1:nx-1).*dx;
T(:,ny,nt)=300+(330-300).*(0:1:nx-1).*dx;
T(1,:,nt)=250+(300-250).*(0:1:ny-1).*dy;
T(nx,:,nt)=190+(330-190).*(0:1:ny-1).*dy;
end
%
for nt=2:1:tmax/dt+1
% crank-nicholson - To compute T(x,y,nt)
nt
n=0;
for x=2:1:nx-1
for y=2:1:nx-1
n=n+1;
if x==2 && y~=2 && y~=(ny-1)
B(n)=T(x,y,nt-1)+...
(r/2)*((T(x+1,y,nt-1)+T(x-1,y,nt-1)-2*T(x,y,nt-1))+(T(x,y+1,nt-1)+T(x,y-1,nt-1)-2*T(x,y,nt-1)))...
+(r/2)*T(x-1,y,nt);
elseif y==2 && x~=2 && x~=(nx-1)
B(n)=T(x,y,nt-1)+...
(r/2)*((T(x+1,y,nt-1)+T(x-1,y,nt-1)-2*T(x,y,nt-1))+(T(x,y+1,nt-1)+T(x,y-1,nt-1)-2*T(x,y,nt-1)))...
+(r/2)*T(x,y-1,nt);
elseif x==(nx-1) && y~=2 && y~=(ny-1)
B(n)=T(x,y,nt-1)+...
(r/2)*((T(x+1,y,nt-1)+T(x-1,y,nt-1)-2*T(x,y,nt-1))+(T(x,y+1,nt-1)+T(x,y-1,nt-1)-2*T(x,y,nt-1)))...
+(r/2)*T(x+1,y,nt);
elseif y==(ny-1) && x~=2 && x~=(nx-1)
B(n)=T(x,y,nt-1)+...
(r/2)*((T(x+1,y,nt-1)+T(x-1,y,nt-1)-2*T(x,y,nt-1))+(T(x,y+1,nt-1)+T(x,y-1,nt-1)-2*T(x,y,nt-1)))...
+(r/2)*T(x,y+1,nt);
elseif x==2 && y==2
B(n)=T(x,y,nt-1)+...
(r/2)*((T(x+1,y,nt-1)+T(x-1,y,nt-1)-2*T(x,y,nt-1))+(T(x,y+1,nt-1)+T(x,y-1,nt-1)-2*T(x,y,nt-1)))...
+(r/2)*(T(x-1,y,nt)+T(x,y-1,nt));
elseif x==2 && y==(ny-1)
B(n)=T(x,y,nt-1)+...
(r/2)*((T(x+1,y,nt-1)+T(x-1,y,nt-1)-2*T(x,y,nt-1))+(T(x,y+1,nt-1)+T(x,y-1,nt-1)-2*T(x,y,nt-1)))...
+(r/2)*(T(x-1,y,nt)+T(x,y+1,nt));
elseif x==(nx-1) && y==2
B(n)=T(x,y,nt-1)+...
(r/2)*((T(x+1,y,nt-1)+T(x-1,y,nt-1)-2*T(x,y,nt-1))+(T(x,y+1,nt-1)+T(x,y-1,nt-1)-2*T(x,y,nt-1)))...
+(r/2)*(T(x+1,y,nt)+T(x,y-1,nt));
elseif x==(nx-1) && y==(ny-1)
B(n)=T(x,y,nt-1)+...
(r/2)*((T(x+1,y,nt-1)+T(x-1,y,nt-1)-2*T(x,y,nt-1))+(T(x,y+1,nt-1)+T(x,y-1,nt-1)-2*T(x,y,nt-1)))...
+(r/2)*(T(x+1,y,nt)+T(x,y+1,nt));
else
B(n)=T(x,y,nt-1)+...
(r/2)*((T(x+1,y,nt-1)+T(x-1,y,nt-1)-2*T(x,y,nt-1))+(T(x,y+1,nt-1)+T(x,y-1,nt-1)-2*T(x,y,nt-1)));
end
end
end
% Calculating T at nt timestep
%TT(:)=A\B';
ain=inv(A);
for i=1:size(ain,1)
for j=1:size(ain,2)
if ain(i,j)<10^-5
ain(i,j)=0;
end
end
end
TT(:)=ain*B';
n=0;
for x=2:1:nx-1
for y=2:1:ny-1
n=n+1;
T(x,y,nt)=TT(n);
end
end
end
%% POST PROCESSING
nt=1+100/dt;
figure;set(gcf, 'Position', get(0, 'Screensize'));
pcolor((0:1:nx-1).*dx,(0:1:ny-1).*dy,T(:,:,nt)');shading flat;axis equal;caxis([190, 330]);
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
