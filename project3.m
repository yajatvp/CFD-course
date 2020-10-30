clc
clear
%% Initialization
nx=201;
dx=1/(nx-1);
stability=0.2; % dt/dx
dt=stability*dx;
tmax=5;
nt=tmax/dt;
u=zeros(nx,nt+1); % For upwind
u1=zeros(nx,nt+1); % For central difference
u2=zeros(nx,nt+1); % For Lax
u3=zeros(nx,nt+1); % For Lax Wendroff
%% Initial condition - Discontinuity
u(1:1+0.2/dx,1:end)=1;
u(2+0.2/dx:end,1:end)=0;
u1(1:1+0.2/dx,1:end)=1;
u1(2+0.2/dx:end,1:end)=0;
u2(1:1+0.2/dx,1:end)=1;
u2(2+0.2/dx:end,1:end)=0;
u3(1:1+0.2/dx,1:end)=1;
u3(2+0.2/dx:end,1:end)=0;
%% FD schemes
for ntt=2:1:nt+1
    for x=2:1:nx-1
        u(x,ntt)=u(x,ntt-1)-(dt/dx)*0.5*((u(x,ntt-1))^2-u(x-1,ntt-1)^2); % Upwind 1st order
        
        u1(x,ntt)=u1(x,ntt-1)-(dt/dx)*0.25*((u1(x+1,ntt-1))^2-u1(x-1,ntt-1)^2);  % Central difference
        
        u2(x,ntt)=0.5*(u2(x+1,ntt-1)+u2(x-1,ntt-1))-(dt/dx)*0.25*((u2(x+1,ntt-1))^2-(u2(x-1,ntt-1))^2); % Lax
        
        
        u3(x,ntt)=u3(x,ntt-1)-((dt/dx)*0.25*((u3(x+1,ntt-1))^2-u3(x-1,ntt-1)^2))...
            +0.125*((dt/dx)^2)*(((u3(x+1,ntt-1)+u3(x,ntt-1))*((u3(x+1,ntt-1))^2-(u3(x,ntt-1))^2))  ...
            - ((u3(x,ntt-1)+u3(x-1,ntt-1))*((u3(x,ntt-1))^2-(u3(x-1,ntt-1))^2))); % Lax Wendroff
        
    end
end

%% POST PROCESSING

for ntt=1:0.2/dt:nt+1
    figure;set(gcf, 'Position', get(0, 'Screensize'));
    plot(0:dx:1,u(:,ntt),'LineWidth',2);grid on;hold on;ylim([-0.3 1.3]);
    % plot(0:dx:1,u1(:,ntt),'LineWidth',2); Central diffeerence, UNSTABLE
    plot(0:dx:1,u2(:,ntt),'LineWidth',2);plot(0:dx:1,u3(:,ntt),'LineWidth',2);
    xlabel('x');ylabel('U(x)');
    title(['U(x) at t =',num2str((ntt-1)*dt),' s']);
    legend('Upwind 1st order','Lax','Lax wendroff');
    w = waitforbuttonpress;
    close
end

%% POST PROCESSING - CENTRAL DIFFERENCE, UNSTABLE
%for ntt=1:10:nt+1
figure;set(gcf, 'Position', get(0, 'Screensize'));
plot(0:dx:1,u1(:,1),'LineWidth',2);grid on;hold on;
plot(0:dx:1,u1(:,11),'LineWidth',2);
plot(0:dx:1,u1(:,21),'LineWidth',2);
xlabel('x');ylabel('U(x)');
title(['U(x)']);
legend('t = 0 s','t = 0.025 s','t = 0.05 s');
%end
