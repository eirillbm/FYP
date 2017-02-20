% initial(sys,x)
close all;
clear all;
clc;

% Define step size and number of steps
h=0.05; % step size
N=10000; % number of steps

% Generate time vector based on this
t = 0:h:(N-1)*h;

% Simple input vector: sinusoid with the same number of element as the time vector
%u = zeros(6,length(t));     % empty array, all inputs = 0
%u(6,:) = u(6,:) + sin(t);   % set linear acceleration in x-direction to sin(t)

% Input matrices from SIMA
[P,V,Ac]=read4(h,N);

u = zeros(6,N);
v = zeros(6,N);
w = zeros(6,N);
for n=1:1:6
    u(n,:)= Ac(n,:);
    v(n,:)= V(n,:);
    w(n,:)= P(n,:);
end

% Simple initial condition of state vector (=0) (or define array)
x0 = zeros(12,1);

% Initial condition from SIMA
x0 = [w(:,1) ; v(:,1)];

% Define arrays to hold the state and output values at each time step
x   = zeros(12,N);
y   = zeros(12,N);

% Set the first value of the state vector to be the initial condition
x(:,1)  = x0;

% Definition of matrices
test_method='non-rot';
switch test_method
    case 'physical'  
        A=calcA2(x);
        B=zeros(12,3);
        B(7,1)=1;
        B(8,2)=1;
        B(9,3)=1;
        
        D=zeros(12,3);
        D(10,1)=1;
        D(11,2)=1;
        D(12,3)=1;
        
        E=D;
    case 'sima'
        A=calcA2(x);
        B=calcB(x);
        
        D=zeros(12,6);
        D(10,1)=1;
        D(11,2)=1;
        D(12,3)=1;
    case 'non-rot'
        A=zeros(12);
        A(1,7)=1;
        A(1,8)=0;
        A(1,9)=0;
        A(2,7)=0;
        A(2,8)=1;
        A(2,9)=0;
        A(3,7)=0;
        A(3,8)=0;
        A(3,9)=1;
        A(4,10)=1;
        A(4,11)=0;
        A(4,12)=0;
        A(5,11)=1;
        A(5,12)=0;
        A(6,11)=0;
        A(6,12)=1;
        
        B=zeros(12,6);
        B(7,1)=1;
        B(8,2)=1;
        B(9,3)=1;
        B(10,4)=1;
        B(11,5)=1;
        B(12,6)=1;
        
        D=zeros(12,6);
        D(10,1)=1;
        D(11,2)=1;
        D(12,3)=1;
end
C=zeros(12,12);
C(1,1)=1;
C(2,2)=1;
C(3,3)=1;
C(4,4)=1;
C(5,5)=1;
C(6,6)=1;
C(7,10)=1;
C(8,11)=1;
C(9,12)=1;

if strcmp(test_method,'non-rot')
    vel = zeros(1,size(x,2));
    acc = zeros(6,size(x,2));
    for i = 1:1:6
        vel(2:end) = diff(P(i,:))/h;
        acc(i,2:end)  = diff(vel)/h;
    end
end

% Loop to update the state and output values
int_method='euler';
for k=1:N-1 % need -1 here so that x(n+1) does not lead to an index>N
    %x(:,k+1)=x(:,k)+h*((A(:,:)*x(:,k)+B(:,:)*acc(:,k))); %one stage computation of Fwd Euler through matrix operations
    %y(:,k+1)=C(:,:)*x(:,k+1)+D(:,:)*acc(:,k+1);
    %if(1==0)
        for i = 1:1:12
            switch int_method
                case 'euler'
                %Euler method: 
                x(i,k+1)=x(i,k)+h*((A(i,:)*x(:,k)+B(i,:)*acc(:,k))); %one stage computation of Fwd Euler through matrix operations
                y(i,k+1)=C(i,:)*x(:,k+1)+D(i,:)*acc(:,k+1);

                case 'RK'
                %Runge-Kutta method:
                F_tx = @(x,acc) (A(i,:)*x(:,k)+B(i,:)*acc(:,k));
                k_1 = F_tx(x,acc);                            
                k_2 = F_tx(x+0.5*h,acc+0.5*h*k_1);
                k_3 = F_tx((x+0.5*h),(acc+0.5*h*k_2));
                k_4 = F_tx((x+h),(acc+k_3*h));

                x(i,k+1) = x(i,k) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  % main equation

                y(i,k+1)=C(i,:)*x(:,k+1)+D(i,:)*acc(:,k+1);       % y is an output based on x and u, so is not integrated

                t(k+1) = t(k)+h;
                                          
%                 case 'NB'
%                 %Newmark Beta method
%                 gamma=0.5;
%                 beta=0.25;
%                 udd=(1-gamma)*u;
%                 x(i,k+1)=x(i,k)+h*udd;

                otherwise
                    warning('Unexpected integration method')
            end
        end
    end
    % Calculate new version of A & B
    if strcmp(test_method,'non-rot')==0
        A=calcA2(x(:,k+1));
    end
    if strcmp(test_method,'sima')
        B=calcB(x(:,k+1));
    end
%end

%Wrap angle in radians to [-pi pi]
% x(4,:) = wrapToPi(x(4,:));
% x(5,:) = wrapToPi(x(5,:));
% x(6,:) = wrapToPi(x(6,:));
% y(4,:) = wrapToPi(y(4,:));
% y(5,:) = wrapToPi(y(5,:));
% y(6,:) = wrapToPi(y(6,:));

%Wrap angle in radians to [-pi/2 pi/2]
%  tmp = mod(x(4,:)+pi/2,pi);
%  x(4,:) = tmp+pi*(x(4,:)>0&tmp==0)-pi/2;
%  tmp = mod(y(4,:)+pi/2,pi);
%  y(4,:) = tmp+pi*(y(4,:)>0&tmp==0)-pi/2;

%Set roll angle to 0
% x(4,:) = zeros(1,N);
% y(4,:) = zeros(1,N);

% Plot output vector
figure;
plot(t,y)
title('Output')
legend('x','y','z','phi','theta','psi','p','q','r','a_x','a_y','a_z')

%% debug: plot states with index

figure
plot(t,x(1,:),'k');
hold on 
% for i = 2:1:6
plot(t,x(2,:),'k--');
plot(t,x(3,:),'k:');
plot(t,x(4,:),'r');
hold on 
% for i = 2:1:6
plot(t,x(5,:),'r--');
plot(t,x(6,:),'r:');
hold off

legend('x','y','z','phi','theta','psi')
title('Positions from MATLAB')
xlabel('Time [s]')
ylabel('Translation [m], rotation [rad]')

figure
plot(t,x(7,:),'k');
hold on 
% for i = 2:1:6
plot(t,x(8,:),'k--');
plot(t,x(9,:),'k:');
plot(t,x(10,:),'r');
hold on 
% for i = 2:1:6
plot(t,x(11,:),'r--');
plot(t,x(12,:),'r:');
hold off

legend('u','v','w','p','q','r')
title('Velocities from MATLAB')
xlabel('Time [s]')
ylabel('Linear velocity [m/s], angular velocity [rad/s]')

%Plot difference between MATLAB and SIMA positions
figure
plot(t,x(1,:)-P(1,:),'k');
hold on 
plot(t,x(2,:)-P(2,:),'k--');
plot(t,x(3,:)-P(3,:),'k:');
plot(t,x(4,:)-P(4,:),'r');
plot(t,x(5,:)-P(5,:),'r--');
plot(t,x(6,:)-P(6,:),'r:');
hold off

legend('x','y','z','phi','theta','psi')
title('Difference in positions between MATLAB and SIMA')
xlabel('Time [s]')
ylabel('Translation [m], rotation [rad]')

%Plot difference between MATLAB and SIMA velocities
figure
plot(t,x(7,:)-V(1,:),'k');
hold on 
plot(t,x(8,:)-V(2,:),'k--');
plot(t,x(9,:)-V(3,:),'k:');
plot(t,x(10,:)-V(4,:),'r');
plot(t,x(11,:)-V(5,:),'r--');
plot(t,x(12,:)-V(6,:),'r:');
hold off

legend('u','v','w','p','q','r')
title('Difference in velocities between MATLAB and SIMA')
xlabel('Time [s]')
ylabel('Linear velocity [m/s], angular velocity [rad/s]')

% Find maximum differences
% xmax = find(max(abs(x(1,:)-P(1,:))) == x(1,:)-P(1,:));
% txmax = t(xmax)
% yxmax = x(1,xmax)-P(1,xmax)
% 
% ymax = find(max(abs(x(2,:)-P(2,:))) == x(2,:)-P(2,:));
% tymax = t(ymax)
% yymax = x(2,ymax)-P(2,ymax)
% 
% zmax = find(max(abs(x(3,:)-P(3,:))) == x(3,:)-P(3,:));
% tzmax = t(zmax)
% yzmax = x(3,zmax)-P(3,zmax)
% 
% phimax = find(max(abs(x(4,:)-P(4,:))) == x(4,:)-P(4,:));
% tphimax = t(phimax)
% yphimax =  x(4,phimax)-P(4,phimax)
% 
% thetamax = find(max(abs(x(5,:)-P(5,:))) == x(5,:)-P(5,:));
% tthetamax = t(thetamax)
% ythetamax =  x(5,thetamax)-P(5,thetamax)
% 
% psimax = find(max(abs(x(6,:)-P(6,:))) == x(6,:)-P(6,:));
% tpsimax = t(psimax)
% ypsimax =  x(6,psimax)-P(6,psimax)
% 
% umax = find(max(abs(x(7,:)-V(1,:))) == x(7,:)-V(1,:));
% tumax = t(umax)
% yumax = x(7,umax)-V(1,umax)
% 
% vmax = find(max(abs(x(8,:)-V(2,:))) == x(8,:)-V(2,:));
% tvmax = t(vmax)
% yvmax = x(8,vmax)-V(2,vmax)
% 
% wmax = find(max(abs(x(9,:)-V(3,:))) == x(9,:)-V(3,:));
% twmax = t(wmax)
% ywmax = x(9,wmax)-V(3,wmax)
% 
% pmax = find(max(abs(x(10,:)-V(4,:))) == x(10,:)-V(4,:));
% tpmax = t(pmax)
% ypmax =  x(10,pmax)-V(4,pmax)
% 
% qmax = find(max(abs(x(11,:)-V(5,:))) == x(11,:)-V(5,:));
% tqmax = t(qmax)
% yqmax =  x(11,qmax)-V(5,ymax)
% 
% rmax = find(max(abs(x(12,:)-V(6,:))) == x(11,:)-V(6,:));
% trmax = t(rmax)
% yrmax =  x(12,rmax)-V(6,rmax)

% %Plot percentage difference between MATLAB and SIMA positions
figure
plot(t,(x(1,:)-P(1,:))./P(1,:),'k');
hold on 
plot(t,(x(2,:)-P(2,:))./P(2,:)*100,'k--');
plot(t,(x(3,:)-P(3,:))./P(3,:)*100,'k:');
plot(t,(x(4,:)-P(4,:))./P(4,:)*100,'r');
plot(t,(x(5,:)-P(5,:))./P(5,:)*100,'r--');
plot(t,(x(6,:)-P(6,:))./P(6,:)*100,'r:');
hold off

legend('x','y','z','phi','theta','psi')
title('Perentage difference in positions between MATLAB and SIMA')
xlabel('Time [s]')
ylabel('Translation [m], rotation [rad]')

% %Plot percentage difference between MATLAB and SIMA velocities
figure
plot(t,(x(7,:)-V(1,:))./V(1,:)*100,'k');
hold on 
plot(t,(x(8,:)-V(2,:))./V(2,:)*100,'k--');
plot(t,(x(9,:)-V(3,:))./V(3,:)*100,'k:');
plot(t,(x(10,:)-V(4,:))./V(4,:)*100,'r');
plot(t,(x(11,:)-V(5,:))./V(5,:)*100,'r--');
plot(t,(x(12,:)-V(6,:))./V(6,:)*100,'r:');
hold off

legend('u','v','w','p','q','r')
title('Percentage difference in velocities between MATLAB and SIMA')
xlabel('Time [s]')
ylabel('Linear velocity [m/s], angular velocity [rad/s]')