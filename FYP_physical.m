%************************************************************************%
%                                                                        %
%   Final Year Project Code                                              %
%   For tests with physical data                                         %
%                                                                        %
%   Version: 1.1                                                         %
%   Author: Eirill Bachmann Mehammer                                     %
%                                                                        %
%************************************************************************%

close all;
clearvars;
clc;

% Define type of Kalman filter
kalman_method='time-varying';   %steady-state or time-varying

%% ------ Kinematic Model------- %%

% Define step size and number of steps
h=0.01/3; % step size
N=68617; % number of steps

% Generate time vector based on this
t = 0:h:(N-1)*h;

% Input matrices from measurements
[Pos,V,Acc]=read_physical(h,N);

u = zeros(3,N);
v = zeros(3,N);
w = zeros(6,N);
for n=1:3
    u(n,:)= Acc(n,:);
    v(n,:)= V(n,:);
    w(n,:)= Pos(n,:);
    w(n+3,:)= Pos(n+3,:);
end

% Initial conditions
x0 = [w(:,1); (w(1:3,2)-w(1:3,1))/h; v(:,1)];
y0 = [w(:,1); v(:,1); u(:,1)];

% Define arrays to hold the state and output values at each time step
x   = zeros(12,N);
y   = zeros(12,N);

% Set the first value of the state vector to be the initial condition
x(:,1) = x0;
y(:,1) = y0;

% Definition of matrices
%A=calcA(x(:,1));
A=zeros(12);
A(1,7)=1;
A(2,8)=1;
A(3,9)=1;
A(4,10)=1;
A(5,11)=1;
A(6,12)=1;

%B=calcB(x(:,1));
B=zeros(12,3);
B(7,1)=1;
B(8,2)=1;
B(9,3)=1;

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

D=zeros(12,3);
D(10,1)=1;
D(11,2)=1;
D(12,3)=1;

G=eye(12);      %'white Gaussian noise' means the covariance matrix is the identity
H=zeros(12);    %Assume process noise does not affect the output

% G=D;
% G(1,1)=0.01;
% G(2,2)=0.01;
% G(3,3)=0.005;
% G(4,1)=0.01;
% G(5,2)=0.01;
% G(6,3)=0.005;
% G(7,1)=0.01;
% G(8,2)=0.01;
% G(9,3)=0.005;
% 
% H=zeros(12,3);

%Check observability of continuous system
Qc=[C; C*A; C*A^2; C*A^3; C*A^4; C*A^5; C*A^6; C*A^7; C*A^8; C*A^9; C*A^10; C*A^11];
Qc_rank=rank(Qc);

%Define continuous state space system
sys=ss(A,[B G],C,[D H]);

%Discretize state space system
sysd = c2d(sys,0.05);        %0.05=step size from SIMA

%Check observability of discrete system
Qd=[sysd.C; sysd.C*sysd.A; sysd.C*sysd.A^2; sysd.C*sysd.A^3; sysd.C*sysd.A^4; sysd.C*sysd.A^5; sysd.C*sysd.A^6; sysd.C*sysd.A^7; sysd.C*sysd.A^8; sysd.C*sysd.A^9; sysd.C*sysd.A^10; sysd.C*sysd.A^11];
Qd_rank=rank(Qd);

%Differentiate positions to get velocities and accelerations
vel = zeros(6,size(x,2));
acc = zeros(6,size(x,2));
for i = 1:6
    vel(i,2:end) = diff(Pos(i,:))/h;
    acc(i,2:end)  = diff(vel(i,:))/h;
end

% % Plot differentiated velocities
% figure;
% plot(t,vel)
% legend('u','v','w','p','q','r')
% title('Differentiated velocities')
% xlabel('Time [s]')
% ylabel('Linear velocity [m/s], angular velocity [rad/s]')
% 
% %Plot difference between differentiated velocities and measured velocities
% figure
% plot(t,vel(4,:)-V(1,:),'r');
% hold on 
% plot(t,vel(5,:)-V(2,:),'r--');
% plot(t,vel(6,:)-V(3,:),'r:');
% hold off
% legend('p','q','r')
% title('Difference between differentiated velocities and measured velocities')
% xlabel('Time [s]')
% ylabel('Angular velocity [rad/s]')
% 
% % Plot differentiated accelerations
% figure;
% plot(t, acc)
% legend('a_x','a_y','a_z','p_{dot}','q_{dot}','r_{dot}')
% title('Differentiated accelerations')
% xlabel('Time [s]')
% ylabel('Linear acceleration [m/s^2], angular acceleration [rad/s^2]')
% 
% %Plot difference between differentiated accelerations and measured accelerations
% figure
% plot(t,acc(1,:)-Acc(1,:),'k');
% hold on 
% plot(t,acc(2,:)-Acc(2,:),'k--');
% plot(t,acc(3,:)-Acc(3,:),'k:');
% hold off
% legend('a_x','a_y','a_z')
% title('Difference between differentiated accelerations and measured accelerations')
% xlabel('Time [s]')
% ylabel('Linear acceleration [m/s^2]')

% Process noise vectors (white, Gaussian, zero-mean)
 uv=zeros(3,N);
 uv2=zeros(3,N);
 uv3=zeros(3,N);
 uv4=zeros(3,N);

 v1=zeros(12,N);
 for j=1:3
     uv(j,:) = awgn(u(j,:),100,'measured');        %SNR=100dB (0.00000001% noise)
     uv2(j,:) = awgn(u(j,:),100,'measured');
     uv3(j,:) = awgn(u(j,:),100,'measured');
     uv4(j,:) = awgn(u(j,:),20,'measured');        %SNR=8.24dB (15% noise)
 end
% xv = awgn(x,20,'measured');        %SNR=20dB (1% noise)
% xv = awgn(x,13.01,'measured');     %SNR=13.01dB (5% noise)
% xv = awgn(x,10,'measured');        %SNR=10dB (10% noise)
% xv = awgn(x,8.24,'measured');      %SNR=8.24dB (15% noise)
%v1(1:6,:) = (uv - acc);
%v1(7:12,:) = (uv - acc);
%v1=sqrt(10^-6)*randn(12,N)*h;

v1(1:3,:) = (uv - u);
v1(4:6,:) = (uv2 - u);
v1(7:9,:) = (uv3 - u);
v1(10:12,:) = (uv4 - u);

figure;
plot(t,v1);
title('Process noise')

% uv=zeros(3,N);
% v1=zeros(3,N);
%  for j=1:3
%      uv(j,:) = awgn(u(j,:),10,'measured');        %SNR=20dB (1% noise)
%  end
% v1(1:3,:) = (uv - u)*h;

% Loop to update the state and output values
for k=1:N-1 % need -1 here so that x(n+1) does not lead to an index>N
      x(:,k+1)=sysd.A*x(:,k)+sysd.B*[u(:,k); v1(:,k)];
      y(:,k+1)=sysd.C*x(:,k+1)+sysd.D*[u(:,k+1); v1(:,k+1)];
      % Calculate new version of A and B
      %A=calcA(x(:,k+1));
      %B=calcB(x(:,k+1));
      %Define new version of continuous state space system
      %sys=ss(A,[B G],C,[D H]);
      %Discretize new state space system
      %sysd = c2d(sys,0.05);        %0.05=step size from SIMA
end

% Plot outputs
figure;
plot(t,Pos(1,:),'k');
hold on 
plot(t,Pos(2,:),'k--');
plot(t,Pos(3,:),'k:');
plot(t,Pos(4,:),'r');
plot(t,Pos(5,:),'r--');
plot(t,Pos(6,:),'r:');
plot(t,V(1,:),'g');
plot(t,V(2,:),'g--');
plot(t,V(3,:),'g:');
plot(t,Acc(1,:),'b');
plot(t,Acc(2,:),'b--');
plot(t,Acc(3,:),'b:');
hold off
title('True outputs')
legend('x','y','z','phi','theta','psi','p','q','r','a_x','a_y','a_z')

figure;
plot(t,y(1,:),'k');
hold on 
plot(t,y(2,:),'k--');
plot(t,y(3,:),'k:');
plot(t,y(4,:),'r');
plot(t,y(5,:),'r--');
plot(t,y(6,:),'r:');
plot(t,y(7,:),'g');
plot(t,y(8,:),'g--');
plot(t,y(9,:),'g:');
plot(t,y(10,:),'b');
plot(t,y(11,:),'b--');
plot(t,y(12,:),'b:');
hold off
title('Outputs from MATLAB')
legend('x','y','z','phi','theta','psi','p','q','r','a_x','a_y','a_z')

%% Debug: plot states with index

figure
plot(t,x(1,:),'k');
hold on 
plot(t,x(2,:),'k--');
plot(t,x(3,:),'k:');
plot(t,x(4,:),'r');
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
plot(t,x(8,:),'k--');
plot(t,x(9,:),'k:');
plot(t,x(10,:),'r');
plot(t,x(11,:),'r--');
plot(t,x(12,:),'r:');
hold off
legend('u','v','w','p','q','r')
title('Velocities from MATLAB')
xlabel('Time [s]')
ylabel('Linear velocity [m/s], angular velocity [rad/s]')

%Plot difference between MATLAB and measured positions
figure
plot(t,x(1,:)-Pos(1,:),'k');
hold on 
plot(t,x(2,:)-Pos(2,:),'k--');
plot(t,x(3,:)-Pos(3,:),'k:');
plot(t,x(4,:)-Pos(4,:),'r');
plot(t,x(5,:)-Pos(5,:),'r--');
plot(t,x(6,:)-Pos(6,:),'r:');
hold off
legend('x','y','z','phi','theta','psi')
title('Difference in positions between MATLAB and measured')
xlabel('Time [s]')
ylabel('Translation [m], rotation [rad]')

%Plot difference between MATLAB and measured velocities
figure
plot(t,x(10,:)-V(1,:),'r');
hold on
plot(t,x(11,:)-V(2,:),'r--');
plot(t,x(12,:)-V(3,:),'r:');
hold off
legend('p','q','r')
title('Difference in velocities between MATLAB and measured')
xlabel('Time [s]')
ylabel('Angular velocity [rad/s]')

%Plot relative difference between MATLAB and measured positions
pos_rel_err=zeros(6,N);
vel_rel_err=zeros(3,N);
for i=1:3
    pos_rel_err(i,:) = deleteoutliers(abs((x(i,:)-Pos(i,:))./Pos(i,:)),0.05,1);
    pos_rel_err(i+3,:) = deleteoutliers(abs((x(i+3,:)-Pos(i+3,:))./Pos(i+3,:)),0.05,1);
    vel_rel_err(i,:) = deleteoutliers(abs((x(i+6,:)-V(i,:))./V(i,:)),0.05,1);
end
figure
plot(t,pos_rel_err);
legend('x','y','z','phi','theta','psi')
title('Relative difference between MATLAB and measured positions')
xlabel('Time [s]')
ylabel('Translation [m], rotation [rad]')

figure
plot(t,vel_rel_err);
legend('p','q','r')
title('Relative difference between MATLAB and measured velocities')
xlabel('Time [s]')
ylabel('Angular velocity [rad/s]')

%% ------ Kalman filter design------- %%
% Output noise vectors (white, Gaussian, zero-mean)
yv=zeros(12,N);
 for j=1:12
     yv(j,:) = awgn(y(j,:),20,'measured');        %SNR=100dB (0.00000001% noise)
 end
% yv = awgn(y,20,'measured');          %SNR=20dB (1% noise)
% yv = awgn(y,13.01,'measured');       %SNR=13.01dB (5% noise)
% yv = awgn(y,10,'measured');          %SNR=10dB (10% noise)
% yv = awgn(y,8.24,'measured');        %SNR=8.24dB (15% noise)
v2 = yv - y;
%v2=sqrt(10^-0)*randn(12,N);
% figure;
% plot(t,v2);
% title('Measurement noise')

%Check reachability
Reach=[sysd.B sysd.A*sysd.B sysd.A^2*sysd.B sysd.A^3*sysd.B sysd.A^4*sysd.B sysd.A^5*sysd.B sysd.A^6*sysd.B sysd.A^7*sysd.B sysd.A^8*sysd.B sysd.A^9*sysd.B sysd.A^10*sysd.B sysd.A^11*sysd.B];
R_rank=rank(Reach);

%Q=diag(var(transpose(v1))); 
%R=diag(var(transpose(v2))); 

Q=cov(transpose(v1));    % Process noise covariance matrix
R=cov(transpose(v2));    % Output noise covariance matrix
%Nn=cov(v1,v2);
Nn=zeros(12);        % Correlation between process and sensor noise is zero 

%Check positivity
e1=[sysd.B(:,4:15) zeros(12);sysd.D(:,4:15) eye(12)];
e2=[Q Nn;Nn' R];
e3=e1';
e=eig(e1*e2*e3);

% Define arrays to hold the Kalman state and output values at each time step
kalmx   = zeros(12,N);
kalmy   = zeros(12,N);

% Initial conditions
kalmx(:,1) = x0;
kalmy(:,1) = y0;

% Measurements
truey = [Pos;V;Acc];

% ----- Design the discrete Kalman filter ----- % 
switch kalman_method
    case 'steady-state'    %Steady-state Kalman filter 
        [kalmf,K,P] = kalman(sysd,Q,R,Nn);

        % Loop to update state and output values of Kalman filter
        for k=1:N-1
            kalmx(:,k+1)=sysd.A*kalmx(:,k)+sysd.B(:,1:3)*u(:,k)+K*(truey(:,k)-kalmy(:,k)); 
            kalmy(:,k+1)=sysd.C*kalmx(:,k+1)+sysd.D(:,1:3)*u(:,k+1);
        end
        
    case 'time-varying'  %Time-varying Kalman filter      
        % Initial condition of P
        P=var(x0);
        
        % Loop to update P, K, kalmx and kalmy at each time step
        for i = 2:N  
            % Check for signal loss
            if any(truey(1:6,i-1))
                R=cov(transpose(v2));    % Output noise covariance matrix
            else
                R=10^20*R;              % Increase value of R in case of signal loss
            end
            P = sysd.A*(P-P*sysd.C'/(R+sysd.C*P*sysd.C'))*sysd.A'+Q;
            K = P*sysd.C'/(sysd.C*P*sysd.C'+R);
            kalmx(:,i) = sysd.A*kalmx(:,i-1)+sysd.B(:,1:3)*u(:,i-1)+K*(truey(:,i-1)-kalmy(:,i-1));     % kalmx[n+1|n]
            kalmy(:,i) = sysd.C*kalmx(:,i)+sysd.D(:,1:3)*u(:,i);        % kalmy[n+1|n]
        end
        
        % Compare K & P from time-varying and steady-state
        [kalmf,Ks,Ps] = kalman(sysd,Q,R); 

    otherwise
        warning('Unexpected Kalman filtering method')
end

%-----Compare the filtered response kalmy with the measured response yv-------%
% Plot estimated output
figure;
plot(t,kalmy(1,:),'k');
hold on 
plot(t,kalmy(2,:),'k--');
plot(t,kalmy(3,:),'k:');
plot(t,kalmy(4,:),'r');
plot(t,kalmy(5,:),'r--');
plot(t,kalmy(6,:),'r:');
plot(t,kalmy(7,:),'g');
plot(t,kalmy(8,:),'g--');
plot(t,kalmy(9,:),'g:');
plot(t,kalmy(10,:),'b');
plot(t,kalmy(11,:),'b--');
plot(t,kalmy(12,:),'b:');
hold off
title('Kalman filter output')
legend('x','y','z','phi','theta','psi','p','q','r','a_x','a_y','a_z')

figure
plot(t,kalmx(1,:),'k');
hold on 
plot(t,kalmx(2,:),'k--');
plot(t,kalmx(3,:),'k:');
plot(t,kalmx(4,:),'r');
plot(t,kalmx(5,:),'r--');
plot(t,kalmx(6,:),'r:');
hold off
legend('x','y','z','phi','theta','psi')
title('Estimated positions')
xlabel('Time [s]')
ylabel('Translation [m], rotation [rad]')

figure
plot(t,kalmx(7,:),'k');
hold on 
plot(t,kalmx(8,:),'k--');
plot(t,kalmx(9,:),'k:');
plot(t,kalmx(10,:),'r');
plot(t,kalmx(11,:),'r--');
plot(t,kalmx(12,:),'r:');
hold off
legend('u','v','w','p','q','r')
title('Estimated velocities')
xlabel('Time [s]')
ylabel('Linear velocity [m/s], angular velocity [rad/s]')

kalmy_rel_err=zeros(12,N);
for i=1:1:12
    figure
    subplot(3,1,1), plot(t,y(i,:),t,kalmy(i,:),'--'),
    ylabel('Output')
    xlabel('No. of samples')
    legend('Measured output','Estimated output')
    subplot(3,1,2), plot(t,y(i,:)-kalmy(i,:),'--'),
    ylabel('Estimation error')
    xlabel('No. of samples')
    
    [kalmy_rel_err(i,:)] = deleteoutliers(abs((kalmy(i,:)-y(i,:))./y(i,:)),0.05,1);
    subplot(3,1,3), plot(t,kalmy_rel_err(i,:),'--');
    ylabel('Estimation relative error')
    xlabel('No. of samples')
    legend('Estimation relative error')
end

kalmx_rel_err=zeros(12,N);
for i=1:1:12
    figure
    subplot(3,1,1), plot(t,x(i,:),t,kalmx(i,:),'--'),
    ylabel('Output')
    xlabel('No. of samples')
    legend('Measured states','Estimated states')
    subplot(3,1,2), plot(t,x(i,:)-kalmx(i,:),'--'),
    ylabel('Estimation error')
    xlabel('No. of samples')
    
    [kalmx_rel_err(i,:)] = deleteoutliers(abs((kalmx(i,:)-x(i,:))./x(i,:)),0.05,1);
    subplot(3,1,3), plot(t,kalmx_rel_err(i,:),'--'),
    ylabel('Relative estimation error')
    xlabel('No. of samples')
end

%RMSD
%rmsd_measy=zeros(1,12);
rmsd_esty=zeros(1,12);
%rmsd_measx=zeros(1,12);
rmsd_estx=zeros(1,12);
for i=1:1:12
    %rmsd_measy(i) = sqrt(1/N*sum((y(i,:)-planty(i,:)-v2(i,:)).^2));
    rmsd_esty(i) = sqrt(1/N*sum((y(i,:)-kalmy(i,:)).^2));
    %rmsd_measx(i) = sqrt(1/N*sum(v1(i,:).^2));
    rmsd_estx(i) = sqrt(1/N*sum((x(i,:)-kalmx(i,:)).^2));
end

% Correlation coefficients
correlation_coefx=zeros(2,24);
correlation_coefy=zeros(2,24);
r_sqrx=zeros(1,12);
r_sqry=zeros(1,12);
for i=1:12
    correlation_coefx(:,2*i-1:2*i) = corrcoef(x(i,:),kalmx(i,:));
    r_sqrx(i) = power(correlation_coefx(1,2*i),2);
    correlation_coefy(:,2*i-1:2*i) = corrcoef(y(i,:),kalmy(i,:));
    r_sqry(i) = power(correlation_coefy(1,2*i),2);
end

% Error covariance before filtering (measurement error)
MeasErr = zeros(12,N);
MeasErrCov = zeros(1,12);
for i =1:12
    MeasErr(i,:) = v2(i,:);
    MeasErrCov(i) = sum(MeasErr(i,:).*MeasErr(i,:))/length(MeasErr(i,:));
end

% Plot measurement error covariance
figure;
plot(MeasErrCov)
hold on
title('Error covariance')

% Error covariance after filtering (estimation error)
EstErr = zeros(12,N);
EstErrCov = zeros(1,12);
for i =1:12
    EstErr(i,:) = y(i,:)-kalmy(i,:);
    EstErrCov(i) = sum(EstErr(i,:).*EstErr(i,:))/length(EstErr(i,:));
end

% Plot estimation error covariance
plot(EstErrCov)
legend('Measurement error covariance','Estimation error covariance')