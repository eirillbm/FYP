%************************************************************************%
%                                                                        %
%   Final Year Project Code                                              %
%   For tests with virtual data                                          %
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
h=0.05; % step size
N=10000; % number of steps

% Generate time vector based on this
t = 0:h:(N-1)*h;

% Input matrices from SIMA
%[Pos,V,Acc]=read_virtual(h,N);
%save('SIMA_data.mat','Pos','V','Acc');
load('SIMA_data.mat');

u = zeros(6,N);
v = zeros(6,N);
w = zeros(6,N);
for n=1:1:6
    u(n,:)= Acc(n,:);
    v(n,:)= V(n,:);
    w(n,:)= Pos(n,:);
end

% Initial condition from SIMA
x0 = [w(:,1) ; v(:,1)];
y0 = [w(:,1) ; v(4,1); v(5,1); v(6,1); u(1,1); u(2,1); u(3,1)];

% Define arrays to hold the state and output values at each time step
x   = zeros(12,N);
y   = zeros(12,N);

% Set the first value of the state vector to be the initial condition
x(:,1)  = x0;
y(:,1) = y0;

% Definition of matrices
A=zeros(12);
A(1,7)=1;
A(2,8)=1;
A(3,9)=1;
A(4,10)=1;
A(5,11)=1;
A(6,12)=1;

B=zeros(12,6);
B(7,1)=1;
B(8,2)=1;
B(9,3)=1;
B(10,4)=1;
B(11,5)=1;
B(12,6)=1;

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

D=zeros(12,6);
D(10,1)=1;
D(11,2)=1;
D(12,3)=1;

%Check observability of continuous system
Qc=[C; C*A; C*A^2; C*A^3; C*A^4; C*A^5; C*A^6; C*A^7; C*A^8; C*A^9; C*A^10; C*A^11];
Qc_rank=rank(Qc);

%Define continuous state space system
sys=ss(A,B,C,D);

%Discretize state space system
sysd = c2d(sys,0.05);        %0.05=step size from SIMA

%Check observability of discrete system
Qd=[sysd.C; sysd.C*sysd.A; sysd.C*sysd.A^2; sysd.C*sysd.A^3; sysd.C*sysd.A^4; sysd.C*sysd.A^5; sysd.C*sysd.A^6; sysd.C*sysd.A^7; sysd.C*sysd.A^8; sysd.C*sysd.A^9; sysd.C*sysd.A^10; sysd.C*sysd.A^11];
Qd_rank=rank(Qd);

%Differentiate positions to get velocities and accelerations
vel = zeros(6,size(x,2));
acc = zeros(6,size(x,2));
for i = 1:1:6
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
% %Plot difference between differentiated velocities and SIMA velocities
% figure
% plot(t,vel(1,:)-V(1,:),'k');
% hold on 
% plot(t,vel(2,:)-V(2,:),'k--');
% plot(t,vel(3,:)-V(3,:),'k:');
% plot(t,vel(4,:)-V(4,:),'r');
% plot(t,vel(5,:)-V(5,:),'r--');
% plot(t,vel(6,:)-V(6,:),'r:');
% hold off
% legend('u','v','w','p','q','r')
% title('Difference between differentiated velocities and SIMA velocities')
% xlabel('Time [s]')
% ylabel('Linear velocity [m/s], angular velocity [rad/s]')
% 
% % Plot differentiated accelerations
% figure;
% plot(t, acc)
% legend('a_x','a_y','a_z','p_{dot}','q_{dot}','r_{dot}')
% title('Differentiated accelerations')
% xlabel('Time [s]')
% ylabel('Linear acceleration [m/s^2], angular acceleration [rad/s^2]')
% 
% %Plot difference between differentiated accelerations and SIMA accelerations
% figure
% plot(t,acc(1,:)-Acc(1,:),'k');
% hold on 
% plot(t,acc(2,:)-Acc(2,:),'k--');
% plot(t,acc(3,:)-Acc(3,:),'k:');
% plot(t,acc(4,:)-Acc(4,:),'r');
% plot(t,acc(5,:)-Acc(5,:),'r--');
% plot(t,acc(6,:)-Acc(6,:),'r:');
% hold off
% legend('a_x','a_y','a_z','p_{dot}','q_{dot}','r_{dot}')
% title('Difference between differentiated accelerations and SIMA accelerations')
% xlabel('Time [s]')
% ylabel('Linear acceleration [m/s^2], angular acceleration [rad/s^2]')

% Loop to update the state and output values
for k=1:N-1 % need -1 here so that x(n+1) does not lead to an index>N
      x(:,k+1)=sysd.A(:,:)*x(:,k)+sysd.B(:,:)*acc(:,k); %one stage computation of Fwd Euler through matrix operations
      y(:,k+1)=sysd.C(:,:)*x(:,k+1)+sysd.D(:,:)*acc(:,k+1);
end

% Plot true output
% figure;
% plot(t,y(1,:),'k');
% hold on 
% plot(t,y(2,:),'k--');
% plot(t,y(3,:),'k:');
% plot(t,y(4,:),'r');
% plot(t,y(5,:),'r--');
% plot(t,y(6,:),'r:');
% plot(t,y(7,:),'g');
% plot(t,y(8,:),'g--');
% plot(t,y(9,:),'g:');
% plot(t,y(10,:),'b');
% plot(t,y(11,:),'b--');
% plot(t,y(12,:),'b:');
% hold off
% title('True output')
% legend('x','y','z','phi','theta','psi','p','q','r','a_x','a_y','a_z')

%% Debug: plot states with index

% figure
% plot(t,x(1,:),'k');
% hold on 
% plot(t,x(2,:),'k--');
% plot(t,x(3,:),'k:');
% plot(t,x(4,:),'r');
% plot(t,x(5,:),'r--');
% plot(t,x(6,:),'r:');
% hold off
% legend('x','y','z','phi','theta','psi')
% title('Positions from MATLAB')
% xlabel('Time [s]')
% ylabel('Translation [m], rotation [rad]')
% 
% figure
% plot(t,x(7,:),'k');
% hold on 
% plot(t,x(8,:),'k--');
% plot(t,x(9,:),'k:');
% plot(t,x(10,:),'r');
% plot(t,x(11,:),'r--');
% plot(t,x(12,:),'r:');
% hold off
% legend('u','v','w','p','q','r')
% title('Velocities from MATLAB')
% xlabel('Time [s]')
% ylabel('Linear velocity [m/s], angular velocity [rad/s]')
% 
% %Plot difference between MATLAB and SIMA positions
% figure
% plot(t,x(1,:)-Pos(1,:),'k');
% hold on 
% plot(t,x(2,:)-Pos(2,:),'k--');
% plot(t,x(3,:)-Pos(3,:),'k:');
% plot(t,x(4,:)-Pos(4,:),'r');
% plot(t,x(5,:)-Pos(5,:),'r--');
% plot(t,x(6,:)-Pos(6,:),'r:');
% hold off
% legend('x','y','z','phi','theta','psi')
% title('Difference in positions between MATLAB and SIMA')
% xlabel('Time [s]')
% ylabel('Translation [m], rotation [rad]')
% 
% %Plot difference between MATLAB and SIMA velocities
% figure
% plot(t,x(7,:)-V(1,:),'k');
% hold on 
% plot(t,x(8,:)-V(2,:),'k--');
% plot(t,x(9,:)-V(3,:),'k:');
% plot(t,x(10,:)-V(4,:),'r');
% plot(t,x(11,:)-V(5,:),'r--');
% plot(t,x(12,:)-V(6,:),'r:');
% hold off
% legend('u','v','w','p','q','r')
% title('Difference in velocities between MATLAB and SIMA')
% xlabel('Time [s]')
% ylabel('Linear velocity [m/s], angular velocity [rad/s]')
% 
% 
% %Plot relative difference between MATLAB and SIMA positions and velocities
% pos_rel_err=zeros(6,N);
% vel_rel_err=zeros(6,N);
% for i=1:6
%     pos_rel_err(i,:) = deleteoutliers(abs((x(i,:)-Pos(i,:))./Pos(i,:)),0.05,1);
%     vel_rel_err(i,:) = deleteoutliers(abs((x(i+6,:)-V(i,:))./V(i,:)),0.05,1);
% end
% figure
% plot(t,pos_rel_err);
% legend('x','y','z','phi','theta','psi')
% title('Relative difference in positions between MATLAB and SIMA')
% xlabel('Time [s]')
% ylabel('Translation [m], rotation [rad]')
% 
% figure
% plot(t,vel_rel_err);
% legend('u','v','w','p','q','r')
% title('Relative difference in velocities between MATLAB and SIMA')
% xlabel('Time [s]')
% ylabel('Linear velocity [m/s], angular velocity [rad/s]')

%% ------ Kalman filter design------- %%
% Process noise vectors (white, Gaussian, zero-mean)
 uv=zeros(6,N);
 uv2=zeros(6,N);

 v1=zeros(12,N);
 for j=1:6
     uv(j,:) = awgn(acc(j,:),100,'measured');        %SNR=10dB (10% noise)
     uv2(j,:) = awgn(acc(j,:),100,'measured');      %SNR=100dB (0.00000001% noise)
 end
%SNR=20dB (1% noise), SNR=16.99dB (2% noise), SNR=13.01dB (5% noise), SNR=8.24dB (15% noise)
v1(1:6,:) = (uv - acc);
v1(7:12,:) = (uv2 - acc);
%v1=sqrt(10^-6)*randn(12,N)*h;
%v1(1:6,:) = (uv - acc)*h;
%v1(7:12,:) = (uv2 - acc)*h;

figure;
plot(t,v1);
title('Process noise')

G=eye(12);      %'white Gaussian noise' means the covariance matrix is the identity
H=zeros(12);    %Assume process noise does not affect the output

%figure;
%plot(t,xv);
%title('States with process noise')

% Output noise vectors (white, Gaussian, zero-mean)
yv=zeros(12,N);
 for j=1:12
     yv(j,:) = awgn(y(j,:),100,'measured');       %SNR=10dB (10% noise) 
 end
%SNR=100dB (0.00000001% noise), SNR=20dB (1% noise), SNR=13.01dB (5% noise), SNR=8.24dB (15% noise)
v2 = (yv - y);
%v2=sqrt(10^-100)*randn(12,N);
figure;
plot(t,v2);
title('Measurement noise')

%Check reachability
Reach=[sysd.B sysd.A*sysd.B sysd.A^2*sysd.B sysd.A^3*sysd.B sysd.A^4*sysd.B sysd.A^5*sysd.B sysd.A^6*sysd.B sysd.A^7*sysd.B sysd.A^8*sysd.B sysd.A^9*sysd.B sysd.A^10*sysd.B sysd.A^11*sysd.B];
R_rank=rank(Reach);

% Specify the plant model. The sample time is set to -1, to mark the model as discrete without specifying a sample time.
plantx = zeros(12,N);       % Define arrays to hold the plant state and output values at each time step
planty = zeros(12,N);
plantx(:,1) = x0;           % Set the first value of the plant state vector to be the initial condition from SIMA
planty(:,1) = y0;

Plant = ss(A,[B G],C,[D H]);

%Discretize noisy state space system
Plantd = c2d(Plant,0.05);        %0.05=step size from SIMA

% Generate plant states and outputs
for k=1:N-1
    plantx(:,k+1)=Plantd.A*plantx(:,k)+Plantd.B*[acc(:,k); v1(:,k)]; 
    planty(:,k+1)=Plantd.C*plantx(:,k+1)+Plantd.D*[acc(:,k+1); v1(:,k+1)]+v2(:,k+1);
end

%Q=diag(var(transpose(v1))); 
%R=diag(var(transpose(v2))); 

Q=cov(transpose(v1));    % Process noise covariance matrix
R=cov(transpose(v2));    % Output noise covariance matrix
%Nn=cov(v1,v2);        % Correlation between process and sensor noise 
Nn=zeros(12);

%Check positivity
e1=[Plantd.B(:,7:18) zeros(12);Plantd.D(:,7:18) eye(12)];
e2=[Q Nn;Nn' R];
e3=e1';
e=eig(e1*e2*e3);

% Define arrays to hold the Kalman state and output values at each time step
kalmx   = zeros(12,N);
kalmy   = zeros(12,N);

% Initial conditions
kalmx(:,1) = x0;
kalmy(:,1) = y0;

%Simulate signal loss
planty(1,2000)=0;
planty(1,4000)=0;
planty(1,6000)=0;
planty(1,8000)=0;
planty(2,2000)=0;
planty(2,4000)=0;
planty(2,6000)=0;
planty(2,8000)=0;
planty(3,2000)=0;
planty(3,4000)=0;
planty(3,6000)=0;
planty(3,8000)=0;
planty(4,2000)=0;
planty(4,4000)=0;
planty(4,6000)=0;
planty(4,8000)=0;
planty(5,2000)=0;
planty(5,4000)=0;
planty(5,6000)=0;
planty(5,8000)=0;
planty(6,2000)=0;
planty(6,4000)=0;
planty(6,6000)=0;
planty(6,8000)=0;
planty(1,2001)=0;
planty(1,4001)=0;
planty(1,6001)=0;
planty(1,8001)=0;
planty(2,2001)=0;
planty(2,4001)=0;
planty(2,6001)=0;
planty(2,8001)=0;
planty(3,2001)=0;
planty(3,4001)=0;
planty(3,6001)=0;
planty(3,8001)=0;
planty(4,2001)=0;
planty(4,4001)=0;
planty(4,6001)=0;
planty(4,8001)=0;
planty(5,2001)=0;
planty(5,4001)=0;
planty(5,6001)=0;
planty(5,8001)=0;
planty(6,2001)=0;
planty(6,4001)=0;
planty(6,6001)=0;
planty(6,8001)=0;
planty(1,2002)=0;
planty(1,4002)=0;
planty(1,6002)=0;
planty(1,8002)=0;
planty(2,2002)=0;
planty(2,4002)=0;
planty(2,6002)=0;
planty(2,8002)=0;
planty(3,2002)=0;
planty(3,4002)=0;
planty(3,6002)=0;
planty(3,8002)=0;
planty(4,2002)=0;
planty(4,4002)=0;
planty(4,6002)=0;
planty(4,8002)=0;
planty(5,2002)=0;
planty(5,4002)=0;
planty(5,6002)=0;
planty(5,8002)=0;
planty(6,2002)=0;
planty(6,4002)=0;
planty(6,6002)=0;
planty(6,8002)=0;
planty(1,2003)=0;
planty(1,4003)=0;
planty(1,6003)=0;
planty(1,8003)=0;
planty(2,2003)=0;
planty(2,4003)=0;
planty(2,6003)=0;
planty(2,8003)=0;
planty(3,2003)=0;
planty(3,4003)=0;
planty(3,6003)=0;
planty(3,8003)=0;
planty(4,2003)=0;
planty(4,4003)=0;
planty(4,6003)=0;
planty(4,8003)=0;
planty(5,2003)=0;
planty(5,4003)=0;
planty(5,6003)=0;
planty(5,8003)=0;
planty(6,2003)=0;
planty(6,4003)=0;
planty(6,6003)=0;
planty(6,8003)=0;
planty(1,2004)=0;
planty(1,4004)=0;
planty(1,6004)=0;
planty(1,8004)=0;
planty(2,2004)=0;
planty(2,4004)=0;
planty(2,6004)=0;
planty(2,8004)=0;
planty(3,2004)=0;
planty(3,4004)=0;
planty(3,6004)=0;
planty(3,8004)=0;
planty(4,2004)=0;
planty(4,4004)=0;
planty(4,6004)=0;
planty(4,8004)=0;
planty(5,2004)=0;
planty(5,4004)=0;
planty(5,6004)=0;
planty(5,8004)=0;
planty(6,2004)=0;
planty(6,4004)=0;
planty(6,6004)=0;
planty(6,8004)=0;
planty(1,2005)=0;
planty(1,4005)=0;
planty(1,6005)=0;
planty(1,8005)=0;
planty(2,2005)=0;
planty(2,4005)=0;
planty(2,6005)=0;
planty(2,8005)=0;
planty(3,2005)=0;
planty(3,4005)=0;
planty(3,6005)=0;
planty(3,8005)=0;
planty(4,2005)=0;
planty(4,4005)=0;
planty(4,6005)=0;
planty(4,8005)=0;
planty(5,2005)=0;
planty(5,4005)=0;
planty(5,6005)=0;
planty(5,8005)=0;
planty(6,2005)=0;
planty(6,4005)=0;
planty(6,6005)=0;
planty(6,8005)=0;
planty(1,2006)=0;
planty(1,4006)=0;
planty(1,6006)=0;
planty(1,8006)=0;
planty(2,2006)=0;
planty(2,4006)=0;
planty(2,6006)=0;
planty(2,8006)=0;
planty(3,2006)=0;
planty(3,4006)=0;
planty(3,6006)=0;
planty(3,8006)=0;
planty(4,2006)=0;
planty(4,4006)=0;
planty(4,6006)=0;
planty(4,8006)=0;
planty(5,2006)=0;
planty(5,4006)=0;
planty(5,6006)=0;
planty(5,8006)=0;
planty(6,2006)=0;
planty(6,4006)=0;
planty(6,6006)=0;
planty(6,8006)=0;
planty(1,2007)=0;
planty(1,4007)=0;
planty(1,6007)=0;
planty(1,8007)=0;
planty(2,2007)=0;
planty(2,4007)=0;
planty(2,6007)=0;
planty(2,8007)=0;
planty(3,2007)=0;
planty(3,4007)=0;
planty(3,6007)=0;
planty(3,8007)=0;
planty(4,2007)=0;
planty(4,4007)=0;
planty(4,6007)=0;
planty(4,8007)=0;
planty(5,2007)=0;
planty(5,4007)=0;
planty(5,6007)=0;
planty(5,8007)=0;
planty(6,2007)=0;
planty(6,4007)=0;
planty(6,6007)=0;
planty(6,8007)=0;
planty(1,2008)=0;
planty(1,4008)=0;
planty(1,6008)=0;
planty(1,8008)=0;
planty(2,2008)=0;
planty(2,4008)=0;
planty(2,6008)=0;
planty(2,8008)=0;
planty(3,2008)=0;
planty(3,4008)=0;
planty(3,6008)=0;
planty(3,8008)=0;
planty(4,2008)=0;
planty(4,4008)=0;
planty(4,6008)=0;
planty(4,8008)=0;
planty(5,2008)=0;
planty(5,4008)=0;
planty(5,6008)=0;
planty(5,8008)=0;
planty(6,2008)=0;
planty(6,4008)=0;
planty(6,6008)=0;
planty(6,8008)=0;
planty(1,2009)=0;
planty(1,4009)=0;
planty(1,6009)=0;
planty(1,8009)=0;
planty(2,2009)=0;
planty(2,4009)=0;
planty(2,6009)=0;
planty(2,8009)=0;
planty(3,2009)=0;
planty(3,4009)=0;
planty(3,6009)=0;
planty(3,8009)=0;
planty(4,2009)=0;
planty(4,4009)=0;
planty(4,6009)=0;
planty(4,8009)=0;
planty(5,2009)=0;
planty(5,4009)=0;
planty(5,6009)=0;
planty(5,8009)=0;
planty(6,2009)=0;
planty(6,4009)=0;
planty(6,6009)=0;
planty(6,8009)=0;

% ----- Design the discrete Kalman filter ----- %
switch kalman_method
    case 'steady-state'    %Steady-state Kalman filter 
        [kalmf,K,P] = kalman(Plantd,Q,R,Nn);
        
        % Loop to update state and output values of Kalman filter
        for k=1:N-1
            kalmx(:,k+1)=sysd.A*kalmx(:,k)+sysd.B*acc(:,k)+K*(planty(:,k)-kalmy(:,k)); 
            kalmy(:,k+1)=sysd.C*kalmx(:,k+1)+sysd.D*acc(:,k+1);
        end
        
    case 'time-varying'  %Time-varying Kalman filter      
        % Initial condition of P
        P=var(x0);
        
        % Loop to update P, K, kalmx and kalmy at each time step
        for i = 2:N    
            % Check for signal loss
            if any(planty(1:6,i-1))
                R=cov(transpose(v2));    % Output noise covariance matrix
            else
            %if planty(1,i)==planty(2,i)==planty(3,i)==planty(4,i)==planty(5,i)==planty(6,i)==0
                R=10^20*R;              % Increase value of R in case of signal loss
            end
            P = sysd.A*(P-P*sysd.C'/(R+sysd.C*P*sysd.C'))*sysd.A'+Q;
            K = P*sysd.C'/(sysd.C*P*sysd.C'+R);
            kalmx(:,i) = sysd.A*kalmx(:,i-1)+sysd.B*acc(:,i-1)+K*(planty(:,i-1)-kalmy(:,i-1));     % kalmx[n+1|n]
            kalmy(:,i) = sysd.C*kalmx(:,i)+sysd.D*acc(:,i);        % kalmy[n+1|n]
        end
        
        % Compare K & P from time-varying and steady-state
        [kalmf,Ks,Ps] = kalman(Plant,Q,R); 

    otherwise
        warning('Unexpected Kalman filtering method')
end
        
%-----Compare the filtered response kalmy with the measured response yv-------%
% Plot measured output
figure;
plot(t,planty(1,:),'k');
hold on 
plot(t,planty(2,:),'k--');
plot(t,planty(3,:),'k:');
plot(t,planty(4,:),'r');
plot(t,planty(5,:),'r--');
plot(t,planty(6,:),'r:');
plot(t,planty(7,:),'g');
plot(t,planty(8,:),'g--');
plot(t,planty(9,:),'g:');
plot(t,planty(10,:),'b');
plot(t,planty(11,:),'b--');
plot(t,planty(12,:),'b:');
hold off
title('Measured output')
legend('x','y','z','phi','theta','psi','p','q','r','a_x','a_y','a_z')

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

% Plot measured states
figure;
plot(t,plantx(1,:),'k');
hold on 
plot(t,plantx(2,:),'k--');
plot(t,plantx(3,:),'k:');
plot(t,plantx(4,:),'r');
plot(t,plantx(5,:),'r--');
plot(t,plantx(6,:),'r:');
plot(t,plantx(7,:),'g');
plot(t,plantx(8,:),'g--');
plot(t,plantx(9,:),'g:');
plot(t,plantx(10,:),'b');
plot(t,plantx(11,:),'b--');
plot(t,plantx(12,:),'b:');
hold off
title('Measured states')
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

% % Compare the true and filtered responses graphically.
% figure
% title('Kalman filter response (output)')
% subplot(6,2,1), plot(t,y(1,:),t,kalmy(1,:),'--'),
% ylabel('Out1')
% subplot(6,2,2), plot(t,y(2,:),t,kalmy(2,:),'--'),
% ylabel('Out2')
% subplot(6,2,3), plot(t,y(3,:),t,kalmy(3,:),'--'),
% ylabel('Out3')
% subplot(6,2,4), plot(t,y(4,:),t,kalmy(4,:),'--'),
% ylabel('Out4')
% subplot(6,2,5), plot(t,y(5,:),t,kalmy(5,:),'--'),
% ylabel('Out5')
% subplot(6,2,6), plot(t,y(6,:),t,kalmy(6,:),'--'),
% ylabel('Out6')
% subplot(6,2,7), plot(t,y(7,:),t,kalmy(7,:),'--'),
% ylabel('Out7')
% subplot(6,2,8), plot(t,y(8,:),t,kalmy(8,:),'--'),
% ylabel('Out8')
% subplot(6,2,9), plot(t,y(9,:),t,kalmy(9,:),'--'),
% ylabel('Out9')
% subplot(6,2,10), plot(t,y(10,:),t,kalmy(10,:),'--'),
% ylabel('Out10')
% subplot(6,2,11), plot(t,y(11,:),t,kalmy(11,:),'--'),
% ylabel('Out11')
% subplot(6,2,12), plot(t,y(12,:),t,kalmy(12,:),'--'),
% ylabel('Out12')
% 
% figure
% title('Kalman filter response (error)')
% subplot(6,2,1), plot(t,y(1,:)-yv(1,:),t,y(1,:)-kalmy(1,:),'--'),
% ylabel('Err1')
% subplot(6,2,2), plot(t,y(2,:)-yv(2,:),t,y(2,:)-kalmy(2,:),'--'),
% ylabel('Err2')
% subplot(6,2,3), plot(t,y(3,:)-yv(3,:),t,y(3,:)-kalmy(3,:),'--'),
% ylabel('Err3')
% subplot(6,2,4), plot(t,y(4,:)-yv(4,:),t,y(4,:)-kalmy(4,:),'--'),
% ylabel('Err4')
% subplot(6,2,5), plot(t,y(5,:)-yv(5,:),t,y(5,:)-kalmy(5,:),'--'),
% ylabel('Err5')
% subplot(6,2,6), plot(t,y(6,:)-yv(6,:),t,y(6,:)-kalmy(6,:),'--'),
% ylabel('Err6')
% subplot(6,2,7), plot(t,y(7,:)-yv(7,:),t,y(7,:)-kalmy(7,:),'--'),
% ylabel('Err7')
% subplot(6,2,8), plot(t,y(8,:)-yv(8,:),t,y(8,:)-kalmy(8,:),'--'),
% ylabel('Err8')
% subplot(6,2,9), plot(t,y(9,:)-yv(9,:),t,y(9,:)-kalmy(9,:),'--'),
% ylabel('Err9')
% subplot(6,2,10), plot(t,y(10,:)-yv(10,:),t,y(10,:)-kalmy(10,:),'--'),
% ylabel('Err10')
% subplot(6,2,11), plot(t,y(11,:)-yv(11,:),t,y(11,:)-kalmy(11,:),'--'),
% ylabel('Err11')
% subplot(6,2,12), plot(t,y(12,:)-yv(12,:),t,y(12,:)-kalmy(12,:),'--'),
% ylabel('Err12')
% 
% figure
% title('Percentage error')
% subplot(6,2,1), plot(t,abs((yv(1,:)-y(1,:))./y(1,:))*100,t,abs((kalmy(1,:)-y(1,:))./y(1,:))*100,'--'),
% ylabel('%Err1')
% subplot(6,2,2), plot(t,abs((yv(2,:)-y(2,:))./y(2,:))*100,t,abs((kalmy(2,:)-y(2,:))./y(2,:))*100,'--'),
% ylabel('%Err2')
% subplot(6,2,3), plot(t,abs((yv(3,:)-y(3,:))./y(3,:))*100,t,abs((kalmy(3,:)-y(3,:))./y(3,:))*100,'--'),
% ylabel('%Err3')
% subplot(6,2,4), plot(t,abs((yv(4,:)-y(4,:))./y(4,:))*100,t,abs((kalmy(4,:)-y(4,:))./y(4,:))*100,'--'),
% ylabel('%Err4')
% subplot(6,2,5), plot(t,abs((yv(5,:)-y(5,:))./y(5,:))*100,t,abs((kalmy(5,:)-y(5,:))./y(5,:))*100,'--'),
% ylabel('%Err5')
% subplot(6,2,6), plot(t,abs((yv(6,:)-y(6,:))./y(6,:))*100,t,abs((kalmy(6,:)-y(6,:))./y(6,:))*100,'--'),
% ylabel('%Err6')
% subplot(6,2,7), plot(t,abs((yv(7,:)-y(7,:))./y(7,:))*100,t,abs((kalmy(7,:)-y(7,:))./y(7,:))*100,'--'),
% ylabel('%Err7')
% subplot(6,2,8), plot(t,abs((yv(8,:)-y(8,:))./y(8,:))*100,t,abs((kalmy(8,:)-y(8,:))./y(8,:))*100,'--'),
% ylabel('%Err8')
% subplot(6,2,9), plot(t,abs((yv(9,:)-y(9,:))./y(9,:))*100,t,abs((kalmy(9,:)-y(9,:))./y(9,:))*100,'--'),
% ylabel('%Err9')
% subplot(6,2,10), plot(t,abs((yv(10,:)-y(10,:))./y(10,:))*100,t,abs((kalmy(10,:)-y(10,:))./y(10,:))*100,'--'),
% ylabel('%Err10')
% subplot(6,2,11), plot(t,abs((yv(11,:)-y(11,:))./y(11,:))*100,t,abs((kalmy(11,:)-y(11,:))./y(11,:))*100,'--'),
% ylabel('%Err11')
% subplot(6,2,12), plot(t,abs((yv(12,:)-y(12,:))./y(12,:))*100,t,abs((kalmy(12,:)-y(12,:))./y(12,:))*100,'--'),
% ylabel('%Err12')
% 

planty_rel_err=zeros(12,N);
kalmy_rel_err=zeros(12,N);

for i=1:1:12
    figure
    subplot(3,1,1), plot(t,y(i,:),t,kalmy(i,:),'--', t, planty(i,:),':'),
    ylabel('Output')
    xlabel('Time [s]')
    legend('Actual output','Estimated output','Measured output')
    subplot(3,1,2), plot(t,y(i,:)-planty(i,:),t,y(i,:)-kalmy(i,:),'--'),
    ylabel('Absoulte error [m]')
    xlabel('Time [s]')
    legend('Measurement error','Estimation error')
    
    [planty_rel_err(i,:)] = deleteoutliers(abs((planty(i,:)-y(i,:))./y(i,:)),0.05,1);
    [kalmy_rel_err(i,:)] = deleteoutliers(abs((kalmy(i,:)-y(i,:))./y(i,:)),0.05,1);
    subplot(3,1,3), plot(t,planty_rel_err(i,:),t,kalmy_rel_err(i,:),'--'),
    ylabel('Relative error')
    xlabel('Time [s]')
    legend('Measurement relative error','Estimation relative error')
end

plantx_rel_err=zeros(12,N);
kalmx_rel_err=zeros(12,N);
for i=1:1:12
    figure
    subplot(3,1,1), plot(t,x(i,:),t,kalmx(i,:),'--', t, x(i,:)+v1(i,:),':'),
    ylabel('Output')
    xlabel('Time [s]')
    legend('Actual states','Estimated states','Measured states')
    subplot(3,1,2), plot(t,v1(i,:),t,x(i,:)-kalmx(i,:),'--'),
    ylabel('Absolute error [m]')
    xlabel('Time [s]')
    legend('Measurement error','Estimation error')
    
    [plantx_rel_err(i,:)] = deleteoutliers(abs(v1(i,:)./x(i,:)),0.05,1);
    [kalmx_rel_err(i,:)] = deleteoutliers(abs((kalmx(i,:)-x(i,:))./x(i,:)),0.05,1);
    subplot(3,1,3), plot(t,plantx_rel_err(i,:),t,kalmx_rel_err(i,:),'--'),
    ylabel('Relative error')
    xlabel('Time [s]')
    legend('Measurement relative error','Estimation relative error')
end

% for i=1:1:12
%     figure
%     subplot(3,1,1), plot(t,x(i,:),t,kalmy(i+12,:),'--', t, x(i,:)+v1(i,:),':'),
%     ylabel('Output')
%     xlabel('No. of samples')
%     legend('Actual states','Estimated states','Measured states')
%     subplot(3,1,2), plot(t,v1(i,:),t,x(i,:)-kalmy(i+12,:),'--'),
%     ylabel('Error')
%     xlabel('No. of samples')
%     legend('Measurement error','Estimation error')
%     subplot(3,1,3), plot(t,abs(v1(i,:)./x(i,:))*100,t,abs((kalmy(i+12,:)-x(i,:))./x(i,:))*100,'--'),
%     ylabel('Percentage error')
%     xlabel('No. of samples')
%     legend('Measurement percentage error','Estimation percentage error')
% end

%RMSD
rmsd_measy=zeros(1,12);
rmsd_esty=zeros(1,12);
rmsd_measx=zeros(1,12);
rmsd_estx=zeros(1,12);
for i=1:1:12
    rmsd_measy(i) = sqrt(1/N*sum((y(i,:)-planty(i,:)).^2));
    rmsd_esty(i) = sqrt(1/N*sum((y(i,:)-kalmy(i,:)).^2));
    rmsd_measx(i) = sqrt(1/N*sum(v1(i,:).^2));
    rmsd_estx(i) = sqrt(1/N*sum((x(i,:)-kalmx(i,:)).^2));
end
% for i=1:1:12
%     rmsd_meas1 = sqrt(sum((y(i,:)-yv(i,:)).^2)/N);
%     rmsd_est1 = sqrt(sum((y(i,:)-kalmy(i,:)).^2)/N);
% end

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