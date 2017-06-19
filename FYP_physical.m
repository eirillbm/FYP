%************************************************************************%
%                                                                        %
%   Final Year Project Code                                              %
%   For tests with physical data                                         %
%                                                                        %
%   Version: 1.1                                                         %
%   Author: Eirill Bachmann Mehammer                                     %
%                                                                        %
%************************************************************************%

%noiseval=[23.01 30 33.01 40 43.01 50];
%for m=1:length(noiseval)
close all;
clearvars;
clc;

% Define type of Kalman estimator
kalman_method='steady-state';   %steady-state or time-varying

%% ------ Kinematic Model ------ %%

% Define step size and number of steps
h=0.01/3; % step size
N=68617; % number of steps

% Generate time vector based on this
t = 0:h:(N-1)*h;

% Input matrices from measurements
[Pos,Vel,Acc,V_est]=read_physical(N);

% Store input data in new vectors to preserve original values
u = zeros(3,N);
v = zeros(3,N);
w = zeros(6,N);
for n=1:3
    u(n,:)= Acc(n,:);
    v(n,:)= Vel(n,:);
    w(n,:)= Pos(n,:);
    w(n+3,:)= Pos(n+3,:);
end

%Differentiate angular velocities to get angular accelerations
acc = zeros(3,N);
for i = 1:3
    acc(i,2:end)  = diff(Vel(i,:))/h;
end

% Definition of matrices
A=zeros(12);
A(1,7)=1;
A(2,8)=1;
A(3,9)=1;
A(4,10)=1;
A(5,11)=1;
A(6,12)=1;

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

%Check observability of continuous system
Obs_c=[C; C*A; C*A^2; C*A^3; C*A^4; C*A^5; C*A^6; C*A^7; C*A^8; C*A^9; C*A^10; C*A^11];
Oc_rank=rank(Obs_c);

%Define continuous state space system
sys=ss(A,[B G],C,[D H]);

%Discretize state space system
sysd = c2d(sys,0.01/3);        %0.01/3=step size

%Check observability and reachability of discrete system
Obs=[sysd.C; sysd.C*sysd.A; sysd.C*sysd.A^2; sysd.C*sysd.A^3; sysd.C*sysd.A^4; sysd.C*sysd.A^5; sysd.C*sysd.A^6; sysd.C*sysd.A^7; sysd.C*sysd.A^8; sysd.C*sysd.A^9; sysd.C*sysd.A^10; sysd.C*sysd.A^11];
Od_rank=rank(Obs);

Reach=[sysd.B sysd.A*sysd.B sysd.A^2*sysd.B sysd.A^3*sysd.B sysd.A^4*sysd.B sysd.A^5*sysd.B sysd.A^6*sysd.B sysd.A^7*sysd.B sysd.A^8*sysd.B sysd.A^9*sysd.B sysd.A^10*sysd.B sysd.A^11*sysd.B];
R_rank=rank(Reach);

% Process noise vectors (white, Gaussian, zero-mean)
uv=zeros(3,N);
uv2=zeros(3,N);
uv3=zeros(3,N);
uv4=zeros(3,N);

v1=zeros(12,N);
 
for j=1:3
    uv(j,:) = awgn(u(j,:),20,'measured');          %SNR=20dB (1% noise)  
    uv2(j,:) = awgn(u(j,:),100,'measured');        %SNR=100dB (0.00000001% noise)
    uv3(j,:) = awgn(u(j,:),100,'measured');        %SNR=100dB (0.00000001% noise)
end

v1(1:3,:) = (uv - u);
v1(4:6,:) = (uv2 - u);
v1(7:9,:) = (uv3 - u);
v1(10:12,:)=acc;

% figure;
% plot(t,v1);
% title('Process noise')

% Initial conditions
x0 = [w(:,1); (w(1:3,2)-w(1:3,1))/h; v(:,1)];
y0 = [w(:,1); v(:,1); u(:,1)];

% Define arrays to hold the state and output values at each time step
x   = zeros(12,N);
y   = zeros(12,N);

% Set the first value of the state vector to be the initial condition
x(:,1) = x0;
y(:,1) = y0;

% Output noise vectors (white, Gaussian, zero-mean)
yv=zeros(12,N);
 for j=1:12
     yv(j,:) = awgn(y(j,:),20,'measured');        %SNR=20dB (1% noise)
 end
 
v2 = yv - y;

% figure;
% plot(t,v2);
% title('Measurement noise')

% Define noise covariance matrices
Q=cov(transpose(v1));    % Process noise covariance matrix
R=cov(transpose(v2));    % Output noise covariance matrix
Nn=zeros(12);        % Correlation between process and sensor noise is zero 

% Loop to update the state and output values
for k=1:N-1 % need -1 here so that x(n+1) does not lead to an index>N
    x(:,k+1)=sysd.A*x(:,k)+sysd.B*[u(:,k); v1(:,k)];
    y(:,k+1)=sysd.C*x(:,k+1)+sysd.D*[u(:,k+1); v1(:,k+1)]+v2(:,k+1);
end

% Check positivity
e1=[sysd.B(:,4:15) zeros(12);sysd.D(:,4:15) eye(12)];
e2=[Q Nn;Nn' R];
e3=e1';
e=eig(e1*e2*e3);

% % Plot true outputs
% figure;
% plot(t,Pos(1,:),'k');
% hold on 
% plot(t,Pos(2,:),'k--');
% plot(t,Pos(3,:),'k:');
% plot(t,Pos(4,:),'r');
% plot(t,Pos(5,:),'r--');
% plot(t,Pos(6,:),'r:');
% plot(t,Vel(1,:),'g');
% plot(t,Vel(2,:),'g--');
% plot(t,Vel(3,:),'g:');
% plot(t,Acc(1,:),'b');
% plot(t,Acc(2,:),'b--');
% plot(t,Acc(3,:),'b:');
% hold off
% title('True outputs')
% legend('x','y','z','phi','theta','psi','p','q','r','a_x','a_y','a_z')
% 
% %Plot outputs from MATLAB
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
% title('Outputs from MATLAB')
% legend('x','y','z','phi','theta','psi','p','q','r','a_x','a_y','a_z')
% 
%% ------ Debug: plot states with index ------ %%
 
% Plot MATLAB positions
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

% Plot difference between MATLAB and measured positions
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
% title('Difference in positions between MATLAB and measured')
% xlabel('Time [s]')
% ylabel('Translation [m], rotation [rad]')
% 
% %Plot difference between MATLAB and measured velocities
% figure
% plot(t,x(10,:)-Vel(1,:),'r');
% hold on
% plot(t,x(11,:)-Vel(2,:),'r--');
% plot(t,x(12,:)-Vel(3,:),'r:');
% hold off
% legend('p','q','r')
% title('Difference in velocities between MATLAB and measured')
% xlabel('Time [s]')
% ylabel('Angular velocity [rad/s]')
% 
% %Plot relative difference between MATLAB and measured positions
% pos_rel_err=zeros(6,N);
% vel_rel_err=zeros(3,N);
% for i=1:3
%     pos_rel_err(i,:) = deleteoutliers(abs((x(i,:)-Pos(i,:))./Pos(i,:)),0.05,1);
%     pos_rel_err(i+3,:) = deleteoutliers(abs((x(i+3,:)-Pos(i+3,:))./Pos(i+3,:)),0.05,1);
%     vel_rel_err(i,:) = deleteoutliers(abs((x(i+6,:)-V(i,:))./V(i,:)),0.05,1);
% end
% figure
% plot(t,pos_rel_err);
% legend('x','y','z','phi','theta','psi')
% title('Relative difference between MATLAB and measured positions')
% xlabel('Time [s]')
% ylabel('Translation [m], rotation [rad]')
% 
% figure
% plot(t,vel_rel_err);
% legend('p','q','r')
% title('Relative difference between MATLAB and measured velocities')
% xlabel('Time [s]')
% ylabel('Angular velocity [rad/s]')

%% ------ Kalman estimator design ------ %%

% Define arrays to hold the Kalman state and output values at each time step
kalmx   = zeros(12,N);
kalmy   = zeros(12,N);

% Initial conditions
kalmx(:,1) = x0;
kalmy(:,1) = y0;

% Measurements
truey = [Pos;Vel;Acc];
truex = [Pos;V_est(1:3,:);Vel];

switch kalman_method
    case 'steady-state'    %Steady-state Kalman estimator 
        [kalmf,K,P] = kalman(sysd,Q,R,Nn);

        % Loop to update state and output values of Kalman estimator
        for k=1:N-1
            kalmx(:,k+1)=sysd.A*kalmx(:,k)+sysd.B(:,1:3)*u(:,k)+K*(truey(:,k)-kalmy(:,k)); 
            kalmy(:,k+1)=sysd.C*kalmx(:,k+1)+sysd.D(:,1:3)*u(:,k+1);
        end
        
    case 'time-varying'  %Time-varying Kalman estimator      
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
        warning('Unexpected Kalman estimation method')
end

%% ------ Compare the measured and estimated states and outputs ------ %%
% % Plot estimated outputs
% figure;
% plot(t,kalmy(1,:),'k');
% hold on 
% plot(t,kalmy(2,:),'k--');
% plot(t,kalmy(3,:),'k:');
% plot(t,kalmy(4,:),'r');
% plot(t,kalmy(5,:),'r--');
% plot(t,kalmy(6,:),'r:');
% plot(t,kalmy(7,:),'g');
% plot(t,kalmy(8,:),'g--');
% plot(t,kalmy(9,:),'g:');
% plot(t,kalmy(10,:),'b');
% plot(t,kalmy(11,:),'b--');
% plot(t,kalmy(12,:),'b:');
% hold off
% title('Kalman estimator output')
% legend('x','y','z','phi','theta','psi','p','q','r','a_x','a_y','a_z')
% 
% % Plot estimated states
% figure
% plot(t,kalmx(1,:),'k');
% hold on 
% plot(t,kalmx(2,:),'k--');
% plot(t,kalmx(3,:),'k:');
% plot(t,kalmx(4,:),'r');
% plot(t,kalmx(5,:),'r--');
% plot(t,kalmx(6,:),'r:');
% hold off
% legend('x','y','z','phi','theta','psi')
% title('Estimated positions')
% xlabel('Time [s]')
% ylabel('Translation [m], rotation [rad]')
% 
% figure
% plot(t,kalmx(7,:),'k');
% hold on 
% plot(t,kalmx(8,:),'k--');
% plot(t,kalmx(9,:),'k:');
% plot(t,kalmx(10,:),'r');
% plot(t,kalmx(11,:),'r--');
% plot(t,kalmx(12,:),'r:');
% hold off
% legend('u','v','w','p','q','r')
% title('Estimated velocities')
% xlabel('Time [s]')
% ylabel('Linear velocity [m/s], angular velocity [rad/s]')

% Compute and plot absolute and reltive errors
kalmy_rel_err=zeros(12,N);
 for i=1:1:12
    figure
    subplot(3,1,1), plot(t,truey(i,:),t,kalmy(i,:),'--'),
    ylabel('Output')
    xlabel('No. of samples')
    legend('Measured output','Estimated output')
    subplot(3,1,2), plot(t,truey(i,:)-kalmy(i,:),'--'),
    ylabel('Estimation error')
    xlabel('No. of samples')
    
    [kalmy_rel_err(i,:)] = deleteoutliers(abs((kalmy(i,:)-truey(i,:))./truey(i,:)),0.05,1);
    subplot(3,1,3), plot(t,kalmy_rel_err(i,:),'--');
    ylabel('Estimation relative error')
    xlabel('No. of samples')
    legend('Estimation relative error')
end

kalmx_rel_err=zeros(12,N);
for i=1:1:12
    figure
    subplot(3,1,1), plot(t,truex(i,:),t,kalmx(i,:),'--'),
    ylabel('Output')
    xlabel('No. of samples')
    legend('Measured states','Estimated states')
    subplot(3,1,2), plot(t,truex(i,:)-kalmx(i,:),'--'),
    ylabel('Estimation error')
    xlabel('No. of samples')

    [kalmx_rel_err(i,:)] = deleteoutliers(abs((kalmx(i,:)-truex(i,:))./truex(i,:)),0.05,1);
    subplot(3,1,3), plot(t,kalmx_rel_err(i,:),'--'),
    ylabel('Relative estimation error')
    xlabel('No. of samples')
end

% Compute RMSD values
rmsd_esty=zeros(1,12);
rmsd_estx=zeros(1,12);
for i=1:1:12
    rmsd_esty(i) = sqrt(1/N*sum((truey(i,:)-kalmy(i,:)).^2));
    rmsd_estx(i) = sqrt(1/N*sum((truex(i,:)-kalmx(i,:)).^2));
end

% Compute correlation coefficients
correlation_coefx=zeros(2,24);
correlation_coefy=zeros(2,24);
r_sqrx=zeros(1,12);
r_sqry=zeros(1,12);
for i=1:12
    correlation_coefx(:,2*i-1:2*i) = corrcoef(truex(i,:),kalmx(i,:));
    r_sqrx(i) = power(correlation_coefx(1,2*i),2);
    correlation_coefy(:,2*i-1:2*i) = corrcoef(truey(i,:),kalmy(i,:));
    r_sqry(i) = power(correlation_coefy(1,2*i),2);
end

% Compute measurement error
MeasErr = zeros(12,N);
MeasErrCov = zeros(1,12);
for i =1:12
    MeasErr(i,:) = v2(i,:);
    MeasErrCov(i) = sum(MeasErr(i,:).*MeasErr(i,:))/length(MeasErr(i,:));
end

% % Plot measurement error covariance
% figure;
% plot(MeasErrCov)
% hold on
% title('Error covariance')

% Compute estimation error
EstErr = zeros(12,N);
EstErrCov = zeros(1,12);
for i =1:12
    EstErr(i,:) = truey(i,:)-kalmy(i,:);
    EstErrCov(i) = sum(EstErr(i,:).*EstErr(i,:))/length(EstErr(i,:));
end

% % Plot estimation error covariance
% plot(EstErrCov)
% legend('Measurement error covariance','Estimation error covariance')

% % Write to file
% fileID = fopen('noiseval.txt','a');
% fprintf(fileID,'%f ',noiseval(m));
% fprintf(fileID,'\n');
% fprintf(fileID,'%f ',rmsd_estx);
% fprintf(fileID,'\n');
% fprintf(fileID,'%f ',r_sqrx);
% fprintf(fileID,'\n');
% fclose(fileID);
% %end