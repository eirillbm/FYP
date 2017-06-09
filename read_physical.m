function[Pos,V,Acc]=read_physical(h,N)

% Read in measured values
load('Test1.mat');
t = data.time';
Acc = data.linearAcc_measured';
V = data.angularVel';
Pos = data.position_predicted';

% % Read in estimated values
% V_est = data.velocity_estimated';
% Pos_est = data.position_estimated';

% Convert angular variables to radians from degrees
V = V*pi/180;
Pos(4:6,:) = Pos(4:6,:)*pi/180;

% Rotate accelerations and angular velocities from local to global reference frame
Rot=zeros(3);
for k=1:N
    Rot=[cos(Pos(6,k))*cos(Pos(5,k)), -sin(Pos(6,k))*cos(Pos(4,k))+cos(Pos(6,k))*sin(Pos(5,k))*sin(Pos(4,k)),sin(Pos(6,k))*sin(Pos(4,k))+cos(Pos(6,k))*cos(Pos(4,k))*sin(Pos(5,k));
    sin(Pos(6,k))*cos(Pos(5,k)),cos(Pos(6,k))*cos(Pos(4,k))+sin(Pos(4,k))*sin(Pos(5,k))*sin(Pos(6,k)),-cos(Pos(6,k))*sin(Pos(4,k))+sin(Pos(5,k))*sin(Pos(6,k))*cos(Pos(4,k));
    -sin(Pos(5,k)),cos(Pos(5,k))*sin(Pos(4,k)),cos(Pos(5,k))*cos(Pos(4,k))];
    Acc(:,k)=Rot*Acc(:,k);
    V(:,k)=Rot*V(:,k);
end

% Subtract gravitational acceleration from acceleration measurements
Acc(3,:)=Acc(3,:)+9.80665;

% Plot input positions
figure;
plot(t,Pos(1:3,:))
title('Input translations')
legend('x','y','z')
xlabel('Time [s]')
ylabel('Translation [m]')

figure;
plot(t,Pos(4:6,:))
title('Input rotations')
legend('phi','theta','psi')
xlabel('Time [s]')
ylabel('Rotation [rad]')

% Plot input velocities
figure;
plot(t, V)
legend('p','q','r')
title('Input velocities')
xlabel('Time [s]')
ylabel('Angular velocity [rad/s]')

% Plot input accelerations
figure;
plot(t, Acc)
legend('a_x','a_y','a_z')
title('Input accelerations')
xlabel('Time [s]')
ylabel('Linear acceleration [m/s^2]')

% % Plot estimated positions
% figure;
% plot(t,Pos_est)
% title('Estimated positions')
% legend('x','y','z','phi','theta','psi')
% xlabel('Time [s]')
% ylabel('Translation [m], rotation [rad]')
% 
% % Plot estimated velocities
% figure;
% plot(t, V_est)
% legend('u','v','w','p','q','r')
% title('Estimated velocities')
% xlabel('Time [s]')
% ylabel('Linear velocity [m/s], angular velocity [rad/s]')




