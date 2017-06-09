function[Pos,V,Acc]=read_virtual(h,N)

% Read in SIMA values

% Read in positions
xtime=xlsread('input.xlsx','A3:FHZ3');
xpos=xlsread('input.xlsx','A5:FHZ5');
ytime=xlsread('input.xlsx','A9:HZQ9');
ypos=xlsread('input.xlsx','A11:HZQ11');
ztime=xlsread('input.xlsx','A15:CNY15');
zpos=xlsread('input.xlsx','A17:CNY17');
xrtime=xlsread('input.xlsx','A21:KKA21');
xrot=xlsread('input.xlsx','A23:KKA23');
yrtime=xlsread('input.xlsx','A27:FZM27');
yrot=xlsread('input.xlsx','A29:FZM29');
zrtime=xlsread('input.xlsx','A33:KOT33');
zrot=xlsread('input.xlsx','A35:KOT35');

% Read in velocities
xvtime=xlsread('input.xlsx','A39:LUR39');
xvel=xlsread('input.xlsx','A41:LUR41');
yvtime=xlsread('input.xlsx','A45:MVH45');
yvel=xlsread('input.xlsx','A47:MVH47');
zvtime=xlsread('input.xlsx','A51:LJF51');
zvel=xlsread('input.xlsx','A53:LJF53');
xvrtime=xlsread('input.xlsx','A57:NGQ57');
xrvel=xlsread('input.xlsx','A59:NGQ59');
yvrtime=xlsread('input.xlsx','A63:NGZ63');
yrvel=xlsread('input.xlsx','A65:NGZ65');
zvrtime=xlsread('input.xlsx','A69:NHZ69');
zrvel=xlsread('input.xlsx','A71:NHZ71');

% Read in accelerations
xatime=xlsread('input.xlsx','A75:NHZ75');
xacc=xlsread('input.xlsx','A77:NHZ77');
yatime=xlsread('input.xlsx','A81:NHL81');
yacc=xlsread('input.xlsx','A83:NHL83');
zatime=xlsread('input.xlsx','A87:MGR87');
zacc=xlsread('input.xlsx','A89:MGR89');
xartime=xlsread('input.xlsx','A93:MQO93');
xracc=xlsread('input.xlsx','A95:MQO95');
yartime=xlsread('input.xlsx','A99:NGA99');
yracc=xlsread('input.xlsx','A101:NGA101');
zartime=xlsread('input.xlsx','A105:NED105');
zracc=xlsread('input.xlsx','A107:NED107');

% Convert angular variables to radians from degrees
xrot    = xrot*pi/180;
yrot    = yrot*pi/180;
zrot    = zrot*pi/180;
xrvel    = xrvel*pi/180;
yrvel    = yrvel*pi/180;
zrvel    = zrvel*pi/180;
xracc    = xracc*pi/180;
yracc    = yracc*pi/180;
zracc    = zracc*pi/180;

% % Subtract gravitational acceleration from acceleration measurements (SIMA uses NWU coordinates)
% xacc     = xacc-9.80665;
% yacc     = yacc-9.80665;
% zacc     = zacc-9.80665;
% xracc    = xracc-9.80665;
% yracc    = yracc-9.80665;
% zracc    = zracc-9.80665;

% Plot input positions from SIMA
% figure;
% plot(xtime,xpos, ytime, ypos, ztime, zpos, xrtime, xrot, yrtime, yrot, zrtime, zrot)
% title('Input')
% legend('x','y','z','phi','theta','psi')

% Generate time vector 
t = 0:h:(N-1)*h;

% Interpolate values
xposi=interp1(xtime, xpos, t, 'linear', 'extrap');
yposi=interp1(ytime, ypos, t, 'linear', 'extrap');
zposi=interp1(ztime, zpos, t, 'linear', 'extrap');
xroti=interp1(xrtime, xrot, t, 'linear', 'extrap');
yroti=interp1(yrtime, yrot, t, 'linear', 'extrap');
zroti=interp1(zrtime, zrot, t, 'linear', 'extrap');

xveli=interp1(xvtime, xvel, t, 'linear', 'extrap');
yveli=interp1(yvtime, yvel, t, 'linear', 'extrap');
zveli=interp1(zvtime, zvel, t, 'linear', 'extrap');
xrveli=interp1(xvrtime, xrvel, t, 'linear', 'extrap');
yrveli=interp1(yvrtime, yrvel, t, 'linear', 'extrap');
zrveli=interp1(zvrtime, zrvel, t, 'linear', 'extrap');

xacci=interp1(xatime, xacc, t, 'linear', 'extrap');
yacci=interp1(yatime, yacc, t, 'linear', 'extrap');
zacci=interp1(zatime, zacc, t, 'linear', 'extrap');
xracci=interp1(xartime, xracc, t, 'linear', 'extrap');
yracci=interp1(yartime, yracc, t, 'linear', 'extrap');
zracci=interp1(zartime, zracc, t, 'linear', 'extrap');

% % Plot SIMA positions
% figure;
% plot(t, xposi, t, yposi, t, zposi, t, xroti, t, yroti, t, zroti)
% legend('x','y','z','phi','theta','psi')
% title('Positions from SIMA')
% xlabel('Time [s]')
% ylabel('Translation [m], rotation [rad]')
% 
% % Plot SIMA velocities
% figure;
% plot(t, xveli, t, yveli, t, zveli, t, xrveli, t, yrveli, t, zrveli)
% legend('u','v','w','p','q','r')
% title('Velocities from SIMA')
% xlabel('Time [s]')
% ylabel('Linear velocity [m/s], angular velocity [rad/s]')
% 
% % Plot SIMA accelerations
% figure;
% plot(t, xacci, t, yacci, t, zacci, t, xracci, t, yracci, t, zracci)
% legend('a_x','a_y','a_z','p_{dot}','q_{dot}','r_{dot}')
% title('Accelerations from SIMA')
% xlabel('Time [s]')
% ylabel('Linear acceleration [m/s^2], angular acceleration [rad/s^2]')

% Transpose the vectors
pos1=xposi';
pos2=yposi';
pos3=zposi';
rot1=xroti';
rot2=yroti';
rot3=zroti';

vel1=xveli';
vel2=yveli';
vel3=zveli';
rvel1=xrveli';
rvel2=yrveli';
rvel3=zrveli';

acc1=xacci';
acc2=yacci';
acc3=zacci';
racc1=xracci';
racc2=yracci';
racc3=zracci';

% Define position, velocity and acceleration matrices
Pos=[pos1,pos2,pos3,rot1,rot2,rot3];
V=[vel1,vel2,vel3,rvel1,rvel2,rvel3];
Acc=[acc1,acc2,acc3,racc1,racc2,racc3];

% Transpose matrices
Pos=Pos';
V=V';
Acc=Acc';




