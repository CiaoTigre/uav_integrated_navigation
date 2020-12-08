%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%   SINS/GPS Intergration Navigation System test version 2.0
%   Indirect kalman filter（反馈校正法）
%
% Reference:
%   - @捷联惯导算法与组合导航原理,严恭敏 P239
% 
% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1


close all
clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load Uav true trajectory data.
addpath UavTrajectorySim;

disp(' ')
disp('Available UAV Truth Trajectory Data Files:')
dir_mat_files = dir('UavTrajectorySim\*.mat');
for nFile=1:length(dir_mat_files)
    fprintf('   %d: %s\n',nFile,dir_mat_files(nFile).name);
end
% nFileChoice = input('Choose a UAV Truth data file (e.g. 1<Enter>): ');
try
%     load(dir_mat_files(nFileChoice).name)
    load(dir_mat_files(1).name)
catch
    error('Selected UAV Truth Trajectory data file (%d) is invalid.\n',nFileChoice);
end



gvar_earth;

% 单次更新中使用的子样数
nn = 2;
% 采样时间
ts = 0.01;
nts = nn*ts;


% 初始姿态、速度、位置
att0 = [0, 0, 90]'*arcdeg;
vn0  = [0, 0, 0]';
pos0 = [34*arcdeg, 108*arcdeg, 100]'; % lattitude, longtitude, height
qbn0 = a2qua(att0);


% 姿态四元数、速度、位置
qbn = qbn0;
vn = vn0;
pos = pos0;

eth = earth(pos, vn);


% *** 添加误差 *** 
% 失准角
phi = [0.1, 0.2, 1]'*arcmin;
qbn = qaddphi(qbn, phi);

% 陀螺零偏，角度随机游走
eb_ref = [0.1, 0.15, 0.2]'*dph;
eb = [0.01, 0.015, 0.02]'*dph;
web = [0.001, 0.001, 0.001]'*dpsh;

% 加计零偏，速度随机游走
db_ref = [800, 900, 1000]'*ug;
db = [80, 90, 100]'*ug;
wdb = [1, 1, 1]'*ugpsHz;

Qk = diag([web', wdb', zeros(1, 9)]')^2*nts;

rk = [[0.1, 0.1, 0.1], [5/Re, 5/Re, 5]]';
Rk = diag(rk)^2;

% 协方差矩阵，x = [phi, delta_vn, delta_p, eb, db]
P0 = diag([[0.1, 0.1, 10]*arcdeg, [1, 1, 1], [10/Re, 10/Re, 10]...
           [0.1, 0.1, 0.1]*dph, [80, 90, 100]*ug]')^2;

% 量测矩阵
Hk = [zeros(6,3), eye(6), zeros(6, 6)];

% Kalman filter initialization
kf = kfinit(Qk, Rk, P0, zeros(15), Hk);

% 与模拟轨迹时长一致
kTime = fix(t_SD/ts);   
err = zeros(kTime, 10);
xkpk = zeros(kTime, 2*kf.n + 1);

pos_ref = zeros(kTime,3);
pos_est = zeros(kTime,3);
pos_gps = zeros(kTime,3);

kk = 1;
t = 0;
for k = 2 : nn : kTime
    t = t + nts;
    
    % 获取模拟轨迹对应的imu输出: 角增量和速度增量（参考值）
    wm(1:nn,:) = imu_SD.wm(k-nn+1:k,:);
    vm(1:nn,:) = imu_SD.vm(k-nn+1:k,:);
    
    % 为IMU参考输出添加噪声
    [wm1, vm1] = imuadderr(wm, vm, eb, web, db, wdb, ts);
    
    % 惯导更新：姿态四元数、速度、位置 
    [qbn, vn, pos, eth] = insupdate(qbn, vn, pos, wm1, vm1, ts);
    
    % 基于模型预测：导航误差系统模型卡尔曼滤波
    kf.Phikk_1 = eye(15) + kfft15(eth, q2mat(qbn), sum(vm1, 1)'/nts)*nts;
    kf = kfupdate(kf);
    
    % 模拟GPS量测数据
    gps = [avp_SD.vn(k,:)'; avp_SD.pos(k,:)'] + rk.*randn(6, 1);
    pos_gps(kk,:) = gps(4:6)';
    % 量测更新 5Hz
    if mod(t, 0.2) < nts
        Zk = [vn', pos']' - gps;
        kf = kfupdate(kf, Zk, 'M');
    end
    
    % Indirect Kalman filter：feedback to IMU (反馈校正法)
    qbn = qdelphi(qbn, kf.Xk(1:3));
    vn  = vn - kf.Xk(4:6);
    pos = pos - kf.Xk(7:9);
    pos_est(kk,:) = pos';
    % 反馈校正：由于反馈项的存在导致卡尔曼滤波的先验估计值始终为零. Ref: 王辰熙
    kf.Xk(1:3) = 0;
    kf.Xk(4:6) = 0;
    kf.Xk(7:9) = 0;
%     kf.Xk(10:12) = 0;
%     kf.Xk(13:15) = 0;
    
        
    % compute the error between estimation & truth data 
    % Note that this 'error' is not the 'state vector' in the Kalman equ. 
    % In indirect kalman filter, the 'state vector' means the error of 
    % the IMU update (respect to True data.)
    qbn_ref = a2qua(avp_SD.att(k,:));
    vn_ref = avp_SD.vn(k,:)';
    pos_ref(kk,:) = avp_SD.pos(k,:);
    err(kk, :) = [qq2phi(qbn, qbn_ref)', (vn - vn_ref)', (pos - pos_ref(kk,:)')', t];
    xkpk(kk, :) = [kf.Xk', diag(kf.Pk)', t]';
    
    kk = kk + 1;
    
%     % 程序运行时显示当前进度
%     if mod(t, 50) == 0
% %         disp(fix(t));
%         disp('...');
%     end

end

% 为了让err有足够的空间，在初始化时我们将其长度设置为len。由于采用多子样算
% 法或者别的某些缘故，err通常“装不满”，该操作便是为了把多余的0拿掉。
err(kk:end, :) = [];
xkpk(kk:end, :) = [];
pos_ref(kk:end,:) = [];
pos_est(kk:end,:) = [];
pos_gps(kk:end,:) = [];
tt = err(:, end);

%% 以下是绘图程序
figure;
subplot(3,3,[1,4]);  
% 横轴是经度Lontitu，纵轴是纬度Latitude
plot(pos_gps(:,2)/arcdeg,pos_gps(:,1)/arcdeg, 'dg','LineWidth',0.1); hold on;
plot(pos_est(:,2)/arcdeg,pos_est(:,1)/arcdeg, 'r','LineWidth',4); hold on;
plot(pos_ref(:,2)/arcdeg,pos_ref(:,1)/arcdeg, 'b','LineWidth',1); hold on;
plot(pos_ref(1,2)/arcdeg,pos_ref(1,1)/arcdeg, 'oc','LineWidth',4);
% axis equal;
grid on;
xlabel('\it\lambda\rm /(\circ)');
ylabel('\itL\rm /(\circ)');
legend('GPS meas.','Est. pos.','True pos.', 'Start');
title('UAV Position')

subplot(3,3,2);
plot(tt, err(:, 1:2)/arcdeg);
grid on;
axis tight;
xlabel('t/s');
ylabel('\it\phi\rm/(\circ)');
legend('\it\phi\rm_E', '\it\phi\rm_N');
title('Pitch & Roll Est. error')

subplot(3,3,3);
plot(tt, err(:, 3)/arcdeg);
grid on;
axis tight;
% ylim([-10,10])
xlabel('t/s');
ylabel('\it\phi\rm_U\rm/(\circ)');
legend('\it\phi\rm_U');
title('Yaw Est. error')

subplot(3,3,5);
plot(tt, err(:, 4:6));
grid on;
axis tight;
xlabel('t/s');
ylabel('\delta\itv^n\rm/(m.s^{-1})');
legend('\delta\itv\rm_E', '\delta\itv\rm_N', '\delta\itv\rm_U');
title('Velocity Est. error')

subplot(3,3,6);
plot(tt, [err(:, 7)*Re, err(:, 8)*Re*cos(pos(1)), err(:, 9)]);
grid on;
axis tight;
ylim([-10,10]); 
xlabel('t/s');
ylabel('\delta\itp\rm/m');
legend('\delta\itL', '\delta\it\lambda', '\delta\ith');
title('Position Est. error')


subplot(3,3,7);  
plot(tt,pos_gps(:,3), ':g','LineWidth',0.1); hold on;
plot(tt,pos_est(:,3), 'r','LineWidth',2); hold on;
plot(tt,pos_ref(:,3), 'b','LineWidth',1); 
grid on;
axis tight;
xlabel('t/s');
ylabel('\ith\rm /(m)');
legend('GPS meas.','Est. Alt.','True Alt.');
title('UAV Altitude')

subplot(3,3,8);
plot(tt, xkpk(:, 10:12)/dph);
grid on;
axis tight;
xlabel('t/s');
ylabel('\it\epsilon\rm/(\circ.h^{-1})');
legend('\it\epsilon_x', '\it\epsilon_y', '\it\epsilon_z');
title('Gyro biases')

subplot(3,3,9);
plot(tt, xkpk(:, 13:15)/ug);
grid on;
axis tight;
xlabel('t/s');
ylabel('\it\nabla\rm/\mu\itg');
legend('\it\nabla_x', '\it\nabla_y', '\it\nabla_z');
title('Accelerometer biases')

% 均方误差收敛图
spk = sqrt(xkpk(:, 16:end-1 ));

msplot(321, tt, spk(:, 1:2)/arcdeg, '\it\phi\rm/(\circ)');
legend('\it\phi\rm_E', '\it\phi\rm_N');

msplot(322, tt, spk(:, 3)/arcdeg, '\it\phi\rm_U\rm/(\circ)');
legend('\it\phi\rm_U');

msplot(323, tt, spk(:, 4:6), '\delta\itv^n\rm/(m.s^{-1})');
legend('\delta\itv\rm_E', '\delta\itv\rm_N', '\delta\itv\rm_U');

msplot(324, tt, [spk(:, 7)*Re, spk(:, 8)*Re*pos(1), spk(:, 9)],...
       '\delta\itp\rm/m');
legend('\delta\itL', '\delta\it\lambda', '\delta\ith');

msplot(325, tt, spk(:, 10:12)/dph, '\it\epsilon\rm/(\circ.h^{-1})');
legend('\it\epsilon_x', '\it\epsilon_y', '\it\epsilon_z');

msplot(326, tt, spk(:, 13:15)/ug, '\it\nabla\rm/\mu\itg');
legend('\it\nabla_x', '\it\nabla_y', '\it\nabla_z');

% 三维轨迹
figure(3)
plot3(pos_gps(:,2)/arcdeg, pos_gps(:,1)/arcdeg, pos_gps(:,3), ':g','LineWidth',0.1); hold on;
plot3(pos_ref(:,2)/arcdeg, pos_ref(:,1)/arcdeg, pos_ref(:,3),'b','LineWidth',2); hold on;
plot3(pos_est(:,2)/arcdeg, pos_est(:,1)/arcdeg, pos_est(:,3),'r','LineWidth',3); hold on;
plot3(pos_ref(1,2)/arcdeg, pos_ref(1,1)/arcdeg, pos_ref(1,3),'oc','LineWidth',10);
xlabel('\it\lambda\rm /(\circ)');
ylabel('\itL\rm /(\circ)');
zlabel('\ith\rm /(m)')
legend('GPS meas.','True pos.','Est. pos.','start')
grid on;
title('3D Trajectory')
