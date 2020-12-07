%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%   SINS/GPS Intergration Navigation System test version 2.0
%   Indirect kalman filter（输出校正法）
%
% References:
%   - @捷联惯导算法与组合导航原理,严恭敏 P239
% 
% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1


close all
clear
clc

gvar_earth;


% load模拟轨迹数据: avp_SD and imu_SD
load('trajectory_simulator_data.mat');


% 单次更新中使用的子样数
nn = 2;
% 采样时间
ts = 0.01;
nts = nn*ts;


% 初始姿态、速度、位置
att0 = [0, 0, 90]'*arcdeg;
vn0 = [0, 0, 0]';
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

rk = [[0.1, 0.1, 0.1], [10/Re, 10/Re, 10]]';
Rk = diag(rk)^2;

% 协方差矩阵，x = [phi, delta_vn, delta_p, eb, db]
P0 = diag([[0.1, 0.1, 10]*arcdeg, [1, 1, 1], [10/Re, 10/Re, 10]...
           [0.1, 0.1, 0.1]*dph, [80, 90, 100]*ug]')^2;

% 量测矩阵
Hk = [zeros(6,3), eye(6), zeros(6, 6)];

% Kalman filter initialization
kf = kfinit(Qk, Rk, P0, zeros(15), Hk);

% 与模拟轨迹时长一致
kTime = fix(100/ts);   
err = zeros(kTime, 10);
xkpk = zeros(kTime, 2*kf.n + 1);

kk = 1;
t = 0;
for k = 2 : nn : kTime
    t = t + nts;
    
    % 获取模拟轨迹对应的imu输出: 角增量和速度增量（参考值）
    wm(1:nn,:) = imu_SD.wm(k-nn+1:k,:);
    vm(1:nn,:) = imu_SD.vm(k-nn+1:k,:);
    
    % 为imu参考输出添加噪声
    [wm1, vm1] = imuadderr(wm, vm, eb, web, db, wdb, ts);
    
    % 惯导更新：姿态四元数、速度、位置 
    [qbn, vn, pos, eth] = insupdate(qbn, vn, pos, wm1, vm1, ts);
    
    % 基于模型预测：导航误差系统模型卡尔曼滤波
    kf.Phikk_1 = eye(15) + kfft15(eth, q2mat(qbn), sum(vm1, 1)'/nts)*nts;
    kf = kfupdate(kf);
    
    % 量测更新：模拟GPS量测数据 5Hz
    if mod(t, 0.2) < nts
        gps = [avp_SD.vn(k,:)'; avp_SD.pos(k,:)'] + rk.*randn(6, 1);
        Zk = [vn', pos']' - gps;
        kf = kfupdate(kf, Zk, 'M');
    end
    
    % 间接滤波：输出矫正法（只更新导航输出，未对惯导系统进行校正）
    qbn_nav = qdelphi(qbn, kf.Xk(1:3));
    vn_nav  = vn - kf.Xk(4:6);
    pos_nav = pos - kf.Xk(7:9);

    % compute the err between the navigation output(estimation) & truth data
    qbn_ref = a2qua(avp_SD.att(k,:));
    vn_ref = avp_SD.vn(k,:)';
    pos_ref = avp_SD.pos(k,:)';
    err(kk, :) = [qq2phi(qbn_nav, qbn_ref)', (vn_nav - vn_ref)', (pos_nav - pos_ref)', t];
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
tt = err(:, end);


%% 绘图
msplot(321, tt, err(:, 1:2)/arcdeg, '\it\phi\rm/(\circ)');  % 度
legend('\it\phi\rm_E', '\it\phi\rm_N');

msplot(322, tt, err(:, 3)/arcdeg, '\it\phi\rm_U\rm/(\circ)'); % 度
legend('\it\phi\rm_U');

msplot(323, tt, err(:, 4:6), '\delta\itv^n\rm/(m.s^{-1})');
legend('\delta\itv\rm_E', '\delta\itv\rm_N', '\delta\itv\rm_U');

msplot(324, tt, [err(:, 7)*Re, err(:, 8)*Re*cos(pos(1)), err(:, 9)],...
    '\delta\itp\rm/m');
legend('\delta\itL', '\delta\it\lambda', '\delta\ith');

msplot(325, tt, xkpk(:, 10:12)/dph, '\it\epsilon\rm/(\circ.h^{-1})');
legend('\it\epsilon_x', '\it\epsilon_y', '\it\epsilon_z');

msplot(326, tt, xkpk(:, 13:15)/ug, '\it\nabla\rm/\mu\itg');
legend('\it\nabla_x', '\it\nabla_y', '\it\nabla_z');

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
