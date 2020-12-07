%% **************************************************************
%名称 inertial navigation test version 1.0
%功能：纯惯导解算程序
%程序@ 捷联惯导系统与组合导航原理 P236
%
%

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
close all
clear all
gvar_earth;

% 惯导更新时使用的多子样数目和采样时间
nn = 2;
ts = 0.1;
nts = nn*ts;

% 姿态、速度、位置的初始值
att = [1, 1, 30]'*arcdeg;
vn = [0, 0, 0]';
pos = [34*arcdeg, 108*arcdeg, 100]';
qbn = a2qua(att);

eth = earth(pos, vn);

% 生成模拟数据
wm = qmulv(qconj(qbn), eth.wien)*ts;
vm = qmulv(qconj(qbn), -eth.gn)*ts;
% 仿真静态imu数据
wm = repmat(wm', nn, 1);
vm = repmat(vm', nn, 1);

% 失准角
phi = [0.1, 0.2, 3]'*arcmin;
% 含有误差的姿态四元数
qbn = qaddphi(qbn, phi);

% 仿真时长
len = fix(3600/ts);

% 纯惯导解算获得的姿态、速度、位置信息。这里为其预分配存储空间
avp = zeros(len, 10);

kk = 1;
t = 0;
%% pure inertial navigation
for k = 1 : nn : len
    % 当前时间
    t = t + nts;
    
    % 纯惯导更新
    [qbn, vn, pos] = insupdate(qbn, vn, pos, wm, vm, ts);
    
    % 保存计算结果
    avp(kk, :) = [q2att(qbn)', vn', pos', t];
    
    kk = kk + 1;
    
    % 在程序运行时显示当前进度
    if mod(t, 100) < nts
        disp(fix(t));
    end
    
end
%% 绘图
avp(kk:end, :) = [];

tt = avp(:, end);

msplot(221, tt, avp(:, 1:2)/arcdeg, 'Att/(\circ)');
legend('\it\theta', '\it\gamma');

msplot(222, tt, avp(:, 3)/arcdeg, '\psi /\circ');

msplot(223, tt ,avp(:, 4:6), 'Vel /m.s^{-1}');
legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U');

msplot(224, tt, deltapos(avp(:, 7:9)), '\DeltaPos /m');
legend('\Delta\itL', '\Delta\it\lambda', '\Delta\ith');
