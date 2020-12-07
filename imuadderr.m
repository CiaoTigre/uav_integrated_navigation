function [ wm, vm ] = imuadderr( wm, vm, eb, web, db, wdb, ts )
% 名称：imu add error
% 功能：给真实imu输出添加误差，模拟实际中带有量测噪声的imu输出。
% 程序@ 捷联惯导算法与组合导航原理 P235
%
% Inputs:
%       wm: n×3 矩阵，每一行对应一个角增量：[theta_x, theta_y, theta_z] rad
%       vm: n×3 矩阵，每一行对应一个速度增量 m/s
%       注意！wm和vm如果由四列，最后一列是对应的时间
%       eb: 陀螺仪常值漂移，即陀螺零偏
%       web:角度随机游走误差
%       db: 加速度计常值偏置
%       wdb: 速度随机游走误差
%       ts: 步长 s
% Outputs:
%       wm: 添加了误差的角增量
%       vm: 添加了误差的速度增量

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% 数据组数
[m, n] = size(wm);
% 步长开平方
sts = sqrt(ts);

% 是否保存了时间戳
switch n
    case 3
        % 列数为3，没有保存时间
        wm = wm + [ts*eb(1) + sts*web(1)*randn(m, 1), ...
                   ts*eb(2) + sts*web(2)*randn(m, 1), ...
                   ts*eb(3) + sts*web(3)*randn(m, 1)];
        
        vm = vm + [ts*db(1) + sts*wdb(1)*randn(m, 1), ...
                   ts*db(2) + sts*wdb(2)*randn(m, 1), ...
                   ts*db(3) + sts*wdb(3)*randn(m, 1)];
        
    case 4
        % 列数为4，最后一列保存时间
        wm(:, 1:3) = wm(:, 1:3) + [ts*eb(1) + sts*web(1)*randn(m, 1), ...
                                   ts*eb(2) + sts*web(2)*randn(m, 1), ...
                                   ts*eb(3) + sts*web(3)*randn(m, 1)];
        
        vm(:, 1:3) = vm(:, 1:3) + [ts*db(1) + sts*wdb(1)*randn(m, 1), ...
                                   ts*db(2) + sts*wdb(2)*randn(m, 1), ...
                                   ts*db(3) + sts*wdb(3)*randn(m, 1)];
    otherwise
        disp('imuadderr.m 输入参数形式有误！');
end

end
