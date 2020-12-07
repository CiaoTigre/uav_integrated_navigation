%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 名称：Opitiaml Based Alignment (GPS velocity aided) version 2.0
% 功能：
%
%
% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
close all
clear all

% 全局变量
gvar_earth;

% 下载数据
% 包含avp_SD.(att, pos, vn) 和imu_SD.(acc, gry)
load('trajectory_simulator_data.mat');

% 数据步长
ts_imu = 0.01;
% 采用多子样算法时使用的子样数
num_subsample = 2;
nts = num_subsample*ts_imu;

% 用msr=ref+err模拟imu输出
imu_err = imuerrorset('selfdefine');
[imu_msr.gyr, imu_msr.acc] = imuadderr(imu_SD.wb, imu_SD.fb, ...
              imu_err.eb, imu_err.web, imu_err.db, imu_err.wdb, ts_imu);

% 时间变量
init_time = 0;
k_init = fix(init_time/ts_imu) + 1;
total_time = 100;

% 状态变量初值
pos_init_ref = avp_SD.pos(k_init, :)';
vn_init_ref = avp_SD.vn(k_init, :)';
att_init_ref = avp_SD.att(k_init, :)';
Cbn0_ref = a2mat(att_init_ref);

% 变量在前一时刻值,第二项是未来可能会设置的初值误差
pos_prv = pos_init_ref + [0, 0, 0]';
vn_prv = vn_init_ref + [0, 0, 0]';
vn_init = vn_init_ref + [0, 0, 0]';
att_prv = att_init_ref + [5, 5, 100]'*arcdeg;
eth_prv = earth(pos_prv, vn_prv);
K_prv = zeros(4, 4);

% 时变姿态矩阵初值
Cntn0_prv = eye(3);
Cbtb0_prv = eye(3);
alpha_sigma = [0, 0, 0]';
beta_prv = [0, 0, 0]';
beta = [0, 0, 0]';
qbtb0 = [1, 0, 0, 0]';
% 为动态变量分配内存
att_cmpt_sv = zeros(length(avp_SD.att), 3);
att0_cmpt_sv = zeros(length(avp_SD.att), 3);
att_err_sv = zeros(length(avp_SD.att), 3);

% 设置循环中一些必要的计数变量
% 获得imu量测信息的次数
Lp_msr = 0;
% 进行更新的次数
Lp_update = 0;
% 当前时刻(自从导航起始时刻而非数据起始时刻)
crt = ts_imu;
%%
% 设计思路是这样的：导航中，每过ts_imu时间imu就会输出一组角增量和速度增量
% 信息。在采用多子样算法(以二子样为例)时，每获得2组imu量测值，导航计算机进行
% 一次惯导更新。Lp的含义就是imu输出信息的次数。
while (crt + init_time) < total_time
    % 进入循环，计数器+1
    Lp_msr = Lp_msr+1;
    
    % 当前时刻数据在avp_SD和imu_SD中对应的编号 
    k = round(crt/ts_imu) + k_init;
    
    % 每获得足够子样数的量测就进行一次惯导更新,否则什么也不做
    if mod(Lp_msr, num_subsample) == 0
        % 进行惯导更新的次数
        Lp_update = Lp_update+1;
        
        % 读入imu量测信息
        % 实际情况下导航计算机是每检测到一次imu输出就会保存一次。但是为
        % 编程方便，这里改为了每隔nts导航计算计就一下子读出num_subsample
        % 次采样，而不是一次读入并保存一个数据
        wm = imu_msr.gyr(k-num_subsample+1 : k, :);
        vm = imu_msr.acc(k-num_subsample+1 : k, :);
        
        % 补偿圆锥/划桨误差
        [phim, dvbm] = cnscl(wm, vm);
        
        % 目前除了计算姿态之外该程序暂不进行惯导更新，pos和vn都用GPS数据
        % 模拟GPS数据(暂时不加误差)
        gps.pos = avp_SD.pos(k, :)';
        gps.vn = avp_SD.vn(k, :)';
        pos = gps.pos;
        vn = gps.vn; 
        
        % 更新Cbtb0和Cntn0
        Cbtb0 = Cbtb0_prv*rv2m(phim);
        Cntn0 = Cntn0_prv*rv2m(eth_prv.winn*nts);
        
        % *** 计算alpha(n0)和beta(b0) ***
        % alpha(n0)
        % 1. 目前用vn0_ref计算alpha, 暂时不考虑滑动窗口。
        % 2. 用gps输出的位置计算wien和gn
        alpha_sigma = alpha_sigma + ...
            Cntn0_prv*(cross(eth_prv.wien, vn_prv) - eth_prv.gn)*nts;
        alpha = Cntn0*vn - vn_init + alpha_sigma;
        % beta(b0)
        beta = beta_prv + Cbtb0_prv*dvbm;
        
        % QUEST法求qbn0
        [ qbn0, K_prv ] = QUEST( beta, alpha, K_prv );
        
        % 
        Cbn0 = q2mat(qbn0); 
        Cbn = Cntn0'*Cbn0*Cbtb0;
        att_cmpt_sv(Lp_update, :) = m2att(Cbn)';
        att_err_sv(Lp_update, :) = atterrnorml(m2att(Cbn)' ...
                                        - avp_SD.att(k, :));
        att0_cmpt_sv(Lp_update, :) = q2att(qbn0)';
        
        % 计算一些必要角速度
        eth = earth(pos, vn);
        
        % 设置下次循环的初值
        Cbtb0_prv = Cbtb0;
        Cntn0_prv = Cntn0;
        beta_prv = beta;
        eth_prv = eth;
        pos_prv = pos;
        vn_prv = vn;
    else
        % do noting
    end
        
    % 下次进入循环的时间
    crt = crt + ts_imu;
    
end

% 删除空余空间
att_cmpt_sv(Lp_update+1 : end, :) = [];
att0_cmpt_sv(Lp_update+1 : end, :) = [];
att_err_sv(Lp_update+1 : end, :) = [];
%% 绘图
Time_axis = (1:1:Lp_update)*nts;
length_imu = (1:1:length(imu_SD.fb))*ts_imu;

% 计算姿态(Cbn0)
msplot(311, Time_axis, att0_cmpt_sv(:, 1)*deg, '时间 /s', 'pitch');
title('Cbn0');
msplot(312, Time_axis, att0_cmpt_sv(:, 2)*deg, '时间 /s', 'roll');
msplot(313, Time_axis, att0_cmpt_sv(:, 3)*deg, '时间 /s', 'yaw');

% 计算姿态(Cbn)
msplot(311, Time_axis, att_cmpt_sv(:, 1)*deg, '时间 /s', 'pitch');
title('change of Cbn');
msplot(312, Time_axis, att_cmpt_sv(:, 2)*deg, '时间 /s', 'roll');
msplot(313, Time_axis, att_cmpt_sv(:, 3)*deg, '时间 /s', 'yaw');

% 姿态误差
msplot(311, Time_axis, att_err_sv(:, 1)*deg, '时间 /s', 'pitch');
title('attitude error');
msplot(312, Time_axis, att_err_sv(:, 2)*deg, '时间 /s', 'roll');
msplot(313, Time_axis, att_err_sv(:, 3)*deg, '时间 /s', 'yaw');

% imu输出模拟值
msplot(311, length_imu, imu_SD.fb(:, 1), '时间 /s', 'fb(1)');
title('accelerometer inpute');
msplot(312, length_imu, imu_SD.fb(:, 2), '时间 /s', 'fb(2)');
msplot(313, length_imu, imu_SD.fb(:, 3), '时间 /s', 'fb(3)');

% imu输出模拟值
msplot(311, length_imu, imu_SD.wb(:, 1), '时间 /s', 'wb(1)');
title('gyroscope inpute');
msplot(312, length_imu, imu_SD.wb(:, 2), '时间 /s', 'wb(2)');
msplot(313, length_imu, imu_SD.wb(:, 3), '时间 /s', 'wb(3)');
