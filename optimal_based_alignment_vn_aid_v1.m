%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 名称：Opitiaml Based Alignment (GPS velocity aided) version 1.0
% 功能：
%
% 缩写的含义：
% msr: measurement; 
% ref: reference; 特指真实值
% sv: save; 
% prv: previous; 对应变量在上次解算后的值或在本次解算开始时刻的值
% opt: optimal

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% 全局变量
close all
gvar_earth;

% ** 下载数据 **
% 数据包含 
% 1. trajectory_ref.(pos, vn, att) 
% 2. imu.ref.(acc, gyr)
% 轨迹仿真程序生成的数据，trajectory的长度会比imu长1。trajectory首个数据为
% t0时刻的avp_ref, imu首个数据为理想imu在t1时刻输出的(t0, t1]内角度和速
% 度增量。
load('trj_and_imu.mat');

% ** 与测量数据相关的量 **
% imu输出频率 (Hz)
f_imu = 100;
ts_imu = 1/f_imu;
% gps输出频率 (Hz)
f_gps = f_imu;
ts_gps = 1/f_gps;
% 总数据长度 (个)
num_imu_data = size(imu_ref.acc, 1);

% *** 给imu数据添加误差 ***
imu_err = imuerrorset('selfdefine');
% add imu error
[imu_msr.gyr, imu_msr.acc] = imuadderr(imu_ref.gyr, imu_ref.acc, ...
                imu_err.eb, imu_err.web, imu_err.db, imu_err.wdb, ts_imu);

% ** 设置全局变量 **
% 第一个"读入"的数据所对应的时刻以及在imu_ref和trajectory_ref中的位置 (s)
% 可以认为导航系统是在第start_time秒启动的
start_time = 0;
first_data_index = round(start_time/ts_imu) + 1;
% 对准时间长度 (s)
total_alignment_time = 92;

% 更新算法子样数 (个)
num_subsample = 2;
nts = num_subsample*ts_imu;

% 初始位置、速度、姿态
pos0 = trajectory_ref.pos(first_data_index, :)';
vn0 = trajectory_ref.vn(first_data_index, :)';
att0 = trajectory_ref.att(first_data_index, :)';

att0_ref = att0;
Cbn0_ref = a2mat(att0);
q0_ref = a2qua(att0);

% *** 一些变量在上一时刻的值 ***
% 速度、位置、姿态 
% 目前还没有准备对初始时刻值添加误差
pos_prv = pos0;
vn_prv = vn0;
att_prv = att0;

% eth保存惯导解算需要的变量，例如wien, winn等。
eth_prv = earth(pos_prv, vn_prv);

% bt系和nt系相对惯性空间变化量
Cntn0_prv = eye(3);
Cbtb0_prv = eye(3);
qbtb0 = [1, 0, 0, 0]';

% 上次更新过后的alpha和beta
alpha_prv = zeros(3, 1);
beta_prv = zeros(3, 1);
% alpha中除了累加项外还有随时间变化的项，因此单独保存alpha中的累加项是值得
% 的,不妨特地设置一个保存alpha中累加项的变量。
alpha_sigma = zeros(3, 1);
K_prv = zeros(4, 4);

% 分配动态变量存储空间
phi_sv = zeros(round(total_alignment_time/ts_imu), 3);
phi0_sv = zeros(round(total_alignment_time/ts_imu), 3);

% 计数变量
% 当前时刻 s
current = 0;
% 当前循环时第i次循环
i = 0;
% index
k = 0;
%%
% 我们的终极目标是让算法能够在实际中运行，实际中情况大概是这样的：
% 1. t=0，惯导开机
% 2. t=ts_imu，imu在此时输出第一组速度增量vm(1)和角增量wm(1)
% 3. t=nts,imu输出凑足了执行一次惯导更新所需要的子样数，导航计算机立即执行
%    惯导更新程序并输出导航结果p(t),vn(t),att(t)
% 4. 导航计算机继续等待imu测得n组测量信息再进行惯导更新
while current < total_alignment_time  
    % 每过nts秒，imu就会测量出num_subsample组采样结果。获得足够imu输出后立
    % 即在current时刻（亦称：当前时刻）进行导航更新
    current = current + nts;
    
    % 当前时刻(current)在trajectory_ref中对应的数据编号
    k = round(current/ts_imu) + first_data_index;
    % 惯导解算循环次数
    i = i+1;
    
    % current时刻，载体的位置、速度、姿态参考值
    pos_ref = trajectory_ref.pos(k, :)';
    vn_ref = trajectory_ref.vn(k, :)';
    att_ref = trajectory_ref.att(k, :)';
    qbn_ref = a2qua(att_ref);
    
    % 用参考速度模拟当前时刻GPS速度,其中位置信息暂时不加误差
    vn_gps = vn_ref + 0.*randn(3, 1);
    pos_gps = pos_ref;
    
    vn_gps_sv(i, :) = vn_gps;
    % 用gps输出的速度作为组合导航系统给出的速度
    vn = vn_gps;
    pos = pos_gps;
    
    % 从imu_ref中读出(current - nts, current]这段时间内imu输出的n组量测信息
    wm = imu_msr.gyr(k-num_subsample+1 : k, :);
    vm = imu_msr.acc(k-num_subsample+1 : k, :);
    
    wm_sv(2*i-1:2*i, :) = wm;
    vm_sv(2*i-1:2*i, :) = vm;
    % 步长圆锥/划桨误差
    [phim, dvbm] = cnscl(wm, vm);
    
    % 更新Cbtb0和Cntn0
    % 用winn(tm-1)计算Cntmntm-1
    Cntn0 = Cntn0_prv*rv2m(eth_prv.winn*nts);
    % 用经过补偿得到的(current - nts, current]这段时间内的角增量算Cbtb0
    % Cbtb0 = Cbtb0_prv*rv2m(phim);
    qbtb0 = qmul(qbtb0, rv2q(phim));
    qbtb0 = qnormlz(qbtb0);
    Cbtb0 = q2mat(qbtb0);
    
    % *** 计算alpha(n0)和beta(b0) ***
    % alpha(n0)
    % 1. 目前用vn0_ref计算alpha, 暂时不考虑滑动窗口。
    % 2. 用gps输出的位置计算wien和gn
    alpha_sigma = alpha_sigma + ... 
                  Cntn0_prv*(cross(eth_prv.wien, vn_prv) - eth_prv.gn)*nts;
    alpha = Cntn0*vn - vn0 + alpha_sigma;
    % beta(b0)
    beta = beta_prv + Cbtb0_prv*dvbm;
    
    % QUEST 方法计算qbn0
    [ qbn0, K ] = QUEST( beta, alpha, K_prv );
    
    % 姿态解算
    Cbn0 = q2mat(qbn0);
    Cbn = Cntn0'*Cbn0*Cbtb0;
    att = m2att(Cbn);
    phi0_sv(i, :) = atterrnorml(q2att(qbn0) - att0_ref)*deg;
    phi_sv(i, :) = atterrnorml(att - att_ref)*deg;
    
    % 计算导航解算时所需要的相关参数
    eth = earth(pos, vn);
    
    % 将本次更新后的导航参数作为下次更新初值
    pos_prv = pos;
    vn_prv = vn;
    att_prv = att;
    eth_prv = eth;
    Cntn0_prv = Cntn0;
    Cbtb0_prv = Cbtb0;
    alpha_prv = alpha;
    beta_prv = beta;
    K_prv = K;
    
end

% 注意：如果end前面加了空格就会出错!
phi_sv(i+1:end, :) = [];
phi0_sv(i+1:end, :) = [];
%% 绘图
time_axis = (1:1:i)*nts;
time_axis_imu = (1:1:2*i)*ts_imu;

% Cbn误差
msplot(311, time_axis, phi_sv(:, 1), 'pitch error / \circ');
msplot(312, time_axis, phi_sv(:, 2), 'roll error / \circ');
msplot(313, time_axis, phi_sv(:, 3), 'yaw error / \circ');

%Cbn0 误差
msplot(311, time_axis, phi0_sv(:, 1), 'pitch error / \circ');
msplot(312, time_axis, phi0_sv(:, 2), 'roll error / \circ');
msplot(313, time_axis, phi0_sv(:, 3), 'yaw error / \circ');

%imu输出
msplot(311, time_axis_imu, wm_sv(:, 1), 'wibbx / rad/s');
msplot(312, time_axis_imu, wm_sv(:, 2), 'wibby / rad/s');
msplot(313, time_axis_imu, wm_sv(:, 3), 'wibbz / rad/s');

%imu输出
msplot(311, time_axis_imu, vm_sv(:, 1), 'fbx / m/s^2');
msplot(312, time_axis_imu, vm_sv(:, 2), 'fby / m/s^2');
msplot(313, time_axis_imu, vm_sv(:, 3), 'fbz / m/s^2');
