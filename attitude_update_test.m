%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 名称：attitude update test
% Function: 姿态更新算法测试程序
% 
% 
% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
clear all
close all
% 计算和保存一些全局变量
gvar_earth;
pos = [46*arcdeg, 126*arcdeg, 100]';
eth = earth(pos, [0, 0, 0]');

% 读入数据
load('angle_motion_data');
% 数据长度
[data_roll, data_column] = size(att_SD);
% 总时间
total_time_of_data = att_SD(end, 4);
% 多子样算法采用的子样数
% 注意！调试时发现：若想"声明"一个数为整形，最好调用round(四舍五入)而不是
% fix(朝零进发)，以免出现时间戳与数据对应异常的问题。
num_sample = round(2);
% imu输出对应步长
ts_imu = 0.01;

% 初始时刻,e.g.设置init_time=20s,采用二子样算法输出的首个姿态将是t=20.02s
% 的姿态
init_time = 0;
% t=init_time的数据在att_SD中的index
k_init = round(init_time/ts_imu) + 1;

% 初始姿态
Cbn0_ref = a2mat(att_SD(k_init, 1:3));

Cbn0 = Cbn0_ref;
Cntn0 = eye(3);
Cbtb0 = eye(3);

% 为动态变量分别存储空间
att_cmpt_sv1 = zeros(data_roll, 3);
att_ref_sv = zeros(data_roll, 3);
att_err_sv = zeros(data_roll, 3);
idx_sv = zeros(data_roll, 3);
% 循环次数
Lp = 0;
Lp_subsample = 0;

% 当前时刻
crt = ts_imu;
i = 0;
%% 
while round((crt + init_time)/ts_imu) <= round(120/ts_imu)
    % 第Lp次循环
    Lp = Lp+1;
    
    % crt时刻数据所对应的编号
    i = round(crt/ts_imu) + k_init;
    idx_sv(Lp, :) = [crt + init_time, i, round((crt + init_time)/ts_imu) - i];
    
    % 用不同方法更新姿态 1.二子样法更新四元数 2.四阶龙格库塔法
    % *** 二子样 ***
    if mod(round(Lp), num_sample) == 0  
        % i是num_sample的整数倍，获得了足够的采样，进行更新
        Lp_subsample = Lp_subsample+1;
        % imu增量形式的输出用上下限时刻速率取平均
        vm = 0.5*ts_imu.*(imu_SD.fb_ref(i-num_sample+1 : i, 1:3)...
            + imu_SD.fb_ref(i-num_sample : i-1, 1:3));
        
        wm = 0.5*ts_imu.*(imu_SD.wb_ref(i-num_sample+1 : i, 1:3)...
            + imu_SD.wb_ref(i-num_sample : i-1, 1:3));
        
        % 补偿圆锥/划桨误差
        [ phim, dvbm ] = cnscl(wm, vm);
        
        % 计算姿态变化量
        Cntn0 = Cntn0*rv2m(eth.winn*num_sample*ts_imu);
        Cbtb0 = Cbtb0*rv2m(phim);
        Cbn = Cntn0'*Cbn0*Cbtb0;
        
        att_cmpt_sv1(Lp_subsample, :) = m2att(Cbn)*deg;
        att_ref_sv(Lp_subsample, :) = att_SD(i, 1:3)*deg;
        att_err_sv(Lp_subsample, :) = att_cmpt_sv1(Lp_subsample, :) -...
            att_ref_sv(Lp_subsample, :);
    else
        % do nothing
    end
    
    % *** 四阶龙格库塔法 ***
    
    
    % 断点位置
    stophere = 1;
    
    % 下次更新的时刻
    crt = crt + ts_imu;
end
crt = crt - ts_imu;

% 清除未使用空间
idx_sv(Lp+1:end, :) = [];
att_cmpt_sv1(Lp_subsample+1:end, :) = [];
att_ref_sv(Lp_subsample+1:end, :) = [];
att_err_sv(Lp_subsample+1:end, :) = [];

%% 绘图
Time_axis_subsample = (1:1:Lp_subsample)*num_sample*ts_imu;
num_fig = 1;
% 姿态计算值
figure(num_fig)
subplot(311)
plot(Time_axis_subsample, att_cmpt_sv1(:, 1));
xlabel('时间 /s'); ylabel('pitch /deg'); hold on; 
plot(Time_axis_subsample, att_ref_sv(:, 1));
legend('计算值','参考值');
subplot(312)
plot(Time_axis_subsample, att_cmpt_sv1(:, 2));
xlabel('时间 /s'); ylabel('pitch /deg'); hold on; 
plot(Time_axis_subsample, att_ref_sv(:, 2));
subplot(313)
plot(Time_axis_subsample, att_cmpt_sv1(:, 3));
xlabel('时间 /s'); ylabel('pitch /deg'); hold on; 
plot(Time_axis_subsample, att_ref_sv(:, 3));


% 姿态误差
msplot(311, Time_axis_subsample, att_err_sv(:, 1), '时间 /s', 'pitch err');
msplot(312, Time_axis_subsample, att_err_sv(:, 2), '时间 /s', 'roll err');
msplot(313, Time_axis_subsample, att_err_sv(:, 3), '时间 /s', 'yaw err');
