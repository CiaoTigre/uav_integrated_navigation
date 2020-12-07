%% **************************************************************
% 名称：angle motion and imu output simulation (version 1.0)
% Function：用户预设姿态，该程序模拟出姿态曲线所对应imu输出
%   
% Inputs:
%       ang_motion: 用户设置的角运动
%       ang_motion.p(r,y)的每一行保存一个正弦运动对应的参数
%       theta = A*sin( w*t + phi) + k;
%       ang_motion.p = [A1(幅值), w1(角速度), phi1(初始相位), k1(中心值)
%                 A2(幅值), w2(角速度), phi2(初始相位), k2(中心值)]
%       Time: [total_time, ts, nn] 
% Outputs:
%       att : 与设置的角运动所对应的姿态
%       imu : 包含 imu.fb_ref, imu.wb_ref, imu.fb_msr, imu.wb_msr

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1



%%
close all;
% 全局变量设置
gvar_earth;
pos = [46*arcdeg, 126*arcdeg, 100]';
eth = earth(pos, [0, 0, 0]');

total_time = 300;
ts = 0.01;
nn = 2;

% 设置角运动1
ang_motion.p = [2.5*arcdeg, 10*arcdeg, 0*arcdeg, 0*arcdeg];
ang_motion.r = [1.5*arcdeg, 10*arcdeg, 0*arcdeg, 0*arcdeg];
ang_motion.y = [10*arcdeg, 12*arcdeg, 0*arcdeg, 0*arcdeg
                6*arcdeg, 10*arcdeg, 0*arcdeg, 0*arcdeg
                8*arcdeg, 5*arcdeg, 0*arcdeg, 0*arcdeg];

% % 设置角运动2
% ang_motion.p = [0*arcdeg, 0*arcdeg, 0*arcdeg, 0*arcdeg];
% ang_motion.r = [0*arcdeg, 0*arcdeg, 0*arcdeg, 0*arcdeg];
% ang_motion.y = [10*arcdeg, 5*arcdeg, 0*arcdeg, 0*arcdeg];

% 设置传感器误差
% imu_err.eb = [0.1, 0.1, 0.1]'*dph;
% imu_err.web = [0.01, 0.01, 0.01]'*dpsh;
% imu_err.db = [500, 500, 500]'*ug;
% imu_err.wdb = [10, 10, 10]'*ugpsHz;
imu_err = imuerrorset('yanP239');

% 为动态变量分配存储空间
att_SD = zeros(fix(total_time/ts) + 10, 4);
imu_SD.fb_ref = zeros(fix(total_time/ts) + 10, 4);
imu_SD.wb_ref = zeros(fix(total_time/ts) + 10, 4);
%% simulation
% 当前时刻,注意crt从0开始计时,att_SD中首元素为t=0时刻的载体姿态
crt = 0;
% 循环次数
i = 0;
while crt < total_time
    
    % 循环次数+1
    i = i+1;
    
    % 分别提取三个姿态角运动信息
    % 计算pitch(crt), roll(crt), yaw(crt)
    [pitch, wnb_p] = caculate_theta( ang_motion.p, crt );
    [roll, wnb_r] = caculate_theta( ang_motion.r, crt );
    [yaw, wnb_y] = caculate_theta( ang_motion.y, crt );
    
    % 姿态
    Cbn = a2mat([pitch, roll, yaw]');
    
    % wnbnx, 注意该角速度并不完全投影于n系，这里只是为了方便才这样表示
    wnbnx = [wnb_p, wnb_r, wnb_y]';
    
    % 计算wibb, Cnxb中的nx表示(n1, n2, b)这三个坐标系
    Cnxb = [ cos(roll), 0, -cos(pitch)*sin(roll)
                    0,  1,            sin(pitch)
            sin(roll),  0,  cos(pitch)*cos(roll)];
    wibb_ref = Cnxb*wnbnx + Cbn'*eth.winn;
    
    % 载体无线运动，加计输出应为b系下重力加速度
    fb_ref = Cbn'*eth.gn;
    
    % 保存结果
    att_SD(i, :) = [pitch, roll, yaw, crt];
    imu_SD.fb_ref(i, :) = [fb_ref', crt];
    imu_SD.wb_ref(i, :) = [wibb_ref', crt];
    
    % 增加一个步长的时间，为下次计算做准备
    crt = crt + ts;
    
end
% 删除存储变量中的空位
att_SD(i+1:end, :) = [];
imu_SD.fb_ref(i+1:end, :) = [];
imu_SD.wb_ref(i+1:end, :) = [];

% imu_ref添加误差获得imu_msr
[ imu_SD.fb_msr, imu_SD.wb_msr ] = imuadderr( imu_SD.fb_ref, imu_SD.wb_ref, ...
     imu_err.eb, imu_err.web, imu_err.db, imu_err.wdb, ts );

% 数据保存
save('angle_motion_data', 'att_SD', 'imu_SD');
%% 绘图
% 时间轴
Time_axis = (1 : 1 : i)*ts;

% 姿态角
msplot(311, Time_axis, att_SD(:, 1)*deg, '时间 /s', 'pitch /\circ');
title('姿态角变化');
msplot(312, Time_axis, att_SD(:, 2)*deg, '时间 /s', 'roll /\circ');
msplot(313, Time_axis, att_SD(:, 3)*deg, '时间 /s', 'yaw /\circ');

% imu.fb_ref
msplot(321, Time_axis, imu_SD.fb_ref(:, 1), '时间 /s', 'fbx m/s^2');
title('加计标准输出');
msplot(323, Time_axis, imu_SD.fb_ref(:, 2), '时间 /s', 'fby m/s^2');
msplot(325, Time_axis, imu_SD.fb_ref(:, 3), '时间 /s', 'fbz m/s^2');
% imu.wb_ref
msplot(322, Time_axis, imu_SD.wb_ref(:, 1)*deg, '时间 /s', 'wibb x rad/s');
title('陀螺标准输出');
msplot(324, Time_axis, imu_SD.wb_ref(:, 2)*deg, '时间 /s', 'wibb y rad/s');
msplot(326, Time_axis, imu_SD.wb_ref(:, 3)*deg, '时间 /s', 'wibb z rad/s');

% imu fb_msr
msplot(321, Time_axis, imu_SD.fb_msr(:, 1), '时间 /s', 'fbx m/s^2');
title('加计测量值');
msplot(323, Time_axis, imu_SD.fb_msr(:, 2), '时间 /s', 'fby m/s^2');
msplot(325, Time_axis, imu_SD.fb_msr(:, 3), '时间 /s', 'fbz m/s^2');
% imu.wb_msr
msplot(322, Time_axis, imu_SD.wb_msr(:, 1)*deg, '时间 /s', 'wibb x rad/s');
title('陀螺测量值');
msplot(324, Time_axis, imu_SD.wb_msr(:, 2)*deg, '时间 /s', 'wibb y rad/s');
msplot(326, Time_axis, imu_SD.wb_msr(:, 3)*deg, '时间 /s', 'wibb z rad/s');