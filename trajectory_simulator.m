%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%   Trajectory simulator
%
% Outputs:
%   - avp_SD.(pos,vn,att):  Reference(true) data of pos, vel & att.
%   - imu_SD.(wm,vm):       Increment of angular & vel.
%
% References:
%   - @捷联惯导算法与组合导航原理,严恭敏 P229
% 
% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1



close all
clear
clc



% 加载全局变量
gvar_earth;
% 步长
ts = 0.01;

% 初始姿态、速度、位置
att0 = [0, 0, 90]'*arcdeg;
vn0 = [0, 0, 0]';
pos0 = [34*arcdeg, 108*arcdeg, 100]'; % lattitude, longtitude, height

% 预设轨迹
%   w_pitch   w_roll   w_yaw   vb_y    time
wat = [0,       0,      0,      0,      10      %静止
       0,       0,      0,      1,      10      %加速
       0,       0,      0,      0,      10      %匀速
       5,       0,      0,      0,      4       %抬头
       0,       0,      0,      0,      10      %匀速
       -5,      0,      0,      0,      4       %低头
       0,       0,      0,      0,      10      %匀速
       0,       10,     0,      0,      1       %横滚
       0,       0,      9,      0,      10      %转弯
       0,       -10,    0,      0,      1       %横滚
       0,       0,      0,      0,      10      %匀速
       0,       0,      0,      -1,     10      %减速
       0,       0,      0,      0,      10    ];%静止

   
% % 预设轨迹
% %   w_pitch   w_roll   w_yaw   vb_y    time
% wat = [0,       0,      0,      0,      10      %静止
%        0,       0,      0,      1,      10      %加速
%        0,       0,      0,      0,      10      %匀速
%        5,       0,      0,      0,      4       %抬头
%        0,       0,      0,      0,      10      %匀速
%        -5,      0,      0,      0,      4       %低头
%        0,       0,      0,      0,      10      %匀速
%        0,       10,     0,      0,      1       %横滚
%        0,       0,      9,      0,      10      %转弯
%        0,       -10,    0,      0,      1       %横滚
%        0,       0,      0,      0,      10      %匀速
%        0,       0,      9,      0,      10      %匀速
%        0,       10,      0,      0,      1      %匀速
%        0,       0,      0,      -1,     10      %减速
%        0,       0,      0,      0,      10    ];%静止
   

% 把用deg/s表示的角速度转化成rad/s
wat(:, 1:3) = wat(:, 1:3)*pi/180;

[avp_SD.att, avp_SD.vn, avp_SD.pos] = trjprofile(att0, vn0, pos0, wat, ts);
% 这里用的avp2imu即为书上的av2imu函数 @P229
[imu_SD.wm, imu_SD.vm] = avp2imu(avp_SD.att, avp_SD.vn, avp_SD.pos, ts);

% 在imu数据中"插入" t=0时的imu输出,即(-ts,0]这段时间内的角增量和速度增量。
% 这样做使得imu_SD与avp_SD拥有相同的长度，便于数据的调用。
imu_SD.wm = [[0, 0, 0]; imu_SD.wm];
imu_SD.vm = [[0, 0, 0]; imu_SD.vm];

% 保存数据
save('trajectory_simulator_data.mat', 'avp_SD', 'imu_SD');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  作图
tt = (0 : length(avp_SD.att) - 1 )'*ts;

msplot(221, tt, avp_SD.att/arcdeg, 'Att/(\circ)');
legend('\it\theta', '\it\gamma', '\it\psi')
title('Attitutde');

msplot(222, tt, avp_SD.vn, 'Vel /m.s^{-1}');
legend('\itv\rm_E', '\itv\rm_N', '\itv\rm_U')
title('Velocity');

msplot(223, tt, deltapos(avp_SD.pos), '\DeltaPos /m');
legend('\Delta\itL', '\Delta\it\lambda', '\Delta\ith')
title('Relative position');

% 横轴是经度，纵轴是纬度
msplot(224, avp_SD.pos(:, 2)/arcdeg, avp_SD.pos(:, 1)/arcdeg, ...
    '\it\lambda\rm /(\circ)', '\itL\rm /(\circ)');    hold on
plot(avp_SD.pos(1,2)/arcdeg, avp_SD.pos(1, 1)/arcdeg, 'ro');
legend('trajectory','start position')
title('Absolute position');


% imu输出信息作图
msplot(121, tt, imu_SD.wm/ts/arcdeg, ...
    '\it\omega^b_{ib}\rm /(\circ.s^{-1})');
legend('\it\omega^b_{ibx}', '\it\omega^b_{iby}', '\it\omega^b_{ibz}');
title('Angular rate');

msplot(122, tt, imu_SD.vm/ts, '\itf^b\rm_{sf}/(m.s^{-2})');
legend('\itf^b\rm_{sf\itx}', '\itf^b\rm_{sf\ity}', '\itf^b\rm_{sf\itz}');
title('Acceleration');
