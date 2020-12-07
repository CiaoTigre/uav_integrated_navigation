function [ qbn, vn, pos, eth ] = insupdate( qbn, vn, pos, wm, vm, ts )
% 名称：Inertial Navigation System UPDATE 
% 功能：纯惯导解算更新
% 原理@ 捷联惯导算法与组合导航原理 P80 (姿态)、P82-86(速度)、P87(位置)
% 程序@ P235

% 在(tm-1, tm]这段时间内，imu输出n组数据，考虑到实际中二子样更新算法应用最
% 为广泛，这里就拿二子样的情况举例：
% 1.在 tm-0.5 时刻，imu输出(tm-1, tm-0.5]这段时间内载体的角增量wm(1, :)和速度
% 增量vm(1, :)
% 2.在 tm 时刻，imu输出(tm-0.5, tm]这段时间内载体的角增量wm(2, :)和速度
% 增量vm(2, :)
% 同时已知tm-1时刻的导航参数:qbn(tm-1),vn(tm-1),pos(tm-1)
% 利用上述信息，计算tm时刻载体的qbn、vn、pos
%
% Inputs:
%       qbn: tm-1 时刻载体姿态
%        vn: tm-1 时刻载体速度
%       pos: tm-1 时刻位置
%        wm: n×3 矩阵，每一行中保存一个时刻的角增量 (rad)
%        vm: n×3 矩阵，每一行中保存一个时刻的速度增量 (m/s)
%        ts: imu 输出对应步长

% Outputs:
%       qbn: tm 时刻姿态
%        vn: tm 时刻速度
%       pos: tm 时刻位置
%       eth: 含有了与tm-1时刻地球相关参数的结构体

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% nn：单次惯导解算中用到的imu数据子样数
nn = size(wm, 1);
% nts：单次惯导解算对应步长
nts = nn*ts;

% 圆锥/划桨误差补偿
[phim, dvbm] = cnscl(wm, vm);

% 地球相关参数计算
eth = earth(pos, vn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 速度更新
% 原理@ P82-86， vn1 的计算公是式(4.1.50)。
% dvbm 的准确含义是：
% 在原始imu输出(raw measurement)基础上，补偿了旋转误差和划桨误差的，比力项
% (f_sf_bt)产生的速度增量在btm-1系下的投影。
% 之后简称为：“比力速度增量”。
vn1 = vn + rv2m(-eth.winn*nts/2)*qmulv(qbn, dvbm) + eth.gcc*nts;
% 更新位置时采用[tm-1, tm]这段时间内的平均速度
vn = (vn + vn1)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 位置更新
pos = pos + [vn(2)/eth.RMh, vn(1)/eth.clRNh, vn(3)]'*nts;
% tm-1 时刻的速度
vn = vn1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 姿态更新：
% phim 的准确含义是：
% 与“由btm-1系转动到btm系”这个转动所对应的等效旋转矢量。
% 它是对(raw measurement)输出相加，并进行圆锥误差补偿所获得的结果
qbn = qmul(rv2q(-eth.winn*nts), qmul(qbn, rv2q(phim)));
qbn = qnormlz(qbn);

end
