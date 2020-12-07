function [ eth ] = earth( pos, vn )
% 名称：EARTH
% 功能：
% 程序@ 捷联惯导算法与组合导航原理 P235

% Inputs:
%       pos: [longitude, latitude, high]    (rad, rad, m)
%       vn: [ve, vn, vu]
% Outputs:
%       eth: a structure which contains some variables that we need in
%       navigation computation.

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% 将输入参数调整为列向量
if size(pos, 1) == 1
    pos = pos';
end
if size(vn, 1) == 1
    vn = vn';
end

global Re ff wie g0

% 第一偏心率
ee = sqrt(2*ff - ff^2);
e2 = ee^2;

% 与纬度有关的三角函数
eth.sl = sin(pos(1)); 
eth.cl = cos(pos(1));
eth.tl = eth.sl/eth.cl;
eth.sl2 = eth.sl*eth.sl;
sl4 = eth.sl2*eth.sl2;

sq = 1 - e2*eth.sl2;
sq2 = sqrt(sq);

eth.RM = Re*(1 - e2)/sq/sq2 ;
eth.RMh = eth.RM + pos(3);
eth.RN = Re/sq2;
eth.RNh = eth.RN + pos(3);
eth.clRNh = eth.cl*eth.RNh;

% 这里在使用上下角标表示物理量的时候直接采用熟悉的方式，即投影坐标系
% 在上，相对运动的两个坐标系在下。
% eth.win^n = eth.wie^n + eth.wen^n  => eth.winn = eth.wien + eth.wenn;
eth.vn = vn;

eth.wien = wie*[0, eth.cl, eth.sl]';
eth.wenn = [-vn(2)/eth.RMh, vn(1)/eth.RNh, vn(1)/eth.RNh*eth.tl]';
eth.winn = eth.wien + eth.wenn;
% wienn = wien + winn ,表示有害加速度项括号里角速度之和
eth.wienn = eth.wien + eth.winn;

% grs80 重力模型
gLh = g0*(1 + 5.27094e-3*eth.sl2 + 2.32718e-5*sl4) - 3.086e-6*pos(3);

eth.gn = [0, 0, -gLh]';
% 比力方程中重力加速度与哥式力加速度之和
eth.gcc = eth.gn - cross(eth.wienn, vn);

end
