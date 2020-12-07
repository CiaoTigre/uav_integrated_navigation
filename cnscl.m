function [ phim, dvbm ] = cnscl( wm, vm )
% 名称：coning motion and sculling motion compensate
% Function: 圆锥划桨误差补偿(多子样法)
% 程序@ 捷联惯导算法与组合导航原理 P234
% 原理@ P36(圆锥误差补偿)、P87(划桨误差补偿)
%________________________________________________________________________
% Inputs:
%       wm: n×3矩阵，每行中的三个元素为单子样时间内的角增量，n为子样数目 （rad）
%       vm: 速度增量
% Output:
%       phim: 经过误差补偿得到的角增量 (3×1 vector)
%       dvbm: 经过误差补偿得到的速度增量 (3×1 vector)

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% 2-6 子样补偿系数。注意！与 P38 表2.6.3相比，此处调换了矩阵中元素的顺序，
% 不需要再取相反的顺序。
cs = [ [   2,    0,    0,    0,     0]/3
       [   9,   27,    0,    0,     0]/20
       [  54,   92,  214,    0,     0]/105
       [ 250,  525,  650, 1375,     0]/504
       [2315, 4558, 7296, 7834, 15797]/4620 ];

% 子样直接求和
wmm = sum(wm, 1);
vmm = sum(vm, 1);

% 为修正量分配空间
dphim = zeros(1, 3);
scullm = zeros(1, 3);

% 子样数
n = size(wm, 1);
if n > 1        % 为真则代表可以用多子样法进行补偿
    
    % csw即 cs*w  
    %       cs是 1×n-1行向量，
    %       w是 n-1×3的矩阵，相当于把wm的最后一行拿掉了
    %       csw是 1×3行向量，表示P36 式(2.6.24) 中的Sigma(ki*theta_mi)
    
    % 目的是提取theta1*theta2中的前一项
    % 二子样法的系数在系数矩阵cs的第一行(n-1)保存，有一个修正项(1: n-1)，依次类推。
    % e.g. 三子样的情况：
    % [k1, k2]*[theta1_1, theta1_2, theta1_3
    %           theta2_1, theta2_1, theta2_1];
    csw = cs(n-1, 1:n-1)*wm(1:n-1, :);
    csv = cs(n-1, 1:n-1)*vm(1:n-1, :);
    
    % 圆锥误差补偿 (前后采样时刻角增量的耦合)
    dphim = cross(csw, wm(n, :));
    % 划桨误差补偿 (前后采样时刻角增量与速度增量之间的耦合)
    scullm = cross(csw, vm(n, :)) + cross(csv, wm(n, :));
    
end

% 补偿：
phim = (wmm + dphim )';
dvbm = (vmm + 0.5*cross(wmm, vmm) + scullm)';

end
