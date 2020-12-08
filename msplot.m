function [  ] = msplot( mnp, x, y, xstr, ystr )
%% **************************************************************
% 名称：
% 功能：作图
% 程序@ 捷联惯导系统与组合导航原理 P236
%
% Inputs:
%       mnp: subplot(mnp) = subplot(m, n, p)
%       x, y: x、y轴所对应的数据
%       xstr, ystr: x、y轴的注释

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% 完全按照书上copy
% 如果是第一幅小图，则新建一个figure
% subplot(mnp) = subplot(m, n, p),mnp必须为三位数
if mod(mnp, 10) == 1
    figure;
end

subplot(mnp);
plot(x, y);
grid on;
axis tight;

% x轴和y轴的注释
if nargin == 4
    ystr = xstr;
    xstr = 't/s';
end

xlabel(xstr);   ylabel(ystr);

end
