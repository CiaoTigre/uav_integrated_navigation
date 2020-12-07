function [ vn ] = qmulv( qbn, vb )
% 名称：Quaternion MULtiply Vector
% 功能：用四元数对三维向量进行坐标变换
% 程序@ 捷联惯导系统与组合导航原理 P234
%
% INputs:
%       q: qbn
%       v: vb
% Outputs:
%       vo: vn

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
qi = [0, vb']';
qo = qmul(qmul(qbn, qi), qconj(qbn));

vn = qo(2:4);
end

