function [ vn ] = qmulv( qbn, vb )
% ���ƣ�Quaternion MULtiply Vector
% ���ܣ�����Ԫ������ά������������任
% ����@ �����ߵ�ϵͳ����ϵ���ԭ�� P234
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

