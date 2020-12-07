function [ eth ] = earth( pos, vn )
% ���ƣ�EARTH
% ���ܣ�
% ����@ �����ߵ��㷨����ϵ���ԭ�� P235

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
% �������������Ϊ������
if size(pos, 1) == 1
    pos = pos';
end
if size(vn, 1) == 1
    vn = vn';
end

global Re ff wie g0

% ��һƫ����
ee = sqrt(2*ff - ff^2);
e2 = ee^2;

% ��γ���йص����Ǻ���
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

% ������ʹ�����½Ǳ��ʾ��������ʱ��ֱ�Ӳ�����Ϥ�ķ�ʽ����ͶӰ����ϵ
% ���ϣ�����˶�����������ϵ���¡�
% eth.win^n = eth.wie^n + eth.wen^n  => eth.winn = eth.wien + eth.wenn;
eth.vn = vn;

eth.wien = wie*[0, eth.cl, eth.sl]';
eth.wenn = [-vn(2)/eth.RMh, vn(1)/eth.RNh, vn(1)/eth.RNh*eth.tl]';
eth.winn = eth.wien + eth.wenn;
% wienn = wien + winn ,��ʾ�к����ٶ�����������ٶ�֮��
eth.wienn = eth.wien + eth.winn;

% grs80 ����ģ��
gLh = g0*(1 + 5.27094e-3*eth.sl2 + 2.32718e-5*sl4) - 3.086e-6*pos(3);

eth.gn = [0, 0, -gLh]';
% �����������������ٶ����ʽ�����ٶ�֮��
eth.gcc = eth.gn - cross(eth.wienn, vn);

end
