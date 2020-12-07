function [ m ] = rv2m( rv )
%% **************************************************************
% ���ƣ�Rotation Vector to direction cosine Matrix
% ���ܣ��õ�Ч��תʸ�����㷽�����Ҿ���
%       e.g.: rv = PHIibi ��ʾiϵ��ת����uibi(uibb)ת����phi������bϵ�غ�
% Inputs:
%       rv: 3��1 vector, ��ģ��ʾת���ǶȵĴ�С��rv/norm(rv) Ϊת����������
%          ����ϵ�µ����ꡣ 
% Outputs:
%       m: Cbi

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% ����@ �����ߵ��㷨����ϵ���ԭ���� P233
% ��Ч��תʸ����ģ����
nm2 = rv'*rv ;

if nm2 < 1.e-8
    a = 1 - nm2*(1/6 - nm2/120);
    b = 0.5 - nm2*(1/24 - nm2/720);
else
    nm = sqrt(nm2);
    a = sin(nm)/nm;
    b = (1 - cos(nm))/nm2;
end

VX = skew(rv);
m = eye(3) + a*VX + b*VX*VX;
end
