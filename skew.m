function [ skew_matrix ] =skew( theta )
% ���ƣ�Skew matrix
% ���ܣ�3��1����theta��Ӧ�ķ��Գƾ���
%
% Inputs:
%       theta: (rad*vector) 3��1������ͨ���ǽ��������Ч��תʸ��
% Outputs:
%       skew_matrix: (1) ���Գƾ���

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
skew_matrix=[       0, -theta(3),  theta(2)
             theta(3),         0, -theta(1)
            -theta(2),  theta(1),        0];
end

