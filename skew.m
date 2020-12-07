function [ skew_matrix ] =skew( theta )
% 名称：Skew matrix
% 功能：3×1向量theta对应的反对称矩阵
%
% Inputs:
%       theta: (rad*vector) 3×1向量，通常是角增量或等效旋转矢量
% Outputs:
%       skew_matrix: (1) 反对称矩阵

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
skew_matrix=[       0, -theta(3),  theta(2)
             theta(3),         0, -theta(1)
            -theta(2),  theta(1),        0];
end

