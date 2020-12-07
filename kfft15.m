function [ Ft ] = kfft15( eth, Cbn, fb )
% 名称：Ft matrix in state model for Kalman filter (15 dimension)
% 功能：SINS误差转移矩阵
% 程序@ 捷联惯导算法与组合导航原理 P237
%
% Inputs:
%       eth: 
%       Cbn: 
%       fb: 
% Outputs:
%       Ft: 

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
global g0
% 调整输入形式为列向量
if size(fb, 1) == 1
    fb = fb';
end

O33 = zeros(3, 3);

tl = eth.tl;
secl = 1/eth.cl;

f_RMh = 1/eth.RMh;
f_RNh = 1/eth.RNh;
f_clRNh = 1/eth.clRNh;
f_RMh2 = f_RMh*f_RMh;
f_RNh2 = f_RNh*f_RNh;

vE_clRNh = eth.vn(1)*f_clRNh;
vE_RNh2 = eth.vn(1)*f_RNh2;
vN_RMh2 = eth.vn(2)*f_RMh2;

%%
% P92 (4.2.33a)
Mp1 = [           0,    0,  0
       -eth.wien(3),    0,  0
        eth.wien(2),    0,  0];

% P92 (4.2.33c)
Mp2 = [            0,   0,      vN_RMh2
                   0,   0,     -vE_RNh2
       vE_clRNh*secl,   0,  -vE_RNh2*tl];

% P93 (4.2.38_1)
Maa = skew(-eth.winn);

% P92 (4.2.33b)
Mav = [       0,  -f_RMh,   0
          f_RNh,       0,   0
       f_RNh*tl,       0,   0];

% P93 (4.2.38_2)
Map = Mp1 + Mp2;

% ***P93 (4.2.40 1-3)
Mva = skew(Cbn*fb);
Mvv = skew(eth.vn)*Mav - skew(eth.wienn);

Mvp = skew(eth.vn)*(2*Mp1 + Mp2);
% 对M3的补偿项
scl = eth.sl*eth.cl;
Mvp(3, 1) = Mvp(3, 1) - g0*(5.27094e-3*2*scl + 2.32718e-5*4*eth.sl2*scl);
Mvp(3, 3) = Mvp(3, 3) + 3.086e-6;
% ***P93 (4.2.40 1-3) end

% P94 (4.2.42a)
Mpv = [      0,  f_RMh,  0
       f_clRNh,      0,  0
             0,      0,  1];

% P94 (4.2.42b)
Mpp = [          0,   0,      -vN_RMh2
       vE_clRNh*tl,   0, -vE_RNh2*secl
                 0,   0,             0];

%     phi  dvn  dpos   eb   db
Ft = [Maa, Mav, Map, -Cbn, O33
      Mva, Mvv, Mvp,  O33, Cbn
      O33, Mpv, Mpp,  O33, O33
      zeros(6, 15)];

end
