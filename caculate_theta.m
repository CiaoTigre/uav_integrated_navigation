function [ theta, omega ] = caculate_theta( ang_motion, crt )
% 名称：Caculate theta(angule) form the seted angule motion
% 功能：给定预设角运动形式，计算current时刻的角度
% 注意：目前还没有添加规范化姿态角的功能，因此尽量不要把摇摆幅值设的过大
%
% Inputs:
%       ang_motion: 预设角运动，每一行对应着一个形如：
%                   theta = A*sin( w*t + phi) + k的运动，载体角运动由若干个
%                   正弦运动叠加而成
%       crt: 当前时间
% Output:
%       theta: 跟预设角运动所对应的姿态 (rad)
%       omega: 与预设角运动所对应的当前时刻角速度 (rad/s)

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% i表示 ang_motion 中的第i行
i = 1;

% 如果ang_motion中存在第i行
while i <= size(ang_motion, 1)
    
    if i == 1
        theta = 0;  % rad
        omega = 0;  % rad/s
    end
    
    % ang_motion.p = [A1(幅值), w1(角速度), phi1(初始相位), k1(中心值)
    %                 A2(幅值), w2(角速度), phi2(初始相位), k2(中心值)
    %                                    ...                         ]
    % 注意:为了编程的简单和简洁，ang_motion中的变量单位为rad或rad/s
    A = ang_motion(i, 1);
    w = ang_motion(i, 2);
    phi = ang_motion(i, 3);
    k = ang_motion(i, 4);
    
    % 姿态角
    theta = A*sin( w*crt + phi) + k + theta;
    
    % wnbnx    nx分别是(n n1), (n1 n2), (n1 b)
    omega = w*A*cos(w*crt + phi) + omega;
    
    % i+1, 准备检查ang_motion是否存在下一行
    i = i+1;
    
end

end
