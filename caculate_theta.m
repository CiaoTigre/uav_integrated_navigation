function [ theta, omega ] = caculate_theta( ang_motion, crt )
% ���ƣ�Caculate theta(angule) form the seted angule motion
% ���ܣ�����Ԥ����˶���ʽ������currentʱ�̵ĽǶ�
% ע�⣺Ŀǰ��û����ӹ淶����̬�ǵĹ��ܣ���˾�����Ҫ��ҡ�ڷ�ֵ��Ĺ���
%
% Inputs:
%       ang_motion: Ԥ����˶���ÿһ�ж�Ӧ��һ�����磺
%                   theta = A*sin( w*t + phi) + k���˶���������˶������ɸ�
%                   �����˶����Ӷ���
%       crt: ��ǰʱ��
% Output:
%       theta: ��Ԥ����˶�����Ӧ����̬ (rad)
%       omega: ��Ԥ����˶�����Ӧ�ĵ�ǰʱ�̽��ٶ� (rad/s)

% Author: Kun Gan, Tongji University
% Email : ciaotigre@126.com
% Date  : 2020/12/1
%%
% i��ʾ ang_motion �еĵ�i��
i = 1;

% ���ang_motion�д��ڵ�i��
while i <= size(ang_motion, 1)
    
    if i == 1
        theta = 0;  % rad
        omega = 0;  % rad/s
    end
    
    % ang_motion.p = [A1(��ֵ), w1(���ٶ�), phi1(��ʼ��λ), k1(����ֵ)
    %                 A2(��ֵ), w2(���ٶ�), phi2(��ʼ��λ), k2(����ֵ)
    %                                    ...                         ]
    % ע��:Ϊ�˱�̵ļ򵥺ͼ�࣬ang_motion�еı�����λΪrad��rad/s
    A = ang_motion(i, 1);
    w = ang_motion(i, 2);
    phi = ang_motion(i, 3);
    k = ang_motion(i, 4);
    
    % ��̬��
    theta = A*sin( w*crt + phi) + k + theta;
    
    % wnbnx    nx�ֱ���(n n1), (n1 n2), (n1 b)
    omega = w*A*cos(w*crt + phi) + omega;
    
    % i+1, ׼�����ang_motion�Ƿ������һ��
    i = i+1;
    
end

end
