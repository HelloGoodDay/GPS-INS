format long;
load data.mat;

t = data0(:,1);
lon = 114.0248136140;
lat = 30.541093;
h = 76.3604;

% 整段数据平均后算一个姿态角
data = data0(1:12000, 2:7);
%data = data0(45000:50000, 2:7);
%[data, nepoch] = fitting(data); % 移除Y方向的偏移，改正到0时刻即15060.00
%[data, dt] = state_check(data); % 去除粗差
%[angle_whole] = StaticAlig(data, dt, lat, lon, h);
data = sum(data);
dt = 12000;
[angle_whole] = Alig2(data, dt, lat, lon, h);

% 每秒平均值计算并绘图
t_start = t(1); it_s = 1;
t_end = t(1);   it_e = 1;
n = 1;
for it = 1:a(1)
    if(it == 1)
        t_start = 1;
    else
        t_start = t_end;
        it_s = it_e;
    end
    if(it == a(1))
        it_e = a(1);
        t_end = t(it_e);
    elseif(abs(mod(t(it), 6)) < 0.005)
        it_e = it;
        t_end = t(it);
    else
        continue;
    end
    data = sum(data0(it_s:it_e-1, 2:7));
    dt = t(it_e) - t(it_s);
    if(abs(dt) < 1e-5) 
        continue;
    end
    angle = StaticAlig(data, dt, lat, lon, h);
    ts(n) = t_start;
    angles(n,:) = angle;
    n = n + 1;
end
paint( ts, angles, 'per second', 'second');
clear ts angles n;
