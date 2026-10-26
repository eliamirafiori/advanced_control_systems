robot = importrobot('00001_PRR_zzy.urdf');
showdetails(robot)

figure;
config = homeConfiguration(robot);
show(robot,config);
