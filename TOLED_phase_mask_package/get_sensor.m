function [sensor] = get_sensor(params)

sensor.M = params.sensorRes;
sensor.x = (0:1:sensor.M-1) - floor(sensor.M/2);
t = params.delta2 / params.delta1;
sensor.pixel2angle = @(pid) 90 - 180/pi*acos(pid*params.sensorPitch*t/params.f);
% sensor.angle2pixel = @(theta) round(cos(pi/180*(90-theta))*M*params.delta1/params.scale/lambda);
sensor.angle2pixel = @(theta) round(cos(pi/180*(90-theta))/params.sensorPitch/t*params.f);
