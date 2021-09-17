%   Copyright 2021 Kyoto University
%   Author: Anthony Beaucamp
%   Last modified: 2021/09/17

function [X,Y,Z] = Ruled(diameter,width,length)

    % Generate Square
    xs = linspace(-length/2,-length/2,31);
    ys = linspace(-width/2,width/2,31);    
    zs = linspace(width/2,width/2,31);    
    
    % Generate Circle
    ac = linspace(-pi/4,pi/4,30);
    xc = linspace(length/2,length/2,30);
    yc = 0.5*diameter*sin(ac);
    zc = 0.5*diameter*cos(ac);
    
    % Generate Grid
    X = zeros(30);
    Y = zeros(30);
    Z = zeros(30);
    for i = [1:30],
        X(i,:) = linspace(xs(i),xc(i),30);
        Y(i,:) = linspace(ys(i),yc(i),30);
        Z(i,:) = linspace(zs(i),zc(i),30);
    end;
    