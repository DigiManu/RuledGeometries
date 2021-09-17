%   Copyright 2021 Kyoto University
%   Author: Anthony Beaucamp
%   Last modified: 2021/09/17

%   This script provides an example of cutting wire trajectory calculation
%   for a ruled surface expressed as faces and vertices.

% Create a ruled triangular mesh
[X,Y,Z] = Ruled(15,15,30);
FV = surf2patch(X,Y,Z);

% Compute curvature of triangular mesh data
[Cmean,Cgaussian,Dir1,Dir2,Lambda1,Lambda2] = Curvature(FV);
C.Kmin = reshape(Lambda1,30,30);
C.Kmax = reshape(-Lambda2,30,30);
C.KminU = reshape(Dir2(:,1),30,30);
C.KminV = reshape(Dir2(:,2),30,30);
C.KmaxU = reshape(Dir1(:,1),30,30);
C.KmaxV = reshape(Dir1(:,2),30,30);

% Force min. curvature to be negative
pos = (C.KminU > 0);
C.KminU(pos) = -C.KminU(pos);
C.KminV(pos) = -C.KminV(pos);
C.KmaxU(pos) = -C.KmaxU(pos);
C.KmaxV(pos) = -C.KmaxV(pos);

% Find zero curvature direction angles
K = @(th) (C.Kmax .* cos(th).^2 + C.Kmin .* sin(th).^2);
guess = zeros(size(C.Kmin));
theta = lsqnonlin(K, guess, guess-pi, guess+pi, optimoptions('lsqnonlin','Display','off','TolFun',1e-20));

% Fix angles (to generate continuous solution)
theta(theta<0) = - theta(theta<0);
theta(theta>pi/2) = pi - theta(theta>pi/2);

% Compute zero curvature direction vectors
ZeroU1 =  cos(-theta) .* C.KmaxU - sin(-theta) .* C.KmaxV;
ZeroV1 =  sin(-theta) .* C.KmaxU + cos(-theta) .* C.KmaxV;
ZeroU2 = -cos(theta)  .* C.KmaxU + sin(theta) .* C.KmaxV;
ZeroV2 = -sin(theta)  .* C.KmaxU - cos(theta) .* C.KmaxV;

% Display curvature analysis
figure; 
subplot(2,2,1); surf(C.Kmin); title('Min Curvature');
subplot(2,2,2); surf(C.Kmax); title('Max Curvature');
subplot(2,1,2); surf(X,Y,Z, 'facealpha', 0.33); shading interp; daspect([1 1 1]);
hold on; quiver3(X,Y,Z,C.KminU,C.KminV,zeros(size(C.KminU))); 
hold on; quiver3(X,Y,Z,C.KmaxU,C.KmaxV,zeros(size(C.KminU)),'color','r'); 
title('Min (blue) / Max (red) Curvature Direction');
view(0,90);

% Display zero curvature direction
figure; 
subplot(2,2,1); surf(X,Y,Z, 'facealpha', 0.33); shading interp; daspect([1 1 1]);
hold on; quiver3(X,Y,Z,ZeroU1,ZeroV1,zeros(size(C.KminU))); 
title('Zero Curvature Direction 1');
subplot(2,2,2); surf(X,Y,Z, 'facealpha', 0.33); shading interp; daspect([1 1 1]);
hold on; quiver3(X,Y,Z,ZeroU2,ZeroV2,zeros(size(C.KminU))); 
title('Zero Curvature Direction 2');

% Display ruled trajectories
subplot(2,1,2); 
surf(X,Y,Z, 'facealpha', 0.33); 
shading interp; daspect([1 1 1]);
title('Wire Trajectory');
view(0,90);

% Prepare scattered interpolants
fZ  = scatteredInterpolant(X(:),Y(:),Z(:));
fU1 = scatteredInterpolant(X(:),Y(:),ZeroU1(:));
fV1 = scatteredInterpolant(X(:),Y(:),ZeroV1(:));
fU2 = scatteredInterpolant(X(:),Y(:),ZeroU2(:));
fV2 = scatteredInterpolant(X(:),Y(:),ZeroV2(:));

% Starting Point Loop
for j = [2:30]
    % Initialize
    sX1 = []; sY1 = [];
    sX2 = []; sY2 = [];
    sX1(1) = X(j,1);
    sY1(1) = Y(j,1);
    sX2(1) = X(j,1);
    sY2(1) = Y(j,1);

    % Trajectory Loop
    while (sX1(end) < 15 && ~isnan(sX1(end)))
        sX1(end+1) = sX1(end) + fU1(sX1(end),sY1(end));
        sY1(end+1) = sY1(end) + fV1(sX1(end),sY1(end));
    end;
      
    while (sX2(end) < 15 && ~isnan(sX2(end)))
        sX2(end+1) = sX2(end) + fU2(sX2(end),sY2(end));
        sY2(end+1) = sY2(end) + fV2(sX2(end),sY2(end));
    end;
    
    % Compute Z coords
    sZ1 = fZ(sX1,sY1);
    sZ2 = fZ(sX2,sY2);

    % Remove NaN values
    mask1 = ~isnan(sX1) & ~isnan(sY1) & ~isnan(sZ1);
    mask2 = ~isnan(sX2) & ~isnan(sY2) & ~isnan(sZ2);
    sX1 = sX1(mask1); sY1 = sY1(mask1); sZ1 = sZ1(mask1);
    sX2 = sX2(mask2); sY2 = sY2(mask2); sZ2 = sZ2(mask2);

    % Compute vectors & angles
    vX1 = diff(sX1); vY1 = diff(sY1); 
    vX2 = diff(sX2); vY2 = diff(sY2);
    theta1 = atan2(vY1,vX1);
    theta2 = atan2(vY2,vX2);

    % Compute STD of curve direction
    StdDev1 = sqrt(mean((theta1-mean(theta1)).^2));
    StdDev2 = sqrt(mean((theta2-mean(theta2)).^2));

    % Show straighter trajectory
    if StdDev1 < StdDev2
      % Fit line (to improve accuracy of numerical approx.)
      sY1 = polyval(polyfit(sX1,sY1,1), sX1);
      sZ1 = polyval(polyfit(sX1,sZ1,1), sX1);      
      hold on; plot3(sX1, sY1, sZ1, 'b');
    else
      % Fit line (to improve accuracy of numerical approx.)
      sY2 = polyval(polyfit(sX2,sY2,1), sX2);
      sZ2 = polyval(polyfit(sX2,sZ2,1), sX2);  
      hold on; plot3(sX2, sY2, sZ2, 'b');
    end;
 end;
 