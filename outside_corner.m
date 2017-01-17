% Calculates velocity field for corner flow in the outside corner of
% subduction. x0,y0 are the origin of dipping slab (xslab,litho).
% x1,y1 are the coordinate of the point to calculate the velocities for.
% t0 is the slabdip and U is the subduction velocity.

function [ux,uy] = outside_corner(x,y,t0,U);

A = -U.*t0./(t0 + sin(t0));
B = 0;
C = 2.*U.*cos(t0./2).^2./(t0 + sin(t0));
D = -U.*sin(t0)./(t0 + sin(t0));

r = sqrt(x.^2 + y.^2);
thet = acos(x./sqrt(x.^2 + y.^2));
%psi  = r.*(A.*sin(thet) + B.*cos(thet) + C.*thet.*sin(thet) + D.*thet.*cos(thet));

ur = A.*cos(thet) - B.*sin(thet) + C.*(sin(thet) + thet.*cos(thet)) + D.*(cos(thet) - thet.*sin(thet));
uth = -(A.*sin(thet) + B.*cos(thet) + C.*thet.*sin(thet) + D.*thet.*cos(thet));

ux = ur.*cos(thet) - uth.*sin(thet);
uy = ur.*sin(thet) + uth.*cos(thet);


