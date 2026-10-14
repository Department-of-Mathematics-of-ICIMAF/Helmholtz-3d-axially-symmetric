function [xx,yy,zz,u] = cylinder_xyz(varargin)
%CYLINDER Generate cylinder.
%   [X,Y,Z] = CYLINDER(R,N) forms the unit cylinder based on the generator
%   curve in the vector R. Vector R contains the radius at equally
%   spaced points along the unit height of the cylinder. The cylinder
%   has N points around the circumference. SURF(X,Y,Z) displays the
%   cylinder.
%
%   [X,Y,Z] = CYLINDER(R), and [X,Y,Z] = CYLINDER default to N = 20
%   and R = [1 1].
%
%   Omitting output arguments causes the cylinder to be displayed with
%   a SURF command and no outputs to be returned.
%
%   CYLINDER(AX,...) plots into AX instead of GCA.
%
%   See also SPHERE, CYLINDER.
 

% Parse possible Axes input
narginchk(0,5);
[cax,args,nargs] = axescheck(varargin{:});

n = 20;
r = [1 1]';
if nargs > 0, r = args{1}; end
if nargs > 1, n = args{2}; end
if nargs > 2, z_value = args{3}; end
if nargs > 3
    u_value = args{4};
    u = u_value*ones(1,n+1);
end

    
r = r(:); % Make sure r is a vector.
m = length(r); if m==1, r = [r;r]; m = 2; end

if nargs > 4
    theta_end = args{5};
    theta = linspace(0,theta_end,n+1);
else
    theta = (0:n)/n*2*pi;
end

sintheta = sin(theta);

x = r * cos(theta);
y = r * sintheta;
if nargs > 2
    z = z_value*ones(1,n+1);
else
    z = (0:m-1)'/(m-1) * ones(1,n+1);
end

if nargout == 0
    cax = newplot(cax);
    surf(x,y,z,'parent',cax)
else
    xx = x; yy = y; zz = z;
end
