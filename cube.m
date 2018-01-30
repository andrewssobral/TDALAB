function FV=cube(p0,r)
% create a cube for patch.
% p0, r: center and radius of the cube.
% patch(FV) shows the cube
if numel(r)==1
    r=r(ones(3,1));
end
x=[p0(1)-r(1),p0(1)+r(1)];
y=[p0(2)-r(2),p0(2)+r(2)];
z=[p0(3)-r(3),p0(3)+r(3)];

FV.vertices=[x([1 2 2 1 1 2 2 1]) ;
    y([1 1 2 2 1 1 2 2]) ;
    z([1 1 1 1 2 2 2 2])]';

FV.faces=[1 2 3 4
    5 6 7 8
    1 2 6 5
    4 3 7 8
    1 5 8 4
    2 6 7 3];

FV.facevertexcdata=rand(length(FV.vertices),3);
FV.facecolor='interp';