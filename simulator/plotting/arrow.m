function [x, y] =  arrow(base, dir, scalex, scaley)
%{

    Takes in position of centre of base of triangle and direction
    and scale in x and y and outputs coordinates for corners of 
    triangle. If scalex = scaley then this is an equilateral triangle.

%}

    dir = dir/norm(dir);
    r = sqrt(3)/2;

    x1 = (base(1) + dir(1)*r*scalex);
    y1 = (base(2) + dir(2)*r*scaley);

    x2 = (base(1) + dir(2)/2*scalex);
    y2 = (base(2) - dir(1)/2*scaley);

    x3 = (base(1) - dir(2)/2*scalex);
    y3 = (base(2) + dir(1)/2*scaley);

    x = [x1, x2, x3];
    y = [y1, y2,y3];

end