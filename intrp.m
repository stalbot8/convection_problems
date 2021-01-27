function ym = intrp(x1,xm,x2,y1,y2)

ym = y1 + (xm - x1)*(y2 - y1)/(x2 - x1);

end