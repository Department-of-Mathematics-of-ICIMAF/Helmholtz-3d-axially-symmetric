function[curvenurbs]=abline(a,b)

knotscurve=[0 0 0 1 1 1];
coefscurve(1:4,1)=[a(1) a(2) 0 1];
coefscurve(1:4,2)=[(a(1)+b(1))/2 (a(2)+b(2))/2 0 1];
coefscurve(1:4,3)=[b(1) b(2) 0 1];
curvenurbs = nrbmak(coefscurve,knotscurve);