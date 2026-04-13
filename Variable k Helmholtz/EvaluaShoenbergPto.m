function[cxy,alphaxy,rhoxy]=EvaluaShoenbergPto(x,y,knots,ccoef,alphacoef,rhocoef)

[n,m]=size(x);
cxy=zeros(n,m);
alphaxy=zeros(n,m);
rhoxy=zeros(n,m);

for i=1:n
  for j=1:m
   % Evaluamos los B-splines B_i^k(y) en los ptos y(i,:) de la fila i de y 
   Bmat = spcol(knots,3,y(i,j));
   % Evaluamos el spline Vc en los ptos y(:,j), Vc=Bmat*ccoef;
   cxy(i,j)=Bmat*ccoef;
   alphaxy(i,j)=Bmat*alphacoef;
   rhoxy(i,j)=Bmat*rhocoef;
  end
end