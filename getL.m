%compute and form an L matrix
function L = getL(IP,x,y,X,Y,Z,R,f,n)
	% n = length(x)   % where n is the number of data points. =len(y)=len(X)=len(Y)
	L = zeros(n*2,1);
	for m=1:n
		a=m*2; %index navigations
		L(a-1) = -(x(m)*( R(3,1)*(X(m)-IP(1)) + R(3,2)*(Y(m)-IP(2)) + R(3,3)*(Z(m)-IP(3)) ) + f * ( R(1,1)*(X(m)-IP(1)) + R(1,2)*(Y(m)-IP(2)) + R(1,3)*(Z(m)-IP(3)) ));
    L(a) = -(y(m)*( R(3,1)*(X(m)-IP(1)) + R(3,2)*(Y(m)-IP(2)) + R(3,3)*(Z(m)-IP(3)) ) + f * ( R(2,1)*(X(m)-IP(1)) + R(2,2)*(Y(m)-IP(2)) + R(2,3)*(Z(m)-IP(3)) ));
		
end