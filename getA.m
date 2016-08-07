function A = getA(IP,x,y,X,Y,Z,R,f,drw,drp,drk,n)
  format long
%%%%%%% input parameters to func getAMatrix: are: 
    %Initial Parameters IP
    %data values as matrices x,y,X,Y
    %Rotational matrix R
    %focal length f, with converted units (m)
    %differential elements of the rotations drw, drp, drk

	%calculate differentials that make the elements of the A matrix

	%%F1 derived coefficients
  % n = length(x) % where n is the number of data points. =len(y)=len(X)=len(Y)
  A = zeros(n*2,6);
	for m=1:n
		%index navigations
		a = m*2;
        %%elements from differentials of F1
		A(a-1,1) = ( -x(m)*R(3,1) - f*R(1,1) );
		A(a-1,2) = ( -x(m)*R(3,2) - f*R(1,2) );
		A(a-1,3) = ( -x(m)*R(3,3) - f*R(1,3) );
        
		A(a-1,4) = x(m)*( drw(3,2)*(Y(m)-IP(2)) + drw(3,3)*(Z(m)-IP(3)) ) + f*( drw(1,2)*(Y(m)-IP(2)) + drw(1,3)*(Z(m)-IP(3)) );  %XO = IP(1) YO = IP(2) ZO = IP(3)

		A(a-1,5) = x(m)*( drp(3,1)*(X(m)-IP(1)) + drp(3,2)*(Y(m)-IP(2)) + drp(3,3)*(Z(m)-IP(3)) ) + f*( drp(1,1)*(X(m)-IP(1)) + drp(1,2)*(Y(m)-IP(2)) + drp(1,3)*(Z(m)-IP(3)) );

		A(a-1,6) = f*( drk(1,1)*(X(m)-IP(1)) + drk(1,2)*(Y(m)-IP(2)) + drk(1,3)*(Z(m)-IP(3)) );

		%%elements from differentials of F2
		A(a,1) = ( -y(m)*R(3,1) - f*R(2,1) );
		A(a,2) = ( -y(m)*R(3,2) - f*R(2,2) );
		A(a,3) = ( -y(m)*R(3,3) - f*R(2,3) );
        
		A(a,4) = y(m)*( drw(3,2)*(Y(m)-IP(2)) + drw(3,3)*(Z(m)-IP(3)) ) + f*( drw(2,2)*(Y(m)-IP(2)) + drw(2,3)*(Z(m)-IP(3)) );
        
		A(a,5) = y(m)*( drp(3,1)*(X(m)-IP(1)) + drp(3,2)*(Y(m)-IP(2)) + drp(3,3)*(Z(m)-IP(3)) ) + f*( drp(2,1)*(X(m)-IP(1)) + drp(2,2)*(Y(m)-IP(2)) + drp(2,3)*(Z(m)-IP(3)) );
		
        A(a,6) = f*( drk(2,1)*(X(m)-IP(1)) + drk(2,2)*(Y(m)-IP(2)) + drk(2,3)*(Z(m)-IP(3)) );
end
