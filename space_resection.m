clc

% DATA MATRIX

x= [
	-0.064097
	0.064017
	-0.023087
	0.064353
];

y=[
	0.072924
	0.077183
	-0.039191
	-0.037129
];

X=[
	2294.96
	2274.05
	2736.69
	2734.89
];

Y=[
	4735.29
	5263.62
	4928.94
	5269.75
];

Z=[
	1236.46
	1242.56
	1313.25
	1274.90
];

%focal length of camera
f = 0.11205;

%initial parameters
XO=2597.46; YO=4997.81; ZO=1701.76; wo=0;  po=0;  ko=degtorad(89.6387);

IP = [
	XO
	YO
	ZO
	wo
	po
	ko
];

%Rotational elements rij
r11 = cos(po)*cos(ko);  
r12 = cos(wo)*sin(ko) + sin(wo)*sin(po)*cos(ko);
r13 = sin(wo)*sin(ko) - cos(wo)*sin(po)*cos(ko);

r21 = -cos(po)*sin(ko);
r22 = cos(wo)*cos(ko) + cos(wo)*sin(po)*sin(ko);
r23 = sin(wo)*cos(ko) + cos(wo)*sin(po)*sin(ko);

r31 = sin(po);
r32 = -sin(wo)*cos(po);
r33 = cos(wo)*cos(po);

R = [
	r11 r12 r13
	r21 r22 r23
	r31 r32 r33
];

	% consider using elements below for further code:

	a1 = cos(wo);
	a2 = cos(po);
	a3 = cos(ko);  %cosines of rotational

	b1 = sin(wo);
	b2 = sin(po);
	b3 = sin(ko);  %sines of rotational

%substitutional computations, F1o, F2o

f1O = x(1,1)*( r31*(X(1,1)-XO) + r32*(Y(1,1)-YO) + r33*(Z(1,1)-ZO) ) + f * ( r11*(X(1,1)-XO) + r12*(Y(1,1)-YO) + r13*(Z(1,1)-ZO) );

f2O = y(1,1)*( r31*(X(1,1)-XO) + r32*(Y(1,1)-YO) + r33*(Z(1,1)-ZO) ) + f * ( r21*(X(1,1)-XO) + r22*(Y(1,1)-YO) + r23*(Z(1,1)-ZO) );

%static L matrix 
% L = [
% 	f1O
% 	f2O
% ]

%differentials of rotational elements
dr11w =  0;                  dr11p = -b2*a3;       dr11k = -a2*b3;
dr12w = -b1*b3 + a1*b2*a3;   dr12p =  b1*a2*a3;    dr12k =  a1*a3 - b1*b2*b3;
dr13w =  a1*b3 + b1*b2*a3;   dr13p = -a1*a2*a3;    dr13k =  b1*a3 + a1*b2*b3;

dr21w =  0;                  dr21p =  b2*b3;       dr21k = -a2*a3;
dr22w = -b1*a3 - a1*b2*b3;   dr22p = -b1*a2*b3;    dr22k = -a1*b3 - b1*b2*a3;
dr23w =  a1*a3 - b1*b2*b3;   dr23p =  a1*a2*b3;    dr23k = -b1*b3 + a1*b2*a3;

dr31w =  0;                  dr31p =  a2;          dr31k = 0;
dr32w = -a1*a2;              dr32p =  b1*b2;       dr32k = 0;
dr33w = -b1*a2;              dr33p = -a1*b2;       dr33k = 0;

drw = [
	dr11w dr12w dr13w
	dr21w dr22w dr23w
	dr31w dr32w dr33w
];

drp = [
	dr11p dr12p dr13p
	dr21p dr22p dr23p
	dr31p dr32p dr33p
];

drk = [
	dr11k dr12k dr13k
	dr21k dr22k dr23k
	dr31k dr32k dr33k
];

%number of data points n :
n = length(x); %  = len(y) =len(X) = len(Y) =len(Z)!

%%%iterative solution using the three data points AB and C to solve for the six unknowns
iter = 0; %iteration number #	
DXI = ones(6,1);

while DXI~=zeros(6,1) && iter<=1e2  %loop to convergence using DXI=ZEROS  comparison
    Au = getUnique(IP,x,y,X,Y,Z,R,f,drw,drp,drk,n)
    L=getL(IP,x,y,X,Y,Z,R,f,3);
	DXI = (Au'*L);
	IP = IP-DXI		
	iter=iter+1;
end	

% for m=1:1    % 50 iterations
% 	A=getA(IP,x,y,X,Y,Z,R,f,drw,drp,drk,n);   %call A matrix Func
% 	L=getL(IP,x,y,X,Y,Z,R,f,n);  %call L matrix Func
%     Au = getUnique(IP,x,y,X,Y,Z,R,f,drw,drp,drk,n);
% 	%DX = getdx(A,L);  %call dx matrix Func,, Least Square Adjustment
% 	%IP = IP+DX;
% end