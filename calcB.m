function B=calcB(x)

%Update the B matrix
B=zeros(12,6);
B(7,1)=cos(x(6))*cos(x(5));
B(7,2)=sin(x(6))*cos(x(5));
B(7,3)=-sin(x(5));
B(8,1)=-sin(x(6))*cos(x(4))+cos(x(6))*sin(x(5))*sin(x(4));
B(8,2)=cos(x(6))*cos(x(4))+sin(x(4))*sin(x(5))*sin(x(6));
B(8,3)=cos(x(5))*sin(x(4));
B(9,1)=sin(x(6))*sin(x(4))+cos(x(6))*cos(x(4))*sin(x(5));
B(9,2)=-cos(x(6))*sin(x(4))+sin(x(5))*sin(x(6))*cos(x(4));
B(9,3)=cos(x(5))*cos(x(4));
B(10,4)=1;
B(10,6)=-sin(x(5));
B(11,5)=cos(x(4));
B(11,6)=cos(x(5))*sin(x(4));
B(12,5)=-sin(x(4));
B(12,6)=cos(x(5))*cos(x(4));