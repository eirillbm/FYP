function A=calcA2(x)

%Update the A matrix
A=zeros(12);
A(1,7)=cos(x(6))*cos(x(5));
A(1,8)=-sin(x(6))*cos(x(4))+cos(x(6))*sin(x(5))*sin(x(4));
A(1,9)=sin(x(6))*sin(x(4))+cos(x(6))*cos(x(4))*sin(x(5));
A(2,7)=sin(x(6))*cos(x(5));
A(2,8)=cos(x(6))*cos(x(4))+sin(x(4))*sin(x(5))*sin(x(6));
A(2,9)=-cos(x(6))*sin(x(4))+sin(x(5))*sin(x(6))*cos(x(4));
A(3,7)=-sin(x(5));
A(3,8)=cos(x(5))*sin(x(4));
A(3,9)=cos(x(5))*cos(x(4));
A(4,10)=1;
A(4,11)=sin(x(4))*tan(x(5));
A(4,12)=cos(x(4))*tan(x(5));
A(5,11)=cos(x(4));
A(5,12)=-sin(x(4));
A(6,11)=sin(x(4))/cos(x(5));
A(6,12)=cos(x(4))/cos(x(5));