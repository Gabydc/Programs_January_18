 function [a,b,z]=matrixf(x,y,l,s0,s)


m=x*y*l-2*x;

n=x;
 I = speye(n,n);
 E = sparse(2:n,1:n-1,1,n,n);
 D = -s0*(E+E'-2*I);
 a1 = kron(D,I)+kron(I,D);

for i=x*(y-1)+1:x*y
 a1(i,i)=3*s0+s;
 a1(i,i+y)=-s;
end

a11(1:x*y,x+1:x*y+2*x)=a1;
 h=y;
 for i=1:x
 a1((i-1)*h+1,:)=a1((i-1)*h+1,:)./2;
 a1((i-1)*h+1,(i-1)*h+2)=a1((i-1)*h+1,(i-1)*h+2).*2;
 a1((i-1)*h+y,:)=a1((i-1)*h+y,:)./2;
 a1((i-1)*h+y,(i-1)*h+x-1)=a1((i-1)*h+y,(i-1)*h+x-1).*2;  
 end

   a=a1(x+1:x*y,x+1:x*y+x);
   z(1:size(a,1),1)=1;
 



for i=x+1:2*x
 a11(i-y,i)=3*s+s0; 
 a11(i-y,i-y)=-s;
 if (i~=x+1) & (i~=2*x)
   a11(i-y,i+1)=-s; 
   a11(i-y,i-1)=-s; 
 end
 
end


 h=y;
 for i=1:x
 a11((i-1)*h+1,:)=a11((i-1)*h+1,:)./2;
 
  a11((i-1)*h+1,(i-1)*h+y+2)=a11((i-1)*h+1,(i-1)*h+y+2)*2;
  
  a11((i-1)*h+y,:)=a11((i-1)*h+y,:)./2;

  a11((i-1)*h+y,(i-1)*h+(2*y-1))=a11((i-1)*h+y,(i-1)*h+2*y-1)*2; 
  
 end
 
  a11(x+1-y,x+2)=-s;
 a11(2*x-y,2*x-1)=-s; 

 I = speye(n,n);
 E = sparse(2:n,1:n-1,1,n,n);
 D = -s*(E+E'-2*I);
 a2 = kron(D,I)+kron(I,D);
 
for i=x*(y-1):x*y

 a2(i,i+x)=-s;

end

a22(:,x+1:x+size(a2,2))=a2;

for i=x+1:2*x

 a22(i-x,i-x)=-s;
end



 h=y;
 for i=1:x
 a22((i-1)*h+1,:)=a22((i-1)*h+1,:)./2;
 a22((i-1)*h+1,(i-1)*h+x+2)=a22((i-1)*h+1,(i-1)*h+x+2).*2;
 a22((i-1)*h+x,:)=a22((i-1)*h+x,:)./2;
 a22((i-1)*h+x,(i-1)*h+2*y-1)=a22((i-1)*h+x,(i-1)*h+2*y-1).*2;  
 end

sa1=size(a11);
sa2=size(a22);

for i=2:l
hx=size(a,1);
hy=size(a,2);

    if mod(i,2)==1
        a(hx+1:hx+sa1(1),hy+1-2*y:hy+sa1(2)-2*y)=a11;   
        z(hx+1:hx+sa1(1),i)=1;
        
    else
         a(hx+1:hx+sa2(1),hy-2*y+1:hy-2*y+sa2(2))=a22;
         z(hx+1:hx+sa2(1),i)=1;

    end
end
% spy(z)
z=sparse(z(1:size(a,1)-x,:));
 a=a(1:size(a,1)-x,1:size(a,2)-2*y);
a=a/((x-1)^(-2));
b=zeros(length(a),1);
b(1,1)=s0/(2*(x-1)^(-2));
b(x,1)=s0/(2*(x-1)^(-2));
b(2:x-1,1)=s0/((x-1)^(-2));

