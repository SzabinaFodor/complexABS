function [x_solution, iflag] = S3rr( A, b )
%solve the linear system A*x=b by algorithms of orthogonally scaled ABS class
% v(i)=A*p(i)
% z(i)=e(i)
% w(i)=e(i)
%inputs:
%   A,b: the matrix and column vector in the system: A * x = b,
%        A = [a_1,...,a_m]';
%outputs:
%    x_solution: the solution of the system A*x = b;
%    flag:
%        0: succed to find a solution of A*x=b;
%       -i: fail to solved the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    global H
    global x
	global H
	global x
	keps=4*eps;
	iflag = 0;
    
	[n,m] = size(A);
	assert( n == m, 'A is not a square matrix!' );
	assert( n == m, 'b must have as many elements as the columns of A!' );


	for k = 1 : m
        r = A * x - b;
        z=r;
        p = H'*z; % z(i)=r(i)
        w=r';   %w=r(i)
    
    
        v=A*p;
        ptp = v' * v;
    
        HATv=H*(A'*v);
        if norm(HATv) > keps
            p=p/ptp;
            x=x-(v'*r)*p;
            wH=w*H;
            ww=wH/(wH*HATv);
            H = H - HATv * ww;  
        elseif norm(HATv) < keps && norm(ptp)< keps
			iflag=iflag+1;
		else
			iflag = -k;
			break;
        end
    end
    x_solution = x;
end