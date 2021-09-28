function doit
    addpath('abs');
    addpath('testproblems');
    

    global x;
    global H;
 
    [A, m, n]=loadTestProblem('qc324.mat');
    y = generatecomplex(n);
	b=A*y;
     
    disp('complex ABS orthogonally scaled subclass');
    H = eye(n);
	x = zeros(n,1);

    [solution, iflag]=S3rr(A,b);
    norm(A*solution-b)
end

function [A, m, n]=loadTestProblem(name)
    S = load(name);
    A = S.Problem.A;
    [m,n]=size(A);
end