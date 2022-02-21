function x=Solve(A,b)
% See ichol, pcg and minres:
% https://es.mathworks.com/help/matlab/ref/ichol.html
% https://es.mathworks.com/help/matlab/ref/minres.html
% https://es.mathworks.com/help/matlab/ref/pcg.html
% MUMPS, PETSc
n=size(A,1);
mIter=ceil(size(A,1)/5);
Tol=1e-12;
fl=0;
if n<30000
    x=A\b;
else
    alpha = max(sum(abs(A),2)./diag(A))-2;
    L = ichol(A, struct('type','ict','droptol',1e-3,'diagcomp',alpha));
    [x,fl] = pcg(A,b,Tol,mIter,L,L'); % 1200 it, 17 sec
end
if fl>0 || norm(A*x-b)/norm(b)>1e-10
    warning('No convergence in pcg. Trying minres.')
   [x,fl]=minres(A,b,Tol,mIter); 
end
if fl>0
    warning('Error in iterative solver. Trying direct. May be memory consuming.')
    x=A\b;
end
end

function Test(A,b,opt)
switch(opt)
    case(-1) % Direct Matlab with reordering. Output for 33600 and 67200
        tic;
        p=symrcm(A);
        x=A(p,p)\b;
        x(p)=x;
        t=toc; % 2 sec, 0.2GB / 36 sec, 3.2 GB
    case(0) % Direct Matlab no reordering
        tic;
        %p=symrcm(A);
        %x=A(p,p)\b;
        %x(p)=x;
        x=A\b;
        t=toc; % 2 sec, 0.2GB / 37 sec, 3.1 GB
    case(1) % Conjug gradient with no prec -> Slow convergence      
        tic;
        [x,fl,rr,it,rv] = pcg(A,b,Tol,mIter);   % 1500 it, 16 sec / 5160, 139
        t1=toc;
    case(2) % Conjug gradient with prec michol -> Nonpositive pivot
        tic;
        opts.type = 'nofill';
        opts.michol = 'on';
        L = ichol(A,opts); 
        [x,fl2,rr2,it2,rv2] = pcg(A,b,Tol,mIter,L,L');        
        t2=toc;
    case(3) % Better than no prec, but still slow.
        tic;
        alpha = max(sum(abs(A),2)./diag(A))-2;
        L = ichol(A, struct('type','ict','droptol',1e-3,'diagcomp',alpha));
        [x,fl3,rr3,it3,rv3] = pcg(A,b,Tol,mIter,L,L'); % 1200 it, 17 sec / 3917, 154
        t3=toc;
    case(4) % minres. Linear convergence
        tic;
        [x,fl4,rr4,it4,rv4]=minres(A,b,Tol,mIter); % 1460 it, 15 sec / 4655, 129, fl4=3
        t4=toc;
    case(5) % pcg with diagonal precondicioner.
        tic;
        L=diag(diag(A));
        [x,fl5,rr5,it5,rv5]=pcg(A,b,Tol,mIter,L); % / 4549 141
        t5=toc;
    case(6) % gmres. Long and memory consuming
        tic;
        [x,fl6,rr6,it6,rv6]=gmres(A,b,n/10,Tol,mIter); % 1460 it, 15 sec / 4655, 129, fl4=3
        t6=toc;
end
end