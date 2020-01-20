function [mol, FORCE, Etot] = NLCG(mol,NLCG_opts,int_cord)

% Amartya S. Banerjee and Chao Yang
% Computational Research Division, Lawrence Berkeley National Lab



natoms = sum(mol.natoms);

imax = NLCG_opts.outiter_max ;
jmax = NLCG_opts.ineriter_max ;
n = NLCG_opts.n ;
nlcg_tol = NLCG_opts.NLCG_tol1;
tol1 = nlcg_tol * natoms * 3 ;

tol2 = NLCG_opts.NLCG_tol2 ;
sigma_0 = NLCG_opts.step_size0 ;

xyzlist = mol.xyzlist;
x = xyzlist(:);
%{
x1 = repmat(xyzlist(1,:),natoms,1);
x = xyzlist - x1;
if natoms>=2
    x(2,1) = norm(x(1,:)-x(2,:));
    x(2,2) = 0;
    x(2,3) = 0;
    if natoms>=3
        x(3,1) = (x(3,:)*x(2,:)')/x(2,1);
        x(3,2) = sqrt(x(3,:)*x(3,:)'-x(3,1)^2);
        x(3,3) = 0;
        if natoms>=4
            x(4:end,1) = x(4:end,:)*x(2,:)' / x(2,1);
            x(4:end,2) = (x(3,:)*x(4:end,:)'-x(3,1)*x(4:end,1)) / x(3,2);
            x(4:end,3) = sqrt(diag(x(4:end,:)*x(4:end,:)')'-x(4:end,1)^2-x(4:end,2).^2);
        end
    end
end
x = x(:);
%}     
    
ksopts = setksopt;
ksopts = setksopt(ksopts,'maxscfiter',100,'scftol',1e-6,'betamix',0.08,'cgtol',1e-10,'maxcgiter',20);


fcalls = 0;

i = 0;
k = 0;

fprintf('\n \n Going for force evaluation type 1 (outside i loop). \n\n');
[g,~] = ksfg(x,mol,ksopts,int_cord); ; fcalls = fcalls + 1;
fprintf('\n \n Force evaluation type 1, norm(g) = %-16.10e.', norm(g));
fprintf('\n Maximum force component = %-16.10e.\n\n', max(abs(g)));

r = -g ;
M = ones(size(x)) ; % Need to set REAL PRECONDITIONER
s = r./M ;
d = s ;
delta_new = dot(r,d) ;
% delta_0 = delta_new ;

if (imax == 0)
    fprintf('\n \n Returning after single shot calculation ... \n\n');
    v = x;
    return;
end

while i<imax && delta_new>tol1
    
    j = 0 ;
    delta_d = dot(d,d) ;
    alph = -sigma_0 ;
    
    fprintf('\n \n Going for force evaluation type 2 (inside i loop, before j loop).');
    fprintf('\n Iteration i = %d.\n\n', i);
    [g,info] = ksfg(x+sigma_0*d,mol,ksopts,int_cord); fcalls = fcalls + 1;
    fprintf('\n \n Force evaluation type 2, norm(g) = %-16.10e \n\n', norm(g));
    
    eta_prev = dot(g,d) ;
    
    while (j < jmax) && (alph*alph*delta_d > tol2)
        
        fprintf('\n \n Going in for force evaluation type 3 (inside j-loop).');
        fprintf('\n Iteration i = %d, j = %d.\n\n', i, j);
        
        
        [g,~] = ksfg(x,mol,ksopts,int_cord); fcalls = fcalls + 1;
        fprintf('\n \n Force evaluation type 3, norm(g) = %-16.10e \n\n', norm(g)); 
        
        eta = dot(g,d) ;
        alph = alph*(eta/(eta_prev-eta)) ;
        x = x + alph*d ;
        
        eta_prev = eta ;
        j = j + 1 ;
    end
    
    
    % fprintf('\n No. of inner CG iterations = %d. \n', j);
    fprintf('\n \n Going in for force evaluation type 4 (inside i loop, after j-loop).');
    fprintf('\n Iteration i = %d.\n\n', i);
    
    [g,~] = ksfg(x+sigma_0*d,mol,ksopts,int_cord); fcalls = fcalls + 1;
    fprintf('\n \n Force evaluation type 4, norm(g) = %-16.10e \n\n', norm(g));
    
    r = -g;
    delta_old = delta_new ;
    delta_mid = dot(r,s) ;
    
    M = ones(size(x)) ; % Need to set REAL PRECONDITIONER
    
    s = r./M ;
    delta_new = dot(r,s) ;
    bet = (delta_new-delta_mid)/delta_old ;
    k = k + 1 ;
    
    Etot(i+1)=info.Etot;
    FORCE{i+1} = reshape(g, [natoms 3]);
    
    if (k==n || bet <=0)
        d = s ;
        k = 0 ;
    else
        d = s + bet*d ;
    end
    
    i = i + 1;
    
    fprintf('\n\n Outer CG iteration no. %d completed. \n',i);
    fprintf('\n No. of function (SCF) evaluation calls so far = %d. \n', ...
    fcalls);
    fprintf('\n Current force norm = %-16.10e.\n',norm(g));
    fprintf('\n Maximum force component = %-16.10e.\n\n', max(abs(g)));
    
end %while

fprintf('\n Total no. of outer NLCG iterations = %d. ', i);
fprintf('\n Total no. of function (SCF) evaluation calls = %d. \n', ...
    fcalls);
fprintf('\n Current force norm = %-16.10e.',norm(g));
fprintf('\n Force vector = ');
display(g);

fprintf('\n Maximum force component = %-16.10e. \n\n', max(abs(g)));
xyz =  reshape(x,[natoms 3]);
mol = set(mol,'xyzlist',xyz);
[mol,~, ~, info] = scf4m(mol, ksopts,int_cord);
FORCE{i+1} =  mol.xyzforce;
Etot(i+1)=info.Etot;



function [g,info] = ksfg(x, mol, ksopts,int_cord)
natoms = sum(mol.natoms);
xyz = reshape(x, [natoms 3]);
mol = set(mol,'xyzlist',xyz);

% Use suitable solver here
% [mol,H, X, info] = scf(mol, ksopts);
[mol,~, ~, info] = scf4m(mol, ksopts,int_cord);

%KSMolViewer(mol);

%f = info.Etotvec(end);
xyzforces = mol.xyzforce;
g = -xyzforces(:);

% fprintf('norm(g) = %11.3e\n', norm(g));
%xyz = mol.xyzlist
%pause
end

end %function
