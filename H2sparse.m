%% H2 optimal sparsity promoting code
 % This function takes as input the system parameters, the control
 % specifications, and a set of sparsity tuning gamma values. For each
 % value of gamma, it returns the sparse feedback matrix F, and the
 % corresponding number of non-zero elements nnz and performance J.

%% ADMM main function

function solpath = H2sparse(A,B1,B2,Q,R,gam)

% Initialization (Set parameters for the optimization)
 options = struct('method','wl1','reweightedIter',20, ...
 'gamval',gam,'rho',100,'maxiter',100,'blksize',[1 1]);

% Data preprocessing
rho     = options.rho;
blksize = options.blksize;

% set the number of reweighted scheme for 
% the weighted l1 or the weighted sum of Frobenius norms
if strcmp(options.method, 'wl1') || strcmp(options.method, 'blkwl1')
    reweighted_Max_Iter = options.reweightedIter;
else
    reweighted_Max_Iter = 1;
end

% size of the input matrix
[n,m] = size(B2);

% preallocate memory for output variables
solpath.F    = zeros(m,n);
solpath.nnz  = 0;
solpath.J    = 0;
solpath.gam  = 0;

% use the centralized gain as the initial condition for ADMM
F = lqr(A,B2,Q,R);
G = F;
Lambda = zeros(m,n);

% weight matrix for weighted l1
if strcmp(options.method, 'wl1')
    W = ones(size(F));
end

% weight matrix for sum of Frobenius norm
% sub_mat_size is numbers of rows and columns of partitioned block submatrices
if strcmp(options.method, 'blkwl1')
    sub_mat_size = size(F)./blksize;
    W = ones(sub_mat_size);
end

% absolute and relative tolerances for the stopping criterion of ADMM
eps_abs = 1.e-4;
eps_rel = 1.e-2;

% stopping criterion tolerance for the Anderson-Moore method and Newton's
% method
tolAM = 1.e-2; 

ADMM_Max_Iter = options.maxiter;

% control of display messages
quiet = 1;

% Solve the sparsity-promoting optimal control problem for each value of
% gamma

% for k = 1:length(gamval)
    
    % gam = gamval(k);     
    
    for reweightedstep = 1 : reweighted_Max_Iter

        % Solve the minimization problem using ADMM
        for ADMMstep = 1 : ADMM_Max_Iter

            % ========================================================
            % F-minimization step using Anderson-Moore method 
            % ========================================================
            U = G - Lambda/rho;
            F = Fmin(A,B1,B2,Q,R,U,rho,F,tolAM);
%             F = Fmin_L(A,B1,B2,Q,R,U,rho,F,tolAM);
            % ========================================================
                                                
            % ========================================================
            % G-minimization step 
            % ========================================================
            V = F + Lambda/rho;
                    % shrinkage operator
                    Gnew = shrinkage(V,W,gam,rho);
            if ~isreal(Gnew)
                error('The solution to G-minimization step is not real!')
            end
                        
            % dual residual
            resG = norm(G - Gnew,'fro');
               G = Gnew;
            
            % primal residual
            resFG = norm(F - G,'fro');                        
            % ===========================================================        

            % ==================== update the dual variable ===============
            Lambda = Lambda + rho * ( F - G );
            % =============================================================

            % stoppin criterion for ADMM
            % evaluate the primal epsilon and the dual epsilon
            eps_pri  = sqrt(n*m) * eps_abs + eps_rel * max(norm(F,'fro'), norm(G,'fro'));
            eps_dual = sqrt(n*m) * eps_abs + eps_rel * norm(Lambda,'fro');

            if  (resFG < eps_pri)  &&  (rho*resG < eps_dual)
                break;
            end    
            
            if ~quiet
                disp([num2str(ADMMstep),'   ',num2str(gam,'%6.1E'),'   ',num2str(resFG,'%6.1E'),'    ',num2str(resG,'%6.1E')])
            end
            
        end
        
        if (ADMMstep == ADMM_Max_Iter) && (~quiet)
            disp('Maximum number of ADMM steps reached!')
            disp(['The primal residual is ', num2str(resFG,'%10.4E')]);
            disp(['The dual residual is ', num2str(rho*resG,'%10.4E')]);
        end        
        
        if max( real( eig( A - B2*G ) ) ) < 0
            F = G;
        else
            if ~quiet
                disp(['Gamma value ',num2str(gam,'%6.1E'),' may be too large.'])
            end
        end
        
        % update the weight matrix W for weighted l1 norm
        eps = 1.e-3;
        Wnew= 1./(abs(F) + eps);           
        if norm(Wnew - W)/norm(W) < 1.e-2
            if ~quiet
                disp(['Re-weighted scheme converges in ', num2str(reweightedstep),' steps.'])
            end
            break;
        else
            W = Wnew;
        end
        
    end   %Re-weighted loop ends here
 % record the feedback gain F, H2 norm, and the number of nonzero
        % entries in F
        solpath.gam     = gam;
        solpath.F       = F;        
        solpath.J       = trace(B1' * lyap((A - B2*F)', Q + F'*R*F) * B1);
        solpath.nnz = nnz(F);
        fprintf('%4d\t\t%5.4f\t\t', solpath.nnz, solpath.J);
    % Computations for a fixed gamma value end here.
% end             



