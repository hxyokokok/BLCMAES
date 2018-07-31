% @Author: hxy
% @Date:   2017-11-15 11:48:12
% @Last Modified by:   hxy
% @Last Modified time: 2018-03-06 20:00:16

% i --- injection correction
% e --- elite preservation
% c --- constrained handling
% s --- advanced stopping criteria

% 2017/11/15 modify some termination criteria

function ins = BLCMAES(BI)
global ulFunctionEvaluations;
global llFunctionEvaluations;
ulFunctionEvaluations = 0;
llFunctionEvaluations = 0;

elite = [];
record = [];
CMA = initCMAES(BI);
maxIter = ceil(BI.UmaxFEs/CMA.lambda);
ImprIter = ceil(BI.UmaxImprFEs/CMA.lambda);
for iter = 1 : maxIter
    % sampling
    for i = 1 : CMA.lambda
        U = CMA.xmean + CMA.sigma * randn(1,BI.dim) .* CMA.D * CMA.B';
        flag = U > BI.xrange(2,:);
        U(flag) = (CMA.xmean(flag) + BI.xrange(2,flag))/2;
        flag = U < BI.xrange(1,:);
        U(flag) = (CMA.xmean(flag) + BI.xrange(1,flag))/2;
        POP(i).UX = U(1:BI.u_dim);
        POP(i).LX = U(BI.u_dim+1:end);
    end
    % lower level search
    for i = 1 : CMA.lambda
        [POP(i).LX,POP(i).LF,POP(i).LC,POP(i).RF] = lowerLevelSearch(POP(i).UX,CMA,BI);
        [POP(i).UF,POP(i).UC] = UL_evaluate(POP(i).UX,POP(i).LX,BI.fn);
        POP(i).UFEs = ulFunctionEvaluations;
        POP(i).LFEs = llFunctionEvaluations;
    end
    % fitness assignment && refinement
    POP = assignUpperFitness(POP,BI);

    % find current best and optimal solution
    RF_idx = find([POP.RF]);
    if isempty(RF_idx) RF_idx = 1:length(POP); end
    [~,bestIdx] = max([POP(RF_idx).fit]);
    bestIdx = RF_idx(bestIdx);
    bestIndv = POP(bestIdx);

    % replace the elite
    if upperLevelComparator(bestIndv,elite,BI) || rand > 0.5
    	bestIndv = refine(bestIndv,CMA,BI);
    	POP(bestIdx) = bestIndv;
    	if upperLevelComparator(bestIndv,elite,BI)
    		elite = bestIndv;
    	end
    else
    	elite = refine(elite,CMA,BI);
    end
    % re-assign fitness
    POP = assignUpperFitness(POP,BI);
    % [elite.UF elite.LF CMA.sigma]
%     [ CMA.xmean CMA.sigma]
    % [elite.UF elite.LF elite.UC elite.LC CMA.sigma ]

    % recording
    elite.UFEs = ulFunctionEvaluations;
    elite.LFEs = llFunctionEvaluations;
    record = [record;elite];

    %% termination check
    reachTheTarget_b = abs(bestIndv.UF-BI.u_fopt) < BI.u_ftol;
    reachTheTarget_e = abs(elite.UF-BI.u_fopt) < BI.u_ftol;
    reachTheMaxFEs = ulFunctionEvaluations >= BI.UmaxFEs;
    reachTheFlat = iter>ImprIter && abs(record(iter).UF-record(iter-ImprIter+1).UF)/(abs(record(1).UF)+abs(record(iter).UF)) < BI.u_ftol; 
    reachTheFlat = reachTheFlat && abs(record(iter).UF-record(iter-ImprIter+1).UF) < BI.u_ftol;
    
    if reachTheTarget_b || reachTheTarget_e || reachTheMaxFEs || reachTheFlat
        if reachTheTarget_b
            elite = bestIndv;
        end
        break;
    end
    CMA = updateCMAES(CMA,POP,BI);
end

ins.UF = elite.UF;
ins.LF = elite.LF;
ins.UX = elite.UX;
ins.LX = elite.LX;
ins.UFEs = ulFunctionEvaluations;
ins.LFEs = llFunctionEvaluations;
ins.record = record;

function cfit_ = combineConstraintWithFitness(fit_,c)
if all(c==0)
	cfit_ = fit_;
else
    cfit_ = -c;
    isFeasible = c == 0;
    if any(isFeasible)
        cfit_(isFeasible) = fit_(isFeasible);
        cfit_(~isFeasible) = min(fit_(isFeasible)) - c(~isFeasible);
    end
end

function POP = assignLowerFitness(POP)
fit_ = combineConstraintWithFitness([POP.LF],[POP.LC]);
for i = 1 : length(POP)
    POP(i).fit = fit_(i);
end

function POP = assignUpperFitness(POP,BI)
CV = [POP.UC];
if ~BI.isLowerLevelConstraintsIncludedInUpperLevel
    CV = CV + [POP.LC];
end
fit_ = combineConstraintWithFitness([POP.UF],CV);
for i = 1 : length(POP)
	POP(i).fit = fit_(i);
	% ????? how to handle RF
	% according to latest experimental results, incorporating RF here seems to be meaningless
end

function Q = refine(P,CMA,BI)
Q = P;
[Q.LX,Q.LF,Q.LC,Q.RF] = lowerLevelSearch(Q.UX,CMA,BI);
if lowerLevelComparator(Q,P)
    Q.RF = max(Q.RF,P.RF);
    [Q.UF,Q.UC] = UL_evaluate(Q.UX,Q.LX,BI.fn);
else
	Q = P;
end

function isNoWorseThan = upperLevelComparator(P,Q,BI)
if isempty(Q)
    isNoWorseThan = true;
else
    tmp = assignUpperFitness([P Q],BI);
    isNoWorseThan = tmp(1).fit >= tmp(2).fit;
end

function isNoWorseThan = lowerLevelComparator(P,Q)
if isempty(Q)
    isNoWorseThan = true;
else
    tmp = assignLowerFitness([P Q]);
    isNoWorseThan = tmp(1).fit >= tmp(2).fit;
end

function [bestLX,bestLF,bestLC,bestRF] = lowerLevelSearch(xu,CMA,BI)
sigma0 = 1;
LCMA.xmean = CMA.xmean(BI.u_dim+1:end);
LCMA.sigma = sigma0;
LCMA.C = CMA.C(BI.u_dim+1:end,BI.u_dim+1:end) * CMA.sigma^2;
LCMA.pc = CMA.pc(BI.u_dim+1:end) * CMA.sigma;
LCMA.ps = zeros(1,BI.l_dim);
lambda = 4+floor(3*log(BI.l_dim));
mu = floor(lambda/2);
weights = log(mu+1/2)-log(1:mu);
weights = weights/sum(weights);
mueff=sum(weights)^2/sum(weights.^2);
cc = (4+mueff/BI.l_dim) / (BI.l_dim+4 + 2*mueff/BI.l_dim);
cs = (mueff+2) / (BI.l_dim+mueff+5);
c1 = 2 / ((BI.l_dim+1.3)^2+mueff);
cmu = min(1-c1, 2 * (mueff-2+1/mueff) / ((BI.l_dim+2)^2+mueff));
damps = 1 + 2*max(0, sqrt((mueff-1)/(BI.l_dim+1))-1) + cs;
chiN=BI.l_dim^0.5*(1-1/(4*BI.l_dim)+1/(21*BI.l_dim^2));
[LCMA.B,LCMA.D] = eig(LCMA.C);
LCMA.D = sqrt(diag(LCMA.D))';
LCMA.invsqrtC = LCMA.B * diag(LCMA.D.^-1) * LCMA.B';
cy = sqrt(BI.l_dim)+2*BI.l_dim/(BI.l_dim+2);
bestIndv = [];
bestRF = false;

maxIter = ceil(BI.LmaxFEs / lambda);
ImprIter = ceil(BI.LmaxImprFEs / lambda);
record = zeros(1,maxIter);
for iter = 1 : maxIter
	for i = 1 : lambda
	    Q(i).LX = LCMA.xmean + LCMA.sigma * randn(1,BI.l_dim) .* LCMA.D * LCMA.B';
	    flag = Q(i).LX > BI.l_ub;
	    Q(i).LX(flag) = (LCMA.xmean(flag) + BI.l_ub(flag))/2;
	    flag = Q(i).LX < BI.l_lb;
	    Q(i).LX(flag) = (LCMA.xmean(flag) + BI.l_lb(flag))/2;
	    [Q(i).LF,Q(i).LC] = LL_evaluate(xu,Q(i).LX,BI.fn);
	end
	Q = assignLowerFitness(Q);
	[~, Xindex] = sort(-[Q.fit]);
	xold = LCMA.xmean;
	X = cell2mat(arrayfun(@(q)q.LX,Q','UniformOutput',false));
	Y = bsxfun(@minus,X(Xindex(1:mu),:),xold) / LCMA.sigma;
	Y = bsxfun(@times,Y,min(1,cy./sqrt(sum((Y*LCMA.invsqrtC').^2,2))));
	delta_xmean = weights * Y;
	LCMA.xmean = LCMA.xmean + delta_xmean * LCMA.sigma;
	C_mu = Y' * diag(weights) * Y;
	LCMA.ps = (1-cs)*LCMA.ps + sqrt(cs*(2-cs)*mueff) * delta_xmean * LCMA.invsqrtC;
	LCMA.pc = (1-cc)*LCMA.pc + sqrt(cc*(2-cc)*mueff) * delta_xmean;
	LCMA.C = (1-c1-cmu) * LCMA.C + c1 * (LCMA.pc'*LCMA.pc) + cmu * C_mu; 
	delta_sigma = (cs/damps)*(norm(LCMA.ps)/chiN - 1);
	LCMA.sigma = LCMA.sigma * exp(min(CMA.delta_sigma_max,delta_sigma));
	LCMA.C = triu(LCMA.C) + triu(LCMA.C,1)';
	[LCMA.B,LCMA.D] = eig(LCMA.C);
	LCMA.D = sqrt(diag(LCMA.D))';
    LCMA = repairCMA(LCMA);
	LCMA.invsqrtC = LCMA.B * diag(LCMA.D.^-1) * LCMA.B';

	% elite-preservation right????
	if lowerLevelComparator(Q(Xindex(1)),bestIndv)
		bestIndv = Q(Xindex(1));
	end
    record(iter) = bestIndv.LF;

    if (iter>ImprIter && abs(record(iter)-record(iter-ImprIter+1))/(abs(record(1))+abs(record(iter))) < 1e-4) ...
		|| (iter > ImprIter && abs(record(iter) - record(iter-ImprIter+1)) < 10*BI.l_ftol) ...
        || LCMA.sigma / sigma0 < 1e-2 ...
        || LCMA.sigma / sigma0 > 1e2
    	bestRF = true;
        break; 
    end
end
bestLX = bestIndv.LX;
bestLF = bestIndv.LF;
bestLC = bestIndv.LC;

function CMA = initCMAES(BI)
CMA.lambda = 4+floor(3*log(BI.dim));
% CMA.sigma = min(BI.xrange(2,:)-BI.xrange(1,:))/2;
CMA.sigma = 0.3*median(BI.xrange(2,:)-BI.xrange(1,:));
CMA.mu = floor(CMA.lambda/2);
CMA.weights = log(CMA.mu+1/2)-log(1:CMA.mu);
CMA.weights = CMA.weights/sum(CMA.weights);
CMA.mueff=sum(CMA.weights)^2/sum(CMA.weights.^2);
CMA.cc = (4+CMA.mueff/BI.dim) / (BI.dim+4 + 2*CMA.mueff/BI.dim);
CMA.cs = (CMA.mueff+2) / (BI.dim+CMA.mueff+5);
CMA.c1 = 2 / ((BI.dim+1.3)^2+CMA.mueff);
CMA.cmu = min(1-CMA.c1, 2 * (CMA.mueff-2+1/CMA.mueff) / ((BI.dim+2)^2+CMA.mueff));
CMA.damps = 1 + 2*max(0, sqrt((CMA.mueff-1)/(BI.dim+1))-1) + CMA.cs;
CMA.chiN=BI.dim^0.5*(1-1/(4*BI.dim)+1/(21*BI.dim^2));
CMA.pc = zeros(1,BI.dim);
CMA.ps = zeros(1,BI.dim);
CMA.B = eye(BI.dim);
CMA.D = ones(1,BI.dim);
CMA.C = CMA.B * diag(CMA.D.^2) * CMA.B';
CMA.invsqrtC = CMA.B * diag(CMA.D.^-1) * CMA.B';
CMA.xmean = (BI.xrange(2,:)-BI.xrange(1,:)).*rand(1,BI.dim)+BI.xrange(1,:);
CMA.cy = sqrt(BI.dim)+2*BI.dim/(BI.dim+2);
CMA.delta_sigma_max = 1;

function CMA = updateCMAES(CMA,POP,BI)
[~, Xindex] = sort(-[POP.fit]);
xold = CMA.xmean;
X = cell2mat(arrayfun(@(p)[p.UX p.LX],POP(Xindex(1:CMA.mu))','UniformOutput',false));
Y = bsxfun(@minus,X,xold)/CMA.sigma;
Y = bsxfun(@times,Y,min(1,CMA.cy./sqrt(sum((Y*CMA.invsqrtC').^2,2))));
delta_xmean = CMA.weights * Y;
CMA.xmean = CMA.xmean + delta_xmean * CMA.sigma;
C_mu = Y' * diag(CMA.weights) * Y;
CMA.ps = (1-CMA.cs)*CMA.ps + sqrt(CMA.cs*(2-CMA.cs)*CMA.mueff) * delta_xmean * CMA.invsqrtC;
CMA.pc = (1-CMA.cc)*CMA.pc + sqrt(CMA.cc*(2-CMA.cc)*CMA.mueff) * delta_xmean;
CMA.C = (1-CMA.c1-CMA.cmu) * CMA.C + CMA.c1 * (CMA.pc'*CMA.pc) + CMA.cmu * C_mu; 
delta_sigma = (CMA.cs/CMA.damps)*(norm(CMA.ps)/CMA.chiN - 1);
CMA.sigma = CMA.sigma * exp(min(delta_sigma,CMA.delta_sigma_max));
CMA.C = triu(CMA.C) + triu(CMA.C,1)';
[CMA.B,CMA.D] = eig(CMA.C);
CMA.D = sqrt(diag(CMA.D))';
CMA = repairCMA(CMA);
CMA.invsqrtC = CMA.B * diag(CMA.D.^-1) * CMA.B';

function [F,C] = LL_evaluate(xu,xl,fn)
[F,~,C] = llTestProblem(xl,fn,xu);
C = sum(max(0,C));

function [F,C] = UL_evaluate(UPop,LPOP,fn)
[F,~,C] = ulTestProblem(UPop, LPOP, fn);
C = sum(max(0,C));

function Model = repairCMA(Model)
dim = length(Model.D);
% limit condition of C to 1e14
if any(Model.D<=0)
    Model.D(Model.D<0) = 0;
    tmp = max(Model.D)/1e7;
    Model.C = Model.C + tmp * eye(dim);
    Model.D = Model.D + tmp * ones(1,dim);
end
if max(Model.D) > 1e7 * min(Model.D)
    tmp = max(Model.D)/1e7 - min(Model.D);
    Model.C = Model.C + tmp * eye(dim);
    Model.D = Model.D + tmp * ones(1,dim);
end
% rescale sigma
if Model.sigma > 1e7 * max(Model.D)
    fac = Model.sigma / max(Model.D);
    Model.sigma = Model.sigma / fac;
    Model.D = Model.D * fac;
    Model.pc = Model.pc * fac;
    Model.C = Model.C * fac^2;
end