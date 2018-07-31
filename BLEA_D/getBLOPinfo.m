% @Author: hxy
% @Date:   2017-10-16 18:09:13
% @Last Modified by:   hxy
% @Last Modified time: 2018-03-08 11:56:29
function BI = getBLOPinfo(benchmark,fno,dim)
% BI.UmaxFEs = 2000;
% BI.LmaxFEs = 500;
% BI.UmaxImprFEs = 400;
% BI.LmaxImprFEs = 100;

if nargin > 1
	fno = fno(:);
	if length(fno)>1
		BI = [];
		for i = fno'
	        if nargin == 3
	            BI = [BI;getBLOPinfo(benchmark,i,dim)];
	        else
	            BI = [BI;getBLOPinfo(benchmark,i)];
	        end
		end
		return;
	end
end
if strcmp(benchmark,'TP')
	fn = sprintf('tp%d',fno);
	assert(ismember(fno,1:10));
	ulStoppingCriteria = 1e-4;
	llStoppingCriteria = 1e-5;

	%% check where the constraints are from
	switch fno
		case {2,3,4,5,6,7,8}
			% this happens when the two levels have the same constraints
			% thus the duplications should be removed
			isLowerLevelConstraintsIncludedInUpperLevel = true;
		otherwise
			isLowerLevelConstraintsIncludedInUpperLevel = false;
	end
	%% population size
	if fno == 10
	    ulPopSize=200;
	    llPopSize=200;
	else
	    ulPopSize=50;
	    llPopSize=50;
	end
	% dimensionality
	switch fno
		case {1,2,3,5,6,7,8}
			ulDim = 2; llDim = 2;
		case 4
			ulDim = 2; llDim = 3;
		case 9
			ulDim = 5; llDim = 5;
		case 10
			ulDim = 10; llDim = 10;
	end
	%% boundaries and number of constraints
	switch fno
		case 1
			ulDimMin = [-30 -30];
			ulDimMax = [30 15];
			llDimMin = 0*ones(1,llDim);
			llDimMax = 10*ones(1,llDim);
			ul_ieqcon_num = 2;
			ll_ieqcon_num = 0;
		case 2
			ulDimMin = zeros(1,ulDim);
			ulDimMax = 50*ones(1,ulDim);
			llDimMin = -10*ones(1,llDim);
			llDimMax = 20*ones(1,llDim);
			ul_ieqcon_num = 3;
			ll_ieqcon_num = 2;
		case 3
			ulDimMin = zeros(1,ulDim);
			ulDimMax = 10*ones(1,ulDim);
			llDimMin = zeros(1,llDim);
			llDimMax = 10*ones(1,llDim);
			ul_ieqcon_num = 3;
			ll_ieqcon_num = 2;
		case 4
			ulDimMin = zeros(1,ulDim);
			ulDimMax = 1*ones(1,ulDim);
			llDimMin = zeros(1,llDim);
			llDimMax = 1*ones(1,llDim);
			ul_ieqcon_num = 3;
			ll_ieqcon_num = 3;
		case 5
			ulDimMin = zeros(1,ulDim);
			ulDimMax = 10*ones(1,ulDim);
			llDimMin = zeros(1,llDim);
			llDimMax = 10*ones(1,llDim);
			ul_ieqcon_num = 2;
			ll_ieqcon_num = 2;
		case 6
			ulDimMin = zeros(1,ulDim);
			ulDimMax = 2*ones(1,ulDim);
			llDimMin = zeros(1,llDim);
			llDimMax = 2*ones(1,llDim);
			ul_ieqcon_num = 4;
			ll_ieqcon_num = 4;
		case 7
			ulDimMin = zeros(1,ulDim);
			ulDimMax = 10*ones(1,ulDim);
			llDimMin = zeros(1,llDim);
			llDimMax = [1 10];
			ul_ieqcon_num = 4;
			ll_ieqcon_num = 2;
		case 8
			ulDimMin = zeros(1,ulDim);
			ulDimMax = 50*ones(1,ulDim);
			llDimMin = -10*ones(1,llDim);
			llDimMax = 20*ones(1,llDim);
			ul_ieqcon_num = 3;
			ll_ieqcon_num = 2;
		case 9
			ulDimMin = -1*ones(1,ulDim);
			ulDimMax = 1*ones(1,ulDim);
			llDimMin = -pi*ones(1,llDim);
			llDimMax = pi*ones(1,llDim);
			ul_ieqcon_num = 0;
			ll_ieqcon_num = 0;
		case 10
			ulDimMin = -1*ones(1,ulDim);
			ulDimMax = 1*ones(1,ulDim);
			llDimMin = -pi*ones(1,llDim);
			llDimMax = pi*ones(1,llDim);
			ul_ieqcon_num = 0;
			ll_ieqcon_num = 0;
	end
	ul_eqcon_num = 0;
	ll_eqcon_num = 0;

	%% optimal values
	% from BLEAQ2
	% ulBestKnownFunctionValue = [225, 0, -18.6787, -29.2, -3.6, -1.2091, -1.96, 0, 0, 0];
    % llBestKnownFunctionValue = [100, 100, -1.0156, 3.2, -2.0, 7.6145, 1.96, 100, 1.0, 1.0];
    % my results
%     does TP3 has better results?
	ulBestKnownFunctionValue = [225, 0, -18.6787, -29.2, -3.6, -1.20987, -1.96146, 0, 0, 0];
    llBestKnownFunctionValue = [100, 100, -1.0156, 3.2, -2.0, 7.61728, 1.96146, 100, 1.0, 1.0];

    ulOpt = -ulBestKnownFunctionValue(fno);
    llOpt = -llBestKnownFunctionValue(fno);
end
if strcmp(benchmark,'SMD')
	fn = sprintf('smd%d',fno);
	assert(ismember(fno,1:12));
	ulStoppingCriteria = 1e-4;
	llStoppingCriteria = 1e-5;
	eps = 0.00001;
    if nargin == 3
    	if dim == 5
    		ulDim = 2;
    		llDim = 3;
		elseif dim == 10 || dim == 20
            ulDim = dim / 2;
            llDim = dim / 2;
        else
            error('no supported dimensionality');
        end
        % if dim >= 5 && mod(dim,5) == 0
        %     ulDim = dim * 2 / 5;
        %     llDim = dim - ulDim;
        % else
        %     error('no supported dimensionality');
        % end
    else
        ulDim=2; llDim=3;
    end
    %% check the constraints
	isLowerLevelConstraintsIncludedInUpperLevel = false;
    %% population size
    ulPopSize=20;
    llPopSize=30;
    %% decision variable decomposition
    switch fno
    	case 6
			r = floor(ulDim/2);
			p = ulDim - r;
			q = floor((llDim - r)/2 - eps);
			s = ceil((llDim - r)/2 + eps);
		case {1,2,3,4,5,7,8,9,10,11,12}
			r = floor(ulDim/2); 
			p = ulDim - r; 
			q = llDim - r;
    end
    %% optimum function values
	[~,~,ulOpt,llOpt]=getOptimalSolutionSMD(ulDim,llDim,sprintf('smd%d',fno));
    %% boundaries
	switch fno
		case 1 %ok
			ulDimMin = -5*ones(1,ulDim);                    
			ulDimMax = 10*ones(1,ulDim);                    
			llDimMin = [-5*ones(1,q) -pi/2*ones(1,r)+eps];  
			llDimMax = [10*ones(1,q)  pi/2*ones(1,r)-eps];  
		case 2 %ok
			ulDimMin = -5*ones(1,ulDim); 
			ulDimMax = [10*ones(1,p) 1*ones(1,r)];
			llDimMin = [-5*ones(1,q) eps*ones(1,r)];
			llDimMax = [10*ones(1,q) exp(1)*ones(1,r)];
		case 3 %ok
			ulDimMin = -5*ones(1,ulDim);
			ulDimMax = 10*ones(1,ulDim);
			llDimMin = [-5*ones(1,q) -pi/2*ones(1,r)+eps];
			llDimMax = [10*ones(1,q)  pi/2*ones(1,r)-eps];
		case 4 %ok
			ulDimMin = [-5*ones(1,p) -1*ones(1,r)];
			ulDimMax = [10*ones(1,p) 1*ones(1,r)];
			llDimMin = [-5*ones(1,q) zeros(1,r)];
			llDimMax = [10*ones(1,q) exp(1)*ones(1,r)];
		case 5 %ok
			ulDimMin = -5*ones(1,ulDim);
			ulDimMax = 10*ones(1,ulDim);
			llDimMin = [-5*ones(1,q) -5*ones(1,r)];
			llDimMax = [10*ones(1,q) 10*ones(1,r)];
		case 6 %ok
			ulDimMin = -5*ones(1,ulDim);
			ulDimMax = 10*ones(1,ulDim);
			llDimMin = -5*ones(1,llDim);
			llDimMax = 10*ones(1,llDim);
		case 7 %ok
			ulDimMin = -5*ones(1,ulDim);
			ulDimMax = [10*ones(1,p) 1*ones(1,r)];
			llDimMin = [-5*ones(1,q) eps*ones(1,r)];
			llDimMax = [10*ones(1,q) exp(1)*ones(1,r)];
		case 8 %ok
			ulDimMin = -5*ones(1,ulDim);
			ulDimMax = 10*ones(1,ulDim);
			llDimMin = [-5*ones(1,q) -5*ones(1,r)];
			llDimMax = [10*ones(1,q) 10*ones(1,r)];
		case 9 %ok
			ulDimMin = -5*ones(1,ulDim);
			ulDimMax = [10*ones(1,p) 1*ones(1,r)];
			llDimMin = [-5*ones(1,q) -1+eps*ones(1,r)];
			llDimMax = [10*ones(1,q) -1+exp(1)*ones(1,r)];
		case 10	%ok		
			ulDimMin = -5*ones(1,ulDim);
			ulDimMax = 10*ones(1,ulDim);
			llDimMin = [-5*ones(1,q) -pi/2*ones(1,r)+eps];
			llDimMax = [10*ones(1,q)  pi/2*ones(1,r)-eps];
		case 11 %ok
			ulDimMin = [-5*ones(1,p) -1*ones(1,r)];
			ulDimMax = [10*ones(1,p) 1*ones(1,r)];
			llDimMin = [-5*ones(1,q) 1/exp(1)*ones(1,r)];
			llDimMax = [10*ones(1,q) exp(1)*ones(1,r)];
		case 12 %???? inconsistent with the paper
			ulDimMin = [-5*ones(1,p) -1*ones(1,r)];
			ulDimMax = [10*ones(1,p) 1*ones(1,r)];
			llDimMin = [-5*ones(1,q) -pi/4*ones(1,r)+eps];
			llDimMax = [10*ones(1,q)  pi/4*ones(1,r)-eps];
	end
	%% number of constraints
	switch fno
		case {1,2,3,4,5,6,7,8}
			ul_ieqcon_num = 0;
			ll_ieqcon_num = 0;
		case 9
			ul_ieqcon_num = 1;
			ll_ieqcon_num = 1;
		case 10
			ul_ieqcon_num = p+r;
			ll_ieqcon_num = q;
		case 11
			ul_ieqcon_num = r;
			ll_ieqcon_num = 1;
		case 12
			ul_ieqcon_num = p+2*r;
			ll_ieqcon_num = q+1;
	end
	ul_eqcon_num = 0;
	ll_eqcon_num = 0;
end
if strcmp(benchmark,'DecisionMaking')
	fn = benchmark;
	ulStoppingCriteria = 1e-4;
	llStoppingCriteria = 1e-5;

	isLowerLevelConstraintsIncludedInUpperLevel = false;

	%% population size
    ulPopSize=50; llPopSize=50;
	% dimensionality
	ulDim = 3; llDim = 3;
	%% boundaries and number of constraints
	ulDimMin = [0 0 0];
	ulDimMax = [250 250 1];
	llDimMin = [0 0 0];
	llDimMax = [70 70 70];
	ul_ieqcon_num = 4;
	ll_ieqcon_num = 5;
	ul_eqcon_num = 0;
	ll_eqcon_num = 0;

	%% optimal values
	ulBestKnownFunctionValue = 0;
    llBestKnownFunctionValue = 0;

    ulOpt = -ulBestKnownFunctionValue(1);
    llOpt = -llBestKnownFunctionValue(1);
end
if strcmp(benchmark,'GoldMining')
	fn = benchmark;
	ulStoppingCriteria = 1e-4;
	llStoppingCriteria = 1e-5;

	isLowerLevelConstraintsIncludedInUpperLevel = false;

	%% population size
    ulPopSize=50; llPopSize=50;
	% dimensionality
	ulDim = 2; llDim = 1;
	%% boundaries and number of constraints
	ulDimMin = [0 0];
	ulDimMax = [100 1];
	llDimMin = 0;
	llDimMax = 100;
	ul_ieqcon_num = 0;
	ll_ieqcon_num = 1;
	ul_eqcon_num = 0;
	ll_eqcon_num = 0;

	%% optimal values
	ulBestKnownFunctionValue = 0;
    llBestKnownFunctionValue = 0;

    ulOpt = -ulBestKnownFunctionValue(1);
    llOpt = -llBestKnownFunctionValue(1);
end
BI.u_ieqcon_num = ul_ieqcon_num;
BI.l_ieqcon_num = ll_ieqcon_num;
BI.u_eqcon_num = ul_eqcon_num;
BI.l_eqcon_num = ll_eqcon_num;

BI.u_dim = ulDim;
BI.l_dim = llDim;
BI.dim = ulDim + llDim;

BI.u_lb = ulDimMin;
BI.u_ub = ulDimMax;
BI.l_lb = llDimMin;
BI.l_ub = llDimMax;
BI.xrange = [ulDimMin llDimMin;ulDimMax llDimMax];

BI.u_ftol = 1e-6;
BI.l_ftol = 1e-6;

% from the original BLEAQ2 source code
BI.ulStoppingCriteria = ulStoppingCriteria;
BI.llStoppingCriteria = llStoppingCriteria;

BI.u_fopt = ulOpt;
BI.l_fopt = llOpt;

BI.fn = fn;

% from BLEAQ2
BI.u_N = ulPopSize;
BI.l_N = llPopSize;
BI.u_maxGen = 2000;
BI.l_maxGen = 2000;

BI.isLowerLevelConstraintsIncludedInUpperLevel = isLowerLevelConstraintsIncludedInUpperLevel;