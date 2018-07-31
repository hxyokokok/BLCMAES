% cd D:\soft\BLEAQ2\bleaq2-master 
% cd D:\soft\NBLEA
% cd D:\soft\BLEAQ
cd 'D:\workspace\Computational Intelligence'

BI = getBLOPinfo('SMD',10,5);
% BI = getBLOPinfo('TP',5);
% BI = getBLOPinfo('GoldMining');
% BI = getBLOPinfo('DecisionMaking');
% GoldMining
% 10+10
% BI.UmaxFEs = 5000;
% BI.UmaxImprFEs = 750;
% BI.LmaxFEs = 500;
% BI.LmaxImprFEs = 50;

% 5+5
% BI.UmaxFEs = 3500;
% BI.UmaxImprFEs = 500;
% BI.LmaxFEs = 350;
% BI.LmaxImprFEs = 35;

% 2+3
BI.UmaxFEs = 2500;
BI.UmaxImprFEs = 350;
BI.LmaxFEs = 250;
BI.LmaxImprFEs = 25;

% 2+1
% BI.UmaxFEs = 1500;
% BI.UmaxImprFEs = 70;
% BI.LmaxFEs = 150;
% BI.LmaxImprFEs = 15;

BI.u_N = 50;
BI.l_N = 50;

% disp(BI);

rng_settings = rng;
% rng(rng_settings);

tic;

if isempty(strfind(pwd,'Intelligence'))
    ins = wapper(BI);
else
    ins = BLCMAESv4(BI);
end



ins.runTime = toc;
fprintf('UF = %.6f, LF = %.6f\n', ins.UF, ins.LF);
fprintf('UFAcc = %g, LFAcc = %g, UFEs = %d, LFEs = %d\n', abs(ins.UF-BI.u_fopt), abs(ins.LF-BI.l_fopt), ins.UFEs, ins.LFEs);

function ins = wapper(BI)
[ins.UF,ins.LF,ins.UX,ins.LX,ins.UFEs,ins.LFEs,ins.record] = ulSearch(BI.fn, ...
BI.u_N, BI.u_maxGen, ...
BI.u_dim, BI.u_lb, BI.u_ub, ...
BI.l_N, BI.l_maxGen, ...
BI.l_dim, BI.l_lb, BI.l_ub, ...
BI.ulStoppingCriteria, BI.llStoppingCriteria);
end


