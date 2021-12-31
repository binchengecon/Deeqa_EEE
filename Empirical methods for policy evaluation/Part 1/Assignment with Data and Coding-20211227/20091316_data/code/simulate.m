%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SIMULATE.M: generate simulated datasets using POB algorithm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRELIMINARIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path(path,'pakes_ostrovsky_berry')
path(path,'..\external')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET DEFAULT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D.pstay = 0;                  % probability change in z this period = change in z last period
                              % (used to control degree of serial correlation; pstay=0 is random
                              % walk)
D.niter = 100000;             % number of periods to simulate
D.beta=.9;                    % discount factor
D.a=.3;                       % parameter in the entry fee distribution
D.k0=1/D.a;                   % minimum entry fee
D.sigma=.75;                  % parameter in the exit fee distribution, default is 0.75
D.J=1;                        % number of potential entrants
D.zmax=45;                    % maximum z
D.zmin=1;                     % minimum z
D.nmax=1;                     % maximum number of firms in the industry
D.gmax=3;                     % number of different possible growth rates
D.fc=3;                       % fixed cost of being alive each period
D.comp=0.2;                     % differentiation factor; comp should be between 0 and 1 and indicates
                              % the share of the market that is contested; a firm is always a
                              % monopolist in a (1-comp) share of its market
D.name='';                    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE OPTIONS FOR INDIVIDUAL RUNS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OPT = [];

% 1 firm, random walk
new = D;
new.name = '1firm_rw';
OPT = [OPT;new];

% 2 firm, random walk
new = D;
new.name = '2firm_rw';
new.nmax = 2;
OPT = [OPT;new];

% 1 firm, serial correlation
new = D;
new.name = '1firm_sc';
new.pstay = 0.8;
OPT = [OPT;new];

% 2 firm, serial correlation
new = D;
new.name = '2firm_sc';
new.nmax = 2;
new.pstay = 0.8;
OPT = [OPT;new];

save ..\output\simulations\options.mat OPT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN SIMULATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:length(OPT)
    
    load ..\output\simulations\options.mat OPT;
    struct2var(OPT(ii));
    disp(['Starting: ' name])
    pause(1)
    run d0212;
    run g0212;
    eval(['save ..\output\simulations\' name '.txt ts -ASCII']);
    nn = ts(:,2);
    sum(nn==0)
    sum(nn==1)
    sum(nn==2)
    
    disp(['Completed: ' name])

    pause(1)
    clear all;

end

delete equilibrium.mat
delete intermediate.mat


