% Find an equilibrium

% Last period's growth rate is an additional a state variable
% continuous distribution of exit fees (as Exponential)

% Define transition matrix based on parameter pstay
p = (1-pstay)/3;
g_trans = [(pstay + p)  p           p;...
           p            (pstay+p)   p;...
           p            p (         pstay+p)];

%Compute one-period profits
pi=zeros(nmax, zmax);
for i=1:nmax
    for j=1:zmax
        pi(i,j)= profit(i,j,fc,comp);
    end
end


%We hardcode 3 different growth rates, -1,0,1 corresponding to 1..3
%vc, ve have dimensions nmax*zmax*gmax (gmax=5)

vc(:,:,1)=pi/(1-beta); %initial "guess" of the continuation value
vc(:,:,2)=pi/(1-beta);
vc(:,:,3)=pi/(1-beta);
ve=vc; %--- ve(n,z) in the code is equal to VE(n-1,z) in the paper

% ve is entrant continuation value - value of being in conditioning on actually entering
% it is _not_ equal to vc because of conditioning on different information

ve2=ve-.0011;
vc2=vc-.0011;
% we use these to check for convergence


% main loop
time=cputime;
nit = 0;
while (max(max(max(abs(vc-vc2))))>.001) & (nit <5000)
    nit = nit+1;
    Checker=vc-vc2;
    diff=[max(max(max(Checker))) min(min(min(Checker))) nit/500]

    save intermediate
    vc=vc2;
    ve=ve2;
    
    % calculate b_x(x, n-1)---call it bi[ncumbent]
    % bi(a,b,c,d) is the probability that (a-1) firms exit given that
    % there are b firms (including self) and z,g=c,d
    bi=zeros(nmax, nmax, zmax, gmax);
    for i=0:nmax-1
        for j=1:nmax
            for z=1:zmax
                for g=1:gmax
                    if i<=j-1
                        F=1-exp(-sigma*vc(j,z,g));
                        % F is the probability that a given incumbent stays
                        bi(i+1, j,z,g)=nchoosek(j-1,i)*F^(j-1-i)*(1-F)^(i);
                    end
                end
            end
        end
    end
    
    % calculate b_x(x, n)---call it be[ntrant]
    % be(a,b,c) is the probability that (a-1) firms exit given that
    % there are (b-1) firms and z=c
    be=zeros(nmax, nmax, zmax,gmax);
    for i=0:nmax-1
        for j=1:nmax-1
            for z=1:zmax
                for g=1:gmax
                    if i<=j
                        F=1-exp(-sigma*vc(j,z,g));
                        be(i+1, j+1,z,g)=nchoosek(j,i)*F^(j-i)*(1-F)^(i);
                    end
                end
            end
        end
    end
    
    be(1,1,:,:)=ones(zmax,gmax); %if there are 0 guys, 0 of them exit with probability 1
    
    
    % pbentr(a,b,c) is the probability (from the point of view of an incumbent)
    % that there are (a-1) entrants given that there are (b) incumbents
    % and z=c
    pbentr=zeros(J+1, nmax, zmax,gmax);
    for n=1:nmax-1
        for z=1:zmax
            for g=1:gmax
                rr=max(0,(beta*ve(n+1,z,g)-k0));
                F=1-exp(-a*rr)*(1+a*rr);
                %F is the probability that a given potential entrant enters
                for j=0:J
                    pbentr(j+1,n,z,g)=nchoosek(J,j)*F^(j)*(1-F)^(J-j);
                end
            end
        end
    end
    
    pbentr(1,nmax,:,:)=ones(zmax,gmax); %if there are n incumbents, nobody can enter
    
    
    % pbentrC(a,b,c) is the probability (from the point of view of a potential entrant)
    % that there are (a-1) other entrants given that there are (b-1) incumbents
    % and z=c
    pbentrC=zeros(J, nmax,zmax,gmax);
    for n=0:nmax-1
        for z=1:zmax
            for g=1:gmax
                rr=max(0,(beta*ve(n+1,z,g)-k0));
                F=1-exp(-a*rr)*(1+a*rr);
                for j=0:(J-1)
                    pbentrC(j+1,n+1,z,g)=nchoosek(J-1, j)*F^j*(1-F)^(J-1-j);
                end
            end
        end
    end       
    
    
    % compute the new value of vc2 based on vc, bi(vc, ve), pbentr (vc, ve)
    % for each new#firms, new z, continuation value is
    % pi(new#firms, new z) + {beta*vc if stay| expected exit fee if exit}
    vc2=zeros(nmax, zmax,gmax);
    for n=1:nmax
        for z=1:zmax
            for g=1:gmax
                for dx=0:(n-1)
                    for de=0:J
                        transitionrow = g_trans(g,:);
                        dzl=-1;
                        if z == zmin
                            dzl = 0;
                            transitionrow = transitionrow /sum(transitionrow(2:3));
                        end;
                        dzh=1;
                        if z == zmax
                            dzh = 0;
                            transitionrow = transitionrow /sum(transitionrow(1:2));
                        end;
                        for dz=dzl:dzh
                            nnew=min(n+de-dx, nmax);
                            znew=z+dz;
                            gnew=dz+2;
                            pdz=transitionrow(dz+2);
                            alpha=1-exp(-sigma*(vc(nnew,znew,gnew))); %prob. that the firm stays next period
                            vc2(n,z,g)=vc2(n,z,g) + ...
                                (pi(nnew,znew) + ...
                                alpha*beta*vc(nnew, znew, gnew) + ...
                                (1-alpha)*beta*(vc(nnew,znew,gnew)+1/sigma)) * ...
                                pbentr(de+1, n, z, g)*bi(dx+1, n, z, g)*pdz;
                        end
                    end
                end
            end
        end
    end
    
    ve2=zeros(nmax, zmax,gmax);
    for n=0:nmax-1
        for z=1:zmax
            for g=1:gmax
                for dx=0:n
                    for de=1:J
                        transitionrow = g_trans(g,:);
                        dzl=-1;
                        if z == zmin
                            dzl = 0;
                            transitionrow = transitionrow /sum(transitionrow(2:3));
                        end;
                        dzh=1;
                        if z == zmax
                            dzh = 0;
                            transitionrow = transitionrow /sum(transitionrow(1:2));
                        end;
                        for dz=dzl:dzh
                            nnew=min(n+de-dx, nmax);
                            znew=z+dz;
                            gew=dz+2;
                            pdz=transitionrow(dz+2);
                            alpha=1-exp(-sigma*(vc(nnew,znew,gnew)));
                            ve2(n+1,z,g)=ve2(n+1,z,g) + ...
                                (pi(nnew,znew) + ...
                                alpha*beta*vc(nnew, znew, gnew) + ...
                                (1-alpha)*beta*(vc(nnew,znew,gnew)+1/sigma)) * ...    
                                pbentrC(de, n+1, z, g)*bi(dx+1, n+1, z, g)*pdz;
                        end                                                        
                    end
                end
            end
        end
    end
    
    % make sure the algorithm did not hit +/- infinity somewhere
    if max(isnan(vc2(:)))>0 | max(isnan(ve2(:)))>0
        error('NaN encountered in value function')
    end
end



'Done'
nit
save equilibrium
