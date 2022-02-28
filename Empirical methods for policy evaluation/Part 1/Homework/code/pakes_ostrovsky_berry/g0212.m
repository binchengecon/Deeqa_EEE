% generate time series data
% columns: z, n, #exits, #entrants
% initial z, n are chosen randomly

% niter=100000;   % set externally

% this line randomly resets matlab's random number generator
% should be commented for "reproducible" simulations
rand('seed',61382)

tseries=zeros(niter,4);

zz=ceil(rand*zmax);
nn=floor(rand*(nmax+1));
gg=2;

% pbentrR(a,b,c) is the probability (from the point of view of an econometrician)
% that there are (a-1) entrants given that there are (b-1) incumbents
% and z=c


pbentrR=zeros(J+1, nmax+1, zmax, gmax);
for n=0:nmax-1
    for z=1:zmax
        for g=1:gmax
            rr=max(0,(beta*ve(n+1,z,g)-k0));
            F=1-exp(-a*rr)*(1+a*rr);
            for j=0:J
                pbentrR(j+1,n+1,z,g)=nchoosek(J,j)*F^j*(1-F)^(J-j);
            end
        end
    end
end

pbentrR(1,nmax+1,:,:)=ones(zmax,gmax); %if there are n guys, nobody can enter

% bR(a,b,c) is the probability that (a-1) firms exit given that
% there are (b-1) firms and z=c, from the point of view of an econometrician

bR=zeros(nmax+1, nmax+1, zmax, gmax);
for i=0:nmax
    for j=1:nmax
        for z=1:zmax
            for g=1:gmax
                if i<=j
                    F=1-exp(-sigma*vc(j,z,g));
                    bR(i+1, j+1,z,g)=nchoosek(j,i)*F^(j-i)*(1-F)^(i);
                end
            end
        end
    end
end

bR(1,1,:,:)=ones(zmax,gmax); %if there are 0 guys, 0 of them exit with probability 1

cer=cumsum(pbentrR);
cxr=cumsum(bR);

entrydraw=rand;
k=0;
while (cer(k+1,nn+1,zz,gg) < entrydraw) 
    k=k+1;
end;
nentr=k;

exitdraw=rand;
k=0;
while (cxr(k+1,nn+1,zz,gg) < exitdraw) 
    k=k+1;
end;
nexit=k;

%tseries=zeros(niter+100000,4);
tseries=zeros(niter+100000,4+nmax);
tseries(1,1)=zz;
tseries(1,2)=nn;
tseries(1,3)=nexit;
tseries(1,4)=nentr;

for i = 1:nmax
    tseries(1,4+i) = pi(i,zz);
end;


%we will later discard the first 100K observations
for i=2:(niter+100000)
    lastg = gg;
    transitionrow = g_trans(lastg,:);
    transitioncdf = cumsum(transitionrow);
    newgg = 1;
    gdraw = rand;
    while (gdraw > transitioncdf(newgg))
        newgg = newgg + 1;
    end;
    % check for boundary conditions
    if zz == zmin & newgg == 1
        newgg = 2;
    elseif zz == zmax & newgg == 3
        newgg = 2;
    end;
    gg = newgg;
        
%    gg = 3 + dzmin + floor(zdraw * (dzmax - dzmin + 1));
    
    zz = zz + gg - 2;
    nn=min(nn-nexit+nentr, nmax);
    
    entrydraw=rand;
    k=0;
    while (cer(k+1,nn+1,zz,gg) < entrydraw) 
        k=k+1;
    end;
    nentr=k;
    
    exitdraw=rand;
    k=0;
    while (cxr(k+1,nn+1,zz,gg) < exitdraw) 
        k=k+1;
    end;
    nexit=k;
    
    tseries(i,1)=zz;
    tseries(i,2)=nn;
    tseries(i,3)=nexit;
    tseries(i,4)=nentr;
    for j = 1:nmax
        tseries(i,4+j) = pi(j,zz);
    end;
end

ts=tseries(100001:(niter+100000),:);