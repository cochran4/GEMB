function P = CalculatePower

% Effect sizes
delta = 2*10^(-4);

P(:,1) = SubCalculatePower(1000,delta,delta);
P(:,2) = SubCalculatePower(1000,delta/2,3*delta/2);
P(:,3) = SubCalculatePower(10000,delta,delta);
P(:,4) = SubCalculatePower(10000,delta/2,3*delta/2);
writetable(array2table(P'),'Power_Raw.xlsx');
function P = SubCalculatePower(k,effect1,effect2)

% Joint density function 
nu    = 10;
joint = @(x,y) exp( log( ncfpdf(finv(1-x,nu-1,k-nu),nu-1,k-nu,effect1*k) ) - ...
                    log( fpdf(finv(1-x,nu-1,k-nu),nu-1,k-nu) ) + ...
                    log( ncfpdf(finv(1-y,nu-1,k-nu),nu-1,k-nu,effect2*k) ) - ...
                    log( fpdf(finv(1-y,nu-1,k-nu),nu-1,k-nu) ) );
% joint = @(x,y) exp( (2*norminv(1-x/2)-effect1*sqrt(n))*effect1*sqrt(n) + ...
%                     (2*norminv(1-y/2)-effect2*sqrt(n))*effect2*sqrt(n) );

% Tolerance
tol   = 10^(-6);
err   = 1;
m     = 10^6;
pcrit = 0.0253;
S     = zeros(3,1);
P     = zeros(3,1);
it    = 1;
minit = 10; 

% Estimate area using monte carlo
while it < minit && err > tol
    disp(it)
    % Generate random numbers 
    u = rand(m,1);
    v = rand(m,1);
    
    
    % Singleton
    rg1  = rand(m,1) < 1/(2-pcrit); 
    xy   =   rg1 .*[u*pcrit,v] + ...
          (1-rg1).*[u*(1-pcrit)+pcrit,pcrit*v];
    S(1) = 0.05*mean( joint( xy(:,1), xy(:,2) ) );
    
    % Gene set
    xy       = sqrt(0.1)*[u,v];
    ix       = xy(:,1) + xy(:,2) >= sqrt(0.1);                
    xy(ix,:) = sqrt(0.1) - xy(ix,[2 1]);  
    S(2)     = 0.05*mean( joint( xy(:,1), xy(:,2) ) );
    
    % Weighted gene set
    xy       = sqrt(0.1/3)*[3*u,v];
    ix       =  xy(:,1) + 3*xy(:,2) >= 3*sqrt(0.1/3);
    xy(ix,:) = [3*sqrt(0.1/3)-3*xy(ix,2) sqrt(0.1/3)-xy(ix,1)/3];  
    S(3)     = 0.05*mean( joint( xy(:,1), xy(:,2) ) );
    
    % Calculate error, update P value
    Pnew = ((it-1)/it)*P + S/it;
    err = max( abs( Pnew - P ) );
    P   = Pnew;
    it  = it + 1;
    
end

