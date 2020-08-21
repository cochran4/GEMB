function SerotoninSensitivity

% Initialization
V5ht0   = velocitySerotonin;
y0_5ht  = [0.14; 0.86; 20.6; 2.26; 0.5; 21.45; 0.000768; 5.26; 144.9];

% Indices for modulated parameters (pathway,parameter #)
ix  = [1,4; 2,3; 2,6; 3,2; 4,2; 5,2; 6,2; 14,1; 7,2; 8,2; 9,1];

% Initialize parameters
np    = size(ix,1); 
nruns = 1000;
T     = [0,100]; % 200 simulated hours


% Build design matrix
rng(0)
X     = 1 + lhsnorm(zeros(np,1),eye(np)*0.05^2,nruns);

% Run through design matrix
for i=1:nruns
    % Display counter
    if mod(i-1,100)==0, disp(i), end
    
    % Rescale parameters
    V5ht    = V5ht0;
    for j=1:np
        V5ht(ix(j,1)).param(ix(j,2)) = V5ht(ix(j,1)).param(ix(j,2))*X(i,j);
    end
    
    % Run simulation
    hdl    = @(t,y)Monoamine(t,y,V5ht);
    [t,Y]=ode23t(hdl,T,y0_5ht);
    
    % Collect output of extracellular dopamine
    e5ht(i,1) = Y(end,7);
end

% Calculate PRCC
rho           = partialcorr( [X,e5ht],'Type','Spearman');
serotonin_prcc = rho(1:end-1,end);

% Save weights
save('serotonin_prcc','serotonin_prcc')


% Turn into weights
weights     = zeros(11,1);
weights(1)  = abs(serotonin_prcc(1))/3;               % TPH1
weights(2)  = abs(serotonin_prcc(1))/3;               % TPH2
weights(3)  = mean( abs(serotonin_prcc([2,3])) );     % QDPR
weights(4)  = abs(serotonin_prcc(4));                 % SLC7A5
weights(4)  = abs(serotonin_prcc(5));                 % DDC
weights(6)  = abs(serotonin_prcc(6))/2;               % SLC18A1
weights(7)  = abs(serotonin_prcc(6))/2;               % SLC18A2
weights(8)  = mean( abs(serotonin_prcc([7,8])) );     % SLC6A4
weights(9)  = mean( abs(serotonin_prcc([9,10])) )/2;  % MAOA
weights(10) = mean( abs(serotonin_prcc([9,10])) )/2;  % MAOB
weights(11) = 0.5*abs(serotonin_prcc(1))/3 + ...
              0.5*abs(serotonin_prcc(11));            % HTR1A
tmp.Gene   = {'TPH1', 'TPH2', 'QDPR', 'SLC7A5', 'DDC',  ...
              'SLC18A1', 'SLC18A2', 'SLC6A4', ...
              'MAOA', 'MAOB', 'HTR1A'}';
tmp.Weight = weights/sum(weights);

% Save weights
writetable( struct2table(tmp), 'Weights_Serotonin.xlsx');