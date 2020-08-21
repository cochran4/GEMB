function v = velocitySerotonin
%--------------------------------------------------------------------------
% Velocity information for Serotonin 
%--------------------------------------------------------------------------

% Tryptophan hydroxylase
v(1).idx   = [1,1;2,-1; 3,-1; 4,1]; % Substance index and positive/negative contribution
v(1).fname = 'v_TPH(y(2),y(3),y(7),v(1).param)';
v(1).param = [40,20,1000,400];

% Dihydropteridine reductase
% Note:
%   Concentrations of NADP and NADPH were not in paper
%   Concentrations from dopamine did not match rates in figure
%   Contrations 365 and 5000 matched rates in figure
%   Not sure where these number came from, since they were from earler
%   version
v(2).idx   = [1, -1;2,1];
v(2).fname = 'v_4Var(y(1),365,y(2),5000,v(2).param)';  
v(2).param = [100 75 5000 10 75 3];

% Neutral amino acid transporter
% Note:
%   Parameters from text, NOT Table
%   Parameters in table were from dopamine paper
v(3).idx   = [3,1];
v(3).fname = 'v_1Var(97,v(3).param)';
v(3).param = [330, 700]; 

% Aromatic amino acid decarboyxlase
v(4).idx   = [4, -1; 5, 1];
v(4).fname = 'v_1Var(y(4),v(4).param)';
v(4).param = [160 400];

% Vesicular monoamine transporter
% Note:
%   Changed last parameter (kout) in order to match rates in steady state
%   figure; kout in table is 40
%   Parameters in table were from dopamine paper
v(5).idx   = [5, -1; 6, 1];
v(5).fname = 'v_1Var(y(5),v(5).param(1:2))-v(5).param(3)*y(6)';
v(5).param = [0.198 3500  115.88]; 

% Fluox and Serotonin transporter
v(6).idx   = [5, 1; 7, -1];
v(6).fname = 'fluox(t)*v_1Var(y(7),v(6).param)';
v(6).param = [0.17 4700];

% Catabolism of c5ht
v(7).idx   = [5,-1; 8, 1];
v(7).fname  = 'v_1Var(y(5),v(7).param)';
v(7).param = [95,1000];

% Catabolism of e5ht
v(8).idx   = [7,-1; 8, 1];
v(8).fname = 'v_1Var(y(7),v(8).param)';
v(8).param = [95,1000];

% Release of v5ht out of vesciles
v(9).idx   = [6,-1;7,1];
v(9).fname = 'v(9).param*release(y(7))*fire(t)*y(6)';
v(9).param = 1;

% Catabolism of 5hiaa
v(10).idx   = [8,-1];
v(10).fname = 'v(10).param*y(8)';
v(10).param = 1;

% Catabolism of trp-pool
% Note:
%   Parameter was chosen to match steady state figure
%   Parameter is 0.2 in table, which is the same value for dopamine
v(11).idx   = [9,-1];
v(11).fname = 'v(11).param*y(9)';
v(11).param = 0.8; 

% Catabolism of trp
% Note:
%   Parameter was chosen to match steady state figure
%   Parameter is 0.2 in table, which is the same value for dopamine
v(12).idx   = [3,-1];
v(12).fname = 'v(12).param*y(3)';
v(12).param = 1.8; 

% Flow in and out of trp-pool
% Note:
%   Parameters was chosen to match steady state figure 
%   Parameters were constrained to be 10 x the other
%   Parameters were 6 and 0.6 in table, which is same as dopamine
v(13).idx   = [3,-1; 9,1];
v(13).fname = 'v(13).param(1)*y(3) - v(13).param(2)*y(9)';
v(13).param = [19,1.9]; 

% Removal of external Serotonin
v(14).idx   = [7,-1];
v(14).fname = 'v(14).param*y(7)';
v(14).param = 400;
