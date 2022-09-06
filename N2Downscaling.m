%                   N2Downscaling.m
%-----------------------------------------------
%           Dr. Joakim Munkhammar, PhD 2022
% 
% This function generates minute resolution global,
% beam and diffuse clear-sky indices from hour 
% resolution global clear-sky index. This is based
% on the Markov-chain mixture distribition model for
% multiple time-series generation at Github: 
% github.com/JoakimMunkhammar/MultiComponentMarkov
%
% The specific downscaling method was published 
% in proceedings of (Solar Integration Workshop 2022)
% *Link*
%
% GlobalIn has to be Global Horizontal Irradiance
% clear-sky index on format [1,N]. Output is Global,
% Beam and Diffuse clear-sky index components on
% minute resolution for the hourly inputs.
% 
% For use, a choice of transition matrices and parameters
% have to be done below by commenting out those which
% will not be used.

function [Global,Beam,Diffuse] = N2Downscaling(GlobalIn)

% NaN-warning on input data
if sum(isnan(GlobalIn)>0)
    disp('Input data contains NaN')
end

% Load transition matrices and parameters
% Comment out those not used, either NREL
% or SMHI transition matrix data should be used.
% Default is NREL setting.

% NREL matrices and parameters
load('PNREL1.mat')
load('PNREL2.mat')
load('PNREL3.mat')
load('PNREL4.mat')
MinBHI = [0 0 0 0];
MaxBHI = [1.0069 0.8603 0.8675 0.8696];
MinDHI = [0 0.1129 0.0713 0.0749];
MaxDHI = [0.8785 0.6880 0.7803 0.5638];

% SMHI matrices and parameters
%load('PSMHI1.mat')
%load('PSMHI2.mat')
%load('PSMHI3.mat')
%load('PSMHI4.mat')
%MinBHI = [0 0 0 0];
%MaxBHI = [0.7036 0.9333 0.9297 1.0069];
%MinDHI = [0 0 0 0];
%MaxDHI = [1.1788 1.0105 0.8232 0.8785];

% The total number of input points
M = size(GlobalIn,2);

% Initiang the level variable
GlobalLevel = zeros(M,1);

% Assigning level values to input values
GlobalLevel(find(GlobalIn<0.5)) = 1;
GlobalLevel(intersect(find(GlobalIn>=0.5),find(GlobalIn<0.7))) = 2;
GlobalLevel(intersect(find(GlobalIn>=0.7),find(GlobalIn<0.9))) = 3;
GlobalLevel(find(GlobalIn>=0.9)) = 4;

N=30; % N number of states
T=60; % Pre-set number of minutes
Prob=0.0000000001; % P_b, base probability in transition matrix 
NewState = zeros(T*M,2); % Initiate the vector
Statei = floor(29*rand(2,1))+1; % Initial condition of the state
for k=1:M
  
    % Determine which transition matrix to use
    if GlobalLevel(k)==1
        P = P1;
    end
    if GlobalLevel(k)==2
        P = P2;
    end
    if GlobalLevel(k)==3
        P = P3;
    end
    if GlobalLevel(k)==4
        P = P4;
    end

    % The loop for creating double output time-series
    % NewDist=(Beam,Diffuse) as output
    for t=1:T 
        Pnew = squeeze(P(Statei(1),Statei(2),:,:))./sum(sum(squeeze(P(Statei(1),Statei(2),:,:))));
        Pnew2 = reshape(Pnew,[1,N^2]);
        Prow = zeros(N^2,1);
        for i=1:N^2 % Make a CDF of the transition matrix
            Prow(i+1) = sum(Pnew2(1:i));
        end
        position = find(Prow(:)<rand(1),1,'last'); % Sample from the CDF
        positionmat = zeros(N^2,1);
        positionmat(position) = 1;
        positionmat = reshape(positionmat,[N,N]);
        [Statei(1),Statei(2)] = find(positionmat==1); % Express the sample as states   
        NewState(t+(k-1)*60,:) = Statei(:); % NewState is the new state at time t for each component
        Level(t+(k-1)*60) = GlobalLevel(k);
    end % Below: add the random uniform number (from each emission probability)
end

% Generating the noise on the output
% Returning minute resolution Global, Beam and Diffuse clear-sky index
% components
Beam = MaxBHI(Level)'.*(NewState(:,1)-1)/N+(1/N)*(MaxBHI(Level)'-MinBHI(Level)').*rand(T*M,1);
Diffuse = MaxDHI(Level)'.*(NewState(:,2)-1)/N+(1/N)*(MaxDHI(Level)'-MinDHI(Level)').*rand(T*M,1);
Global = Beam+Diffuse;