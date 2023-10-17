 function [Data] = getForceAndForceCursor(Data,setup)

% Returns data structure containing the forceCursor and the rawForces. Also
%   returns variables ftCal and ftCalManual, which are automatically and 
%   manually calculated raw forces.  This allows user to check that forces
%   are being calculated correctly. 

% User specifies "setup" - the setup specifies the bias, F/T calibration
%   matrix, and sample period. 
%
%   Calib Exp - setup used for calibration in 2019. Here, the F/T was
%   inside a plastic box mounted to 80/20 with a metal handbar attached.
%   There was no plastic handguard

%   EL - setup used for EL experiments prior to 2021. Here, the F/T was
%   inside a plastic box mounted to 80/20 with a metal handbar attached.
%   There was often a plastic handguard used
%
%   shoulder_bar - setup for posture shoulder device with metal handbar
%   mounted vertically
%
%   shoulder_posture_bci - setup for earl's posture/bci experiments in the 
%   shoulder posture device (7/15/21)
%

%% Identify force transducer analog channels
    algCh = Data(1).Definitions.analogChannelNames;
    chName = {'Force X','Force Y','Force Z','Torque X','Torque Y','Torque Z'};
    chIdx = NaN(size(chName));
    for i = 1:length(chName)
        if ~sum(strcmpi(algCh,chName{i}))
            error('Force transducer channel %s not found',chName{i})
        end
        chIdx(i) = find(strcmpi(algCh,chName{i}),1);
    end

%% Load setup specific data
    switch setup
        case 'Calib Exp'
            bias = [0.122480627877372,1.17789258248465,-0.945744447980608,-0.726117640904018,0.831427143641881,-0.832979437691825];
            samplePeriod = 3;
            cal = load('FTcalMatrix.mat');
        case 'EL'
            bias = [0.114942589393029,1.23510163292518,-0.864382826585036,-0.846398620605469,0.750161611703726,-0.811235567533053];
            samplePeriod = 3;
            cal = load('FTcalMatrix.mat');
        case 'shoulder_bar'
            bias = [0.329106,1.22282,-1.05621,-0.884235,0.664902,-0.923956];
            samplePeriod = 50;
            cal = load('FTcalMatrix20210519.mat');
        case 'shoulder_posture_bci'
            bias = [0.330926238415288,1.989506536745558,9.999389648437500,-9.999694824218750,9.999389648437500,-0.557138432521446];
            samplePeriod = 50;
            cal = load('FTcalMatrix20210519.mat');
    end
    
%% Calculate force cursor position from analog data
    cal = cal.cal;
    for i = 1:length(Data)
        %Load Analog Data and forceTransformation
        ftTransform = Data(i).Definitions.forceTransformation;
        ftAlg = Data(i).TrialData.analogData(:,chIdx);
        
        %Get Actual Time
        time = getActualTime(Data(i).Definitions.analogChannelNames,Data(i).TrialData.analogData);
        downsampledTime = time(1:samplePeriod:size(time,1));
        
        %Remove bias, select channels
        ftCal = cal.ftCalMat*(ftAlg-ones(size(ftAlg,1),1)*ftTransform.offset)';
        allFTOutput = ftCal';
        ftCal = ftCal(ftTransform.chSelect==1,:);

        %Get Raw forces and Torques - no rotation applied.
        absoluteForcesAndTorques = cal.ftCalMat*(ftAlg-ones(size(ftAlg,1),1)*ftTransform.offset)';
        absoluteForcesAndTorques = [downsampledTime'; absoluteForcesAndTorques(:,1:samplePeriod:end)]';
        
        %Calculate transformation matrices.
        thX = calcRotMat(ftTransform.rotation(1),1);
        thY = calcRotMat(-ftTransform.rotation(2),2);
        thZ = calcRotMat(ftTransform.rotation(3),3);
        scaleMat = diag([ftTransform.scaling 1]);
        transMat = [[eye(3) ftTransform.translation']; 0 0 0 1];

        %Apply transformations. Rotate about Z, Y, X, Scale, Translate.
        ftc = transMat*scaleMat*thX*thY*thZ*[ftCal; ones(1,size(ftCal,2))];
        ftCalRotated = thX*thY*thZ*[ftCal; ones(1,size(ftCal,2))]; %These are the exact forces from which the forceCursor is calculated

        %Downsample
        ftc = [downsampledTime'; ftc(1:3,1:samplePeriod:end)]';
        ftCal = [1:samplePeriod:size(ftCal,2); ftCal(1:3,1:samplePeriod:end)]'; 
        ftCalRotated = [1:samplePeriod:size(ftCalRotated,2); ftCalRotated(1:3,1:samplePeriod:end)]';

        %Calculate rawForces Manually
        ftAlg = ftAlg-bias;
        ftCalManual = cal.ftCalMat*(ftAlg)';
%         %Alan's Data
%         rotMatX = getRotMat(-2.56535,1);
%         rotMatY = getRotMat(0.00655,2);
%         rotMatZ = getRotMat(22.19945,3);
% %         %Calib Data
% %         rotMatX = getRotMat(0.0406232896147512,1);
% %         rotMatY = getRotMat(-0.00162400925882222,2);
% %         rotMatZ = getRotMat(20.4955596261533,3);
% 
%         %ZYX Fixed Angles
%         thX = calcRotMat(2.56535,1);
%         thY = calcRotMat(-0.00655,2);
%         thZ = calcRotMat(-22.19945,3);
% 
%         %NEED TO CORRECT TIME IF POSSIBLE - TO-DO 10/24/20
%         
        ftCalManualIDTest = [1:samplePeriod:size(ftCalManual,2); thX*thY*thZ*[ftCalManual(1:3,1:samplePeriod:end); ones(1,size(ftCalManual(1:3,1:samplePeriod:end),2))]]'; %downsample
%         ftCalManual = [1:samplePeriod:size(ftCalManual,2); rotMatX*rotMatY*rotMatZ*ftCalManual(1:3,1:samplePeriod:end)]'; %calculate every 3rd sample - These are the forces with X and Y rotations included
%         ftCalManual = [downsampledTime'; rotMatX*rotMatY*rotMatZ*ftCalManual(1:3,1:samplePeriod:end)]'; %downsample - These are the forces with X and Y rotations included
        allFTOutput = allFTOutput(1:samplePeriod:end,:);
        %Add total force
        ftCalManual = ftCalManualIDTest;
        %Uncomment line below (manual ID test for testing purposes only)
%         ftCalManual = ftCalRotated;
        ftCalManual(:,1) = downsampledTime';
        ftCalManual(:,5) = sqrt(ftCalManual(:,2).^2 + ftCalManual(:,3).^2 + ftCalManual(:,4).^2);
        
        Data(i).TrialData.Marker.forceCursor = ftc;
        Data(i).TrialData.Marker.rawForces = ftCalManual;
        Data(i).TrialData.Marker.allFTOutput = allFTOutput;
        Data(i).TrialData.Marker.absoluteForcesAndTorques = absoluteForcesAndTorques;
    end

%% Function for calculating rotation matrices (based on clockwise +?)
    function rotMat = calcRotMat(theta,dim)
        %Calculate rotation matrices about arbitrary dimension
        rotMat = eye(4);
        rotDim = mod(dim+[0 1],3)+1; %X: [2 3], Y: [3 1], Z: [1 2]
        rotMat(rotDim,rotDim) = cosd(theta);
        rotMat(rotDim(2),rotDim(1)) = -sind(theta);
        rotMat(rotDim(1),rotDim(2)) = sind(theta);
    end


    function R = getRotMat(theta,dim)
        if dim == 1
            R = [1 0 0; 0 cosd(theta) -sind(theta); 0 sind(theta) cosd(theta)];
        elseif dim == 2
            R = [cosd(theta) 0 sind(theta); 0 1 0; -sind(theta) 0 cosd(theta)];
        elseif dim == 3
            R = [cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; 0 0 1];
        end
    end

end