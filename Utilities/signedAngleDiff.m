function [diff] = signedAngleDiff(angle1Vec,angle2Vec)
%Computes signed CCW angle [-180,180] traversed moving from angle 1 to angle 2.
%Input angles are specified in [0,360] deg 

if size(angle1Vec,1) > 1
    error('inputs must be row vectors')
end

numAngles = length(angle1Vec);

if length(angle2Vec) ~= numAngles
    error('inputs must have equal length')
end

diff = NaN(1,numAngles);

for i = 1:numAngles
    angle1 = angle1Vec(i); angle2 = angle2Vec(i);
    %Compute counterclockwise
    if angle2 > angle1
        counterClock = angle2-angle1;
    else
        counterClock = 360-(angle1-angle2);
    end

    %Compute clockwise
    if angle2 > angle1
        clock = 360-(angle2-angle1);
    else
        clock = angle1-angle2;
    end

    %Choose shorter
    if clock < counterClock
        diff(i) = -clock;
    else
        diff(i) = counterClock;
    end
end

end