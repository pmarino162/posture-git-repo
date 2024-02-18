function [pcmap,tcmap,rainbow] = getColorMaps(numPostures)
    %Colormaps used for posture paper
    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\purpleAndGreen2.mat')
   
    pcmap = orli;
    tcmap = purpleAndGreen;
    rainbow = customRainbow;
    
    switch numPostures
        case 2
            pcmap = pcmap([1,5],:);
        case 3
            pcmap = pcmap([1,3,5],:);
        case 4 %This covers the elbow/shoulder dataset, for which posture IDs are: 1=I30 3=E30 4=F30 5=eE30
            pcmap = pcmap([5,2,4,1,3],:);
            %I30 (1) - dark blue
            %E30 (3) - light blue
            %F30 (4) - red
            %eE30(5) - yellow
        %No case for 4, because postures are labeled 1,3,4,5 in that exp
        case {7,6}
            pcmap = vertcat(pcmap,rainbow(4:5,:));
    end
    
end