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
        %No case for 4, because postures are labeled 1,2,4,5 in that exp
        case {7,6}
            pcmap = vertcat(pcmap,rainbow(4:5,:));
    end
    
end