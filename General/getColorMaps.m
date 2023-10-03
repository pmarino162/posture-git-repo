function [pcmap,tcmap,rainbow] = getColorMaps()
    %Colormaps used for posture paper
    
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\orli.mat')
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\customRainbow.mat')
    load('C:\Users\pmari\OneDrive - University of Pittsburgh\Documents\Matlab Repository\marino\Palettes\purpleAndGreen2.mat')
   
    pcmap = orli;
    tcmap = purpleAndGreen;
    rainbow = customRainbow;
    
end