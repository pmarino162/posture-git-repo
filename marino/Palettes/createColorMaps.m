clear;

%% Red orange yellow green 
    red = [199, 0, 0];
    orange = [255, 98, 0];
    yellow = [255, 195, 0];
    green = [17, 148, 0];
    
    column1 = [interp1([0,1],[red(1),orange(1)],linspace(0,1,25)),interp1([0,1],[orange(1),yellow(1)],linspace(0,1,25)),interp1([0,1],[yellow(1),green(1)],linspace(0,1,25))]';
    column2 = [interp1([0,1],[red(2),orange(2)],linspace(0,1,25)),interp1([0,1],[orange(2),yellow(2)],linspace(0,1,25)),interp1([0,1],[yellow(2),green(2)],linspace(0,1,25))]';
    column3 = [interp1([0,1],[red(3),orange(3)],linspace(0,1,25)),interp1([0,1],[orange(3),yellow(3)],linspace(0,1,25)),interp1([0,1],[yellow(3),green(3)],linspace(0,1,25))]';
    
    royg = [column1,column2,column3];
    royg = royg./255;
    
    clearvars red orange yellow green column1 column2 column3 
   % save('royg.mat')
 
%% Orli's colormap
    red = [186,61,22];
    orange = [222,119,48];
    %yellow = [213,216,110];
    yellow = [213,200,10]; %made it darker
    blue = [71,137,191];
    darkblue = [0,80,117];
    colors = vertcat(red,orange,yellow,blue,darkblue);
    colors = colors./255;
    figure
    for i = 1:size(colors,1)
        plot(i,1,'.','MarkerSize',200,'Color',colors(i,:));
        hold on 
    end
    orli = colors;
    clearvars i red orange yellow blue darkblue colors
    save('orli.mat')
    
%% Haystacks: Sunset 
    pink = [205,148,146];
    orange = [243,147,109];
    yellow = [251,179,92];
    green = [105,94,73];
    blue = [69,104,181];
    peri = [138,158,199];
    brown = [172,116,95];
    colors = vertcat(pink,orange,yellow,green,blue,peri,brown);
    colors = colors./255;
    figure
    for i = 1:size(colors,1)
        plot(i,1,'.','MarkerSize',200,'Color',colors(i,:));
        hold on
        
    end
    
    %% Haystacks: End of summer
	purple = [105,119,154];
    pink = [184,111,94];
    orange = [239,165,114];
    yellow = [240,194,117];
    green = [107,129,76];
    blue = [83,118,142];
    brown = [134,106,67];
    colors = vertcat(pink,orange,yellow,green,blue,purple,brown);
    colors = colors./255;
    figure
    for i = 1:size(colors,1)
        plot(i,1,'.','MarkerSize',200,'Color',colors(i,:));
        hold on
        
    end
    
    %% Haystacks: End of summer 2
	purple = [103,124,173];
    pink = [164,92,83];
    orange = [210,133,92];
    yellow = [202,169,72];
    green = [91,136,107];
    brown = [131,90,57];
    gray = [123,121,121];
    colors = vertcat(pink,orange,yellow,green,purple,brown,gray);
    colors = colors./255;
    figure
    for i = 1:size(colors,1)
        plot(i,1,'.','MarkerSize',200,'Color',colors(i,:));
        hold on
        
    end
    
    haystacks = colors;
    clearvars i pink orange yellow green purple brown gray colors
    save('haystacks.mat')
%% ColorBrewer Rainbow
cbRainbow = [...
    158,1,66;...
    213,62,79;...
    244,109,67;...
    253,174,97;...
    254,224,139;...
    255,255,191;...
    230,245,152;...
    171,221,164;...
    102,194,165;...
    50,136,189;...
    94,79,162];

cbRainbow = cbRainbow./255;
save('cbRainbow.mat')

%% ColorBrewer Rainbow8
clear
cbRainbow8 = [...
    213,62,79;...
    244,109,67;...
    253,174,97;...
    254,224,139;...
    230,245,152;...
    171,221,164;...
    102,194,165;...
    50,136,189];

cbRainbow8 = cbRainbow8./255;
save('cbRainbow8.mat')

%% ColorBrewer Red Yellow Blue
clear
cbRYB = [...
    165,0,38;...
    215,48,39;...
    244,109,67;...
    253,174,97;...
    254,224,144;...
    255,255,191;...
    224,243,248;...
    171,217,233;...
    116,173,209;...
    69,117,180;...
    49,54,149];

cbRYB = cbRYB./255;
save('cbRYB.mat')

%% ColorBrewer Red Blue
clear
cbRB = [...
    03,0,31;...
    178,24,43;...
    214,96,77;...
    244,165,130;...
    253,219,199;...
    247,247,247;...
    209,229,240;...
    146,197,222;...
    67,147,195;...
    33,102,172;...
    5,48,97];

cbRB = cbRB./255;
save('cbRB.mat')

%% ColorBrewer Pink Green
clear
cbPG = [...
    142,1,82;...
    197,27,125;...
    222,119,174;...
    241,182,218;...
    253,224,239;...
    247,247,247;...
    230,245,208;...
    184,225,134;...
    127,188,65;...
    77,146,33;...
    39,100,25];
    
cbPG = cbPG./255;  
save('cbPG.mat')

%% Custom Rainbow
clear
customRainbow = [...
    199, 0, 0;...  
    255, 98, 0;...
    255, 195, 0;...
    17, 148, 0;...
    0, 161, 179;...
    0, 71, 153;...
    83, 0, 184;...
    179, 0, 152];
customRainbow = customRainbow./255;  

    figure
    for i = 1:size(customRainbow,1)
        plot(i,1,'.','MarkerSize',200,'Color',customRainbow(i,:));
        hold on
        
    end

%save('customRainbow.mat')

%% Custom Rainbow No Green
clear
customRainbowNoGreen = [...
    230, 0, 0;...  
    255, 128, 0;...
    240, 208, 0;...
    0, 224, 202;...
    0, 176, 224;...
    0, 67, 224;...
    153, 0, 224;...
    224, 0, 131];
customRainbowNoGreen = customRainbowNoGreen./255;  
save('customRainbowNoGreen.mat')

%% Custom Summer
clear
customSummer = [horzcat(interp1([0,1],[0,240],linspace(0,1,5))',interp1([0,1],[128,208],linspace(0,1,5))',interp1([0,1],[96,0],linspace(0,1,5))')];
customSummer = customSummer./255;  
save('customSummer.mat')

%% Custom Blue
clear
customBlue = [horzcat(interp1([0,1],[0,0],linspace(0,1,5))',interp1([0,1],[67,176],linspace(0,1,5))',interp1([0,1],[224,224],linspace(0,1,5))')];
customBlue = customBlue./255;  
save('customBlue.mat')

%% Rainbow2d
clear
R=[1 0;
   1 0];
G=[1 1
   0 0];
B=[0 0
   0 1];
R = interp2(R,8);
G = interp2(G,8);
B = interp2(B,8);
I = uint8(230*cat(3,R,G,B));
image(I)
pixThird = round(size(I,1)/3);
pixHalf = round(size(I,1)/2);
rainbow2d = [...
%Top row    
I(1,1,1), I(1,1,2), I(1,1,3);
I(1,pixThird,1), I(1,pixThird,2), I(1,pixThird,3);
I(1,2*pixThird,1), I(1,2*pixThird,2), I(1,2*pixThird,3);
I(1,end,1), I(1,end,2), I(1,end,3);
%Middle posture
I(pixThird,pixHalf,1), I(pixThird,pixHalf,2), I(pixThird,pixHalf,3);
%Bottom row
I(2*pixThird,pixThird,1), I(2*pixThird,pixThird,2), I(2*pixThird,pixThird,3);
I(2*pixThird,2*pixThird,1), I(2*pixThird,2*pixThird,2), I(2*pixThird,2*pixThird,3);
];

figure
hold on
for i = 1:4
   plot(i,3,'.','MarkerSize',20,'Color',rainbow2d(i,:)) 
end
    plot(2.5,2,'.','MarkerSize',20,'Color',rainbow2d(5,:)) 
for i = 6:7
    plot(i-4,1,'.','MarkerSize',20,'Color',rainbow2d(i,:)) 
end
clearvars B G R i I pixHalf pixThird
%save('rainbow2d.mat')

%% Green and purple colormap for targets
clear;
    yellow = [255, 195, 0]./255;
    green = [17, 148, 0]./255;
    darkpurple = [83, 0, 184]./255;
    purple = [179, 0, 152]./255;
    lightblue = [0, 161, 179]./255;
    salmon = [1 .549 .412];
    darkblue = [0, 71, 153]./255;
    
    numSteps = 4;
    
    column = 1;
    column1 = [interp1([0,1],[green(column),lightblue(column)],linspace(0,1,2)),...
        interp1([0,1],[lightblue(column),purple(column)],linspace(1/2,1,2)),...
        interp1([0,1],[purple(column),salmon(column)],linspace(1/2,1,2)),...
        interp1([0,1],[salmon(column),green(column)],linspace(1/4,1/2,2))]';
    column = 2;
    column2 = [interp1([0,1],[green(column),lightblue(column)],linspace(0,1,2)),...
        interp1([0,1],[lightblue(column),purple(column)],linspace(1/2,1,2)),...
        interp1([0,1],[purple(column),salmon(column)],linspace(1/2,1,2)),...
        interp1([0,1],[salmon(column),green(column)],linspace(1/4,1/2,2))]';
    column = 3;
    column3 = [interp1([0,1],[green(column),lightblue(column)],linspace(0,1,2)),...
        interp1([0,1],[lightblue(column),purple(column)],linspace(1/2,1,2)),...
        interp1([0,1],[purple(column),salmon(column)],linspace(1/2,1,2)),...
        interp1([0,1],[salmon(column),green(column)],linspace(1/4,1/2,2))]';
    
%         interp1([0,1],[purple(1),salmon(1)],linspace(0,1,numSteps))]';
%     column2 = [interp1([0,1],[green(2),purple(2)],linspace(0,1,numSteps)),...
%         interp1([0,1],[purple(2),salmon(2)],linspace(0,1,numSteps))]';
%     column3 = [interp1([0,1],[green(3),purple(3)],linspace(0,1,numSteps)),...
%         interp1([0,1],[purple(3),salmon(3)],linspace(0,1,numSteps))]';
%    
%     numSteps = 8;
%    column1 = [interp1([0,1],[green(1),purple(1)],linspace(0,1,numSteps))]';
%     column2 = [interp1([0,1],[green(2),purple(2)],linspace(0,1,numSteps))]';
%     column3 = [interp1([0,1],[green(3),purple(3)],linspace(0,1,numSteps))]';
%     
    
    purpleAndGreen = [column1,column2,column3];
    



    figure
    for i = 1:size(purpleAndGreen,1)
        plot(i,1,'.','MarkerSize',200,'Color',purpleAndGreen(i,:));
        hold on
        
    end
    clearvars  i numSteps purple yellow green darkpurple lightblue column1 column2 column3 darkblue salmon column
   % save('purpleAndGreen.mat')
   
%% Purple and green 2
clear;

    %Even brighter
%     green = [25, 210, 0]./255;
%     purple = [220, 0, 185]./255;
    
    
    green = [25, 210, 1]./255;
    purple = [247, 1, 155]./255;
    
%     %Purple to green over all 8 targets
%     numSteps = 8;
%     column1 = [interp1([0,1],[green(1),purple(1)],linspace(0,1,numSteps))]';
%     column2 = [interp1([0,1],[green(2),purple(2)],linspace(0,1,numSteps))]';
%     column3 = [interp1([0,1],[green(3),purple(3)],linspace(0,1,numSteps))]';
    

    %Purple to green over first 5 targets
    numSteps = 5;
    column1 = [interp1([0,1],[green(1),purple(1)],linspace(0,1,numSteps))]';
    column2 = [interp1([0,1],[green(2),purple(2)],linspace(0,1,numSteps))]';
    column3 = [interp1([0,1],[green(3),purple(3)],linspace(0,1,numSteps))]';
    column1 = vertcat(column1,flip(column1(2:4)));
    column2 = vertcat(column2,flip(column2(2:4)));
    column3 = vertcat(column3,flip(column3(2:4)));
    
    
    purpleAndGreen = [column1,column2,column3];
    
    figure
    for i = 1:size(purpleAndGreen,1)
        plot(i,1,'.','MarkerSize',200,'Color',purpleAndGreen(i,:));
        hold on
        
    end
   clearvars  i numSteps purple  green  lightblue column1 column2 column3   
   save('purpleAndGreen2.mat')
