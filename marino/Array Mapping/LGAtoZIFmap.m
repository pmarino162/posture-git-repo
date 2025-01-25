tdtCh_LGA = [35 40 43 50 67:69 70 74 76 77 81 82 85 87 90 96];
tdtCh_LGA = [21 24 26 60];
[tdt2adapterMap_LGA, adapter2cereportMap_LGA, ~,~,~,~] = getMaps('earlRH_LGAOmnetics');

[tdt2adapterMap_ZIF, adapter2cereportMap_ZIF, ~,~,~,~] = getMaps('earlRHM1');

tdtCh_ZIF = nan(1,length(tdtCh_LGA));

for i = 1:length(tdtCh_ZIF)
    cur_tdtCh_LGA = tdtCh_LGA(i);
    LGA_adapterPin = tdt2adapterMap_LGA{find([tdt2adapterMap_LGA{:,1}]==cur_tdtCh_LGA),2};
    cereportInd = find(contains(adapter2cereportMap_LGA(:,1),LGA_adapterPin));
    cereportPin = adapter2cereportMap_LGA{cereportInd,2};
    ZIF_adapterPin = adapter2cereportMap_ZIF{contains(adapter2cereportMap_ZIF(:,2),cereportPin),1};
    tdtCh_ZIF(i) = tdt2adapterMap_ZIF{contains(tdt2adapterMap_ZIF(:,2),ZIF_adapterPin),1};
end

tdtCh_ZIF