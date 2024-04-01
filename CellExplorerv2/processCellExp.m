load('~/Downloads/normalized_waveforms.mat');
load('~/Downloads/narrow_acg.mat');
load('~/Downloads/CellTypes.mat');
load ~/code/CellExplorer/+groundTruthData/groundTruth_cell_metrics.cellinfo.mat
load ~/code/CellExplorer/+groundTruthData/groundTruth_cell_list.mat

X = normalizedWaveforms';

fN = fieldnames(CellTypes);

for f=1:length(fN);
    for g=setdiff(1:length(fN),f)
        fprintf('\n %d, %d',f,g);
        intersect(CellTypes.(fN{f}), CellTypes.(fN{g}))
    end
end

CellTypes.VGAT = setdiff(CellTypes.VGAT, CellTypes.PV);
newCellTypes = CellTypes;

%%
CellTypeNames = char(size(X,1),4);

for f=1:length(fN)
    d = min(length(fN{f}),4)
    fNc = fN{f}
    indices = CellTypes.(fNc);
    CellTypeNames(indices,[1:d]) = repmat(strtrim(fNc(1:d)),length(indices),1);
end
Y = narrow_acg';
% save('-mat','finalWaveforms.mat','X','Y','CellTypeNames');

%%
acgN = gct.acg.narrow_normalized;
acg = gct.acg.narrow;
acgl = gct.acg.log10_rate;
acgW = gct.acg.wide;


save('-mat','CellExplorer_ACG.mat','acgN','acg','acgW','acgl');

%%
nonAllenNeurons = setdiff(1:size(gt.dataTable,1), ...
    strmatch('Allen', gt.dataTable(:,6)))
validNeurons = zeros(1,length(gt.dataTable(:,1)));
validNeurons(nonAllenNeurons) = 1;
gct = groundTruth_cell_metrics;
features = cat(2, gct.acg_asymptote', ...
                        gct.acg_c',...
                        gct.cv2',...
                        gct.troughToPeak',...
                        gct.troughtoPeakDerivative',...
                        gct.acg_d',...
                        gct.acg_h',...
                        gct.ab_ratio',...
                        gct.acg_refrac',...
                        gct.acg_tau_burst',...
                        gct.acg_tau_decay',...
                        gct.acg_tau_rise',...
                        gct.thetaModulationIndex');
                    
features = cat(2, gct.acg_asymptote', ...
                        gct.acg_c',...
                        gct.cv2',...
                        gct.troughToPeak',...
                        gct.troughtoPeakDerivative',...
                        gct.acg_d',...
                        gct.acg_h',...
                        gct.ab_ratio',...
                        gct.acg_refrac',...
                        gct.acg_tau_burst',...
                        gct.acg_tau_decay',...
                        gct.acg_tau_rise');
                    
areaIds = zeros(1,length(gt.dataTable(:,1)));
visMatch = strmatch('VIS', gt.dataTable(:,5));
hipMatch = strmatch('HIP',gt.dataTable(:,5));
ca1Match = strmatch('CA1',gt.dataTable(:,5));
ca3Match = strmatch('CA3',gt.dataTable(:,5));

areaIds(visMatch) = 1;
areaIds(hipMatch) = 2;
areaIds(ca1Match) = 1;
areaIds(ca3Match) = 1;
                       
 source = {};    
 area = {};
for n=1:length(validNeurons)
    if validNeurons(n) == 0
        source{n} = 'Allen';
        sourceFull{n} = ['Allen' CellTypeNames(n,:)];
    else
        source{n} = 'NonAllen';
        sourceFull{n} = ['NonAllen' CellTypeNames(n,:)];

    end
    
    Temp = gt.dataTable(n,5);
    if areaIds(n) == 1
        
        area{n} = Temp{1};
    elseif areaIds(n) == 2
        area{n} = Temp{1};
    else
        area{n} = 'OTHER';
    end
end
    
save('-mat','features_cellExp_new_Nov2023.mat','features', ...
    'validNeurons', 'nonAllenNeurons','source','sourceFull','area');