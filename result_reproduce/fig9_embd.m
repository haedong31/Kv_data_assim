%% Figure 9 Low-dimensional embedding
clc
clearvars
close all

load('embd_cell_names.mat')

exp_num = "exp51";
base_dir = fullfile(pwd,strcat("calib_",exp_num));

% pktof
pktof = readtable(fullfile(base_dir,"pktof.csv"));
pidx = unique(pktof.param);

nrow = size(pktof,1)/length(pidx);
pktof2 = array2table(NaN(nrow,length(pidx)+2));
pktof2.Properties.VariableNames(1:2) = {'Group','File'};
pktof2.Group = pktof(pktof.param==pidx(1),:).Group;
pktof2.File = pktof(pktof.param==pidx(1),:).file;

for i=1:length(pidx)
    pktof2.Properties.VariableNames(i+2) = {strcat('P',num2str(pidx(i)))};
    pktof2(:,i+2)= pktof(pktof.param==pidx(i),'value');
end

% pkslow1
pkslow1 = readtable(fullfile(base_dir,"pkslow1.csv"));
pidx = unique(pkslow1.param);

nrow = size(pkslow1,1)/length(pidx);
pkslow12 = array2table(NaN(nrow,length(pidx)+1));
pkslow12.Properties.VariableNames(1) = {'Group'};
pkslow12.Group = pkslow1(pkslow1.param==pidx(1),:).Group;

for i=1:length(pidx)
    pkslow12.Properties.VariableNames(i+1) = {strcat('P',num2str(pidx(i)))};
    pkslow12(:,i+1)= pkslow1(pkslow1.param==pidx(i),'value');
end

% pkslow2
pkslow2 = readtable(fullfile(base_dir,"pkslow2.csv"));
pidx = unique(pkslow2.param);

nrow = size(pkslow2,1)/length(pidx);
g = categorical(pkslow2(pkslow2.param==pidx(1),:).Group);
v1 = pkslow2(pkslow2.param==pidx(1),:).value;
v2 = pkslow2(pkslow2.param==pidx(2),:).value;
pkslow22 = table(g,v1,v2);


% pkss
pkss = readtable(fullfile(base_dir,"pkss.csv"));
pidx = unique(pkss.param);

pkss2 = array2table(NaN(nrow,length(pidx)+1));
pkss2.Properties.VariableNames(1) = {'Group'};
pkss2.Group = pkss(pkss.param==pidx(1),:).Group;

for i=1:length(pidx)
    pkss2.Properties.VariableNames(i+1) = {strcat('P',num2str(pidx(i)))};
    pkss2(:,i+1)= pkss(pkss.param==pidx(i),'value');
end

c = NaN(size(pktof2,1),3);
idx1 = pktof2.Group=="WT";
idx2 = pktof2.Group=="Mgat1KO";

c(idx1,:) = repmat([0 0 1],sum(idx1),1);
c(idx2,:) = repmat([1 0 0],sum(idx2),1);

pmx = [table2array(pktof2(:,3:end)),table2array(pkslow12(:,2:end)),table2array(pkslow22(:,2:end)),table2array(pkss2(:,2:end))];

rng(7981)
embd = tsne(normalize(pmx),'Algorithm','exact','Distance','minkowski','NumDimensions',3);

det(cov(embd(idx1,:)))
det(cov(embd(idx2,:)))

% scatter plot (embedding)
fig = figure('Color','w','Position',[100,100,600,500]);
orient(fig,'landscape')
s1 = scatter3(embd(idx1,1),embd(idx1,2),embd(idx1,3),'o','CData',c(pktof2.Group=="WT",:));
s1.SizeData = 90;
s1.LineWidth = 1.5;
hold on
s2 = scatter3(embd(idx2,1),embd(idx2,2),embd(idx2,3),'^','CData',c(pktof2.Group=="Mgat1KO",:));
s2.SizeData = 90;
s2.LineWidth = 1.5;
hold off
xlabel("Embedding Dim1")
ylabel("Embedding Dim2")
zlabel("Embedding Dim3")
axis tight
grid off
legend(["WT","MGAT1KO"],'location','best')
set(gca,'FontWeight','bold','LineWidth',1.5)

%% Quantify the variability
pmx_wt = embd(idx1,:);
pmx_ko = embd(idx2,:);
cov_wt = cov(pmx_wt);
cov_ko = cov(pmx_ko);
det(cov_wt)
det(cov_ko)