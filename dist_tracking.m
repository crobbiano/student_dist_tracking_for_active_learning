%% New script for identifying and splitting students into groups based on
%  their # of standard deviations away from the mean at the end of the
%  sophomore year.  This should be done at the GPA level.
%
%  Want to then track those groups of students in each of the junior year
%  ECE courses and watch how their distributions change from the sophomore
%  year to the junior year.
%
%  Then compare the changes in these distributions from sophomore year to
%  junior year for both the students that are pre-KI and post-KI and identify
%  if any differences have occured due to the introduction of the KI.

clear all
clc
addpath('histogram_distance')
%% Set Up
load('..\correlational_studies\scripts\rawdata_processing\grades_postKI.mat')

FileName = 'grades_postKI.csv';
FilePath = '..\correlational_studies\scripts\rawdata_processing\';

File = [FilePath '\' FileName];

data = importdata(File, ',');

Classes = data.textdata;
Classes = Classes(numNonClasses+1:end);
ki_status = find(ismember(data.textdata, 'POSTKI'));%6;
num_junior_classes = sum(ismember(data.textdata, 'ECE3'));

%%
Grades         = data.data(:,numNonClasses+1:end);
Gpas           = data.data(:,numNonClasses);
gradesGpas     = mean(Grades,2);
[m,n]          = size(Grades);

%% First do the pre-KI then do post-KI
% column 6 is pre/post ki indicator. 1==post
all(~isnan(Grades(:,end-6:end)))
preKI_Grades = Grades(data.data(:,ki_status)<1,:);
preKI_Gpas = Gpas(data.data(:,ki_status)<1);

% FIXME - don't use nanmean
JuniorGrades = preKI_Grades(:,end-5:end);
JuniorGpas   = nanmean(JuniorGrades,2);
[mj,nj]         = size(JuniorGrades);
pdj = fitdist(JuniorGpas,'Normal');

% FIXME - don't use nanmean
nonJuniorGrades = preKI_Grades(:,end-10:end-6);
nonJuniorGpas   = nanmean(nonJuniorGrades,2);
nonJuniorGpas_pre_KI   = nonJuniorGpas;
[m2,n2]         = size(Grades);
pd = fitdist(nonJuniorGpas,'Normal');

lowest_performance_range = pd.mu - (.5*pd.std + 3*pd.std);
lower_performance_range = pd.mu - (.5*pd.std + 2*pd.std);
low_performance_range = pd.mu - (.5*pd.std + 1*pd.std);
average_performance_range = pd.mu - (.5*pd.std + 0*pd.std);
high_performance_range = pd.mu + (.5*pd.std + 0*pd.std);
higher_performance_range = pd.mu + (.5*pd.std + 1*pd.std);
highest_performance_range = pd.mu + (.5*pd.std + 2*pd.std);

% lowester_performers_mask = (nonJuniorGpas < lowest_performance_range);
lowest_performers_mask = (nonJuniorGpas < lower_performance_range);
lower_performers_mask = (nonJuniorGpas >= lower_performance_range).*(nonJuniorGpas < low_performance_range);
low_performers_mask = (nonJuniorGpas >= low_performance_range).*(nonJuniorGpas < average_performance_range);
average_performers_mask = (nonJuniorGpas >= average_performance_range).*(nonJuniorGpas < high_performance_range);
high_performers_mask = (nonJuniorGpas >= high_performance_range).*(nonJuniorGpas < higher_performance_range);
higher_performers_mask = (nonJuniorGpas >= higher_performance_range).*(nonJuniorGpas < highest_performance_range);

% Sanity check
% if sum(lowest_performers_mask.*lower_performers_mask.*low_performers_mask.*average_performers_mask.*high_performers_mask.*higher_performers_mask) ~= 0
%     error('Overlapping groups')
% elseif sum(lowest_performers_mask)+sum(lower_performers_mask)+sum(low_performers_mask)+sum(average_performers_mask)+sum(high_performers_mask)+sum(higher_performers_mask) ~= mj
%     error('Missing students from groups')
% end

% FIXME - might want to add in anon_id to track these easier, but we should
% be able to make this work as is
lowest_performers = [Grades(lowest_performers_mask>0   , :) , nonJuniorGpas(lowest_performers_mask>0)  , JuniorGpas(lowest_performers_mask>0)  ] ;
lower_performers = [Grades(lower_performers_mask>0     , :) , nonJuniorGpas(lower_performers_mask>0)   , JuniorGpas(lower_performers_mask>0)   ] ;
low_performers = [Grades(low_performers_mask>0         , :) , nonJuniorGpas(low_performers_mask>0)     , JuniorGpas(low_performers_mask>0)     ] ;
average_performers = [Grades(average_performers_mask>0 , :) , nonJuniorGpas(average_performers_mask>0) , JuniorGpas(average_performers_mask>0) ] ;
high_performers = [Grades(high_performers_mask>0       , :) , nonJuniorGpas(high_performers_mask>0)    , JuniorGpas(high_performers_mask>0)    ] ;
higher_performers = [Grades(higher_performers_mask>0   , :) , nonJuniorGpas(higher_performers_mask>0)  , JuniorGpas(higher_performers_mask>0)  ] ;
%% Plot the individual group distributions - Fake normal dists
figure(2); clf;
doFakeFit = 0;
if doFakeFit
    subplot(4, 2, [1 3]); hold on
    
    mu1=mean(lowest_performers(:   , end-1));
    s1=std(lowest_performers(:   , end-1));
    x1=linspace(mu1-4*s1,mu1+4*s1,200);
    pdfx1=1/sqrt(2*pi)/s1*exp(-(x1-mu1).^2/(2*s1^2));
    pd1 = makedist('Normal', mu1, s1);
    pdfx11 = pdf(pd1, x1);
    p1 = plot(x1,pdfx11, 'k--*');
    
    mu2=mean(lower_performers(:   , end-1));
    s2=std(lower_performers(:   , end-1));
    x2=linspace(mu2-4*s2,mu2+4*s2,200);
    pdfx2=2/sqrt(2*pi)/s2*exp(-(x2-mu2).^2/(2*s2^2));
    pdfx22 = pdf('Normal', x2, mu2, s2);
    p2 = plot(x2,pdfx22, 'b--');
    
    mu3=mean(low_performers(:   , end-1));
    s3=std(low_performers(:   , end-1));
    x3=linspace(mu3-4*s3,mu3+4*s3,200);
    pdfx3=3/sqrt(2*pi)/s3*exp(-(x3-mu3).^2/(2*s3^2));
    pdfx33 = pdf('Normal', x3, mu3, s3);
    p3 = plot(x3,pdfx33, 'r--');
    
    mu4=mean(average_performers(:   , end-1));
    s4=std(average_performers(:   , end-1));
    x4=linspace(mu4-4*s4,mu4+4*s4,200);
    pdfx4=4/sqrt(2*pi)/s4*exp(-(x4-mu4).^2/(2*s4^2));
    pdfx44 = pdf('Normal', x4, mu4, s4);
    p4 = plot(x4,pdfx44, 'g--');
    
    mu5=mean(high_performers(:   , end-1));
    s5=std(high_performers(:   , end-1));
    x5=linspace(mu5-4*s5,mu5+4*s5,200);
    pdfx5=5/sqrt(2*pi)/s5*exp(-(x5-mu5).^2/(2*s5^2));
    pdfx55 = pdf('Normal', x5, mu5, s5);
    p5 = plot(x5,pdfx55, 'm--');
    
    mu6=mean(higher_performers(:   , end-1));
    s6=std(higher_performers(:   , end-1));
    x6=linspace(mu6-4*s6,mu6+4*s6,200);
    pdfx6=6/sqrt(2*pi)/s6*exp(-(x6-mu6).^2/(2*s6^2));
    pdfx66 = pdf('Normal', x6, mu6, s6);
    p6 = plot(x6,pdfx66, 'c--');
    
    mu7=mean(nonJuniorGpas);
    s7=std(nonJuniorGpas);
    x7=linspace(mu7-4*s7,mu7+4*s7,200);
    pdfx7=7/sqrt(2*pi)/s7*exp(-(x7-mu7).^2/(2*s7^2));
    pdfx77 = pdf('Normal', x7, mu7, s7);
    p7 = plot(x7,pdfx77, 'k--');
    
    mu8=mean(lowest_performers(:   , end));
    s8=std(lowest_performers(:   , end));
    x8=linspace(mu8-4*s8,mu8+4*s8,200);
    pdfx8=8/sqrt(2*pi)/s8*exp(-(x8-mu8).^2/(2*s8^2));
    pdfx88 = pdf('Normal', x8, mu8, s8);
    p8 = plot(x8,pdfx88, 'k');
    
    mu9=mean(lower_performers(:   , end));
    s9=std(lower_performers(:   , end));
    x9=linspace(mu9-4*s9,mu9+4*s9,200);
    pdfx9=9/sqrt(2*pi)/s9*exp(-(x9-mu9).^2/(2*s9^2));
    pdfx99 = pdf('Normal', x9, mu9, s9);
    p9 = plot(x9,pdfx99, 'b');
    
    mu10=mean(low_performers(:   , end));
    s10=std(low_performers(:   , end));
    x10=linspace(mu10-4*s10,mu10+4*s10,200);
    pdfx10=10/sqrt(2*pi)/s10*exp(-(x10-mu10).^2/(2*s10^2));
    pdfx1010 = pdf('Normal', x10, mu10, s10);
    p10 = plot(x10,pdfx1010, 'r');
    
    mu11=mean(average_performers(:   , end));
    s11=std(average_performers(:   , end));
    x11=linspace(mu11-4*s11,mu11+4*s11,200);
    pdfx11=11/sqrt(2*pi)/s11*exp(-(x11-mu11).^2/(2*s11^2));
    pdfx1111 = pdf('Normal', x11, mu11, s11);
    p11 = plot(x11,pdfx1111, 'g');
    
    mu12=mean(high_performers(:   , end));
    s12=std(high_performers(:   , end));
    x12=linspace(mu12-4*s12,mu12+4*s12,200);
    pdfx12=12/sqrt(2*pi)/s12*exp(-(x12-mu12).^2/(2*s12^2));
    pdfx1212 = pdf('Normal', x12, mu12, s12);
    p12 = plot(x12,pdfx1212, 'm');
    
    mu13=mean(higher_performers(:   , end));
    s13=std(higher_performers(:   , end));
    x13=linspace(mu13-4*s13,mu13+4*s13,200);
    pdfx13=13/sqrt(2*pi)/s13*exp(-(x13-mu13).^2/(2*s13^2));
    pdfx1313 = pdf('Normal', x13, mu13, s13);
    p13 = plot(x13,pdfx1313, 'c');
    
    mu14=mean(JuniorGpas);
    s14=std(JuniorGpas);
    x14=linspace(mu14-4*s14,mu14+4*s14,200);
    pdfx14=14/sqrt(2*pi)/s14*exp(-(x14-mu14).^2/(2*s14^2));
    pdfx1414 = pdf('Normal', x14, mu14, s14);
    p14 = plot(x14,pdfx1414, 'k');
    
    legend([p2, p3, p4, p5, p6, p7],['Lower N=' num2str(sum(lower_performers_mask))],['Low N=' num2str(sum(low_performers_mask))],...
        ['Average N=' num2str(sum(average_performers_mask))],['High N=' num2str(sum(high_performers_mask))],...
        ['Higher N=' num2str(sum(higher_performers_mask))], 'Overall' )
    
    xlim([1 4.5])
    ylim([0 3.3])
    grid minor
    title('PRE-KI: PDFs of group performances, dashed is prior to 3rd year, solid is 3rd year')
end
%% New plotting things

%figure(4); clf;
%subplot(2,1,1)
subplot(4, 2, [5]); hold on
hold on
% h11 = histogram(lowest_performers(:  , end-1) , 'BinWidth' , .1);
[h11n, h11edges] = histcounts(lowest_performers(:  , end-1) , 'BinWidth' , .1);
% h12 = histogram(lower_performers(:   , end-1) , 'BinWidth' , .1);
[h12n, h12edges] = histcounts(lower_performers(:  , end-1) , 'BinWidth' , .1);
% h13 = histogram(low_performers(:     , end-1) , 'BinWidth' , .1);
[h13n, h13edges] = histcounts(low_performers(:  , end-1) , 'BinWidth' , .1);
% h14 = histogram(average_performers(: , end-1) , 'BinWidth' , .1);
[h14n, h14edges] = histcounts(average_performers(:  , end-1) , 'BinWidth' , .1);
% h15 = histogram(high_performers(:    , end-1) , 'BinWidth' , .1);
[h15n, h15edges] = histcounts(high_performers(:  , end-1) , 'BinWidth' , .1);
% h16 = histogram(higher_performers(:  , end-1) , 'BinWidth' , .1);
[h16n, h16edges] = histcounts(higher_performers(:  , end-1) , 'BinWidth' , .1);


% edgemax = max([h11edges, h12edges, h13edges, h14edges, h15edges, h16edges]);
% edgemin = min([h11edges, h12edges, h13edges, h14edges, h15edges, h16edges]);
% edgespacing = mean(diff(h11edges));
edgemax = 4; edgemin = 0; edgespacing = .1;
edgesnew = edgemin:edgespacing:edgemax;
newMat = zeros(6, length(edgesnew));
tol = .002;
for idxx = 1:length(h11edges)-1
    fidx = find(abs(edgesnew - h11edges(idxx))<tol);
    if fidx
        newMat(1, fidx) = h11n(idxx);
    end
end
for idxx = 1:length(h12edges)-1
    fidx = find(abs(edgesnew - h12edges(idxx))<tol);
    if fidx
        newMat(2, fidx) = h12n(idxx);
    end
end
for idxx = 1:length(h13edges)-1
    fidx = find(abs(edgesnew - h13edges(idxx))<tol);
    if fidx
        newMat(3, fidx) = h13n(idxx);
    end
end
for idxx = 1:length(h14edges)-1
    fidx = find(abs(edgesnew - h14edges(idxx))<tol);
    if fidx
        newMat(4, fidx) = h14n(idxx);
    end
end
for idxx = 1:length(h15edges)-1
    fidx = find(abs(edgesnew - h15edges(idxx))<tol);
    if fidx
        newMat(5, fidx) = h15n(idxx);
    end
end
for idxx = 1:length(h16edges)-1
    fidx = find(abs(edgesnew - h16edges(idxx))<tol);
    if fidx
        newMat(6, fidx) = h16n(idxx);
    end
end
bar(edgesnew, newMat', 1, 'stacked')
% xlim([.5 4.1])
ylims = ylim;
ylims(2) = ylims(2) + 5;
ylim([ylims])
%ylim([0 60])
grid minor
title('PRE-KI end of sophomore year')


%subplot(2,1,2)
%hold on
subplot(4, 2, [7]); hold on
% h21 = histogram(lowest_performers(:  , end) , 'BinWidth' , .1);
[h21n, h21edges] = histcounts(lowest_performers(:  , end) , 'BinWidth' , .1);
% h22 = histogram(lower_performers(:   , end) , 'BinWidth' , .1);
[h22n, h22edges] = histcounts(lower_performers(:  , end) , 'BinWidth' , .1);
% h23 = histogram(low_performers(:     , end) , 'BinWidth' , .1);
[h23n, h23edges] = histcounts(low_performers(:  , end) , 'BinWidth' , .1);
% h24 = histogram(average_performers(: , end) , 'BinWidth' , .1);
[h24n, h24edges] = histcounts(average_performers(:  , end) , 'BinWidth' , .1);
% h25 = histogram(high_performers(:    , end) , 'BinWidth' , .1);
[h25n, h25edges] = histcounts(high_performers(:  , end) , 'BinWidth' , .1);
% h26 = histogram(higher_performers(:  , end) , 'BinWidth' , .1);
[h26n, h26edges] = histcounts(higher_performers(:  , end) , 'BinWidth' , .1);

edgemax = 4; edgemin = 0; edgespacing = .1;
edgesnew = edgemin:edgespacing:edgemax;
newMat = zeros(6, length(edgesnew));
tol = .002;
for idxx = 1:length(h21edges)-1
    fidx = find(abs(edgesnew - h21edges(idxx))<tol);
    if fidx
        newMat(1, fidx) = h21n(idxx);
    end
end
for idxx = 1:length(h22edges)-1
    fidx = find(abs(edgesnew - h22edges(idxx))<tol);
    if fidx
        newMat(2, fidx) = h22n(idxx);
    end
end
for idxx = 1:length(h23edges)-1
    fidx = find(abs(edgesnew - h23edges(idxx))<tol);
    if fidx
        newMat(3, fidx) = h23n(idxx);
    end
end
for idxx = 1:length(h24edges)-1
    fidx = find(abs(edgesnew - h24edges(idxx))<tol);
    if fidx
        newMat(4, fidx) = h24n(idxx);
    end
end
for idxx = 1:length(h25edges)-1
    fidx = find(abs(edgesnew - h25edges(idxx))<tol);
    if fidx
        newMat(5, fidx) = h25n(idxx);
    end
end
for idxx = 1:length(h26edges)-1
    fidx = find(abs(edgesnew - h26edges(idxx))<tol);
    if fidx
        newMat(6, fidx) = h26n(idxx);
    end
end
bar(edgesnew, newMat', 1, 'stacked')
% xlim([.5 4.1])
ylim([ylims])
grid minor
title('PRE-KI end of junior year')

% figure
% b = bar(data);
% b.FaceColor = 'flat';
% b.CData(specificColumn,:) = [1 0 0];
% set(gca, 'xTick', x)
% return
%% First do the pre-KI then do post-KI
postKI_Grades = Grades(data.data(:,6)>0,:);
postKI_Gpas = Gpas(data.data(:,6)>0);

% FIXME - don't use nanmean
JuniorGrades = postKI_Grades(:,end-6:end);
JuniorGpas   = nanmean(JuniorGrades,2);
[mj,nj]         = size(JuniorGrades);
pdj = fitdist(JuniorGpas,'Normal');

% FIXME - don't use nanmean
nonJuniorGrades = postKI_Grades(:,end-10:end-7);
nonJuniorGpas   = nanmean(nonJuniorGrades,2);
nonJuniorGpas_post_KI   = nonJuniorGpas;
[m2,n2]         = size(Grades);
pd = fitdist(nonJuniorGpas,'Normal');

lowest_performance_range = pd.mu - (.5*pd.std + 3*pd.std);
lower_performance_range = pd.mu - (.5*pd.std + 2*pd.std);
low_performance_range = pd.mu - (.5*pd.std + 1*pd.std);
average_performance_range = pd.mu - (.5*pd.std + 0*pd.std);
high_performance_range = pd.mu + (.5*pd.std + 0*pd.std);
higher_performance_range = pd.mu + (.5*pd.std + 1*pd.std);
highest_performance_range = pd.mu + (.5*pd.std + 2*pd.std);

% lowester_performers_mask = (nonJuniorGpas < lowest_performance_range);
lowest_performers_mask = (nonJuniorGpas < lower_performance_range);
lower_performers_mask = (nonJuniorGpas >= lower_performance_range).*(nonJuniorGpas < low_performance_range);
low_performers_mask = (nonJuniorGpas >= low_performance_range).*(nonJuniorGpas < average_performance_range);
average_performers_mask = (nonJuniorGpas >= average_performance_range).*(nonJuniorGpas < high_performance_range);
high_performers_mask = (nonJuniorGpas >= high_performance_range).*(nonJuniorGpas < higher_performance_range);
higher_performers_mask = (nonJuniorGpas >= higher_performance_range).*(nonJuniorGpas < highest_performance_range);

% Sanity check
% if sum(lowest_performers_mask.*lower_performers_mask.*low_performers_mask.*average_performers_mask.*high_performers_mask.*higher_performers_mask) ~= 0
%     error('Overlapping groups')
% elseif sum(lowest_performers_mask)+sum(lower_performers_mask)+sum(low_performers_mask)+sum(average_performers_mask)+sum(high_performers_mask)+sum(higher_performers_mask) ~= mj
%     error('Missing students from groups')
% end

% FIXME - might want to add in anon_id to track these easier, but we should
% be able to make this work as is
lowest_performers = [Grades(lowest_performers_mask>0   , :) , nonJuniorGpas(lowest_performers_mask>0)  , JuniorGpas(lowest_performers_mask>0)  ] ;
lower_performers = [Grades(lower_performers_mask>0     , :) , nonJuniorGpas(lower_performers_mask>0)   , JuniorGpas(lower_performers_mask>0)   ] ;
low_performers = [Grades(low_performers_mask>0         , :) , nonJuniorGpas(low_performers_mask>0)     , JuniorGpas(low_performers_mask>0)     ] ;
average_performers = [Grades(average_performers_mask>0 , :) , nonJuniorGpas(average_performers_mask>0) , JuniorGpas(average_performers_mask>0) ] ;
high_performers = [Grades(high_performers_mask>0       , :) , nonJuniorGpas(high_performers_mask>0)    , JuniorGpas(high_performers_mask>0)    ] ;
higher_performers = [Grades(higher_performers_mask>0   , :) , nonJuniorGpas(higher_performers_mask>0)  , JuniorGpas(higher_performers_mask>0)  ] ;
%% Plot the individual group distributions
%figure(3); clf; hold on
if doFakeFit
    subplot(4, 2, [2 4]); hold on
    
    mu1=mean(lowest_performers(:   , end-1));
    s1=std(lowest_performers(:   , end-1));
    x1=linspace(mu1-4*s1,mu1+4*s1,200);
    pdfx1=1/sqrt(2*pi)/s1*exp(-(x1-mu1).^2/(2*s1^2));
    pd1 = makedist('Normal', mu1, s1);
    pdfx11 = pdf(pd1, x1);
    p1 = plot(x1,pdfx11, 'k--*');
    
    mu2=mean(lower_performers(:   , end-1));
    s2=std(lower_performers(:   , end-1));
    x2=linspace(mu2-4*s2,mu2+4*s2,200);
    pdfx2=2/sqrt(2*pi)/s2*exp(-(x2-mu2).^2/(2*s2^2));
    pdfx22 = pdf('Normal', x2, mu2, s2);
    p2 = plot(x2,pdfx22, 'b--');
    
    mu3=mean(low_performers(:   , end-1));
    s3=std(low_performers(:   , end-1));
    x3=linspace(mu3-4*s3,mu3+4*s3,200);
    pdfx3=3/sqrt(2*pi)/s3*exp(-(x3-mu3).^2/(2*s3^2));
    pdfx33 = pdf('Normal', x3, mu3, s3);
    p3 = plot(x3,pdfx33, 'r--');
    
    mu4=mean(average_performers(:   , end-1));
    s4=std(average_performers(:   , end-1));
    x4=linspace(mu4-4*s4,mu4+4*s4,200);
    pdfx4=4/sqrt(2*pi)/s4*exp(-(x4-mu4).^2/(2*s4^2));
    pdfx44 = pdf('Normal', x4, mu4, s4);
    p4 = plot(x4,pdfx44, 'g--');
    
    mu5=mean(high_performers(:   , end-1));
    s5=std(high_performers(:   , end-1));
    x5=linspace(mu5-4*s5,mu5+4*s5,200);
    pdfx5=5/sqrt(2*pi)/s5*exp(-(x5-mu5).^2/(2*s5^2));
    pdfx55 = pdf('Normal', x5, mu5, s5);
    p5 = plot(x5,pdfx55, 'm--');
    
    mu6=mean(higher_performers(:   , end-1));
    s6=std(higher_performers(:   , end-1));
    x6=linspace(mu6-4*s6,mu6+4*s6,200);
    pdfx6=6/sqrt(2*pi)/s6*exp(-(x6-mu6).^2/(2*s6^2));
    pdfx66 = pdf('Normal', x6, mu6, s6);
    p6 = plot(x6,pdfx66, 'c--');
    
    mu7=mean(nonJuniorGpas);
    s7=std(nonJuniorGpas);
    x7=linspace(mu7-4*s7,mu7+4*s7,200);
    pdfx7=7/sqrt(2*pi)/s7*exp(-(x7-mu7).^2/(2*s7^2));
    pdfx77 = pdf('Normal', x7, mu7, s7);
    p7 = plot(x7,pdfx77, 'k--');
    
    mu8=mean(lowest_performers(:   , end));
    s8=std(lowest_performers(:   , end));
    x8=linspace(mu8-4*s8,mu8+4*s8,200);
    pdfx8=8/sqrt(2*pi)/s8*exp(-(x8-mu8).^2/(2*s8^2));
    pdfx88 = pdf('Normal', x8, mu8, s8);
    p8 = plot(x8,pdfx88, 'k');
    
    mu9=mean(lower_performers(:   , end));
    s9=std(lower_performers(:   , end));
    x9=linspace(mu9-4*s9,mu9+4*s9,200);
    pdfx9=9/sqrt(2*pi)/s9*exp(-(x9-mu9).^2/(2*s9^2));
    pdfx99 = pdf('Normal', x9, mu9, s9);
    p9 = plot(x9,pdfx99, 'b');
    
    mu10=mean(low_performers(:   , end));
    s10=std(low_performers(:   , end));
    x10=linspace(mu10-4*s10,mu10+4*s10,200);
    pdfx10=10/sqrt(2*pi)/s10*exp(-(x10-mu10).^2/(2*s10^2));
    pdfx1010 = pdf('Normal', x10, mu10, s10);
    p10 = plot(x10,pdfx1010, 'r');
    
    mu11=mean(average_performers(:   , end));
    s11=std(average_performers(:   , end));
    x11=linspace(mu11-4*s11,mu11+4*s11,200);
    pdfx11=11/sqrt(2*pi)/s11*exp(-(x11-mu11).^2/(2*s11^2));
    pdfx1111 = pdf('Normal', x11, mu11, s11);
    p11 = plot(x11,pdfx1111, 'g');
    
    mu12=mean(high_performers(:   , end));
    s12=std(high_performers(:   , end));
    x12=linspace(mu12-4*s12,mu12+4*s12,200);
    pdfx12=12/sqrt(2*pi)/s12*exp(-(x12-mu12).^2/(2*s12^2));
    pdfx1212 = pdf('Normal', x12, mu12, s12);
    p12 = plot(x12,pdfx1212, 'm');
    
    mu13=mean(higher_performers(:   , end));
    s13=std(higher_performers(:   , end));
    x13=linspace(mu13-4*s13,mu13+4*s13,200);
    pdfx13=13/sqrt(2*pi)/s13*exp(-(x13-mu13).^2/(2*s13^2));
    pdfx1313 = pdf('Normal', x13, mu13, s13);
    p13 = plot(x13,pdfx1313, 'c');
    
    mu14=mean(JuniorGpas);
    s14=std(JuniorGpas);
    x14=linspace(mu14-4*s14,mu14+4*s14,200);
    pdfx14=14/sqrt(2*pi)/s14*exp(-(x14-mu14).^2/(2*s14^2));
    pdfx1414 = pdf('Normal', x14, mu14, s14);
    p14 = plot(x14,pdfx1414, 'k');
    
    legend([p2, p3, p4, p5, p6, p7],['Lower N=' num2str(sum(lower_performers_mask))],['Low N=' num2str(sum(low_performers_mask))],...
        ['Average N=' num2str(sum(average_performers_mask))],['High N=' num2str(sum(high_performers_mask))],...
        ['Higher N=' num2str(sum(higher_performers_mask))], 'Overall' )
    
    xlim([1 4.5])
    ylim([0 3.3])
    grid minor
    title('POST-KI: PDFs of group performances, dashed is prior to 3rd year, solid is 3rd year')
end
%% New plotting things
% generate random data set between 1 and 20

%figure(6); clf;
%subplot(2,1,1)
subplot(4, 2, [6]); hold on
hold on
% h31 = histogram(lowest_performers(:  , end-1) , 'BinWidth' , .1);
[h31n, h31edges] = histcounts(lowest_performers(:  , end-1) , 'BinWidth' , .1);
% h32 = histogram(lower_performers(:   , end-1) , 'BinWidth' , .1);
[h32n, h32edges] = histcounts(lower_performers(:  , end-1) , 'BinWidth' , .1);
% h33 = histogram(low_performers(:     , end-1) , 'BinWidth' , .1);
[h33n, h33edges] = histcounts(low_performers(:  , end-1) , 'BinWidth' , .1);
% h34 = histogram(average_performers(: , end-1) , 'BinWidth' , .1);
[h34n, h34edges] = histcounts(average_performers(:  , end-1) , 'BinWidth' , .1);
% h35 = histogram(high_performers(:    , end-1) , 'BinWidth' , .1);
[h35n, h35edges] = histcounts(high_performers(:  , end-1) , 'BinWidth' , .1);
% h36 = histogram(higher_performers(:  , end-1) , 'BinWidth' , .1);
[h36n, h36edges] = histcounts(higher_performers(:  , end-1) , 'BinWidth' , .1);

edgemax = 4; edgemin = 0; edgespacing = .1;
edgesnew = edgemin:edgespacing:edgemax;
newMat = zeros(6, length(edgesnew));
tol = .002;
for idxx = 1:length(h31edges)-1
    fidx = find(abs(edgesnew - h31edges(idxx))<tol);
    if fidx
        newMat(1, fidx) = h31n(idxx);
    end
end
for idxx = 1:length(h32edges)-1
    fidx = find(abs(edgesnew - h32edges(idxx))<tol);
    if fidx
        newMat(2, fidx) = h32n(idxx);
    end
end
for idxx = 1:length(h33edges)-1
    fidx = find(abs(edgesnew - h33edges(idxx))<tol);
    if fidx
        newMat(3, fidx) = h33n(idxx);
    end
end
for idxx = 1:length(h34edges)-1
    fidx = find(abs(edgesnew - h34edges(idxx))<tol);
    if fidx
        newMat(4, fidx) = h34n(idxx);
    end
end
for idxx = 1:length(h35edges)-1
    fidx = find(abs(edgesnew - h35edges(idxx))<tol);
    if fidx
        newMat(5, fidx) = h35n(idxx);
    end
end
for idxx = 1:length(h36edges)-1
    fidx = find(abs(edgesnew - h36edges(idxx))<tol);
    if fidx
        newMat(6, fidx) = h36n(idxx);
    end
end
bar(edgesnew, newMat', 1, 'stacked')
% xlim([.5 4.1])
ylims = ylim;
ylims(2) = ylims(2) + 1;
ylim([ylims])
%ylim([0 60])
grid minor
title('POST-KI end of sophomore year')

%subplot(2,1,2)
%hold on;
subplot(4, 2, [8]); hold on
% h41 = histogram(lowest_performers(:  , end) , 'BinWidth' , .1);
[h41n, h41edges] = histcounts(lowest_performers(:  , end) , 'BinWidth' , .1);
% h42 = histogram(lower_performers(:   , end) , 'BinWidth' , .1);
[h42n, h42edges] = histcounts(lower_performers(:  , end) , 'BinWidth' , .1);
% h43 = histogram(low_performers(:     , end) , 'BinWidth' , .1);
[h43n, h43edges] = histcounts(low_performers(:  , end) , 'BinWidth' , .1);
% h44 = histogram(average_performers(: , end) , 'BinWidth' , .1);
[h44n, h44edges] = histcounts(average_performers(:  , end) , 'BinWidth' , .1);
% h45 = histogram(high_performers(:    , end) , 'BinWidth' , .1);
[h45n, h45edges] = histcounts(high_performers(:  , end) , 'BinWidth' , .1);
% h46 = histogram(higher_performers(:  , end) , 'BinWidth' , .1);
[h46n, h46edges] = histcounts(higher_performers(:  , end) , 'BinWidth' , .1);

edgemax = 4; edgemin = 0; edgespacing = .1;
edgesnew = edgemin:edgespacing:edgemax;
newMat = zeros(6, length(edgesnew));
tol = .002;
for idxx = 1:length(h41edges)-1
    fidx = find(abs(edgesnew - h41edges(idxx))<tol);
    if fidx
        newMat(1, fidx) = h41n(idxx);
    end
end
for idxx = 1:length(h42edges)-1
    fidx = find(abs(edgesnew - h42edges(idxx))<tol);
    if fidx
        newMat(2, fidx) = h42n(idxx);
    end
end
for idxx = 1:length(h43edges)-1
    fidx = find(abs(edgesnew - h43edges(idxx))<tol);
    if fidx
        newMat(3, fidx) = h43n(idxx);
    end
end
for idxx = 1:length(h44edges)-1
    fidx = find(abs(edgesnew - h44edges(idxx))<tol);
    if fidx
        newMat(4, fidx) = h44n(idxx);
    end
end
for idxx = 1:length(h45edges)-1
    fidx = find(abs(edgesnew - h45edges(idxx))<tol);
    if fidx
        newMat(5, fidx) = h45n(idxx);
    end
end
for idxx = 1:length(h46edges)-1
    fidx = find(abs(edgesnew - h46edges(idxx))<tol);
    if fidx
        newMat(6, fidx) = h46n(idxx);
    end
end
bar(edgesnew, newMat', 1, 'stacked')
% xlim([.5 4.1])
ylim([ylims])
grid minor
title('POST-KI end of junior year')


%% JS Divergence calculations
display('Comparing distributions between sophomore and junior years')
[nh11, nh21, nbins] = makeSimilarDists(h11n, h11edges, h21n, h21edges, .1);
[nh12, nh22, nbins] = makeSimilarDists(h12n, h12edges, h22n, h22edges, .1);
[nh13, nh23, nbins] = makeSimilarDists(h13n, h13edges, h23n, h23edges, .1);
[nh14, nh24, nbins] = makeSimilarDists(h14n, h14edges, h24n, h24edges, .1);
[nh15, nh25, nbins] = makeSimilarDists(h15n, h15edges, h25n, h25edges, .1);
[nh16, nh26, nbins] = makeSimilarDists(h16n, h16edges, h26n, h26edges, .1);

js11 = jsdiv(nh11/sum(nh11), nh21/sum(nh21));
js12 = jsdiv(nh12/sum(nh12), nh22/sum(nh22));
js13 = jsdiv(nh13/sum(nh13), nh23/sum(nh23));
js14 = jsdiv(nh14/sum(nh14), nh24/sum(nh24));
js15 = jsdiv(nh15/sum(nh15), nh25/sum(nh25));
js16 = jsdiv(nh16/sum(nh16), nh26/sum(nh26));

display([' Pre KI JS Divs - lowest:' num2str(js11) ...
    ', lower:' num2str(js12) ...
    ', low:' num2str(js13) ...
    ', average:' num2str(js14) ...
    ', high:' num2str(js15) ...
    ', higher:' num2str(js16)])

[nh31, nh41, nbins] = makeSimilarDists(h31n, h31edges, h41n, h41edges, .1);
[nh32, nh42, nbins] = makeSimilarDists(h32n, h32edges, h42n, h42edges, .1);
[nh33, nh43, nbins] = makeSimilarDists(h33n, h33edges, h43n, h43edges, .1);
[nh34, nh44, nbins] = makeSimilarDists(h34n, h34edges, h44n, h44edges, .1);
[nh35, nh45, nbins] = makeSimilarDists(h35n, h35edges, h45n, h45edges, .1);
[nh36, nh46, nbins] = makeSimilarDists(h36n, h36edges, h46n, h46edges, .1);

js21 = jsdiv(nh31/sum(nh31), nh41/sum(nh41));
js22 = jsdiv(nh32/sum(nh32), nh42/sum(nh42));
js23 = jsdiv(nh33/sum(nh33), nh43/sum(nh43));
js24 = jsdiv(nh34/sum(nh34), nh44/sum(nh44));
js25 = jsdiv(nh35/sum(nh35), nh45/sum(nh45));
js26 = jsdiv(nh36/sum(nh36), nh46/sum(nh46));
display(['Post KI JS Divs - lowest:' num2str(js21) ...
    ', lower:' num2str(js22) ...
    ', low:' num2str(js23) ...
    ', average:' num2str(js24) ...
    ', high:' num2str(js25) ...
    ', higher:' num2str(js26)])


display(' ')
display('Comparing distributions between Pre and Post KI')
[nh51, nh61, nbins] = makeSimilarDists(h11n, h11edges, h31n, h31edges, .1);
[nh52, nh62, nbins] = makeSimilarDists(h12n, h12edges, h32n, h32edges, .1);
[nh53, nh63, nbins] = makeSimilarDists(h13n, h13edges, h33n, h33edges, .1);
[nh54, nh64, nbins] = makeSimilarDists(h14n, h14edges, h34n, h34edges, .1);
[nh55, nh65, nbins] = makeSimilarDists(h15n, h15edges, h35n, h35edges, .1);
[nh56, nh66, nbins] = makeSimilarDists(h16n, h16edges, h36n, h36edges, .1);

js31 = jsdiv(nh51/sum(nh51), nh61/sum(nh61));
js32 = jsdiv(nh52/sum(nh52), nh62/sum(nh62));
js33 = jsdiv(nh53/sum(nh53), nh63/sum(nh63));
js34 = jsdiv(nh54/sum(nh54), nh64/sum(nh64));
js35 = jsdiv(nh55/sum(nh55), nh65/sum(nh65));
js36 = jsdiv(nh56/sum(nh56), nh66/sum(nh66));

display(['Sophomore JS Divs - lowest:' num2str(js31) ...
    ', lower:' num2str(js32) ...
    ', low:' num2str(js33) ...
    ', average:' num2str(js34) ...
    ', high:' num2str(js35) ...
    ', higher:' num2str(js36)])

[nh71, nh81, nbins] = makeSimilarDists(h21n, h21edges, h41n, h41edges, .1);
[nh72, nh82, nbins] = makeSimilarDists(h22n, h22edges, h42n, h42edges, .1);
[nh73, nh83, nbins] = makeSimilarDists(h23n, h23edges, h43n, h43edges, .1);
[nh74, nh84, nbins] = makeSimilarDists(h24n, h24edges, h44n, h44edges, .1);
[nh75, nh85, nbins] = makeSimilarDists(h25n, h25edges, h45n, h45edges, .1);
[nh76, nh86, nbins] = makeSimilarDists(h26n, h26edges, h46n, h46edges, .1);

js41 = jsdiv(nh71/sum(nh71), nh81/sum(nh81));
js42 = jsdiv(nh72/sum(nh72), nh82/sum(nh82));
js43 = jsdiv(nh73/sum(nh73), nh83/sum(nh83));
js44 = jsdiv(nh74/sum(nh74), nh84/sum(nh84));
js45 = jsdiv(nh75/sum(nh75), nh85/sum(nh85));
js46 = jsdiv(nh76/sum(nh76), nh86/sum(nh86));
display(['   Junior JS Divs - lowest:' num2str(js41) ...
    ', lower:' num2str(js42) ...
    ', low:' num2str(js43) ...
    ', average:' num2str(js44) ...
    ', high:' num2str(js45) ...
    ', higher:' num2str(js46)])
