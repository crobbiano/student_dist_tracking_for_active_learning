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
% addpath('..\..\')
%% Set Up
load('..\correlational_studies_data\data\grades.mat')

FileName = 'grades.csv';
FilePath = '..\correlational_studies_data\data';

File = [FilePath '\' FileName];

data = importdata(File, ',');

Classes = data.textdata;
Classes = Classes(numNonClasses+1:end-numMathClasses);

%%
Grades         = data.data(:,numNonClasses+1:end-numMathClasses);
Gpas           = data.data(:,numNonClasses);
gradesGpas     = mean(Grades,2);
[m,n]          = size(Grades);

JuniorGrades = data.data(:,numNonClasses+1+5:end-numMathClasses);
JuniorGpas   = mean(JuniorGrades,2);
[mj,nj]         = size(Grades);
pdj = fitdist(JuniorGpas,'Normal');

nonJuniorGrades = data.data(:,numNonClasses+1:end-numMathClasses-6);
nonJuniorGpas   = mean(nonJuniorGrades,2);
[m2,n2]         = size(Grades);
pd = fitdist(nonJuniorGpas,'Normal');

lowest_performance_range = pd.mu - (.5*pd.std + 3*pd.std);
lower_performance_range = pd.mu - (.5*pd.std + 2*pd.std);
low_performance_range = pd.mu - (.5*pd.std + 1*pd.std);
average_performance_range = pd.mu - (.5*pd.std + 0*pd.std);
high_performance_range = pd.mu + (.5*pd.std + 0*pd.std);
higher_performance_range = pd.mu + (.5*pd.std + 1*pd.std);
highest_performance_range = pd.mu + (.5*pd.std + 2*pd.std);

lowest_performers_mask = (nonJuniorGpas >= lowest_performance_range).*(nonJuniorGpas < lower_performance_range);
lower_performers_mask = (nonJuniorGpas >= lower_performance_range).*(nonJuniorGpas < low_performance_range);
low_performers_mask = (nonJuniorGpas >= low_performance_range).*(nonJuniorGpas < average_performance_range);
average_performers_mask = (nonJuniorGpas >= average_performance_range).*(nonJuniorGpas < high_performance_range);
high_performers_mask = (nonJuniorGpas >= high_performance_range).*(nonJuniorGpas < higher_performance_range);
higher_performers_mask = (nonJuniorGpas >= higher_performance_range).*(nonJuniorGpas < highest_performance_range);

% Sanity check
if sum(lowest_performers_mask.*lower_performers_mask.*low_performers_mask.*average_performers_mask.*high_performers_mask.*higher_performers_mask) ~= 0
    error('Overlapping groups')
elseif sum(lowest_performers_mask)+sum(lower_performers_mask)+sum(low_performers_mask)+sum(average_performers_mask)+sum(high_performers_mask)+sum(higher_performers_mask) ~= m
    error('Missing students from groups')
end

% FIXME - might want to add in anon_id to track these easier, but we should
% be able to make this work as is
lowest_performers = [Grades(lowest_performers_mask>0   , :) , nonJuniorGpas(lowest_performers_mask>0)  , JuniorGpas(lowest_performers_mask>0)  ] ;
lower_performers = [Grades(lower_performers_mask>0     , :) , nonJuniorGpas(lower_performers_mask>0)   , JuniorGpas(lower_performers_mask>0)   ] ;
low_performers = [Grades(low_performers_mask>0         , :) , nonJuniorGpas(low_performers_mask>0)     , JuniorGpas(low_performers_mask>0)     ] ;
average_performers = [Grades(average_performers_mask>0 , :) , nonJuniorGpas(average_performers_mask>0) , JuniorGpas(average_performers_mask>0) ] ;
high_performers = [Grades(high_performers_mask>0       , :) , nonJuniorGpas(high_performers_mask>0)    , JuniorGpas(high_performers_mask>0)    ] ;
higher_performers = [Grades(higher_performers_mask>0   , :) , nonJuniorGpas(higher_performers_mask>0)  , JuniorGpas(higher_performers_mask>0)  ] ;

%% Plot the individual group distributions
figure(2); clf; hold on

mu1=mean(lowest_performers(:   , end-1));
s1=std(lowest_performers(:   , end-1));
x1=linspace(mu1-4*s1,mu1+4*s1,200);
pdfx1=1/sqrt(2*pi)/s1*exp(-(x1-mu1).^2/(2*s1^2));
pd1 = makedist('Normal', mu1, s1);
pdfx11 = pdf(pd1, x1);
plot(x1,pdfx11, 'k--*');

mu2=mean(lower_performers(:   , end-1));
s2=std(lower_performers(:   , end-1));
x2=linspace(mu2-4*s2,mu2+4*s2,200);
pdfx2=2/sqrt(2*pi)/s2*exp(-(x2-mu2).^2/(2*s2^2));
pdfx22 = pdf('Normal', x2, mu2, s2);
plot(x2,pdfx22, 'b--');

mu3=mean(low_performers(:   , end-1));
s3=std(low_performers(:   , end-1));
x3=linspace(mu3-4*s3,mu3+4*s3,200);
pdfx3=3/sqrt(2*pi)/s3*exp(-(x3-mu3).^2/(2*s3^2));
pdfx33 = pdf('Normal', x3, mu3, s3);
plot(x3,pdfx33, 'r--');

mu4=mean(average_performers(:   , end-1));
s4=std(average_performers(:   , end-1));
x4=linspace(mu4-4*s4,mu4+4*s4,200);
pdfx4=4/sqrt(2*pi)/s4*exp(-(x4-mu4).^2/(2*s4^2));
pdfx44 = pdf('Normal', x4, mu4, s4);
plot(x4,pdfx44, 'g--');

mu5=mean(high_performers(:   , end-1));
s5=std(high_performers(:   , end-1));
x5=linspace(mu5-4*s5,mu5+4*s5,200);
pdfx5=5/sqrt(2*pi)/s5*exp(-(x5-mu5).^2/(2*s5^2));
pdfx55 = pdf('Normal', x5, mu5, s5);
plot(x5,pdfx55, 'm--');

mu6=mean(higher_performers(:   , end-1));
s6=std(higher_performers(:   , end-1));
x6=linspace(mu6-4*s6,mu6+4*s6,200);
pdfx6=6/sqrt(2*pi)/s6*exp(-(x6-mu6).^2/(2*s6^2));
pdfx66 = pdf('Normal', x6, mu6, s6);
plot(x6,pdfx66, 'c--');

% mu7=mean(highest_performers(:   , end-1));
% s7=std(highest_performers(:   , end-1));
% x7=linspace(mu7-4*s7,mu7+4*s7,200);
% pdfx7=7/sqrt(2*pi)/s7*exp(-(x7-mu7).^2/(2*s7^2));
% plot(x7,pdfx7, 'k--');

mu8=mean(lowest_performers(:   , end));
s8=std(lowest_performers(:   , end));
x8=linspace(mu8-4*s8,mu8+4*s8,200);
pdfx8=8/sqrt(2*pi)/s8*exp(-(x8-mu8).^2/(2*s8^2));
pdfx88 = pdf('Normal', x8, mu8, s8);
plot(x8,pdfx88, 'k');

mu9=mean(lower_performers(:   , end));
s9=std(lower_performers(:   , end));
x9=linspace(mu9-4*s9,mu9+4*s9,200);
pdfx9=9/sqrt(2*pi)/s9*exp(-(x9-mu9).^2/(2*s9^2));
pdfx99 = pdf('Normal', x9, mu9, s9);
plot(x9,pdfx99, 'b');

mu10=mean(low_performers(:   , end));
s10=std(low_performers(:   , end));
x10=linspace(mu10-4*s10,mu10+4*s10,200);
pdfx10=10/sqrt(2*pi)/s10*exp(-(x10-mu10).^2/(2*s10^2));
pdfx1010 = pdf('Normal', x10, mu10, s10);
plot(x10,pdfx1010, 'r');

mu11=mean(average_performers(:   , end));
s11=std(average_performers(:   , end));
x11=linspace(mu11-4*s11,mu11+4*s11,200);
pdfx11=11/sqrt(2*pi)/s11*exp(-(x11-mu11).^2/(2*s11^2));
pdfx1111 = pdf('Normal', x11, mu11, s11);
plot(x11,pdfx1111, 'g');

mu12=mean(high_performers(:   , end));
s12=std(high_performers(:   , end));
x12=linspace(mu12-4*s12,mu12+4*s12,200);
pdfx12=12/sqrt(2*pi)/s12*exp(-(x12-mu12).^2/(2*s12^2));
pdfx1212 = pdf('Normal', x12, mu12, s12);
plot(x12,pdfx1212, 'm');

mu13=mean(higher_performers(:   , end));
s13=std(higher_performers(:   , end));
x13=linspace(mu13-4*s13,mu13+4*s13,200);
pdfx13=13/sqrt(2*pi)/s13*exp(-(x13-mu13).^2/(2*s13^2));
pdfx1313 = pdf('Normal', x13, mu13, s13);
plot(x13,pdfx1313, 'c');

% mu14=mean(highest_performers(:   , end));
% s14=std(highest_performers(:   , end));
% x14=linspace(mu14-4*s14,mu14+4*s14,200);
% pdfx14=14/sqrt(2*pi)/s14*exp(-(x14-mu14).^2/(2*s14^2));
% plot(x14,pdfx14, 'k--');
legend(['Lowest N=' num2str(sum(lowest_performers_mask))],['Lower N=' num2str(sum(lower_performers_mask))],['Low N=' num2str(sum(low_performers_mask))],...
    ['Average N=' num2str(sum(average_performers_mask))],['High N=' num2str(sum(high_performers_mask))],['Higher N=' num2str(sum(higher_performers_mask))])
xlim([1.5 4.5])
grid minor
title('PDFs of group performances, dashed is prior to 3rd year, solid is 3rd year')
%% Old plotting things
% h1  = histfit(lower_performers(:   , end-1) , 5);
% h2  = histfit(low_performers(:     , end-1) , 5);
% h3  = histfit(average_performers(: , end-1) , 4);
% h4  = histfit(high_performers(:    , end-1) , 4);
% h5  = histfit(higher_performers(:  , end-1) , 3);
% h6  = histfit(lower_performers(:   , end) , 5);
% h7  = histfit(low_performers(:     , end) , 5);
% h8  = histfit(average_performers(: , end) , 4);
% h9  = histfit(high_performers(:    , end) , 4);
% h10 = histfit(higher_performers(:  , end) , 3);
% 
% delete(h1(1)); h1(2).Color = 'k';  h1(2).LineStyle = '--';  h1(2).YData = h1(2).YData/max(h1(2).YData)
% delete(h2(1)); h2(2).Color = 'r';  h2(2).LineStyle = '--'
% delete(h3(1)); h3(2).Color = 'c';  h3(2).LineStyle = '--'
% delete(h4(1)); h4(2).Color = 'g';  h4(2).LineStyle = '--'
% delete(h5(1)); h5(2).Color = 'm';  h5(2).LineStyle = '--'
% delete(h6(1)); h6(2).Color = 'k';
% delete(h7(1)); h7(2).Color = 'r';
% delete(h8(1)); h8(2).Color = 'c';
% delete(h9(1)); h9(2).Color = 'g';
% delete(h10(1)); h10(2).Color = 'm';
% xlim([1.5 4.5])
% grid minor