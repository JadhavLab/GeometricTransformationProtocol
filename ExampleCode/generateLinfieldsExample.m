%Example script for generating linearized place fields from position and
%spikes

addpath(''); %Add code directory (See https://github.com/JadhavLab/GeometricTransformation)
dataDir = 'D:\DATA_ANALYSIS\STAR_Protocol_Code\'; %Directory to example data and code provided
cd(dataDir);

%example day
day = 1;

%example epoch
epoch = 1;
load('EXpos01.mat')

%Get coodinates for a wtrack.
coords = getcoord_wtrack_example(pos{day}{epoch}.data);
%Click locations in the following order:
%   
%   1    2    3
%   |    |    |
%   |    |    |
%   |    |    |
%   |    |    |
%   4----5----6

task{day}{epoch}.linearcoord = coords; %create a task file with linear coordinates of W-Track nodes
save([dataDir,'EX','task01'],'task');

%Linearize position
[linpos{day}{epoch}.statematrix,linpos{day}{epoch}.segmenttable, linpos{day}{epoch}.trajwells, ...
    linpos{day}{epoch}.wellSegmentInfo, linpos{day}{epoch}.segmentInfo]...
    = wb_linearizeposition_example(pos, task, [day epoch]);
linpos{day}{epoch}.statematrixfields = ['time startwell endwell segment segheaddir velocity lineardist'];
linpos{day}{epoch}.segmenttablefields = ['trajnum segnum segmentID'];
save([dataDir,'EX','linpos01'],'linpos');

%Extract the trajectories 
[trajbound, rewarded, trajtime, wellstend] ...
    = trackperformance_example([day epoch], [], linpos, pos, ...
    [2 1 3]);
trajinfo{day}{epoch}.trajbound = trajbound;
trajinfo{day}{epoch}.rewarded = rewarded;
trajinfo{day}{epoch}.trajtime = trajtime;
trajinfo{day}{epoch}.wellstend = wellstend;
trajinfo{day}{epoch}.descript...
    = 'trajbound: 1=inbound, 0=outbound; rewarded: 1=correct, 0=error; trajtime: start and end time of traj/lap; wellstend = start and end well of traj/lap';
load('EXspikes01.mat')
save([dataDir,'EX','trajinfo01'],'trajinfo');

%Example cell recorded on tetrode 1
tetrode = 1;
cell = 1;
exampleTraj = 3;
index = [day epoch tetrode cell];
includeperiods = [trajinfo{day}{epoch}.trajtime(7,:)]; %example trajectory time
linfields = calclinfields_example(index, includeperiods, linpos, spikes); %Example trial linearized field (traj3)
figure; 
plot(linfields.trajdata{exampleTraj}(:,5)); %column 5 is occ. norm. firing rate
figure; 
imagesc(linfields.trajdata{exampleTraj}(:,5)) %EXAMPLE - One cell, one trial, one trajectory
ylabel('Position')
xlabel('Cell #')
xticks(1)
title('Example single trial field')

