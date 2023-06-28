function [trajbound, rewarded, trajtime, wellstend] = trackperformance_example(index, excludetimes, linpos, pos, correctorder, varargin)
% SJ: Dec 2014
% In my data, I have already computed turn-around errors with 0.75 distance length of travel arm
% Here, I am adding an output: wellstend, which gives start and end well of trajectory
% Can also include a condition that instead of using a mindist travled on arm (eg. % of length of arm),
% can also use a condition that says "travel within x dist of intersection", but tis will reuire Position input
% Skip for now. mindist input should be based on segmentlength for center arm. 


% function [trajbound, rewarded, time] = trackperformance(index, excludetimes, linpos, correctorder, mindist)
%%assumes animal is performing an alternation task, with a middle well and
%2, alternating outside wells, though it does not have to be a W-shaped
%track.  does not track trajectories, as animal may wander between time
%leaves one well and enters another. 
%
% LINPOS is the output of linearizeposition for one day
% INDEX [day epoch]
% EXCLUDETIMES times to exclude from analysis.  This excludes trials whose
% start time falls in an exclude period
% CORRECTORDER is a vector of 3 well numbers, the middle number corresponds
% to the middle well, the first and third ones correspond to outside wells,
% for instance [2 1 3] would mean the middle well is well 1, the outside
% wells are 2 and 3.  these well numbers are determined in createtaskstruct
% MINDISTARM is the minimum length of an excursion from a well for the excursion
% to be counted as a trajectory 
% varargin:       THRSDISTINT is the threshold distance to an intersection for the excursion
% to be counted as a trajectory 
%
% TRAJBOUND defined by where animal heading, 
%   -1 if animal heading to a well not included in correct order
%   1 if headed to middle well = INBOUND
%   0 if headed to outside well = OUTBOUND
% REWARDED, 0 or 1, 1 if rewarded at well heading to, -1 if not
% TIME time of traj start

 
mindistarm=[]; thrsdistint=[]; 
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'mindistarm'
            mindistarm = varargin{option+1};
        case 'thrsdistint'
            thrsdistint = varargin{option+1};  
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

thrsdistint = 5;
day = index(1);
epoch = index(2);
          
midwell=correctorder(2); %middle well
owell1 =correctorder(1); %1st outside well
owell2 =correctorder(3); %2nd outside well

wellExitEnter = linpos{day}{epoch}.statematrix.wellExitEnter;
startwell = wellExitEnter(:,1); %well leaving from
endwell = wellExitEnter(:,2); %well heading to
wellchange = [1; (find((diff(startwell) | diff(endwell))) + 1)];

% For using thrsdistint
if size(pos{day}{epoch}.data,2)>5
    posn = pos{day}{epoch}.data(:,6:7); %  Normal x,y is in 2:3. Smoothened position in 6:7.
else
    posn = pos{day}{epoch}.data(:,2:3);
end
% Get Well and Intersection positions
wellpos = linpos{day}{epoch}.wellSegmentInfo.wellCoord;
intpos(1,:) = linpos{day}{epoch}.segmentInfo.segmentCoords(1,3:4); % Center arm
intpos(2,:) = linpos{day}{epoch}.segmentInfo.segmentCoords(2,3:4); % Left arm - our view
intpos(3,:) = linpos{day}{epoch}.segmentInfo.segmentCoords(4,3:4); % right arm - our view

% we need to get rid of cases where the trajectory indicates that the animal
% has started and ended at the same well but where the excursion distance is
% less than mindistarm, or thrsdistint is not met. 

% Default is to use Threshold distance from intersection. If mindistarm is given, use that insted
if ~(isempty(thrsdistint) || isempty(mindistarm))
    samewell = find(wellExitEnter(wellchange,1) == wellExitEnter(wellchange,2));
    invalid = [];
    for s = 1:length(samewell)
        % check to find the maximum excursion distance for this set of points
        if (samewell(s) < length(wellchange))
            tind = wellchange(samewell(s)):wellchange(samewell(s)+1);
        else
            % this is the last trajectory, so we go to the end of the data
            tind = wellchange(samewell(s)):size(wellExitEnter,1);
        end
            
        if ~isempty(thrsdistint) 
             % USING THRSDISTINT, threshold distance from intersection which is DEFAULT
            currpos = posn(tind,:); % all relevant positions
            currintpos = intpos(wellExitEnter(wellchange((samewell(s))),1),:);
            alldist = dist(currpos,repmat(currintpos,length(currpos),1));
            if (min(alldist) > thrsdistint) % Never gets close enough to intersection
                invalid(s) = 1;
            else
                invalid(s) = 0;
            end  
        else 
             % If thrsdistint is not available, USING MINDISTARM - excursion from well
             % get the distances to the start well
            d = linpos{day}{epoch}.statematrix.linearDistanceToWells(tind, ...
                wellExitEnter(wellchange((samewell(s))),1));
            if (max(d) < mindistarm)
                invalid(s) = 1; % Excursions dont exceed mindist
            else
                invalid(s) = 0;
            end   
        end
        
        
    end
    % get rid of the invalid trajectories
    wellchange(samewell(find(invalid))) = -1;
    wellchange = wellchange(find(wellchange ~= -1));
end

wellSequence = endwell(wellchange); %seq of end wells visited
strt = [1;wellSequence(1:end-1)];
wellstend = [strt wellSequence]; % Well St-end from wellSequence, with the assumption that first well is always 1=midwell. Also see below

correct = 0; %tracks number correct traj
rewarded = -5*ones(length(wellSequence),1); %tracks if each traj correct (1) or incorrect(0)
trajbound = rewarded;
total = 0;%tracks total number traj counted
lastoutwell = []; %tracks last outward bound well
time = linpos{day}{epoch}.statematrix.time(wellchange);
% we also need a list of end times.  This needs to be checked to make sure it
% is correct for the case of an incomplete last trajectory
endtime = [time(2:end) ; linpos{day}{epoch}.statematrix.time(end)];

%makes a vector: correctSeq of same length as wellSeq that is 0 if traj is
%incorrect, 1 if traj is incorrect
i = 1; %first traj
% SJ: For First Traj: it is always marked correct here.
% Even is endwell is midwell, it is marked correct. In wellstend, I am marking first start well as always 1=midweel.
% So the potential error is here for a turnaround error on first trial 1->1 to be marked correct, but we will ignore that
if ( ((wellSequence(i) == owell1) || (wellSequence(i) == owell2)) || (wellSequence(i) == midwell) ) 
    correct = correct+1; % marked correct regardless of what the endwell is
    tmpR = 1;
    total = total+1;
    tmpT = 1;
    if ( ((wellSequence(i) == owell1) || (wellSequence(i) == owell2)) )
        tmpT= 0; % marked outbound is endwell is one of the outer well
    end
else
    total = total+1;  % invalid trajectory: will never happen
    tmpR = 0;
    tmpT = -1;
end
rewarded(1) = tmpR;
trajbound(1) = tmpT;

for i = 2:length(wellSequence)
    tmpR = -5; tmpT = -5;
    if ( ((wellSequence(i) == owell1) || (wellSequence(i) == owell2)) & (wellSequence(i-1) == midwell) ) %is an outbound traj
        tmpT= 0;
        if i==2  % If i==2 is outbound, its the first outbound and is always correct      
            correct = correct+1; 
            tmpR = 1;
            total = total+1;
        elseif ( (wellSequence(i-2) ~= wellSequence(i)) ) %if last outward bound well is different from current outward bound well    
            correct = correct+1;
            tmpR = 1;
            total = total+1;
        else
            total = total+1;
            tmpR = 0;
        end
    elseif ( ((wellSequence(i) == owell1) || (wellSequence(i) == owell2)) & (wellSequence(i-1) ~= midwell) ) % out to out
        tmpR = 0;
        tmpT = 1;
    elseif ((wellSequence(i) == midwell) & (wellSequence(i-1) == midwell)) % in to in turnaround
	%error
        tmpT = 0;
        tmpR = 0;
    elseif ( (wellSequence(i) == midwell) ) %inbound traj % only remaining case is a correct inbound
        tmpT = 1;
        tmpR = 1;
    else 
        total = total+1;
        tmpR = 0;
        tmpT = -1;
    end

    if ( ((wellSequence(i) ==owell1) || (wellSequence(i) == owell2)) )
        lastoutwell = wellSequence(i); %set last outward well to current outwell for next trial
    end
    rewarded(i) = tmpR;
    trajbound(i) = tmpT;
end

%apply excludetimes
include = find(~isExcluded(time, excludetimes));
trajtime = [time(include) endtime(include)];
rewarded = rewarded(include);
trajbound = trajbound(include);
wellstend = wellstend(include,:);

