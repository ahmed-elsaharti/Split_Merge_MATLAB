%% Ahmed Elsaharti 2019

clear all
clc
clf

%Browse File
[filename, pathname] = uigetfile('*.csv', 'Select a data file');
if isequal(filename,0)
   error('No data file selected')
else
   disp(['File selected ', fullfile(pathname, filename)])
end
data = load(fullfile(pathname, filename));

%remove garbage data 
data = data(1:6:end)';

%find the x and y points
angle = pi/180*linspace(120,-120,length(data))';
x = data.*cos(angle);
y = data.*sin(angle);

%remove all the zero, zero points
I = find( sqrt((x).^2+(y).^2) >100);
xy = [x(I) y(I)];

%Remove outliers
idx=[];
for ii=2:(length(xy)-1)
    x1=xy(ii-1,1);
    y1=xy(ii-1,2);
    x2=xy(ii,1);
    y2=xy(ii,2);
    x3=xy(ii+1,1);
    y3=xy(ii+1,2);
    if sqrt((x2-x1)^2+(y2-y1)^2)<800 && sqrt((x2-x3)^2+(y2-y3)^2)<800
    idx=[idx; ii];
    end
end
xy = [xy(1,:); xy(idx,1) xy(idx,2); xy(length(xy),:)];

%Plot initial datset
figure(1)
plot(xy(:,1), xy(:,2), 'r*')
hold on
axis equal

%Initialize the main dataset into a set
Lines{1}=xy;

%Initialize 'checklist' companion matrix 
checklist=[1];


%-----------Keep going in the check loop as long as not all the datasets in
%-----------'Lines' are marked as checked
%-----------This is to guarantee that regardless of order, as long as a set
%-----------is yet to be check we'll stay in the loop and prevents problems
%-----------due to sloppy logic in my code that would have skipped some datasets
%-----------due to order 
while sum(checklist)~=0
    
    % Run the loop as many times as there are lines right now 
    %(this is effective but the moment you split a line into two this is garbage this
    %is why we rely on the outer loop and the checklist).
    for ii=1:length(Lines)
        
        % Run the algorithm on the current dataset in 'Lines' ONLY if it hasn't
        % been already checked and marked for completion
        if checklist(ii)==1
            
                    %Use line fitting function to fit a line to the dataset of the
                    %current iteration
                    [xyfit,RANDSH]=fitline(Lines{ii});
                    
                    
                    %obtain the equation of the line we fit
                    lineparam = polyfit(xyfit(:,1), xyfit(:,2), 1);
                    b=1;
                    a=-lineparam(1);
                    c=-lineparam(2);
                    d=abs(a*Lines{ii}(:,1)+b*Lines{ii}(:,2)+c)/sqrt(a^2+b^2);
                    
                    %Chop off the ends of the fitted line
                    d(1:round(0.25*length(d)))=0;
                    d(round(0.75*length(d)):end)=0;

                    %Run the split only if the max distance between the points and line
                    %is greater than 20
                    if max(d)>20
                                    %Find the index of the split with reference to the d matrix
                                    %(which is also in the same order as the Lines and the fitted
                                    %line matrices
                                    splitidx=find(d==max(d));
                                    
                                    %Plot the split point
                                    plot(Lines{ii}(splitidx,1),Lines{ii}(splitidx,2),'og','Linewidth',3)
                                    
                                    %Check if the current line being checked is the first one in
                                    %the list, this matters only due to math (1-1=0 and we dont
                                    %want that)
                                    if ii~=1
                                        Lines={Lines{1:ii-1},Lines{ii},Lines{ii:end}};
                                        checklist=[checklist(1:ii-1),1,checklist(ii:end)];
                                        Lines{ii+1}=[Lines{ii}(splitidx:end,:)];
                                        Lines{ii}=[Lines{ii}(1:splitidx,:)];
                                        
                                        %Skip to the next line (iterate ii so that when it goes up
                                        %to the beginning of the loop it goes to ii+2 ie: going to
                                        %the line AFTER the one we just created.
                                        ii=ii+1;
                                    else
                                        Lines={Lines{ii},Lines{ii:end}};
                                        checklist=[1,checklist(ii:end)];             
                                        Lines{ii+1}=[Lines{ii}(splitidx:end,:)];
                                        Lines{ii}=[Lines{ii}(1:splitidx,:)];
                                        
                                        %Skip to the next line (iterate ii so that when it goes up
                                        %to the beginning of the loop it goes to ii+2 ie: going to
                                        %the line AFTER the one we just created.
                                        ii=ii+1;
                                    end
                    
                    %Otherwise mark the line as 'checked' (aka is a good fit and no
                    %longer needs splitting
                    else
                        checklist(ii)=0;
                    end
        end
    end
end

for ii=1:length(Lines)
    
    % Fit a line into each of the subsets
    [xyfit,r,alpha]=fitline(Lines{ii});
        
    % Plot the fitted lines
    plot(xyfit(:,1), xyfit(:,2), 'g-','Linewidth',4)
end




% Merge Attempt
finalr=[];
finalalpha=[];
for ii=1:length(Lines)
    % Find alpha and r for all lines
    [xyfit,r,alpha]=fitline(Lines{ii});
    finalr=[finalr; r];
    finalalpha=[finalalpha; alpha];
end

mergedlines=0;
ii=1;
while ii<length(Lines)
    jj=(ii+1)
        if abs(finalr(ii)-finalr(jj))< 100 && abs(finalalpha(ii)-finalalpha(jj))<0.5
           deltar=finalr(ii)-finalr(jj)
           deltaalpha=finalalpha(ii)-finalalpha(jj)
           disp(['MERGE lines ',num2str(ii),' and ',num2str(jj)])
           Lines{ii}=[Lines{ii};Lines{jj}];
           Lines(jj)=[];
           disp(['Getting New Alpha R matrix'])
           mergedlines=mergedlines+1;
           finalr=[];
           finalalpha=[];
           for mm=1:length(Lines)
                % Find alpha and r for all lines
                [xyfit,r,alpha]=fitline(Lines{mm});
                finalr=[finalr; r];
                finalalpha=[finalalpha; alpha];
           end
           %pause
           ii=0;
           
        end
    ii=ii+1
end

% After splitting is done and all lines have been checked and are optimal
% (aka max error between each point and each line is within the set value)
% Plot fitted lines to the segmented datasets (remember the split algorithm
% only helps segment a large dataset into smaller more related ones)
for ii=1:length(Lines)
    
    % Fit a line into each of the subsets
    [xyfit,r,alpha]=fitline(Lines{ii});
    
    % Plot the fitted lines
    plot(xyfit(:,1), xyfit(:,2), 'b-','Linewidth',2);
end
mergedlines
axis square