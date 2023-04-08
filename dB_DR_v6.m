%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dB_DR_v6: Data reduction script for processing laser ablation boron
% isotope data. The methodology is outlined in:
%   Evans, Gerdes, Coenen, Marschall, & Mueller (2021) Accurate correction 
%   for the matrix interference on laser ablation MC-ICPMS boron isotope
%   measurements in CaCO3 and silicate matrices. JAAS 36:1607.
%
% If this script is useful for your research, please cite the above paper.
%
% Derive carbonate standardised d11B values. Briefly:
% 1) Remove outliers and baseline subtract sample/standard analyses
% 2) Identify the location of the NIST analyses to be used for mass bias
%   and drift correction, and apply these corrections
% 3) Calculate 'raw' (NIST-standardised) d11B values
% 4) Calculate Dd11B for three carbonate standards, and derive the
%   power relationship between ~B/Ca (11/10.035) and Dd11B
% 5) Correct sample data using the derived relationship between Dd11B and
%   the 11/10.035 ratio
%
% This script is designed to be used in conjuction with a log file from the
% laser ablation system; functionality is greatly limited without.
%
% Two key parameters are hard-coded, the applicability of which should be 
% checked on your system before using this script:
% 1) The script expects an 8-column Neptune exp data file. If you
%   additionally export ratio data (or differ in any other way),
%   comment/uncomment lines 88 and 91 to read in (e.g.) an 11-column file.
% 2) Five rows of data are thrown away at the start of the analysis
%   (potential surface contamination) and four lines from the end (possible
%   cell washout). See line 263.

clear variables
% Select files
dirF = uigetdir('','Select directory containing data files');

[FileName,PathName] = ...
    uigetfile('*.exp','Pick a file',dirF,'MultiSelect','on');

% Load files
w = whos('FileName');
if strcmp('char',w.class)==1
    fidi{1} = fopen([PathName FileName],'rt');
else
    fidi = cell(size(FileName));
    for i = 1:size(FileName,2)
        fidi{i} = fopen([PathName FileName{1,i}],'rt');
    end
end

% Data file read in requires Neptune block export with 995 cycles/block and
% all blocks in the same data file, with headers between each block. (Note 
% 995 cycles/block was chosen as write errors can occur with longer blocks)
% If no ratio data are calculated by the Neptune software, the file will
% have seven columns. If ratio data are calculated and written (11 columns),
% comment/uncomment 'tempD' below as appropriate.
%
% In the case that this data file read is unapplicable to your system, the
% data processing part of the script requires the following datasets in the
% following formats:
% 1) EITHER 
%   RawData AND anTime
%   RawData: An n*6 matrix of voltage data with the following columns:
%   cycle number, m/z 9.95 voltage, 10B voltage, m/z 10.035 voltage', m/z
%   10.5 voltage, 11B voltage. Not all of these voltages are used in
%   the calculations below (m/z 9.95 and 10.5 can be replaced by NaNs if
%   this information was not collected). 
%   anTime: An n*1 matrix of time stamps, where n is equal to the number
%   of rows of data. Date number format.
%   [use this if you intend to supply a log file via GeoStarData to 
%   determine the locations of the analyses but no MC-ICPMS data file]
%
%   OR
%
%   bg AND data
%   An (n+1)*1 and n*1 cell array of separated background and
%   sample data respectively, where each cell has the same column format as
%   that described for 'RawData' above.
%   [use this if your data are already separated into background and sample
%   segments and thus do not require comparison to a laser log file (but
%   note that full script functionality nonetheless later requires a log
%   file, see #2)]   
% 2) In addition, 'AnID' and 'AnOver' are metadata from the laser ablation
%   system, which are later required for full functionality. See the
%   function 'GeoStarData' for formatting.

for j = 1:size(fidi,2)
    % use this line for 7 columns of data
    tempD = textscan(fidi{1,j}, '%s%s%s%s%s%s%s%s', 'Delimiter','\t',...
        'HeaderLines',0, 'CollectOutput',1);
    % Alternatively use this for 11 data columns
    %tempD = textscan(fidi{1,j}, '%s%s%s%s%s%s%s%s%s%s%s', 'Delimiter',...
    % '\t','HeaderLines',0, 'CollectOutput',1);
   
    anDate = tempD{1}{9,1};     % isolate date the file was written
    anDate = anDate(16:size(anDate,2));

    RawData{j,1} = NaN(size(tempD{1},1),size(tempD{1},2)-2);
    time = cell(size(tempD{1},1),1);
    b = 0; a = 0;
    % Needs optimising...
    for i = 18:size(tempD{1},1)
        if b~=floor((i-17)/1024)    % if the block number changes
            a = 0;                  % reset the counter
        else
            a = a+1;
            if a>1 && a<996                % 995 rows of data per block
                RawData{j,1}(i,1:6) = [str2double(tempD{1}{i,1}) str2double(tempD{1}{i,3}) ...
                    str2double(tempD{1}{i,4}) str2double(tempD{1}{i,5}) ...
                    str2double(tempD{1}{i,6}) str2double(tempD{1}{i,7})];
                tTime = cellstr(tempD{1}{i,2});
                if isnan(str2double(tTime)) || isempty(str2double(tTime)) % if the data were exported with date format timestamp
                    time{i,1} = [anDate,' ',tTime{1}];
                else
                    time{i,1} = str2double(tTime);
                end
            end
        end
        b = floor((i-17)/1024); % block number
    end
    
    % throw out emptry columns, 9 empty rows if the raw data files include
    % calculated ratios, otherwise there will be 6 empty rows
    if size(RawData{j,1},2)==9
        RawData{j,1} = RawData{j,1}(sum(isnan(RawData{j,1}),2)~=9,:);
    else
        RawData{j,1} = RawData{j,1}(sum(isnan(RawData{j,1}),2)~=6,:);
    end
    timeOut{j,1} = time(~cellfun('isempty',time));
end
clearvars -EXCEPT RawData dirF anDate timeOut PathName FileName


for j = 1:size(RawData,1)   % delete empty rows if block stopped maually
    for i = 1:size(RawData{j,1},1)
        if sum(isnan(RawData{j,1}(i,2:6)))==5
            timeOut{j,1}{i,1} = [];
        end
    end
    RawData{j,1}(isnan(RawData{j,1}(:,2)),:) = [];
    timeOut{j,1} = timeOut{j,1}(~cellfun('isempty',timeOut{j,1}));
end

% Extract time stamps for every row of data
anTime = cell(size(timeOut));
for j = 1:size(timeOut,1)
    if ischar(timeOut{1,1}{1,:}) % if the data were exported with a datestring timestamp
        anTime{j,1} = NaN(size(timeOut{j,1},1),1);   % convert times to date numbers
        add = 0;
        for i = 1:size(timeOut{j,1},1)
            t = timeOut{j,1}{i,1};
            t(1,20)='.';
            t(1,3) = '-';
            t(1,6) = '-';
            if i==1
                anTime{j,1}(i,1) = datenum(t,'dd-mm-yyyy HH:MM:SS.FFF');
            else
                if datenum(t,'dd-mm-yyyy HH:MM:SS.FFF')-anTime{j,1}(i-1,1)<0
                    add = 1;
                end
                anTime{j,1}(i,1) = datenum(t,'dd-mm-yyyy HH:MM:SS.FFF') + add;
            end
        end
        clear a b t add
    else % else, if exported as unix format
        d = datetime(cell2mat(timeOut{j,1}),'ConvertFrom','epochtime','TicksPerSecond',1,'Format','dd-MMM-yyyy HH:mm:ss.SSS');
        anTime{j,1} = datenum(d);
    end
end
clear d

sepData = RawData;  % save separated data in case required later


if size(RawData,1)>1   % stick multiple data files together
    RawData = cell2mat(RawData);
    anTime = cell2mat(anTime);
else
    RawData = RawData{1,1};
    anTime = anTime{1,1};
end
tempsort = [RawData anTime];  % sort datasets according to time order
tempsort = sortrows(tempsort,size(tempsort,2));
RawData = tempsort(:,1:size(tempsort,2)-1);
time = tempsort(:,size(tempsort,2));
clear tempsort

% load in GeoStar log files
haveLog = input('do you have log files to load in? (1/0)');

if haveLog==1
    % The function GeoStarData extracts analysis metadata from GeoStar LA
    % software.
    [AnID,AnOver] = GeoStarData(dirF);
    a = 0;
    for i = 1:size(AnOver,1)    % delete very short analyses (log file issue)
        if AnOver(i,2)-AnOver(i,1)<1e-5
            AnOver(i,:) = NaN;
            AnID(i-a,:) = [];
            a = a+1;
        end
    end
    AnOver(isnan(AnOver(:,1)),:) = [];
    
    % user selects start of first analysis to determine the time offset
    % between the two computers
    close(figure(1))
    figure(1)
    plot(anTime(1:2000),RawData(1:2000,6))
    set(gca,'yscale','log')
    set(gcf,'color','w')
    ylim([1e-4 0.5])
    gi = ginput(1);
    close(figure(1))

    % calculate time offset between PCs
    loc = [(1:1:size(anTime,1))' abs(anTime-gi(1))];
    loc = sortrows(loc,2);
    Toff = anTime(loc(1,1),1)-AnOver(1,1);
    clear gi

    % Plot the analysis locations based on the log file data
    figure(1)
    plot(anTime,RawData(:,6))
    set(gca,'yscale','log')
    set(gcf,'color','w')
    ylim([1e-4 0.5])
    xlabel('time (date number)')
    ylabel('^{11}B (V)')
    hold on
    for i = 1:size(AnOver,1)
        fill([AnOver(i,1) AnOver(i,2) AnOver(i,2) AnOver(i,1)]+Toff,...
            [1e-4 1e-4 0.5 0.5],[195/255 221/255 240/255],'edgecolor','none',...
            'facealpha',0.5)
    end
%    axis([anTime(1) anTime(1000) ylim])
    drawnow
%     check = 1;
%     while check~=0
%         check = input(['Want to check everything looks good? Enter start and end row numbers (max ',num2str(size(RawData,1)),'), enter 0 to move on.']);
%         if max(size(check))>1
%             axis([anTime(check(1)) anTime(check(2)) ylim])
%             drawnow
%         end
%     end
    checkOK = input('Did everything work? Enter 1 if so, or 0 to try manual data selection.');
    close(figure(1))

    if checkOK==1
        bg = cell(size(AnOver,1)+1,1);
        data = cell(size(AnOver,1),1);
        [~, ind1] = min(abs(AnOver(1,1)-anTime(:,1)+Toff));    % 1st background
        bg{1,1} = RawData(4:ind1-4,:);
        avSep = mean(AnOver(2:size(AnOver,1),1) - AnOver(1:size(AnOver,1)-1,1));  % average spearation

        for i = 1:size(AnOver,1)
            if max(anTime(:,1)-Toff) - AnOver(i,2) < 0  % if data collection stopped before LA stopped
                data{i,1} = [];
                bg{i+1,1} = [];
            else
                [~, ind1] = min(abs(AnOver(i,1)-anTime(:,1)+Toff));
                [~, ind2] = min(abs(AnOver(i,2)-anTime(:,1)+Toff));
                %%%%%%%%%%%%%%%%%%%%%%%% adjust here! %%%%%%%%%%%%%%%%%%%
                % choose how much data to throw away at the start and end
                % of the analysis here:
                data{i,1} = RawData(ind1+5:ind2-4,:); %%%%
                if i<size(AnOver,1)
                    [~, ind3] = min(abs(AnOver(i+1,1)-anTime(:,1)+Toff));
                    if AnOver(i+1,1)-AnOver(i,1)>4*avSep   % if the analyses are separated by more than usual
                        bg{i+1,1} = RawData(ind3-30:ind3-10,:);
                    else                
                        bg{i+1,1} = RawData(ind2+10:ind3-10,:);
                    end
                else % for the final background
                    bg{i+1,1} = RawData(ind2+10:ind2+30,:);
                end
            end
        end
        clear ind1 ind2 ind3
    end
end

% If no laser log file is available, or reading the data for the log file 
% did not work, the next section will try to locate analyses without this
% information.
% This is a *preliminary* attempt at locating sample versus background 
% segments of data if no laser log files are available. Doing so is
% complicated by the fact that voltage spikes are often present in the
% background data, resulting from the occasional washout of particles. 
% In order to avoid these spikes being identified as analyses, the user is 
% i) first asked to screen through the data, to select the location of 
%   spikes for deletion
% ii) then asked for a cutoff point above which data are considered to be
%   a sample
% Having selected a cutoff point, the number of analyses that would be
% produced from this choice is displayed, for comparison to the known
% number of analyses performed. A new cutoff point can be selected if
% necessary.
%
% This section is not optimised as in our lab, laser log files are
% virtually always available. It is provided only as a workaround for the
% rare cases that the log file data is lost.
if haveLog~=1 || exist('checkOK','var') && checkOK==0      % if log files aren't available

    raw = RawData;  % copy data to new matrix for editing
    raw(:,1) = (1:1:size(raw,1))';

    % cycle through the data in blocks of 5000 rows
    for j = 1:floor(size(raw,1)/5000)+1
        for i = 1:1000
            figure(1)
            plot(raw(:,1),raw(:,6))
            title('Click twice around the spike to zoom in, then hit enter. Hit enter to move to the next section.')
            set(gca,'yscale','log')
            set(gcf,'color','w')
            ylabel('^{11}B (V)')
            xlabel('integration')
            if j==1
                axis([0 5000 1e-4 1])
            elseif floor(size(raw,1)/5000)+1==j
                axis([j*5000-5050 size(raw,1) 1e-4 1])
            else
                axis([j*5000-5050 j*5000 1e-4 1])
            end
            % select the approximate location of the spike
            gi = ginput();
            if size(gi,1)==0
                break
            end
            % select the precise location of the spike
            close(figure(2))
            figure(2)
            plot(raw(find(raw(:,1)==round(gi(1,1),0)):find(raw(:,1)==round(gi(2,1),0)),1),...
                raw(find(raw(:,1)==round(gi(1,1),0)):find(raw(:,1)==round(gi(2,1),0)),6))
            set(gca,'yscale','log')
            set(gcf,'color','w')
            title('click next to the spike on either side')
            gi2 = ginput(2);
            % delete the relevant portion of the data
            raw(find(raw(:,1)==round(gi2(1,1),0)):find(raw(:,1)==round(gi2(2,1),0)),:) = [];
        end
        close(figure(2))
    end
    close(figure(1))
    raw(raw(:,6)<0,:) = [];

    spikeCutRaw = raw;

    h = 0;
    while h==0
        % choose cutoff point above which data are considered to be a
        % sample
        figure(1)
        plot(raw(:,1),raw(:,6))
        ylabel('^{11}B (V)')
        set(gca,'yscale','log')
        set(gcf,'color','w')
        xlabel('integration')
        axis([xlim 1e-4 1])

        thresh = ginput(1);
        thresh = thresh(2);
        where = NaN(size(raw,1),1);
        for i = 4:size(raw,1)-3
            if raw(i,6)>thresh && raw(i+1,6)>thresh && raw(i+2,6)>thresh && raw(i+3,6)>thresh && raw(i-1,6)<thresh 
                where(i) = 1;
            elseif raw(i,6)>thresh && raw(i-1,6)>thresh && raw(i+1,6)<thresh && raw(i+2,6)<thresh && raw(i+3,6)<thresh
                where(i) = 2;
            end
        end
        % This section attempts to catch missed spikes (see above) and
        % enables spike deletion in these regions without re-running the
        % script in its entirety
        if sum(where==2)~=sum(where==1)
            h = 0;
            disp('you have a problematic spike in this region')

            a = [(1:1:size(where,1))' where NaN(size(where,1),1)];
            a(isnan(a(:,2)),:) = [];
            a(1:size(a,1)-1,3) = a(2:size(a,1),2) - a(1:size(a,1)-1,2);
            for i = 2:size(a,1)
                if abs(a(i-1,2)-a(i,2))~=1
                    break
                end
            end      
            close(figure(2))
            figure(2)
            plot(raw(a(i-1,1):a(i+1),1),raw(a(i-1,1):a(i+1),6))       
            set(gca,'yscale','log')
            set(gcf,'color','w')
            title('click right next to the spike on either side')
            gi2 = ginput(2);

            raw(find(raw(:,1)==round(gi2(1,1),0)):find(raw(:,1)==round(gi2(2,1),0)),:) = [];      

        else
            h = input(['this will produce ', num2str(sum(where==1)) ,' analyses. Does that sound good? Enter any number if so, otherwise zero   ']);
        end
    end
    close(figure(1))
    rawSaved = raw;     % duplicate matrix before deletions
    clear h gi2 

    % Define the locations of the start and end of the analyses
    for i = 6:size(raw,1)-9 
        if where(i)==1
            raw(i-5:i+3,:) = 0; % switch to sig.
        elseif where(i)==2
            raw(i-3:i+4,:) = 0; % switch to bg
        end
    end
    raw(raw==0)= NaN;

    loc = 1;
    for i = 1:size(raw,1)-1
        if isnan(raw(i+1,3)) && ~isnan(raw(i,3))
            loc = [loc ; i];
        elseif isnan(raw(i,3)) && ~isnan(raw(i+1,3))
            loc = [loc ; i+1];
        end
    end

    % Separate the data into sample and background segments based on the
    % above locations
    bg = cell(round(size(loc,1)/4),1);
    data = cell(round(size(loc,1)/4-1),1);
    a = 1;
    b = 1;
    for i = 1:size(loc,1) % sort data into bg/samp cells
        if rem(i+3,4)==0 && i==size(loc,1) % special case, final BG
            bg{a,1} = raw(loc(b,1):size(raw,1),:);    
        elseif rem(i+3,4)==0
            bg{a,1} = raw(loc(b,1):loc(b+1,1),:);
            a = a+1;
        elseif rem(i+1,4)==0
            data{a-1,1} = raw(loc(b+2):loc(b+3),:);
            b = b+4;
        end
    end
    clear a b
end
clear haveLog checkOK

% This loop enables user-defined start and end points for the analyses to
% be chosen. If this is not done, the mean of all data will be used.
% (useful if e.g. samples broke during ablation)
manCut = input('Do you want to manually select start/end points for your samples? Press any number if so, otherwise zero.   ');
if manCut~=0
    for i = 1:size(data,1)  % manually cut samples
        if ~isempty(data{i,1})
            close(figure(6))
            figure(6)
            yyaxis left
            scatter(1:1:size(data{i},1),data{i}(:,6)./data{i}(:,3),'.')
            yL = ylim;
            axis([1 size(data{i},1) yL(1).*0.95 yL(2)*1.05])
            xlabel('integration')
            ylabel('raw ^{11}B/^{10}B')
            yyaxis right
            plot(1:1:size(data{i},1),data{i}(:,6),'-','linewidth',1.5)
            ylabel('^{11}B (V)')
            ylim([0 max(data{i}(:,6))*1.2])
            yL = ylim;
            xL = xlim;
            title('click twice, or press enter to use the whole analysis')
            set(gcf,'color','w')
            text(xL(2)*0.05,yL(2)*0.95,['analysis ',num2str(i),' of ',num2str(size(data,1))])
            if exist('AnID','var')
                text(xL(2)*0.05,yL(2)*0.05,['sample: ',AnID(i,:)])
            end
            [x,y] = ginput(2);
            if ~isempty(x)
                if size(x,1)==1
                    data{i}(1:round(x,0),:) = NaN;
                else
                    data{i} = data{i}(round(x(1),0):round(x(2),0),:);
                end
            end
            clear x y
        end
    end
    clear yL xL
end
close(figure(6))

% If no background data recorded at the end of the sequence, duplicate the
% penultimate background measurement
bg = bg(~any(cellfun('isempty', bg), 2), :);
data = data(~any(cellfun('isempty', data), 2), :);
if size(bg,1)==size(data,1) % dulicate final BG if acq. cut during analysis
    bg{size(bg,1)+1,1} = bg{size(bg,1),1};
end
% If data acquisition stopped but ablation continued, delete missed 
% analysis log files
if exist('AnOver','var')
    if size(data,1)<size(AnOver,1)
        AnOver(size(data,1)+1:size(AnOver,1),:) = [];
        AnID(size(data,1)+1:size(AnID,1),:) = [];
    end
end


%% data processing starts here
tempMeanBG = NaN(size(bg,1),size(bg{1,1},2)-1,2); % outlier subtraction
tempMeanSamp = NaN(size(data,1),size(data{1,1},2)-1,2);
bgOut = cell(size(bg));
dataOut = cell(size(data));
for i = 1:size(tempMeanBG,1) % backgrounds (2SD)
    bgOut{i,1} = bg{i,1};
    for j = 1:size(tempMeanBG,2)
        tempMeanBG(i,j,1) = mean(bg{i,1}(:,j+1),'omitnan');
        tempMeanBG(i,j,2) = 2.*std(bg{i,1}(:,j+1),'omitnan');
        bgOut{i,1}(bg{i,1}(:,j+1)>tempMeanBG(i,j,1)+tempMeanBG(i,j,2) | ...
            bg{i,1}(:,j+1)<tempMeanBG(i,j,1)-tempMeanBG(i,j,2),:) = NaN;
    end
end
for i = 1:size(tempMeanSamp,1)
    dataOut{i,1} = data{i,1};
    for j = 1:size(tempMeanBG,2) % samples (3SD)
        tempMeanSamp(i,j,1) = mean(data{i,1}(:,j+1),'omitnan');
        tempMeanSamp(i,j,2) = 3.*std(data{i,1}(:,j+1),'omitnan');
        dataOut{i,1}(data{i,1}(:,j+1)>tempMeanSamp(i,j,1)+tempMeanSamp(i,j,2) | ...
            data{i,1}(:,j+1)<tempMeanSamp(i,j,1)-tempMeanSamp(i,j,2),:) = NaN;
    end
end

% calculate average background values for each segment
meanBG = NaN(size(bg,1)-1,5);
for i = 1:size(bg,1)-1
    meanBG(i,:) = mean([bgOut{i,1}(:,2:6) ; bgOut{i+1,1}(:,2:6)],'omitnan');
end

% For each sample segmenet, ratio everything to 10B and calculate the mean, 
% 2SD, and 2SE of all ratios
dataR = cell(size(data,1),2);
dataOver = NaN(size(data,1),8);
dataSD = NaN(size(data,1),4);
dataSE = NaN(size(data,1),4);
dataN = NaN(size(data,1),4);
for i = 1:size(data,1)
    dataR{i,1} = dataOut{i,1}(:,2:6) - meanBG(i,:);
    dataR{i,2} = repmat(dataR{i,1}(:,5),1,4)./dataR{i,1}(:,1:4);
    dataOver(i,:) = [mean(dataR{i,2},'omitnan') ...
        mean(dataOut{i,1}(:,6),'omitnan') ...
        meanBG(i,5) mean(dataOut{i,1}(:,3),'omitnan') meanBG(i,2)];
    dataSD(i,:) = 2.*std(dataR{i,2},'omitnan');
    dataSE(i,:) = dataSD(i,:)./sum(~isnan(dataOut{i,1}(:,3)))^0.5;
    dataN(i,:) = sum(~isnan(dataOut{i,1}(:,3)));
    
    % remove outliers from 11/10 ratio (removal above was for 
    % concentration spikes)
    tMean = mean(dataR{i,2}(:,2),'omitnan');
    tSTD = std(dataR{i,2}(:,2),'omitnan');
    dataR{i,2}(dataR{i,2}(:,2) > tMean + 2*tSTD | dataR{i,2}(:,2) < tMean - 2*tSTD,:) = NaN;
    
end
clear tSTD tMean

% Check again for log files if the manual sample location option was taken
% but these are nonetheless available
if exist('AnOver','var')==0
    logIn = input('do you have a log file for sample names? (1/0)');
    if logIn==1
        [AnID,AnOver] = GeoStarData(dirF);
    end
end
clear logIn

% Mass bias correction, drift correction, conversion to NIST951 d11B
if exist('AnOver','var')==1     % if log files were loaded
    
    % This section requires AnID AnOver (from GeoStarData) and 
    % dataOver (above)
    anTypes = unique(AnOver(:,10)); % Unique sets of beam diameters
    % Duplicate key datasets as non-relevant portions for a given set of
    % ablation conditions will later be removed
    AnIDsave = AnID;
    AnOverSave = AnOver;
    dataOverSave = dataOver;
    MBsampOut = cell(size(anTypes,1),1);
    mdlOut = cell(size(anTypes,1),1);
    % Data processing is performed separately for each set of ablation
    % conditions
    for j = 1:size(anTypes,1)
        AnID = AnIDsave;    % read the key variables back in
        AnOver = AnOverSave;
        dataOver = dataOverSave;
        
        % Remove data not collected under the current ablation conditions
        dataOver(AnOver(:,10)~=anTypes(j,1),:) = [];
        AnID(AnOver(:,10)~=anTypes(j,1),:) = [];
        AnOver(AnOver(:,10)~=anTypes(j,1),:) = [];
        
        % Find the location of NIST analyses
        nistLoc = contains(cellstr(AnID),'nist','IgnoreCase',true);

        % Create a matrix for each analysis in the current subset of data
        % that defines which NIST analyses are to be used for mass bias
        % and drift correction
        whichNIST = NaN(size(nistLoc,1),2);
        a = 1;
        whichNIST(1,1) = 1;
        for i = 2:size(whichNIST,1)
            if nistLoc(i-1,1)<nistLoc(i,1)
                a = a+1;
            end
            whichNIST(i,1) = a;
        end
        % NIST is usually measured in duplicate to avoid spurious drift
        % correction, calculate mean values before proceeding
        avNIST = NaN(max(whichNIST(:,1)),8);
        for i = 1:max(whichNIST(:,1))
            avNIST(i,:) = mean(dataOver(whichNIST(:,1)==i & nistLoc==1,:),1);
        end
        NISTrat = 4.041648165; % Absolute 11/10 of NIST SRM612 (GeoReM)
        
        % Mass bias correct all sample data using the average of the
        % bracketing NIST analyses (column 1 of MBsamp) and calculate the
        % mean voltage ratio between the bracketing NISTs and the sample
        % (column 2 of MBsamp) to approximately assess [B]
        MBsamp = NaN(size(dataOver,1),4);
        for i = 1:size(whichNIST,1)-1
            if whichNIST(i,1)==max(whichNIST(:,1))
                    MBsamp(i,1) = dataOver(i,2)...
                        *NISTrat/mean([avNIST(whichNIST(i,1),2) avNIST(whichNIST(i,1)-1,2)],'omitnan');
                    MBsamp(i,2) = dataOver(i,5)...
                        /mean([avNIST(whichNIST(i,1),5) avNIST(whichNIST(i,1)-1,5)],'omitnan');
            else
                MBsamp(i,1) = dataOver(i,2)...
                    *NISTrat/mean([avNIST(whichNIST(i,1),2) avNIST(whichNIST(i,1)+1,2)],'omitnan');
                MBsamp(i,2) = dataOver(i,5)...
                    /mean([avNIST(whichNIST(i,1),5) avNIST(whichNIST(i,1)+1,5)],'omitnan');
            end
            MBsamp(i,3) = (MBsamp(i,1)/4.04367 - 1)*1000;   % raw d11B
        end
        % Calculate the approximate [B] (ppm) 
        MBsamp(:,4) = MBsamp(:,2).*34.3; % 34.3 ppm [B] in NIST SRM612
        
        % Find the locations of the secondary carbonate standards
        jcpLoc = contains(cellstr(AnID),'jcp','IgnoreCase',true); % find JCp positions
        jctLoc = contains(cellstr(AnID),'jct','IgnoreCase',true); % find JCp positions
        macsLoc = contains(cellstr(AnID),'macs','IgnoreCase',true); % find JCp positions
    
        % Filter out nano pellets, if measured
        jcpLocNP = contains(cellstr(AnID),'jcpnp','IgnoreCase',true); % find JCp positions
        jctLocNP = contains(cellstr(AnID),'jctnp','IgnoreCase',true); % find JCp positions
        macsLocNP = contains(cellstr(AnID),'macs3np','IgnoreCase',true); % find JCp positions
        jcpLoc(jcpLocNP==1,1) = 0;
        jctLoc(jctLocNP==1,1) = 0;
        macsLoc(macsLocNP==1,1) = 0;

        % Do not inlcude analyses with a repetition rate less than 4
        jcpLoc(AnOver(:,9)<4,1) = 0;
        jctLoc(AnOver(:,9)<4,1) = 0;
        macsLoc(AnOver(:,9)<4,1) = 0;
        % or with negative 11/10.035 ratios
        jcpLoc(dataOver(:,3)<0,1) = 0;
        jctLoc(dataOver(:,3)<0,1) = 0;
        macsLoc(dataOver(:,3)<0,1) = 0;

        % If no secondary carbonate standards were found, do no further
        % data processing
        if sum(jcpLoc)==0 || sum(jctLoc)==0 || sum(macsLoc)==0
            MBsamp(:,5:6) = NaN;
        else        
            % Else calculate the Dd11B of the measured carbonate standards,
            % use this information to derive a correction line, and apply
            % it to the sample data
            cmap = parula(12);
            close(figure(2))
            figure(2)   % plot offset standard data
            hold on
            set(gcf,'color','w')
            xlabel('11/10.089')
            ylabel('\Delta\delta^{11}B')
            scatter(dataOver(macsLoc,3),MBsamp(macsLoc,3)-0.57,20,cmap(3,:),'filled')
            scatter(dataOver(jcpLoc,3),MBsamp(jcpLoc,3)-24.3,20,cmap(6,:),'filled')
            scatter(dataOver(jctLoc,3),MBsamp(jctLoc,3)-16.3,20,cmap(9,:),'filled')
            xL = xlim;
            axis([xL(1) xL(2) ylim])

            % combine datasets for fitting
            tempX = [dataOver(macsLoc,3) ; dataOver(jcpLoc,3) ; dataOver(jctLoc,3)];
            tempY = [MBsamp(macsLoc,3)-0.57 ; MBsamp(jcpLoc,3)-24.3 ; MBsamp(jctLoc,3)-16.3];
            tempY = tempY(tempX>0);
            tempX = tempX(tempX>0);
            tempXY = [tempX tempY];
            
            try % model fails in some instances where poor data quality means fit cannot be found
                % Calcuate power regression parameters and add to the plot
                modelfun = @(b,x) b(1)+b(2).*x.^b(3);
                beta0 = [50 -1000 -1];
                mdl = fitnlm(tempX,tempY,modelfun,beta0);
                p1 = plot((round(min(tempX),1):0.1:round(max(tempX),1)),...
                    mdl.Coefficients{1,1}+...
                    mdl.Coefficients{2,1}.*(round(min(tempX),1):0.1:round(max(tempX),1)).^...
                    mdl.Coefficients{3,1},'-','color','k');
                [beta,R,J,COVB,MSE] = nlinfit(tempX,tempY,modelfun,beta0);
                [ypred,delta] = nlpredci(modelfun,...
                    (round(min(tempX),1):0.1:round(max(tempX),1)),beta,R,'Covar',COVB,...
                    'MSE',MSE,'SimOpt','off');
                lower = ypred - delta;
                upper = ypred + delta;
                p2 = plot((round(min(tempX),1):0.1:round(max(tempX),1)),...
                    [lower;upper],'--','color','k','linewidth',0.5);
                drawnow
    
                % The carbonate standard dataset frequently contains
                % outliers, especially when analysed at low repetition
                % rates or small beam diameters. Manually throw these out
                % by clicking on the figure, one at a time.
                out = input('want to delete outliers? (1/0)');
                while out~=0
                    % Select offending data point
                    gi = ginput(1);
                    %dist = sum((gi(1,:) - tempXY) .^ 2, 2);
                    dist = sum(([gi(1,1)/10 gi(1,2)] - ...
                        [tempXY(:,1)./10 tempXY(:,2)]) .^ 2, 2);
                    scatter(tempXY(dist == min(dist),1),...
                        tempXY(dist == min(dist),2),60,'r')
                    % Remove it from the calibration dataset
                    tempXY(dist == min(dist),:) = [];
                    delete(p1)
                    delete(p2)
                    % Recalculate model and replot
                    mdl = fitnlm(tempXY(:,1),tempXY(:,2),modelfun,beta0);
                    p1 = plot((round(min(tempX),1):0.1:round(max(tempX),1)),...
                        mdl.Coefficients{1,1}+...
                        mdl.Coefficients{2,1}.*(round(min(tempX),1):0.1:round(max(tempX),1)).^...
                        mdl.Coefficients{3,1},'-','color','k');
                    [beta,R,J,COVB,MSE] = nlinfit(tempXY(:,1),tempXY(:,2),modelfun,beta0);
                    [ypred,delta] = nlpredci(modelfun,...
                        (round(min(tempX),1):0.1:round(max(tempX),1)),beta,R,'Covar',COVB,...
                        'MSE',MSE,'SimOpt','on');
                    lower = ypred - delta;
                    upper = ypred + delta;
                    p2 = plot((round(min(tempX),1):0.1:round(max(tempX),1)),...
                        [lower;upper],'--','color','k','linewidth',0.5);
                    drawnow
    
                    out = input('want to delete more outliers? (1/0)');
                end
                clear tempX tempY tempXY xL a beta R J COVB MSE ypred delta gi beta0 p1 p2 upper lower
                close(figure(2))
                
                % Apply the resulting power function to the sample data
                MBsamp(dataOver(:,3)>0,5) = mdl.Coefficients{1,1} + mdl.Coefficients{2,1}.*...
                    dataOver(dataOver(:,3)>0,3).^mdl.Coefficients{3,1};
                MBsamp(:,6) = MBsamp(:,3) - MBsamp(:,5);
            catch
                close(figure(2))
    
                MBsamp(dataOver(:,3)>0,5) = NaN;
                MBsamp(:,6) = NaN;

                warning('Data could not be fit. Only NIST612-corrected data will be output.')
            end
        end

        MBsampOut{j,1} = MBsamp;

        % read the key variables back in (unnecessary except in the last
        % loop step)
        AnID = AnIDsave;
        AnOver = AnOverSave;
        dataOver = dataOverSave;

    end
    clear dataOverSave AnOverSave AnIDsave out

    % Compile the results from multiple beam diameters back into one matrix
    MBsamp = NaN(size(AnOver,1),6);
    for j = 1:size(MBsampOut,1)
        MBsamp(find(AnOver(:,10)==anTypes(j,1),1):...
            find(AnOver(:,10)==anTypes(j,1),1)+size(MBsampOut{j,1},1)-1,:) = ...
            MBsampOut{j,1};
    end
    % Create data table for output
    dataTab = table(AnID,AnOver(:,9),AnOver(:,10),dataOver(:,1),dataOver(:,2),dataOver(:,3),dataOver(:,4),...
        dataOver(:,5),dataOver(:,6),dataOver(:,7),dataOver(:,8),...
        dataSE(:,1),dataSE(:,3),dataSE(:,4),...
        MBsamp(:,1),MBsamp(:,2),MBsamp(:,3),MBsamp(:,4),MBsamp(:,5),MBsamp(:,6),...
        dataSE(:,2)./dataOver(:,2).*1000,...
        'variablenames',{'sample','repetition rate','spot size','11/9.95','11/10','11/10.089','11/10.5',...
        '11B sig. (V)','11B bg (V)','10B sig. (V)','10B bg (V)',...
        '11/9.95 2SE','11/10.089 2SE','11/10.5 2SE',...
        '11/10 mass bias corrected','sample/NIST 11B V ratio','raw d11B',...
        'approx. [B] (ppm)','d11B correction','calibrated d11B',...
        '11/10 2SE (permil)'});

    saveName = input('Enter filename for save (single quotation marks around text in Matlab!):    ');
    
    % Save regression data and data table
    for j = 1:size(mdlOut,1)
        mdlOutSave = evalc('disp(mdlOut{j,1})');        % print regression to text file
        fID = fopen([saveName '_BCa_regression_' num2str(anTypes(j,1)) 'um_spot_data.txt'],'wt');
        fprintf(fID,'%s',mdlOutSave);
        fclose(fID);
        clear fID beta0
    end
    writetable(dataTab,[saveName '_processed_data'])
    save(saveName)
else
   % Output raw data if no way of IDing samples/standards
   warning('you are on your own from here!...')
   dataTab = table(dataOver(:,1),dataOver(:,2),dataOver(:,3),dataOver(:,4),...
        dataOver(:,5),dataOver(:,6),dataOver(:,7),dataOver(:,8),...
        dataSE(:,1),dataSE(:,3),dataSE(:,4),...
        MBsamp(:,1),MBsamp(:,2),MBsamp(:,3),MBsamp(:,4),MBsamp(:,5),MBsamp(:,6),...
        dataSE(:,2)./dataOver(:,2).*1000,...
        'variablenames',{'11/9.95','11/10','11/10.089','11/10.5',...
        '11B sig. (V)','11B bg (V)','10B sig. (V)','10B bg (V)',...
        '11/9.95 2SE','11/10.089 2SE','11/10.5 2SE',...
        '11/10 mass bias corrected','sample/NIST 11B V ratio','raw d11B',...
        'approx. [B] (ppm)','d11B correction','calibrated d11B',...
        '11/10 2SE (permil)'});
    
    saveName = input('Enter filename for save (single quotation marks around text in Matlab!):    ');
    writetable(dataTab,[saveName '_processed_data'])
    save(saveName)
end


% to do
% - option to split long files into two if standards drift?
% - option to base the carboante standard regression on the average of all
%    data for a given standard in the case of very noisy data?
