function [AnID,AnOver] = GeoStarData(dirName)
% Produces a sequence of analyses from Geostar log files selected by user
% given a starting directory input.
%
% AnID - user sample names entered into GeoStar
% AnOver - ablation conditions, see details below for column definitions
%
% Columns give the following information;
% 1 & 2 - start and end date numbers
% 3 & 4 - start x and y co-ordinates
% 5 & 6 - end x and y co-ordinates (NaN if analysis type 1)
% 7 - analysis type (1 = drill, 2 = track, 3 = track with multiple points
%                       4 = map)
% 8 - scan speed (NaN if analysis type 1)
% 9 - repetition rate
% 10 - spot size - NaN if slit was used
% 11 - He MFC flow rate
% 12 - N2 MFC flow rate
%
% NOTE: analyses shorter than 5 s are deleted, as these probably reflect
% samples that were manually stopped and restarted during analysis. If your
% sequence contains real analyses shorter than 5 s then the value on line
% 459 needs to be changed.


[FileName,PathName] = ...
    uigetfile('*.csv','select all log files',dirName,'MultiSelect','on');
%cd(PathName)


if ischar(FileName)
    fileList{1,1} = [PathName,FileName];
else
    fileList = cell(size(FileName,2),1);
    for i = 1:size(FileName,2)      % dir structure for all files in 
        fileList{i,1} = [PathName,FileName{1,i}];
    end
end



SumLogSum = size(fileList,2);        % Amount of log files
ACSV = cell(1,SumLogSum);
FID = nan(1,size(fileList,2));
for i = 1:size(fileList,1)        % Load all log files
    FID(i) = fopen(fileList{i,1});
    fline = fgets(FID(i));
    AllCSV = char(fline);
    while fline ~= -1
        fline = fgets(FID(i));
        AllCSV = char(AllCSV,fline);
    end
    ACSV{i} = AllCSV;
    if size(ACSV{i},1)<5        % Remove blank Geostar logs
        ACSV{i} = [];
    end
end
ACSV = ACSV(~cellfun('isempty',ACSV));
fclose('all');
if isempty(ACSV)
    AnID = [];
    NumLogData = [];
    close(h)
else
    CombLog = cell(1);                  % Combines log files
    CombLog{1} = ACSV{1};
    for i = 2:size(ACSV,2)
        CombLog{1} = char(CombLog{1},ACSV{i});
    end
    A = CombLog{1};         % Messy...
    clear CombLog
    CombLog = A;
    clear A

    CombLogLog = zeros(size(CombLog,1),1);   % Removes all header rows
    for j = 1:size(CombLog,1)
        if ~isempty(strfind(CombLog(j,:),'Timestamp'))
            CombLogLog(j) = 1;
        end
    end
    CombLog = CombLog(CombLogLog==0,:);     % Removes empty cells

    SepLog = cell(1,size(CombLog,1));       % Comma separated information
    for i = 1:size(CombLog,1)               % to cells
        SepLog{i} = textscan(CombLog(i,:),'%s','Delimiter',',');
    end

    LogDate = NaN(size(SepLog,2),1);    % Data from log files
    Sequence = NaN(size(SepLog,2),1);
    SubPoint = NaN(size(SepLog,2),1);
    VertixNo = NaN(size(SepLog,2),1);
    Xum = NaN(size(SepLog,2),1);
    Yum = NaN(size(SepLog,2),1);
    ScanSpeed = NaN(size(SepLog,2),1);
    LaserRep = NaN(size(SepLog,2),1);
    SpotSize = NaN(size(SepLog,2),1);
    MFC1 = NaN(size(SepLog,2),1);
    MFC2 = NaN(size(SepLog,2),1);
    for i = 1:size(SepLog,2)
        if ~isempty(SepLog{i}{1}) && size(SepLog{i}{1},1)>16
        if ~isempty(SepLog{i}{1}{1})
            LogDate(i) = datenum(SepLog{i}{1}(1,:));
            Sequence(i) = str2double(char(SepLog{i}{1}(2)));
            SubPoint(i) = str2double(char(SepLog{i}{1}(3)));
            VertixNo(i) = str2double(char(SepLog{i}{1}(4)));
            Xum(i) = str2double(char(SepLog{i}{1}(6)));
            Yum(i) = str2double(char(SepLog{i}{1}(7)));
            ScanSpeed(i) = str2double(char(SepLog{i}{1}(10)));
            LaserRep(i) = str2double(char(SepLog{i}{1}(12)));
            if ~isnan(str2double(SepLog{i}{1}(14)))
                SpotSize(i) = str2double(char(SepLog{i}{1}(14)));
                
            else
                % <2 changed to <1 on 30.04.22 - probably a mistake as this
                % should never be >1, but flagging in case it wasn't
                if size(SepLog{i}{1}(14),1)<1   % if nothing recorded
                    SpotSize(i) = 0;
                else
                    mLoc = regexp(SepLog{i}{1}{14},'x');
                    switch mLoc
                        case 3
                            sDim = SepLog{i}{1}{14}(1);
                            if size(SepLog{i}{1}{14},2)-mLoc==3
                                lDim = SepLog{i}{1}{14}(5:6);
                            else
                                lDim = SepLog{i}{1}{14}(5:7);
                            end
                        case 4
                            sDim = SepLog{i}{1}{14}(1:2);
                            if size(SepLog{i}{1}{14},2)-mLoc==3
                                lDim = SepLog{i}{1}{14}(6:7);
                            else
                                lDim = SepLog{i}{1}{14}(6:8);
                            end
                        case 5
                            sDim = SepLog{i}{1}{14}(1:3);
                            if size(SepLog{i}{1}{14},2)-mLoc==3
                                lDim = SepLog{i}{1}{14}(7:8);
                            else
                                lDim = SepLog{i}{1}{14}(7:9);
                            end
                        case 6 % added on 30.04.22 to catch square spots of the FIERCE system. Square spots greater than 100 um will still not be read?
                            sDim = SepLog{i}{1}{14}(1:2);
                            if size(SepLog{i}{1}{14},2)-mLoc==5
                                lDim = SepLog{i}{1}{14}(1:2);
                            else
                                lDim = SepLog{i}{1}{14}(4:6);
                            end
                    end
                    SpotSize(i) = str2double([sDim ('99') lDim]); 
                end
            end
            MFC1(i) = str2double(char(SepLog{i}{1}(16)));
            MFC2(i) = str2double(char(SepLog{i}{1}(17)));
        end
        end
    end
    NumLogData = [LogDate Sequence SubPoint VertixNo Xum Yum ScanSpeed...
        LaserRep SpotSize MFC1 MFC2];

    A = SepLog{1}{1}(5);            % Analysis Sequence
    AnID = char(A{1});          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 2:size(SepLog,2)
        if size(SepLog{i}{1},1)<17
            AnID = AnID;
        elseif ~isempty(SepLog{i}{1}{1}) && isempty(SepLog{i}{1}{5}) &&...
                ~isempty(SepLog{i}{1}{2})
            AnID = char(AnID,'NO LABEL');
        elseif ~isempty(SepLog{i}{1}{1}) && ~isempty(SepLog{i}{1}{5})
            A = SepLog{i}{1}(5);
            AnID = char(AnID,A{1});
        end
    end
end

                                                    
IDnum = ones(size(NumLogData,1),1);
for i = 2:size(NumLogData,1)          % Find each new analysis
    if isnan(NumLogData(i,2))
        IDnum(i) = IDnum(i-1);
    else
        IDnum(i) = IDnum(i-1)+1;
    end
end
IDN = IDnum(length(IDnum));         % Number of analyses

AnSep = cell(1,IDN);    % Separate cell for each analysis
for i = 1:size(NumLogData,1)
    if isempty(AnSep{1,IDnum(i)})
        AnSep{1,IDnum(i)} = NumLogData(i,:);
    else
        AnSep{1,IDnum(i)} = [AnSep{1,IDnum(i)} ; NumLogData(i,:)];
    end
end

% A = 0;                          % Removes gaps from AnID (still required?)
% B = AnID;
% for i = 1:size(AnID,1)
%     if strcmp(AnID(i,1),' ')==1
%         B(i-A,:) = [];
%         A = A+1;
%     end
% end        
% AnID = B;

A = 0;
AnOver = NaN(IDN,10);   % Log data to one matrix and analysis type
for i = 1:size(AnSep,2)
    if size(AnSep{1,i},1)>4    % If the run was not aborted before analysis
        if sum(~isnan(AnSep{1,i}(:,4)))==0    % If it is drill data
            AnOver(i,1) = AnSep{1,i}(4,1);        % Start time
            AnOver(i,2) = AnSep{1,i}(5,1);        % End time
            AnOver(i,3) = AnSep{1,i}(4,5);        % X co-ordinate
            AnOver(i,4) = AnSep{1,i}(4,6);        % Y co-ordinate
            AnOver(i,7) = 1;                    % Analysis type
            AnOver(i,8) = 0;                    % zero scan speed when spot
            AnOver(i,9) = AnSep{1,i}(4,8);        % Repetition rate
            AnOver(i,10) = AnSep{1,i}(4,9);       % Spot size
            if AnSep{1,i}(4,10)~=0
                AnOver(i,11) = AnSep{1,i}(4,10);        % MFC1
                AnOver(i,12) = AnSep{1,i}(4,11);        % MFC2
            else
                AnOver(i,11) = AnSep{1,i}(size(AnSep{1,i},1),10);                
                AnOver(i,12) = AnSep{1,i}(size(AnSep{1,i},1),11); 
            end
        else
            if AnSep{1,i}(1,3)==1
                if sum(~isnan(AnSep{1,i}(:,4)))==2 % If it is a point to point track
                                    % If it contains two Vertix entries (why??)
                    place = 1:1:size(AnSep{1,i},1);
                    place = place(~isnan(AnSep{1,i}(:,4)));
                    place = place(length(place));
                    if datenum(AnSep{1,i}(size(AnSep{1,i},1),1))-...    % If the last two entries have the same time
                            datenum(AnSep{1,i}(size(AnSep{1,i},1)-1,1))<1.7361e-05   % = 1.5 seconds in datenum time                        
                        if datenum(AnSep{1,i}(place,1))-...    % If the last vertix number and row above have the same time
                            datenum(AnSep{1,i}(place-1,1))<1.7361e-05
                            if place==size(AnSep{1,i},1)    % If the last vertix number is on the last line
                                AnOver(i,1) = AnSep{1,i}(place-2,1);        % Start time
                                AnOver(i,2) = AnSep{1,i}(place-1,1);        % End time
                                AnOver(i,3) = AnSep{1,i}(place-2,5);        % X co-ordinate start
                                AnOver(i,4) = AnSep{1,i}(place-2,6);        % Y co-ordinate start
                                AnOver(i,5) = AnSep{1,i}(place-1,5);        % X co-ordinate end
                                AnOver(i,6) = AnSep{1,i}(place-1,6);        % Y co-ordinate end
                                AnOver(i,7) = 2;                    % Analysis type
                                AnOver(i,8) = AnSep{1,i}(place,7);        % Scan speed
                                AnOver(i,9) = AnSep{1,i}(place,8);        % Repetition rate
                                AnOver(i,10) = AnSep{1,i}(place,9);      % Spot size
                                if AnSep{1,i}(place-1,10)~=0
                                    AnOver(i,11) = AnSep{1,i}(place-1,10);        % MFC1
                                    AnOver(i,12) = AnSep{1,i}(place-1,11);        % MFC2
                                else
                                    AnOver(i,11) = AnSep{1,i}(size(AnSep{1,i},1),10);                
                                    AnOver(i,12) = AnSep{1,i}(size(AnSep{1,i},1),11); 
                                end
                            else
                                AnOver(i,1) = AnSep{1,i}(place,1);        % Start time
                                AnOver(i,2) = AnSep{1,i}(place+1,1);        % End time
                                AnOver(i,3) = AnSep{1,i}(place,5);        % X co-ordinate start
                                AnOver(i,4) = AnSep{1,i}(place,6);        % Y co-ordinate start
                                AnOver(i,5) = AnSep{1,i}(place+1,5);        % X co-ordinate end
                                AnOver(i,6) = AnSep{1,i}(place+1,6);        % Y co-ordinate end
                                AnOver(i,7) = 2;                    % Analysis type
                                AnOver(i,8) = AnSep{1,i}(place,7);        % Scan speed
                                AnOver(i,9) = AnSep{1,i}(place,8);        % Repetition rate
                                AnOver(i,10) = AnSep{1,i}(place,9);       % Spot size
                                if AnSep{1,i}(place,10)~=0
                                    AnOver(i,11) = AnSep{1,i}(place,10);        % MFC1
                                    AnOver(i,12) = AnSep{1,i}(place,11);        % MFC2
                                else
                                    AnOver(i,11) = AnSep{1,i}(size(AnSep{1,i},1),10);                
                                    AnOver(i,12) = AnSep{1,i}(size(AnSep{1,i},1),11); 
                                end
                            end                          
                        else
                            AnOver(i,1) = AnSep{1,i}(place-1,1);        % Start time
                            AnOver(i,2) = AnSep{1,i}(place,1);        % End time
                            AnOver(i,3) = AnSep{1,i}(place-1,5);        % X co-ordinate start
                            AnOver(i,4) = AnSep{1,i}(place-1,6);        % Y co-ordinate start
                            AnOver(i,5) = AnSep{1,i}(place,5);        % X co-ordinate end
                            AnOver(i,6) = AnSep{1,i}(place,6);        % Y co-ordinate end
                            AnOver(i,7) = 2;                    % Analysis type
                            AnOver(i,8) = AnSep{1,i}(place,7);        % Scan speed
                            AnOver(i,9) = AnSep{1,i}(place,8);        % Repetition rate
                            AnOver(i,10) = AnSep{1,i}(place,9);       % Spot size
                            if AnSep{1,i}(place,10)~=0
                                AnOver(i,11) = AnSep{1,i}(place,10);        % MFC1
                                AnOver(i,12) = AnSep{1,i}(place,11);        % MFC2
                            else
                                AnOver(i,11) = AnSep{1,i}(size(AnSep{1,i},1),10);                
                                AnOver(i,12) = AnSep{1,i}(size(AnSep{1,i},1),11); 
                            end
                        end
                    elseif datenum(AnSep{1,i}(place,1))-...    % last vertix number and row above have the same time - for last in seq. late 2013 on?
                            datenum(AnSep{1,i}(place-1,1))<1.7361e-05
                        if place~=size(AnSep{1,i},1)
                            AnOver(i,1) = AnSep{1,i}(size(AnSep{i},1)-1,1);        % Start time
                            AnOver(i,2) = AnSep{1,i}(size(AnSep{i},1),1);        % End time
                            AnOver(i,3) = AnSep{1,i}(size(AnSep{i},1)-1,5);        % X co-ordinate start
                            AnOver(i,4) = AnSep{1,i}(size(AnSep{i},1)-1,6);        % Y co-ordinate start
                            AnOver(i,5) = AnSep{1,i}(size(AnSep{i},1),5);        % X co-ordinate end
                            AnOver(i,6) = AnSep{1,i}(size(AnSep{i},1),6);        % Y co-ordinate end
                            AnOver(i,7) = 2;                    % Analysis type
                            AnOver(i,8) = AnSep{1,i}(size(AnSep{i},1)-1,7);        % Scan speed
                            AnOver(i,9) = AnSep{1,i}(size(AnSep{i},1)-1,8);        % Repetition rate
                            AnOver(i,10) = AnSep{1,i}(size(AnSep{i},1)-1,9);      % Spot size                       
                            AnOver(i,11) = AnSep{1,i}(size(AnSep{1,i},1),10);                
                            AnOver(i,12) = AnSep{1,i}(size(AnSep{1,i},1),11);                                 
                        else
                            AnOver(i,1) = AnSep{1,i}(place-1,1);        % Start time
                            AnOver(i,2) = AnSep{1,i}(place,1);        % End time
                            AnOver(i,3) = AnSep{1,i}(place-1,5);        % X co-ordinate start
                            AnOver(i,4) = AnSep{1,i}(place-1,6);        % Y co-ordinate start
                            AnOver(i,5) = AnSep{1,i}(place,5);        % X co-ordinate end
                            AnOver(i,6) = AnSep{1,i}(place,6);        % Y co-ordinate end
                            AnOver(i,7) = 2;                    % Analysis type
                            AnOver(i,8) = AnSep{1,i}(place,7);        % Scan speed
                            AnOver(i,9) = AnSep{1,i}(place,8);        % Repetition rate
                            AnOver(i,10) = AnSep{1,i}(place,9);       % Spot size
                            if AnSep{1,i}(place,10)~=0
                                AnOver(i,11) = AnSep{1,i}(place,10);        % MFC1
                                AnOver(i,12) = AnSep{1,i}(place,11);        % MFC2
                            else
                                AnOver(i,11) = AnSep{1,i}(size(AnSep{1,i},1),10);                
                                AnOver(i,12) = AnSep{1,i}(size(AnSep{1,i},1),11); 
                            end
                        end
                    else
                        if size(AnSep{i},1)-place>1
                            AnOver(i,1) = AnSep{1,i}(place-1,1);        % Start time
                            AnOver(i,2) = AnSep{1,i}(place,1);        % End time
                            AnOver(i,3) = AnSep{1,i}(place-1,5);        % X co-ordinate start
                            AnOver(i,4) = AnSep{1,i}(place-1,6);        % Y co-ordinate start
                            AnOver(i,5) = AnSep{1,i}(place,5);        % X co-ordinate end
                            AnOver(i,6) = AnSep{1,i}(place,6);        % Y co-ordinate end
                            AnOver(i,7) = 2;                    % Analysis type
                            AnOver(i,8) = AnSep{1,i}(place,7);        % Scan speed
                            AnOver(i,9) = AnSep{1,i}(place,8);        % Repetition rate
                            AnOver(i,10) = AnSep{1,i}(place,9);       % Spot size
                            AnOver(i,11) = AnSep{1,i}(size(AnSep{i},1),10);        % MFC1
                            AnOver(i,12) = AnSep{1,i}(size(AnSep{i},1),11);        % MFC2
                        else
                        %%%%% not checked!!!! - is this necessary?
                            AnOver(i,1) = AnSep{1,i}(size(AnSep{i},1)-1,1);        % Start time
                            AnOver(i,2) = AnSep{1,i}(size(AnSep{i},1),1);        % End time
                            AnOver(i,3) = AnSep{1,i}(size(AnSep{i},1)-1,5);        % X co-ordinate start
                            AnOver(i,4) = AnSep{1,i}(size(AnSep{i},1)-1,6);        % Y co-ordinate start
                            AnOver(i,5) = AnSep{1,i}(size(AnSep{i},1),5);        % X co-ordinate end
                            AnOver(i,6) = AnSep{1,i}(size(AnSep{i},1),6);        % Y co-ordinate end
                            AnOver(i,7) = 2;                    % Analysis type
                            AnOver(i,8) = AnSep{1,i}(size(AnSep{i},1)-1,7);        % Scan speed
                            AnOver(i,9) = AnSep{1,i}(size(AnSep{i},1)-1,8);        % Repetition rate
                            AnOver(i,10) = AnSep{1,i}(size(AnSep{i},1)-1,9);       % Spot size
                            AnOver(i,11) = AnSep{1,i}(size(AnSep{i},1),10);        % MFC1
                            AnOver(i,12) = AnSep{1,i}(size(AnSep{i},1),11);        % MFC2
                        end
                    end
                elseif sum(~isnan(AnSep{1,i}(:,4)))==1  % If it contains one Vertix Number
                    AnOver(i,1) = AnSep{1,i}(size(AnSep{i},1)-1,1);        % Start time
                    AnOver(i,2) = AnSep{1,i}(size(AnSep{i},1),1);        % End time
                    AnOver(i,3) = AnSep{1,i}(size(AnSep{i},1)-1,5);        % X co-ordinate start
                    AnOver(i,4) = AnSep{1,i}(size(AnSep{i},1)-1,6);        % Y co-ordinate start
                    AnOver(i,5) = AnSep{1,i}(size(AnSep{i},1),5);        % X co-ordinate end
                    AnOver(i,6) = AnSep{1,i}(size(AnSep{i},1),6);        % Y co-ordinate end
                    AnOver(i,7) = 2;                    % Analysis type
                    AnOver(i,8) = AnSep{1,i}(5,7);        % Scan speed
                    AnOver(i,9) = sum(AnSep{1,i}(:,8));        % Repetition rate
                    AnOver(i,10) = AnSep{1,i}(5,9);       % Spot size
                    if isnan(AnOver(i,8))     % if scan speed was not in log file - bug in ?2010?
                        z = (((AnOver(i,5)-AnOver(i,3))^2+...
                            (AnOver(i,6)-AnOver(i,4))^2)^0.5/1000)/...
                            ((AnOver(i,2)-AnOver(i,1))*60*24);
                        z = roundsd(z,2);   % round to 0.1 mm/min
                        AnOver(i,8) = z*1000/60;
                    end
                    if AnSep{1,i}(5,10)~=0
                        AnOver(i,11) = AnSep{1,i}(5,10);        % MFC1
                        AnOver(i,12) = AnSep{1,i}(5,11);        % MFC2
                    else
                        AnOver(i,11) = AnSep{1,i}(size(AnSep{1,i},1),10);                
                        AnOver(i,12) = AnSep{1,i}(size(AnSep{1,i},1),11); 
                    end
                elseif sum(~isnan(AnSep{1,i}(:,4)))>2   % If it is a track with
                                                        % multiple points
                    if datenum(AnSep{1,i}(size(AnSep{i},1),1))-...
                            datenum(AnSep{1,i}(size(AnSep{i},1)-1,1))<3e-5
                        AnOver(i,1) = AnSep{1,i}(4,1);        % Start time
                        AnOver(i,2) = AnSep{1,i}(size(AnSep{i},1)-1,1);        % End time
                        AnOver(i,3) = AnSep{1,i}(4,5);        % X co-ordinate start
                        AnOver(i,4) = AnSep{1,i}(4,6);        % Y co-ordinate start
                        AnOver(i,5) = AnSep{1,i}(size(AnSep{i},1)-1,5);        % X co-ordinate end
                        AnOver(i,6) = AnSep{1,i}(size(AnSep{i},1)-1,6);        % Y co-ordinate end
                        AnOver(i,7) = 3;                    % Analysis type
                        AnOver(i,8) = AnSep{1,i}(5,7);        % Scan speed
                        AnOver(i,9) = AnSep{1,i}(5,8);        % Repetition rate
                        AnOver(i,10) = AnSep{1,i}(5,9);       % Spot size
                        if AnSep{1,i}(5,10)~=0
                            AnOver(i,11) = AnSep{1,i}(5,10);        % MFC1
                            AnOver(i,12) = AnSep{1,i}(5,11);        % MFC2
                        else
                            AnOver(i,11) = AnSep{1,i}(size(AnSep{1,i},1),10);                
                            AnOver(i,12) = AnSep{1,i}(size(AnSep{1,i},1),11); 
                        end
                    else        % Some files have final points >1 s apart.. (?)
                        AnOver(i,1) = AnSep{1,i}(4,1);        % Start time
                        AnOver(i,2) = AnSep{1,i}(size(AnSep{i},1),1);        % End time
                        AnOver(i,3) = AnSep{1,i}(4,5);        % X co-ordinate start
                        AnOver(i,4) = AnSep{1,i}(4,6);        % Y co-ordinate start
                        AnOver(i,5) = AnSep{1,i}(size(AnSep{i},1)-1,5);        % X co-ordinate end
                        AnOver(i,6) = AnSep{1,i}(size(AnSep{i},1)-1,6);        % Y co-ordinate end
                        AnOver(i,7) = 3;                    % Analysis type
                        AnOver(i,8) = AnSep{1,i}(5,7);        % Scan speed
                        AnOver(i,9) = AnSep{1,i}(5,8);        % Repetition rate
                        AnOver(i,10) = AnSep{1,i}(5,9);       % Spot size
                        if AnOver(i,10)~=0
                            AnOver(i,11) = AnSep{1,i}(5,10);        % MFC1
                            AnOver(i,12) = AnSep{1,i}(5,11);        % MFC2
                        else
                            AnOver(i,11) = AnSep{1,i}(size(AnSep{1,i},1),10);                
                            AnOver(i,12) = AnSep{1,i}(size(AnSep{1,i},1),11); 
                        end
                    end
                end
            elseif AnSep{1,i}(1,3)~=1       % If it is a map
                    AnOver(i,1) = AnSep{1,i}(4,1);        % Start time
                    AnOver(i,2) = AnSep{1,i}(5,1);        % End time
                    AnOver(i,3) = AnSep{1,i}(4,5);        % X co-ordinate start
                    AnOver(i,4) = AnSep{1,i}(4,6);        % Y co-ordinate start
                    AnOver(i,5) = AnSep{1,i}(5,5);        % X co-ordinate end
                    AnOver(i,6) = AnSep{1,i}(5,6);        % Y co-ordinate end
                    AnOver(i,7) = 4;                    % Analysis type
                    AnOver(i,8) = AnSep{1,i}(5,7);        % Scan speed
                    AnOver(i,9) = AnSep{1,i}(5,8);        % Repetition rate
                    AnOver(i,10) = AnSep{1,i}(5,9);       % Spot size
                    if AnOver(i,10)~=0
                        AnOver(i,11) = AnSep{1,i}(5,10);        % MFC1
                        AnOver(i,12) = AnSep{1,i}(5,11);        % MFC2
                    else
                        AnOver(i,11) = AnSep{1,i}(size(AnSep{1,i},1),10);                
                        AnOver(i,12) = AnSep{1,i}(size(AnSep{1,i},1),11); 
                    end
            end
        end
    else
        AnID(i-A,:) = [];   % Removes AnID entry if run was aborted
        A = A+1;
    end
end

AnID(isnan(AnOver(:,2)),:) = [];   %Delete NaN rows
AnOver = AnOver(~isnan(AnOver(:,2)),:);
% note, use column two as stopped analyses sometimes have time data in
% column 1, but not 2

% sort log files, in case more than one read in
sortLog = (1:1:size(AnOver,1))';
sortLog(:,2) = AnOver(:,1);
sortLog = sortrows(sortLog,2);
for i = 1:size(AnID,1)
    AnIDtemp(i,:) = AnID(sortLog(i,1),:);
end
AnOver = sortrows(AnOver,1);
AnID = AnIDtemp;

% a = 0;
% for i = 1:size(AnOver,1)
%     if 3600*24*(AnOver(i-a,2)-AnOver(i-a,1))<5
%         AnOver(i-a,:) = [];
%         AnID(i-a,:) = [];
%         a = a + 1;
%     end
% end

end





% probably unnecessary
function y=roundsd(x,n,method)
%ROUNDSD Round with fixed significant digits
%	ROUNDSD(X,N) rounds the elements of X towards the nearest number with
%	N significant digits.
%
%	ROUNDS(X,N,METHOD) uses following methods for rounding:
%		'round' - nearest (default)
%		'floor' - towards minus infinity
%		'ceil'  - towards infinity
%		'fix'   - towards zero
%
%	Examples:
%		roundsd(0.012345,3) returns 0.0123
%		roundsd(12345,2) returns 12000
%		roundsd(12.345,4,'ceil') returns 12.35
%
%	See also Matlab's function ROUND.
%
%	Author: François Beauducel <beauducel@ipgp.fr>
%	  Institut de Physique du Globe de Paris
%	Acknowledgments: Edward Zechmann, Daniel Armyr
%	Created: 2009-01-16
%	Updated: 2010-03-17

%	Copyright (c) 2010, François Beauducel, covered by BSD License.
%	All rights reserved.
%
%	Redistribution and use in source and binary forms, with or without 
%	modification, are permitted provided that the following conditions are 
%	met:
%
%	   * Redistributions of source code must retain the above copyright 
%	     notice, this list of conditions and the following disclaimer.
%	   * Redistributions in binary form must reproduce the above copyright 
%	     notice, this list of conditions and the following disclaimer in 
%	     the documentation and/or other materials provided with the distribution
%	                           
%	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%	ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%	LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%	CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%	SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%	INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%	CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%	ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%	POSSIBILITY OF SUCH DAMAGE.

error(nargchk(2,3,nargin))

if ~isnumeric(x)
		error('X argument must be numeric.')
end

if ~isnumeric(n) | numel(n) ~= 1 | n < 0 | mod(n,1) ~= 0
	error('N argument must be a scalar positive integer.')
end

opt = {'round','floor','ceil','fix'};

if nargin < 3
	method = opt{1};
else
	if ~ischar(method) | ~ismember(opt,method)
		error('METHOD argument is invalid.')
	end
end

og = 10.^(floor(log10(abs(x)) - n + 1));
y = feval(method,x./og).*og;
y(find(x==0)) = 0;
end