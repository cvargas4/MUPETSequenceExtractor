%% Script will take a CSV input from MUPET
% Principally the goal of the script is to identify where multisyllable
% sequences occur in a .wav of USV recordings. 
% MUPET can take multiple files as input, this script will separate out 
% when different files are occuring based on NaN values from ISIs.
% 
% Once multisyllable sequences are identified they are extracted from the
% .wav file of interest. Each sequence is then store as its own .wav file
% within a folder INSIDE THE CURRENTW ORKING DIRECTORY.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Location of CSV File From MUPET
path = '/Users/cesarvargas/Desktop/MouseExperiments/VocalRecordings/B96_10-28_FE_postsurgery/sonograms/repertoires/CSV';
myUSVFolder = path;
USVcsv = 'B96USVs_N60_syllable_sequence.csv';

%location of .wav file
USVFolder = '/Users/cesarvargas/Desktop/MouseExperiments/VocalRecordings/B96_10-28_FE_postsurgery/ch1';
USVFile = 'T0000001.wav';
USVwavFile = fullfile(USVFolder, USVFile);

%% Import USV Data
%Variable Used for USV Data from MSA (in order of column)
% {'Dataset', 'AnimalID', 
% 'SyllableNumber','SyllableStartTime', 'SyllableEndTime',
% 'RepertoireUnit'}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opts = detectImportOptions('B96USVs_N60_syllable_sequence.csv');
getvaropts(opts,{'dataset', 'repertoireFile', ...
'syllableNumber','syllableStartTime_sec_', 'syllableEndTime_sec_',...
'repertoireUnit_RU_Number'});
opts = setvartype(opts,{'dataset', 'repertoireFile', ...
'syllableNumber','syllableStartTime_sec_', 'syllableEndTime_sec_',...
'repertoireUnit_RU_Number'},...
    {'string', 'string', ...
    'int16','double', 'double', 'int16'});

USVstuff = readtable(USVcsv, opts);
USVstuff.Properties.VariableNames = {'Dataset', 'AnimalID',...
     'SyllableNumber','SyllableStartTime', 'SyllableEndTime', 'RepertoireUnit'};

numUSVs = reshape(1:height(USVstuff),1,[])';
numUSVs = array2table(numUSVs);
USVstuff = [USVstuff, numUSVs];
%Imports table of USV data from MUPET and adds variable names for convenience

clear opts

%% Find number of animals
%Finds number of animals or files in csv file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

animals = {1};
animalnum = 1;

animalID = USVstuff.AnimalID(1);
animals{animalnum,:} = animalID;

for i = 1:height(USVstuff)
    if strcmp(USVstuff.AnimalID(i), animalID) == 0
        animalID = USVstuff.AnimalID(i);
        animalnum = animalnum + 1;
        animals = vertcat(animals, animalID);
    end
end

clear animalID animalnum i ;

%% Find syllable lengths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('duration', 'var') == 1  
    USVstuff.duration = [];
    clear duration;
end

duration = reshape(1:height(USVstuff),1,[])';
sylstart = USVstuff.SyllableStartTime;
sylend = USVstuff.SyllableEndTime;

for i = 1:height(USVstuff)
    duration(i,:) = sylend(i) - sylstart(i);

end

duration = array2table(duration);
USVstuff = [USVstuff, duration];

clear sylend sylstart i;

%% Find and fill in ISIs
%Creates ISI column and calculates
%If ISI is negative, it went into new file/animal therefore is NaN
%last ISI will always be NaN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if re-running, will make sure ISI is overwritten properly
if exist('ISI', 'var') == 1  
    USVstuff.ISI = [];
    clear ISI;
end

ISI = reshape(1:height(USVstuff),1,[])';
sylstart = USVstuff.SyllableStartTime;
sylend = USVstuff.SyllableEndTime;

for i = 1:height(USVstuff)
    if i<height(USVstuff)
    ISI(i,:) = sylstart(i+1) - sylend(i);
        if ISI(i,:) < 0
           ISI(i,:) = NaN;
        end
    else
        ISI(i,:) = NaN;
    end
end
ISI = array2table(ISI);
USVstuff = [USVstuff, ISI];

clear sylend sylstart i;

%% Find Sequences
%From Chabout et al. 2015, if ISI is greater than 250ms its a new sequence
%Castelluci et al. 2018, groups are 100-275ms, bouts are >275ms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NaNloci = zeros(height(USVstuff), 1);
for i = 1:height(USVstuff)
    A = isnan(USVstuff.ISI(i));
    if A == 1
        NaNloci(i,:) = A;
    end
end
NaNloci = find(NaNloci);
clear i ans A

%% Creates a variable to allow for ease of identifying sequences
% The following loops will poulate a table with different values from
% USVstuff
% The meaning of each column of values is detailed in the comments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%finds sequence start times
seqinfo = zeros(height(USVstuff), 6);
% column 1 is when seq stopped, 0 indicates it ended, 1 indicates its still
%   going, 2 indicates a start
% column 2 is time in recording 
% column 3 is ISI
% column 4 is length of syllable
% column 5 is USV number in recording
% column 6 is how long sequence was
% column 7 is single or mulisyllable sequence, 0 is single, 1 is multi

for i = 1:length(animals)
    for j = 1:NaNloci(i,:)
        if USVstuff.ISI(j) <= 0.250 
            seqinfo(j,1) = 1;
            seqinfo(j,2) = USVstuff.SyllableStartTime(j);
            seqinfo(j,3) = USVstuff.ISI(j);
            seqinfo(j,4) = USVstuff.duration(j);
            seqinfo(j,5) = USVstuff.SyllableNumber(j);
        else
            seqinfo(j,1) = 0;
            seqinfo(j,2) = USVstuff.SyllableStartTime(j);
            seqinfo(j,3) = USVstuff.ISI(j);
            seqinfo(j,4) = USVstuff.duration(j);
            seqinfo(j,5) = USVstuff.SyllableNumber(j);
        end
    end
end
clear i j

%Determine when sequences started
for j = 1:length(seqinfo)
    if j>1
        if (seqinfo(j,3) <= 0.250) && (seqinfo((j-1),1) == 0)
            seqinfo(j,1) = 2;
        end
    end
end
clear j

%seqtime is when in recording sequence started
last = 0;
if seqinfo(1,1) == 0 && seqinfo(1,3) >= 0.25
    seqtime = sum(seqinfo(1,2)) + sum(seqinfo(1,3));
    seqinfo(1,6) = seqtime;
    seqinfo(1,6) = 0;
    last = 1;
end

for k = last:length(seqinfo)
    if k >=2
        if seqinfo(k,1) == 0
            if seqinfo(k-1,1) == 1
                for count = 1:last
                    seqtime = minus((seqinfo(k,2)), (seqinfo((last+1),2)));
                    seqinfo(k,6) = seqtime;
                    seqinfo(k,7) = 1;
                        %1 denotes sequences with multiple syllables
                end
            else
                    seqtime = minus((seqinfo(k,2)), (seqinfo((last+1),2)));
                    seqinfo(k,6) = seqtime;
                    seqinfo(k,7) = 0;
                %zero denotes sequences with one syllable
            end
            last = k;
        end
    end
end
clear k last ans count seqtime

%column 8 is sequence start in recording
%ONLY CONSIDERS SEQUENCES WITH MULTIPLE SYLLABLES
%if recording ends on single syllable this will also be ginored
for i = 1:(length(seqinfo)-1)
    if (i == 1) && (seqinfo((i+1),1) == 1)
        seqinfo(i,8) = seqinfo(i,2);
    elseif seqinfo(i,1) == 0 && seqinfo(i+1,1) == 1
        seqinfo(i,8) = seqinfo(i,2);
    end
end
clear i

%% convert recording time into ephys time
% Only do if needed and do it last
% This is done, because when I needed to do these experiments I recorded
% from electrodes several minutes before eliciting and recording USVs
% This section helps with then extracting the correct spike times in the 
% electrophysiological recordings (e.g. for PSTH).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timeadd = 900; %time of ephys in seconds when audio recording was started

for i = 1:NaNloci(1)
    if seqinfo(i,9) ~= 0
        seqinfo(i,9) = seqinfo(i,8) + timeadd;
    end
end
clear i

%% Pick out start and stop times of multi syllable sequences
% Finds the start ands top times of multisyllable sequences
% These are then stored in a new variable that is composed of samples instead
% of time. This is due to how the audioread and spectrogram features work in 
% the next section. 

% MAKE SURE TO CHANGE THE recordingRate vriable for your recordings.
% Jarvis lab standard is 250000. Due to Nyquist, this allows us to analyze
% USVs up to 125000Hz. Probably not great to record at frequencies 
% below 250000. Dodotronic Ultramic 384K records at 384000Hz.

% NaNloci is a variable that was created above. Change this value to decide
% which file/animal you want to analyze.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%column 1 is start
%column 2 is end
%column 3 is duration
startstops = zeros(length(seqinfo),3);
row = 1;
for i = 1:NaNloci(1)
    if seqinfo(i,1) == 2
        startstops(row,1) = seqinfo(i,2);
    elseif seqinfo(i,1) == 0 && seqinfo(i,7) == 1
        startstops(row,2) = seqinfo(i,2);
        startstops(row,3) = seqinfo(i,6);
        row = row +1;
    end
end
clear i row

recordingRate = 250000;
startstops = startstops(any(startstops,2),:);
startstopsFS = round(recordingRate*startstops(:,1:2)); %Recordings sample at 250000Hz
NumofSeq = length(startstops); %number of multisyllable sequences

%% Extract Sequences from Audio File
% Using sample-based sequence start times from above, finds these chunks
% in the audio file of the USV. This is the USV file listed at the top 
% of this script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

USVSeqFolder = pwd;
mkdir USVSequences
USVSeqPath = fullfile(pwd, 'USVSequences');

for i = 1:length(startstopsFS)
    samples = [startstopsFS(i,1),startstopsFS(i,2)];
    [USVseq, Fs] = audioread(USVwavFile, samples);
    index = num2str(i);
    disp(['Saving USVSequence_' index]);
    audiowrite([USVSeqPath '/' 'USVSequence_' index '.wav'], USVseq, Fs);
end
clear i ans index

%% Use this to visualize any of the saved frequencies
% rows of startstopsFS is each file, just change the row in the first line
% to visualize that sequence/file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

samples = [startstopsFS(20,1),startstopsFS(20,2)];
[USVseq, Fs] = audioread(USVwavFile, samples);
figure; 
spectrogram(USVseq,125,100,100,Fs,'yaxis', 'MinThreshold', -100);
% view(3) %for 3D viz use this