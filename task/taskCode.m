%MP 2021

%This is the master script used for presenting behavioral stimuli and
%recording responses in the memory-guided saccade task described in
%Panichello, Jonikaitis, Oh, Zhu, Trepka, & Moore, 2024

%Due to the very large number of dependencies (eg., eyelink), it is not practical to run
%this code on machines outside the experimental rig, but we hope it
%nevertheless makes our method more transparent 

clear;
close all;
sca;
rng('shuffle');
KbName('unifykeynames');
%UDPReceiver()

%% define parameters
opts.monkey = 'Joplin';
opts.saveData = true;
opts.doCalibrate = false;
opts.repeatError = false;

opts.screenDist = 30; %inches
opts.screenDiag = 30; %inches [23 30]

opts.fixDvaX = .7;
opts.fixDvaY = .7;

opts.sampleDVA = .7;
opts.sampleAngles = deg2rad(0:45:315);%(2*pi/8):(2*pi/8):(2*pi); %[0 deg2rad(180)];
opts.sampleEcc = 6;

opts.FixAcquireWindowRadius = 3/180*pi; %dva (in radians)
opts.FixHoldWindowRadius = 2/180*pi; %dva (in radians)
opts.ResponseBarrierRadius = 3.5/180*pi; %dva (in radians)

opts.maxAcquireFixTime = 4; %seconds
opts.fixGraceTime = 0;
opts.delay1range = [.4 .4]; %[1.7 3.5];
opts.sampleTime = .050;
opts.delay2range = [1.4 1.6];
opts.targetDelay = .1;
opts.sampleLuminance = 1;
opts.maxResponseTime = 3;
opts.targetHoldTime = 0.000;
opts.FixColor = [0 0 0]*255;
opts.ITI = 1;
opts.breakTO = 0;
opts.noinitTO = 0;
opts.errorTO = 0.2; %just to give him time to see feedback
opts.backgroundColor = [.5 .5 .5]*255;
opts.pauseColor = [.5 .5 .5]*255;
opts.errorColor = [.5 .5 .5]*255;

opts.pMGS = 1;

opts.angLocStart = 30;
opts.angLocWindow = 100;
opts.angLocTemp = .9;

opts.errorWindow = deg2rad(70); %= +/- 35 degrees
opts.doReward = 1;
opts.nDrops = 3;
opts.rewardPulse = 250; %started at 500
opts.interRewardDelay = 160;

opts.exitKey = 'escape'; %for exiting
opts.pauseKey = 'space'; %for pausing
opts.calKey = 'c'; %for calibration

bhv.opts = opts;
bhv.EyeUsed = 0;

%%

%Encode variables
%nums 0 thru 9 reserved for trial 
% 256 thousands hundreds tens ones 256 AorB ChangeNoChange 256 
opts.Encodes.START_TRIAL = 255;

opts.Encodes.FIXATE_ON = 10;
opts.Encodes.FIXATE_ACQUIRED = 11;
opts.Encodes.CUE_ON = 12;
opts.Encodes.DELAY_START = 13;
opts.Encodes.TARGET_ON = 14;
opts.Encodes.TARGET_ACQUIRED = 15;
opts.Encodes.REWARD_ON = 16;

opts.Encodes.CORRECT_TRIAL = 19;
opts.Encodes.INCORRECT_TRIAL = 21;
opts.Encodes.NO_INIT_TRIAL = 22;
opts.Encodes.FIX_BREAK_TRIAL = 23;


%% Create save file
bhv.monkey = opts.monkey;
bhv.date = datestr(now, 'yymmdd');
bhv.weekday = weekday(datenum(bhv.date, 'yymmdd'));


bhv.filename = sprintf('./data/mgs/%s_%s_00_bhv.mat', bhv.monkey, bhv.date);
count = 0;
while exist(bhv.filename, 'file') && (count < 100)
    count = count + 1;
    bhv.filename = sprintf('./data/mgs/%s_%s_%02.0f_bhv.mat', bhv.monkey, bhv.date, count);
end
if count >= 100, error('Couldn''t find a unique filename for today.'); end

%% initialize ptb
windowPtr = 2;
Screen('Preference', 'SkipSyncTests', 1); %need to troubleshoot
[windowPtr, bhv.screen.rect] = Screen('OpenWindow', windowPtr, opts.backgroundColor);

pauseKey = KbName(opts.pauseKey);
exitKey = KbName(opts.exitKey);
calKey = KbName(opts.calKey);

%% initialize reward digital i/o
d = daq("ni");
addoutput(d,"Dev1","port0/line0","Digital")

%% initialize trig digital i/o
trig = daq("ni");
addoutput(trig,"Dev1","port0/line1","Digital")
addoutput(trig,"Dev1","port0/line2","Digital")
addoutput(trig,"Dev1","port0/line3","Digital")
addoutput(trig,"Dev1","port0/line4","Digital")
addoutput(trig,"Dev1","port0/line5","Digital")
addoutput(trig,"Dev1","port0/line6","Digital")
addoutput(trig,"Dev1","port0/line7","Digital")
addoutput(trig,"Dev1","port1/line1","Digital")

%% key conversions
bhv.func.DVAToIn = @(a) tan(a)*opts.screenDist;
bhv.func.InToDVA = @(x) atan(x/opts.screenDist);
bhv.func.DVAToPIX = @(a, res) tan(a)*opts.screenDist*res;
bhv.func.PIXToDVA = @(x, res) atan(x/opts.screenDist/res);

%% determine screen and stimulus parameters
%Determine screen center (in pixels)
bhv.screen.center = [bhv.screen.rect(3)/2 bhv.screen.rect(4)/2];

%Screen parameters
bhv.screen.width  = bhv.opts.screenDiag * cos( atan( bhv.screen.rect(4) / bhv.screen.rect(3) ) ); %in inches
bhv.screen.height = bhv.opts.screenDiag * sin( atan( bhv.screen.rect(4) / bhv.screen.rect(3) ) ); %in inches
bhv.screen.horizRes = bhv.screen.rect(3)./ bhv.screen.width; %in pix/inch
bhv.screen.vertRes = bhv.screen.rect(4) ./ bhv.screen.height; %in pix/inch

bhv.fix.rect = [0, 0, bhv.func.DVAToPIX(deg2rad(opts.fixDvaX),bhv.screen.horizRes), bhv.func.DVAToPIX(deg2rad(opts.fixDvaY),bhv.screen.vertRes) ];
bhv.fix.rect = bhv.fix.rect + [bhv.screen.rect(3:4) bhv.screen.rect(3:4)]/2 - [bhv.fix.rect(3:4) bhv.fix.rect(3:4)]/2;
[FixRect_x, FixRect_y] = RectCenter(bhv.fix.rect);

bhv.sample.rect = [0, 0, bhv.func.DVAToPIX(deg2rad(opts.sampleDVA),bhv.screen.horizRes), bhv.func.DVAToPIX(deg2rad(opts.sampleDVA),bhv.screen.vertRes) ];


bhv.FixAcquireWindowRadius = [bhv.func.DVAToPIX(opts.FixAcquireWindowRadius, bhv.screen.horizRes); bhv.func.DVAToPIX(opts.FixAcquireWindowRadius, bhv.screen.vertRes)] ;
bhv.FixHoldWindowRadius = [bhv.func.DVAToPIX(opts.FixHoldWindowRadius, bhv.screen.horizRes); bhv.func.DVAToPIX(opts.FixHoldWindowRadius, bhv.screen.vertRes)] ;
bhv.ResponseBarrierRadius = [bhv.func.DVAToPIX(opts.ResponseBarrierRadius, bhv.screen.horizRes); bhv.func.DVAToPIX(opts.ResponseBarrierRadius, bhv.screen.vertRes)] ;

bhv.pauseImage = imread('./stimuli/beach.jpg');
bhv.pauseImageTexture = Screen('makeTexture',windowPtr,bhv.pauseImage);
%% prepare eyelink
el=EyelinkInitDefaults(windowPtr);
el.callback = [];
el.backgroundcolour = opts.backgroundColor;
el.calibrationtargetsize = 2;
dummymode = 0;


% Initialization of the connection with the Eyelink Gazetracker.
% exit program if this fails.
if ~EyelinkInit(dummymode, 1)
    fprintf('Eyelink Init aborted.\n');
    cleanup;  % cleanup function
    return;
end
Eyelink('command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA');

%enable on-line drift correction
Eyelink('command', 'drift_correct_cr_disable = OFF');
Eyelink('command', 'online_dcorr_refposn 512,384');
Eyelink('command', 'online_dcorr_button = ON');

%calibrate
if opts.doCalibrate
    calibrate( bhv,el,d,windowPtr);
end

Eyelink('StartRecording');


%% execute trial loop
finished = false;
iTrial = 0;



while ~finished
    %%ITI
    Eyelink('command','clear_screen 0');
    WaitSecs(opts.ITI);
    [x, y] = ellipse_coords(bhv.FixAcquireWindowRadius(1), bhv.FixAcquireWindowRadius(2), 0, FixRect_x, FixRect_y,'r');
    for i = 2:4:301
        Eyelink('Command',sprintf('draw_line %f %f %f %f 15',x(i-1), y(i-1), x(i), y(i) ));
    end
    [x, y] = ellipse_coords(bhv.ResponseBarrierRadius(1), bhv.ResponseBarrierRadius(2), 0, FixRect_x, FixRect_y,'r');
    for i = 2:4:301
        Eyelink('Command',sprintf('draw_line %f %f %f %f 15',x(i-1), y(i-1), x(i), y(i) ));
    end


    %initialize trial
    tempTrial = struct('stopCondition',nan, 'delay1',nan,'delay2',nan, ...
        'HoldOnset_flipTime',nan,'HoldOnset_onsetTime',nan, ...
        'fixOnset_onsetTime',nan, ...
        'fixDotOffset_flipTime', nan, 'fixDotOffset_onsetTime', nan, ...
        'sampleOnset_onsetTime', nan, 'sampleOffset_onsetTime', nan, ...
        'targetOnset_onsetTime', nan, 'targetAcquiredOnset_onsetTime', nan, ...
        'stimXcenter', nan, 'stimYcenter', nan, ...
        'aquisitionTime', nan, ...
        'releaseTime',nan,'holdDuration',nan,'rt',nan, ...
        'responseCoordinate', nan, 'responseAngle', nan,'angError',nan,'isMGS',nan);
    iTrial = iTrial+1;
    tempTrial.iTrial = iTrial;

    tempTrial.delay1 = rand.*diff(opts.delay1range) + min(opts.delay1range);
    tempTrial.delay2 = rand.*diff(opts.delay2range) + min(opts.delay2range);

    tempTrial.isMGS = opts.pMGS >= rand;
    %     if iTrial == 1
    %         tempTrial.isMGS = false;
    %     else
    %         nCompletedTrials = sum(abs([bhv.Trials.stopCondition])==1);
    %         if mod(nCompletedTrials,20) <= 1
    %             tempTrial.isMGS = false;
    %         elseif mod(nCompletedTrials,20) > 1
    %             tempTrial.isMGS = true;
    %         end
    %     end

%     if iTrial <= opts.angLocStart
%         tempTrial.angleID = datasample(1:numel(opts.sampleAngles),1);
%         tempTrial.sampleAngle = opts.sampleAngles(tempTrial.angleID);
%     else
% 
%         idxStart = iTrial-opts.angLocWindow;
%         if idxStart <= 0
%             idx = 1:(iTrial-1);
%         else
%             idx = idxStart:(iTrial-1);
%         end
% 
%         nAngs = numel(opts.sampleAngles);
%         sc = [bhv.Trials(idx).stopCondition];
%         c = [bhv.Trials(idx).angleID];
%         pCorr = nan(1,nAngs);
% 
%         for iAng = 1:nAngs
%             %pCorr(iAng) = sum(sc(c==iAng)==1)./sum(abs(sc(c==iAng))==1);
%             pCorr(iAng) = sum(sc(c==iAng)==1)./sum(c==iAng); %correct for bf too
%         end
% 
%         pCorr(isnan(pCorr)) = 0;
%         V = 1-pCorr;
%         p = exp(V./opts.angLocTemp) ./ sum(exp(V./opts.angLocTemp));
%         tempTrial.angleID = discretesample(p, 1);
%         tempTrial.sampleAngle = opts.sampleAngles(tempTrial.angleID);
%     end
    if (opts.repeatError) && (iTrial > 5) && (bhv.Trials(end).stopCondition ~= 1) %last trial was error
        tempTrial.angleID = bhv.Trials(end).angleID;
        tempTrial.sampleAngle = opts.sampleAngles(tempTrial.angleID);
    else
        tempTrial.angleID = datasample(1:numel(opts.sampleAngles),1);
        tempTrial.sampleAngle = opts.sampleAngles(tempTrial.angleID);
    end

    tempTrial.sampleEcc = opts.sampleEcc;
    [tempTrial.sample_location(1), tempTrial.sample_location(2)] = pol2cart(tempTrial.sampleAngle, tempTrial.sampleEcc);

    tempStimRect = bhv.sample.rect;
    tempTrial.stimXcenter = bhv.screen.center(1) + bhv.func.DVAToPIX(deg2rad(tempTrial.sample_location(1)),bhv.screen.horizRes);
    tempTrial.stimYcenter = bhv.screen.center(2) + bhv.func.DVAToPIX(deg2rad(tempTrial.sample_location(2)),bhv.screen.vertRes);
    tempStimRect = CenterRectOnPoint(tempStimRect,tempTrial.stimXcenter,tempTrial.stimYcenter);

    [thousand, hundred, ten, one] = grabDigits(iTrial);
    sendTrig(trig, opts.Encodes.START_TRIAL); pause(.010);
    sendTrig(trig, thousand); pause(.010);
    sendTrig(trig, hundred);  pause(.010);
    sendTrig(trig, ten); pause(.010);
    sendTrig(trig, one); pause(.010);
    sendTrig(trig, opts.Encodes.START_TRIAL); pause(.010);
    sendTrig(trig, 0); pause(.010);

    %begin stimulus presentation
    %present fixation at leverpress
    quitTask = false;
    [~,~,keyCode] = KbCheck;
    if keyCode(exitKey)
        quitTask = true;
        break
    elseif keyCode(calKey)
        calibrate( bhv,el,d,windowPtr);
    elseif keyCode(pauseKey)
        %Screen('DrawTexture', windowPtr, bhv.pauseImageTexture);
        %Screen('Flip', windowPtr)
        Screen('FillRect', windowPtr, 255*[0 .5 .5]);
        Screen('Flip', windowPtr);
        disp('task paused');
        WaitSecs(1);
        KbWait;
        Screen('FillRect', windowPtr, opts.backgroundColor);
        Screen('Flip', windowPtr);
    end
    Screen('FillRect', windowPtr, opts.backgroundColor);
    Screen('FillOval', windowPtr, opts.FixColor, bhv.fix.rect);
    sendTrig(trig, opts.Encodes.FIXATE_ON); 
    [tempTrial.HoldOnset_flipTime, tempTrial.HoldOnset_onsetTime] = Screen('Flip', windowPtr);
    sendTrig(trig, 0);

    Eyelink('command',sprintf('draw_cross %d %d 15',FixRect_x(:), FixRect_y(:)));
    %acquire fixation
    fixAcquired = false;
    while (~fixAcquired) && (GetSecs < (tempTrial.HoldOnset_flipTime + opts.maxAcquireFixTime))
        %check fix
        [mx, my] = EyelinkGetEyeSample(el, bhv.EyeUsed);
        if ~isnan(mx) && ~isnan(my)
            fixAcquired = PointInEllipse([mx my], [FixRect_x(:) FixRect_y(:)], [bhv.FixAcquireWindowRadius(1) bhv.FixAcquireWindowRadius(2)]);
        end
    end
    if ~fixAcquired %trial not initiated
        Screen('FillRect', windowPtr, opts.backgroundColor);
        tempTrial.releaseTime = Screen('Flip', windowPtr);
        sendTrig(trig, opts.Encodes.NO_INIT_TRIAL); pause(.010);
        sendTrig(trig, 0);
        tempTrial.stopCondition = -2; %no init
        bhv.Trials(iTrial) = tempTrial;
        Eyelink('command',sprintf('draw_cross %d %d 0',FixRect_x(:), FixRect_y(:)));
        WaitSecs(opts.noinitTO);
        continue;
    end

    sendTrig(trig, opts.Encodes.FIXATE_ACQUIRED); pause(.010);
    sendTrig(trig, 0);
    %hold fixation for delay1
    fixBroken = false;
    tempTrial.fixOnset_onsetTime = GetSecs;
    while (GetSecs < (tempTrial.fixOnset_onsetTime + tempTrial.delay1)) && ~fixBroken
        [mx, my] = EyelinkGetEyeSample(el, bhv.EyeUsed);
        if ~isnan(mx) && ~isnan(my)
            fixBroken = ~PointInEllipse([mx my], [FixRect_x(:) FixRect_y(:)], [bhv.FixAcquireWindowRadius(1) bhv.FixAcquireWindowRadius(2)]);
        end
        if fixBroken
            breakTime = GetSecs;
            while (GetSecs < (tempTrial.fixOnset_onsetTime + tempTrial.delay1)) && (GetSecs<(breakTime+opts.fixGraceTime))
                [mx, my] = EyelinkGetEyeSample(el, bhv.EyeUsed);
                if ~isnan(mx) && ~isnan(my)
                    fixBroken = ~PointInEllipse([mx my], [FixRect_x(:) FixRect_y(:)], [bhv.FixAcquireWindowRadius(1) bhv.FixAcquireWindowRadius(2)]);
                end
                if ~fixBroken
                    break
                end
            end
        end
    end
    if fixBroken
        tempTrial.stopCondition = -4;
        tempTrial.releaseTime = Screen('Flip', windowPtr);
        tempTrial.holdDuration = tempTrial.releaseTime - tempTrial.fixOnset_onsetTime;
        sendTrig(trig, opts.Encodes.FIX_BREAK_TRIAL); pause(.010);
        sendTrig(trig, 0);
        bhv.Trials(iTrial) = tempTrial;
        Eyelink('command',sprintf('draw_cross %d %d 0',FixRect_x(:), FixRect_y(:)));
        WaitSecs(opts.breakTO);
        continue;
    end

    %present stimulus
    fixBroken = false;
    Screen('FillRect', windowPtr, opts.backgroundColor);
    Screen('FillOval', windowPtr, opts.FixColor, bhv.fix.rect);
    Screen('FillOval', windowPtr,([1 1 1]*255)*opts.sampleLuminance + opts.backgroundColor*(1-opts.sampleLuminance),tempStimRect);
    Eyelink('command',sprintf('draw_cross %d %d 15',tempTrial.stimXcenter, tempTrial.stimYcenter));
    sendTrig(trig, opts.Encodes.CUE_ON);
    tempTrial.sampleOnset_onsetTime = Screen('Flip', windowPtr);
    sendTrig(trig, 0);
    while (GetSecs < (tempTrial.sampleOnset_onsetTime + opts.sampleTime)) && ~fixBroken
        [mx, my] = EyelinkGetEyeSample(el, bhv.EyeUsed);
        if ~isnan(mx) && ~isnan(my)
            fixBroken = ~PointInEllipse([mx my], [FixRect_x(:) FixRect_y(:)], [bhv.FixHoldWindowRadius(1) bhv.FixHoldWindowRadius(2)]);
        end
        if fixBroken
            breakTime = GetSecs;
            while (GetSecs < (tempTrial.sampleOnset_onsetTime + opts.sampleTime)) && (GetSecs<(breakTime+opts.fixGraceTime))
                [mx, my] = EyelinkGetEyeSample(el, bhv.EyeUsed);
                if ~isnan(mx) && ~isnan(my)
                    fixBroken = ~PointInEllipse([mx my], [FixRect_x(:) FixRect_y(:)], [bhv.FixAcquireWindowRadius(1) bhv.FixAcquireWindowRadius(2)]);
                end
                if ~fixBroken
                    break
                end
            end
        end
    end
    if fixBroken
        tempTrial.stopCondition = -4;
        tempTrial.releaseTime = Screen('Flip', windowPtr);
        tempTrial.holdDuration = tempTrial.releaseTime - tempTrial.fixOnset_onsetTime;
        sendTrig(trig, opts.Encodes.FIX_BREAK_TRIAL); pause(.010);
        sendTrig(trig, 0);
        bhv.Trials(iTrial) = tempTrial;
        Eyelink('command',sprintf('draw_cross %d %d 0',FixRect_x(:), FixRect_y(:)));
        WaitSecs(opts.breakTO);
        continue;
    end

    %hold fixation for delay2
    fixBroken = false;
    Screen('FillRect', windowPtr, opts.backgroundColor);
    Screen('FillOval', windowPtr, opts.FixColor, bhv.fix.rect);
    Eyelink('command',sprintf('draw_cross %d %d 0',tempTrial.stimXcenter, tempTrial.stimYcenter));
    sendTrig(trig, opts.Encodes.DELAY_START);
    tempTrial.sampleOffset_onsetTime = Screen('Flip', windowPtr);
    sendTrig(trig, 0);
    while (GetSecs < (tempTrial.sampleOffset_onsetTime + tempTrial.delay2)) && ~fixBroken
        [mx, my] = EyelinkGetEyeSample(el, bhv.EyeUsed);
        if ~isnan(mx) && ~isnan(my)
            fixBroken = ~PointInEllipse([mx my], [FixRect_x(:) FixRect_y(:)], [bhv.FixHoldWindowRadius(1) bhv.FixHoldWindowRadius(2)]);
        end
        if fixBroken
            breakTime = GetSecs;
            while (GetSecs < (tempTrial.sampleOffset_onsetTime + tempTrial.delay2)) && (GetSecs<(breakTime+opts.fixGraceTime))
                [mx, my] = EyelinkGetEyeSample(el, bhv.EyeUsed);
                if ~isnan(mx) && ~isnan(my)
                    fixBroken = ~PointInEllipse([mx my], [FixRect_x(:) FixRect_y(:)], [bhv.FixAcquireWindowRadius(1) bhv.FixAcquireWindowRadius(2)]);
                end
                if ~fixBroken
                    break
                end
            end
        end
    end
    if fixBroken
        tempTrial.stopCondition = -4;
        tempTrial.releaseTime = Screen('Flip', windowPtr);
        tempTrial.holdDuration = tempTrial.releaseTime - tempTrial.fixOnset_onsetTime;
        sendTrig(trig, opts.Encodes.FIX_BREAK_TRIAL); pause(.010);
        sendTrig(trig, 0);
        bhv.Trials(iTrial) = tempTrial;
        Eyelink('command',sprintf('draw_cross %d %d 0',FixRect_x(:), FixRect_y(:)));
        WaitSecs(opts.breakTO);
        continue;
    end

    %turn off fixation and wait for fixation break
    fixBroken = false;
    targetShown = false;
    Screen('FillRect', windowPtr, opts.backgroundColor);
    %Screen('FillOval', windowPtr,[1 1 1]*255,tempStimRect);
    sendTrig(trig, opts.Encodes.TARGET_ON); %not really target on
    [tempTrial.fixDotOffset_flipTime, tempTrial.fixDotOffset_onsetTime] = Screen('Flip', windowPtr);
    sendTrig(trig, 0);
    Eyelink('command',sprintf('draw_cross %d %d 0',FixRect_x(:), FixRect_y(:)));
    while (GetSecs < (tempTrial.fixDotOffset_flipTime + opts.maxResponseTime)) && ~fixBroken
        if (GetSecs >= (tempTrial.fixDotOffset_flipTime + opts.targetDelay)) && (~targetShown) && ~tempTrial.isMGS
            Screen('FillOval', windowPtr,[0 0 0]*255,tempStimRect);
            tempTrial.targetOnset_onsetTime = Screen('Flip', windowPtr);
            targetShown = true;
            Eyelink('command',sprintf('draw_cross %d %d 15',tempTrial.stimXcenter, tempTrial.stimYcenter));
        end
        [mx, my] = EyelinkGetEyeSample(el, bhv.EyeUsed);
        if ~isnan(mx) && ~isnan(my)
            fixBroken = ~PointInEllipse([mx my], [FixRect_x(:) FixRect_y(:)], [bhv.ResponseBarrierRadius(1) bhv.ResponseBarrierRadius(2)]);
        end
    end
    if ~fixBroken
        Screen('FillRect', windowPtr, opts.backgroundColor);
        tempTrial.targetAcquiredOnset_onsetTime = NaN;
        Screen('Flip', windowPtr);
        tempTrial.stopCondition = -5; %never left fixation at end of trial
        bhv.Trials(iTrial) = tempTrial;
        Eyelink('command',sprintf('draw_cross %d %d 0',FixRect_x(:), FixRect_y(:)));
        WaitSecs(opts.noinitTO);
        continue;
    end

    %reveal target (if not already) and compute angular aerror
    if ~targetShown
        Screen('FillOval', windowPtr,[0 0 0]*255,tempStimRect);
        tempTrial.targetOnset_onsetTime = Screen('Flip', windowPtr);
        targetShown = true;
        Eyelink('command',sprintf('draw_cross %d %d 15',tempTrial.stimXcenter, tempTrial.stimYcenter));
    end

    tempTrial.responseCoordinate = [mx my] - bhv.screen.center;
    [tempTrial.responseAngle, ~] = cart2pol(tempTrial.responseCoordinate(1),tempTrial.responseCoordinate(2));
    tempTrial.responseAngle = wrapTo2Pi(tempTrial.responseAngle);
    disp([tempTrial.sampleAngle, tempTrial.responseAngle])
    tempTrial.angError = wrapToPi(tempTrial.sampleAngle - tempTrial.responseAngle);

    if abs(tempTrial.angError) <= opts.errorWindow/2
        %if we got this far, give reward
        sendTrig(trig, opts.Encodes.CORRECT_TRIAL ); pause(.010);
        sendTrig(trig, 0);
        tempTrial.stopCondition = 1;
        if opts.doReward

            for cur_drop = 1:opts.nDrops
                %Give reward
                write(d,1);
                WaitSecs(opts.rewardPulse/1000);
                write(d,0);
                WaitSecs(opts.interRewardDelay/1000);
            end

        end
    else
        %angular error was too large
        WaitSecs(opts.errorTO);
        Screen('FillRect', windowPtr, opts.backgroundColor);
        Screen('Flip', windowPtr);
        sendTrig(trig, opts.Encodes.INCORRECT_TRIAL ); pause(.010);
        sendTrig(trig, 0);
        tempTrial.stopCondition = -1;
        bhv.Trials(iTrial) = tempTrial;
        Eyelink('command',sprintf('draw_cross %d %d 0',FixRect_x(:), FixRect_y(:)));
        continue;
    end

    sendTrig(trig, 0);

    %return screen to blank
    Screen('FillRect', windowPtr, opts.backgroundColor);
    Screen('Flip', windowPtr);

    bhv.Trials(iTrial) = tempTrial;


end
%% save off results and close
sca;

if opts.saveData
    bhv.juiceEarned = input('juice earned:');
    bhv.juiceSupp   = input('juice supp:');
    bhv.totalJuice  = bhv.juiceEarned+bhv.juiceSupp;
    bhv.weight      = input('weight:');
    bhv.notes   = input('notes:','s');
    save(bhv.filename,'bhv');
end

