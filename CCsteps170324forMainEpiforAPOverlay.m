134%Script to analyse passive and active properties of neurons (Biophysical
%Fingerprint) for EpimiRNA. Reads .mat files containing time, I and V
%trace, where 17 current steps were applied to a neuron

clear all


%Load data
directory = dir('/Users/gmorris/Desktop/IC/*.mat'); %Load any mat files from within 'Input' folder
results = cell(length(directory)+1,9);
results{1,1} = 'File Name';
results{1,2} = 'RMP';
results{1,3} = 'Ri';
results{1,4} = 'First AP amplitude';
results{1,5} = 'First AP rise time';
results{1,6} = 'First AP decay time';
results{1,7} = 'First AP half-width';
results{1,8} = 'Max Rising Slope';
results{1,9} = 'Max Decay Slope';
results{1,10} = 'Firing threshold';
results{1,11} = 'Threshold current pulse';
results{1,12} = 'Max firing freq';

ThresholdAPStorage = [];

    
for I = 1:length(directory)%does the whole script for each file in the input directory
    
    eval(['load /Users/gmorris/Desktop/IC/' directory(I).name]); %opens each mat file as 2 arrays. T = time,  Y = current (col 1) and Vm (col 2)
    filename=directory(I).name(1:(end-4));
    rawfig = figure;
    tempy=[];
    tempy(:,1) = Y(:,2); %use these 3 lines to switch I and V if necessary
    tempy(:,2) = Y(:,1); %%%
    Y = tempy; %%%
    subplot (2,1,1)
    plot(T,Y(:,1))
    subplot(2,1,2)
    plot(T,Y(:,2)) %these 5 lines plot the raw t vs I and t vs V traces
    SR = round(1/(T(2,1)-T(1,1)));
    
    %Having opened and plotted the data, need to pick out regions of
    %current injection and make them into separate frames
    
    Idy = gradient(Y(:,1)); %creates new trace with the gradient of the I trace, to pick out injections
    Iinj1 = find(Idy<-50,1); %finds the start of the first current injection %%%Changed from -0.05
    Injections=zeros(17,1);
    
    for i=1:17
        Injections(i,1) = Iinj1 + SR*(i-1);
       rawfig;
       subplot(2,1,1)
        hold on
        h = line([T(Injections(i,1)) T(Injections(i,1))] , [min(Y(:,1)) max(Y(:,1))]);
        set(h,'color','r')
    end %This loop calculates the start time of all of the other injections, based on there being 17 injections which are 1s apart
    
    frames = struct; %create a new structure into which all data can be placed, sorted by the injection with which they are associated
    for i = 1:17
        num = num2str(i);
        name = ['frame' num];
        Tstart = Injections(i,1)-(0.1*SR);
        Tend = Injections(i,1)+(0.4*SR);
        frames.(name)(:,1) = T(Tstart:Tend);
        frames.(name)(:,2) = Y(Tstart:Tend,1);
        frames.(name)(:,3) = Y(Tstart:Tend,2);
    end %the loop plots time into col 1, I into 2 and Vm into 3 for 100 ms either side of the pulse
    
        Vtraces = figure; %plots all of the Vms on different subplots
    for i = 1:17
        istring = num2str(i);
        framename= ['frame' istring]; 
        subplot(4,6,i);        
        plot(frames.(framename)(:,1),frames.(framename)(:,3))
        title(framename)
    end

        Itraces = figure; %plots all of the injected currents on different subplots
    for i = 1:17
        istring = num2str(i);
        framename= ['frame' istring]; 
        subplot(4,6,i);
        plot(frames.(framename)(:,1),frames.(framename)(:,2))
        title(framename)
    end
    

    
    
     %Start by detecting%APs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     AP = struct; %Makes a new structure to store the APs in
        i=1;
        i2=1;
        i3=1;   

    for i = 1:17
    DetermineThreshold = -15;
    istring = num2str(i);
    framename= ['frame' istring];%for each frame, go through each point and find peaks above 0
    APframename= ['APsforframe' istring];
    AP.(APframename) = [];
        for i2= 251:length(frames.(framename)(:,1))-250
            if all([frames.(framename)(i2,3)>DetermineThreshold , frames.(framename)(i2,3)== max(frames.(framename)(((i2-(0.001*SR)):(i2+(0.001*SR))),3)) , frames.(framename)(i2,3)>frames.(framename)(i2-1,3) , frames.(framename)(i2,3)>=frames.(framename)(i2+1,3)]) == 1 %if a point matches all criteria.....                            
                AP.(APframename)(i3,1)=frames.(framename)(i2,1); %.....then its time is stored
                figure(Vtraces);
                hold on
                subplot(4,6,i);
                APmarker = line ([frames.(framename)(i2,1) frames.(framename)(i2,1)] , [min(frames.(framename)(:,3)) max(frames.(framename)(:,3))]); %each detected AP is plotted as a red line onto the V trace for visual checking
                set(APmarker,'color','r')
                i3=i3+1;
            end
        end
        i3=1;  
        i2=1;
    end


    %----------------------Passive properties-------------------------------------------------
    %Pick out the relevant frames (assumes that those with no APs were subthreshold I injections)
    i=1;
    i2=1;
    PassiveFrames=[];
    for i = 1:17 %Goes through each frame
        istring = num2str(i);
        APframename = ['APsforframe' istring];
        if isempty(AP.(APframename)) == 1; %If there are no APs detected from that frame then assumes passive
            PassiveFrames(i2,1)=i;
            i2=i2+1;
        end
    end
    
    
    
    
    
    %Resting membrane potential (from all frames)
        figure(Vtraces);
        hold on
        rmpvalues=[];
        %for each frame, the RMP is calculated as the mean of the first 50 ms
        %of data
        for i = 1:17
            istring = num2str(i);
            framename= ['frame' istring];
            temprmp = mean(frames.(framename)(1:(round(0.05*SR)),3));
            subplot(4,6,i);                       
            rmpmarker = line([min(frames.(framename)(:,1)) max(frames.(framename)(:,1))] ,[temprmp temprmp]);
            set(rmpmarker,'color','r','LineStyle','--')
            rmpvalues(i,1)=temprmp;
            i=i+1;
        end
         rmpplot = figure;
         scatter(1:17,rmpvalues)
         xlabel('Frame Number')
         ylabel('Resting Membrane Potential (mV)')
         xlim([0 17])


         output_RMP=mean(rmpvalues)
         figure(rmpplot)
         RMPmeanLine = line([0 17] , [output_RMP output_RMP]);
         set(RMPmeanLine,'color','r','LineStyle','--') %So far good to here with EDR data


        %Input resistance-------------------------------------------
        RiCalc=[];
    
        for i = 1:length(PassiveFrames)
            istring = num2str(PassiveFrames(i,1));
            framename= ['frame' istring];
            temprmp = mean(frames.(framename)(1:(round(0.05*SR)),3));
            if mean(frames.(framename)(0.1*SR:0.4*SR,3)) < temprmp; %the voltage excursion to calculate Ri is now the maximum change in V, before Ih has an effect
                evokedV = min(frames.(framename)(0.1*SR:0.4*SR,3));
            else
                evokedV = max(frames.(framename)(0.1*SR:0.4*SR,3));
            end
            figure(Vtraces)
            subplot(4,6,i);                       
            rmpmarker = line([min(frames.(framename)(:,1)) max(frames.(framename)(:,1))] ,[evokedV evokedV]);
            set(rmpmarker,'color','r','LineStyle','--')
            deltaV=evokedV-temprmp;
            initI = mean(frames.(framename)(1:(round(0.05*SR)),2));
            maxI = mean(frames.(framename)(0.1*SR:0.4*SR,2));
            deltaI = maxI-initI;
            RiCalc(i,1)=deltaV;
            RiCalc(i,2)=deltaI;                    
        end
        
        
        %plot the V wrt I values
        RiFig = figure;
        scatter(RiCalc(:,2),RiCalc(:,1))
        hold on
        xlabel('Change in I (nA)')
        ylabel('Change in V (mV)')
        xlim([min(RiCalc(:,2)) max(RiCalc(:,2))])
        ylim([min(RiCalc(:,1)) max(RiCalc(:,1))])
        line([min(RiCalc(:,2)) max(RiCalc(:,2))] , [0 0])
        line([0 0] , [min(RiCalc(:,1)) max(RiCalc(:,1))]) 

        %Fit the polynomial to V vs I
        RiFit = polyfit(RiCalc(:,2),RiCalc(:,1),1);
        RiFitLine = polyval(RiFit,RiCalc(:,2));
        plot(RiCalc(:,2),RiFitLine)


        grad=gradient(RiFitLine,RiCalc(:,2));
        m=grad(9,1);
       % c=polyval(RiFit,0);
       % x=RiCalc(:,2);
       % plot(x,(m*x)+c,'color','k')

        Ri=m; %input resistance in mV/nA -->M Ohm (10^-3 / 10^-9)
        Ri_Output_MegaOhms = Ri %Assuming I is in nA, this seems to work OK compared with raw data



        %Current threshold to elicit one AP--------------------------


       
     %this loop picks out frames with a single AP   
    Candidates=[]; 
    i=1;  
    i2=1;
    NumberofAPs = zeros(17,1);
     for i = 1:17
        istring = num2str(i);
        APframename= ['APsforframe' istring];    
            NumberofAPs(i,1) = length(AP.(APframename));
            if NumberofAPs(i,1) == 1;
            Candidates(i2,1)=i;
            i2=i2+1;
            end
     end
     
    % if isempty(Candidates)==0 %If candidates were found....
         %...then this loop finds the one which had the smallest I
     %       i=1;
      %       for i=1:length(Candidates)
       %      istring = num2str(Candidates(i,1));
        %     framename= ['frame' istring];
         %    Candidates(i,2)=max(frames.(framename)(:,2))-min(frames.(framename)(:,2));
          %   end
     %output_I_threshold = min(Candidates(:,2))
     %else         
      %   threshold_frame = find(NumberofAPs>0,1);
       %  istring = num2str(threshold_frame);
        % framename = ['frame' istring];
         %output_I_threhold = max(frames.(framename)(:,2))-min(frames.(framename)(:,2))
         %disp('I threshold was detected for multiple APs - No traces elicited only one action potential')
     %end



    
    %Single AP properties at threshold current
        %%%%%Now working on just the threshold frame%%%%%%
    ThresholdFrame = find(NumberofAPs>0,1);
    if isempty(ThresholdFrame) == 0
    TFstring = num2str(ThresholdFrame);
    framename = ['frame' TFstring];
    APframename = ['APsforframe' TFstring];
    APstart = floor((AP.(APframename)(1,1) - 0.003)*SR)
    APend = floor((AP.(APframename)(1,1) + 0.003)*SR)
    isolateAPV=Y(APstart:APend,2);
    isolateAPt=T(APstart:APend,1);
    %ThresholdAPStorage(I,:) = isolateAPV;
    FirstAPPlot = figure;
    plot(isolateAPt,isolateAPV);
    APLine = get(FirstAPPlot,'Children');
    set(APLine,'LineWidth',2)
    hold on
    title('Threshold action potential')
    ylabel('Membrane Potential (mV)')
    xlabel('Time(ms)')
    xlim([(T(APstart,1)) (T(APend,1))])
    ylim([min(isolateAPV)-5 max(isolateAPV)+5])
    dy2= diff(isolateAPV); %gives a measure of gradient in mV/sample rate
    dy2 = (dy2*SR)/1000; %multiplies by SR to measure change in V per second, then divides by 1000 to give change per ms
    maxrise = max(dy2);
    maxdecay = min(dy2);
    
    figure(1)
    hold on
    plot(isolateAPV,'color','g');
    APLine = get(FirstAPPlot,'Children');
    set(APLine,'LineWidth',2)
    hold on

    %Find the firing threshold (mV)
    FT = find(dy2>20,1); %finds the first point where rate of change exceeds 20 mV/ms
    output_FiringThreshold = isolateAPV(FT,1)
    figure(FirstAPPlot)
    FTLine = line([min(isolateAPt) max(isolateAPt)] ,[output_FiringThreshold output_FiringThreshold]);
    set(FTLine,'color','r','LineWidth',1)

    %Find the peak of the AP
    output_AP_Peak = max(isolateAPV)
    PeakLine = line([min(isolateAPt) max(isolateAPt)] ,[output_AP_Peak output_AP_Peak]);
    set(PeakLine,'color','r','LineWidth',1)

    %calculate the amplitude of the AP
    output_AP_amplitude = output_AP_Peak - output_FiringThreshold

    %calculate 10, 50 and 90% amplitudes
    TenpercentVm = output_FiringThreshold + 0.1*output_AP_amplitude;
    Tenpercentline = line([min(isolateAPt) max(isolateAPt)] ,[TenpercentVm TenpercentVm]);
    set(Tenpercentline,'color','r','LineStyle','--')

    FiftypercentVm = output_FiringThreshold + 0.5*output_AP_amplitude;
    Fiftypercentline = line([min(isolateAPt) max(isolateAPt)] ,[FiftypercentVm FiftypercentVm]);
    set(Fiftypercentline,'color','r','LineStyle','--')

    NinetypercentVm = output_FiringThreshold + 0.9*output_AP_amplitude;
    Ninetypercentline = line([min(isolateAPt) max(isolateAPt)] ,[NinetypercentVm NinetypercentVm]);
    set(Ninetypercentline,'color','r','LineStyle','--')

    %Find the times at which the AP crosses each Vm and plot them on the raw
    %trace
    %10 %
    firsttenpercentindex = find(isolateAPV>TenpercentVm,1);
    FirstTenPercentTime = isolateAPt(firsttenpercentindex,1);
    FirstTenLine = line([FirstTenPercentTime FirstTenPercentTime] , [min(isolateAPV)-5 max(isolateAPV)+5]);
    set(FirstTenLine,'color','k')

    lasttenpercentindex = find(isolateAPV>TenpercentVm,1,'last');
    LastTenPercentTime = isolateAPt(lasttenpercentindex,1);
    LastTenLine = line([LastTenPercentTime LastTenPercentTime] , [min(isolateAPV)-5 max(isolateAPV)+5]);
    set(LastTenLine,'color','g')

    %50%
    firstfiftypercentindex = find(isolateAPV>FiftypercentVm,1);
    FirstFiftyPercentTime = isolateAPt(firstfiftypercentindex,1);
    FirstFiftyLine = line([FirstFiftyPercentTime FirstFiftyPercentTime] , [min(isolateAPV)-5 max(isolateAPV)+5]);
    set(FirstFiftyLine,'color','b')

    lastfiftypercentindex = find(isolateAPV>FiftypercentVm,1,'last');
    LastFiftyPercentTime = isolateAPt(lastfiftypercentindex,1);
    LastFiftyLine = line([LastFiftyPercentTime LastFiftyPercentTime] , [min(isolateAPV)-5 max(isolateAPV)+5]);
    set(LastFiftyLine,'color','b')

    %90%
    firstninetypercentindex = find(isolateAPV>NinetypercentVm,1);
    FirstNinetyPercentTime = isolateAPt(firstninetypercentindex,1);
    FirstNinetyLine = line([FirstNinetyPercentTime FirstNinetyPercentTime] , [min(isolateAPV)-5 max(isolateAPV)+5]);
    set(FirstNinetyLine,'color','k')

    lastninetypercentindex = find(isolateAPV>NinetypercentVm,1,'last');
    LastNinetyPercentTime = isolateAPt(lastninetypercentindex,1);
    LastNinetyLine = line([LastNinetyPercentTime LastNinetyPercentTime] , [min(isolateAPV)-5 max(isolateAPV)+5]);
    set(LastNinetyLine,'color','g')

    %Calculate rise and decay times
    output_ten_to_ninety_rise_time = FirstNinetyPercentTime - FirstTenPercentTime
    output_ninety_to_ten_decay_time = LastTenPercentTime - LastNinetyPercentTime

    %Calculate half width
    output_half_width = LastFiftyPercentTime - FirstFiftyPercentTime
    end
  


        %Amplitude of AHP
        %Duration of AHP

    %%%%%%%%%% Frequency vs Current Relationship %%%%%%%%%%%%%%%%%%%%%   
    i=1;
    i2=1;
    I_Freq_plot=[];
    for i = 1:17
        istring = num2str(i);
        APframename = ['APsforframe' istring];
        framename= ['frame' istring];
        if length(AP.(APframename))>1
            TimeFromFirstToLastAP = max(AP.(APframename)) - min(AP.(APframename));
            numberofAPs = length(AP.(APframename));
            AvgTimeDiff= TimeFromFirstToLastAP/(numberofAPs-1);
            APfreq=1/(AvgTimeDiff);
            I_Freq_plot(i2,1)=APfreq;
            I_Freq_plot(i2,2) = max(frames.(framename)(:,2)) - min(frames.(framename)(:,2));
            i2=i2+1;
        end
        i=i+1;

    end
    if isempty(I_Freq_plot)==0         
        IvF = figure;
        scatter(I_Freq_plot(:,2),I_Freq_plot(:,1))
        title('Current vs Frequency Relationship')
        xlabel('Injected Current (pA)')
        ylabel('Frequency of evoked APs (Hz)') %Now good to here with EDR data
    end


    %Accommodation properties with at least 10 APs
    i=1;%index for frames
    i2=1;%index for action potentials
    i3=1;%index for thresholds
    i4=1;%index for amplitudes
    i5=1;%index for 10-90% rise times
    i6=1;%index for 90-10% decays
    i7=1;%index for half widths
    accomm_thresholds = struct ; %structure to store all the thresholds
    accomm_amps = struct; %structure to store AP amplitudes
    accomm_tentoninety = struct; %structure to store ten to ninety% rise times
    accomm_ninetytoten = struct; %structure to store 90-10% decay times
    accomm_halfwidth = struct; %structure to store half widths
    GenT=0:0.0001:0.5; %creates a generic time base for accommodation analysis

    for i = 1:17; %everything inside this loop is done for each frame

        istring = num2str(i);
        framename = ['frame' istring];
        APframename = ['APsforframe' istring];

        if isempty(AP.(APframename)) == 0  %if the frame is not empty........

            for i2 = 1:length(AP.(APframename))%for every AP within the frame.....
               APstart = (find(frames.(framename)(:,1)>= AP.(APframename)(i2,1),1)) - 0.0025*SR;
               APend  = (find(frames.(framename)(:,1)>= AP.(APframename)(i2,1),1)) + 0.0025*SR;             
               isolateAPV=frames.(framename)(APstart:APend,3);
               isolateAPt=frames.(framename)(APstart:APend,1);
               dy2= diff(isolateAPV); %gives a measure of gradient in mV/sample rate
               dy2 = (dy2*SR)/1000;


                %Find the firing threshold (mV)
                FT = find(dy2>20,1); %finds the first point where rate of change exceeds 35 mV/ms
                tempFiringThreshold = isolateAPV(FT,1);
                accomm_thresholds.(framename)(i2,1)=tempFiringThreshold;



                %Find the peak of the AP
                tempAP_Peak = max(isolateAPV);

                %calculate the amplitude of the AP
                tempAP_amplitude = tempAP_Peak - tempFiringThreshold;
                accomm_amps.(framename)(i2,1)=tempAP_amplitude;


                %calculate 10, 50 and 90% amplitudes
                TenpercentVm = tempFiringThreshold + 0.1*tempAP_amplitude;

                FiftypercentVm = tempFiringThreshold + 0.5*tempAP_amplitude;

                NinetypercentVm = tempFiringThreshold + 0.9*tempAP_amplitude;

                %Find the times at which the AP crosses each Vm and plot them on the raw
                %trace
                %10 %
                firsttenpercentindex = find(isolateAPV>TenpercentVm,1);
                FirstTenPercentTime = isolateAPt(firsttenpercentindex,1);

                lasttenpercentindex = find(isolateAPV>TenpercentVm,1,'last');
                LastTenPercentTime = isolateAPt(lasttenpercentindex,1);

                %50%
                firstfiftypercentindex = find(isolateAPV>FiftypercentVm,1);
                FirstFiftyPercentTime = isolateAPt(firstfiftypercentindex,1);

                lastfiftypercentindex = find(isolateAPV>FiftypercentVm,1,'last');
                LastFiftyPercentTime = isolateAPt(lastfiftypercentindex,1);

                %90%
                firstninetypercentindex = find(isolateAPV>NinetypercentVm,1);
                FirstNinetyPercentTime = isolateAPt(firstninetypercentindex,1);

                lastninetypercentindex = find(isolateAPV>NinetypercentVm,1,'last');
                LastNinetyPercentTime = isolateAPt(lastninetypercentindex,1);

                %Calculate rise and decay times
                tempten_to_ninety_rise_time = FirstNinetyPercentTime - FirstTenPercentTime;
                accomm_tentoninety.(framename)(i2,1) = tempten_to_ninety_rise_time;


                tempninety_to_ten_decay_time = LastTenPercentTime - LastNinetyPercentTime;
                accomm_ninetytoten.(framename)(i2,1) = tempninety_to_ten_decay_time;



                %Calculate half width
                temphalf_width = LastFiftyPercentTime - FirstFiftyPercentTime;
                accomm_halfwidth.(framename)(i2,1) = temphalf_width;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Amplitude of AHP
                    %Duration of AHP

                i2=i2+1;%move to the next AP           



            end


        end 
         i=i+1;%move on to next frame
     end

        %Frequency vs current relationship
        %Change in inter-AP freq
        %Change in amplitude
        %Change in half width
        %Change in rise slope
        %Change in afterdepol

    %Save and export figures and numbers; close everything before moving to
    %the next file in the input dir
    %mkdir(['C:/Users/gmorris/Desktop/analysis/Output/' filename '/'])
    %saveas (Vtraces, ['C:/Users/gmorris/Desktop/analysis/Output/' filename '/' filename ' voltage traces figure'] , 'fig')
    %saveas (Itraces, ['C:/Users/gmorris/Desktop/analysis/Output/' filename '/' filename ' current traces figure'] , 'fig')
    %saveas (rmpplot, ['C:/Users/gmorris/Desktop/analysis/Output/' filename '/' filename ' resting membrane potential figure'] , 'fig')
    %saveas (RiFig, ['C:/Users/gmorris/Desktop/analysis/Output/' filename '/' filename ' input resistance figure'] , 'fig')
    %if isempty(I_Freq_plot)==0  
     %   saveas (IvF, ['C:/Users/gmorris/Desktop/analysis/Output/' filename '/' filename ' current vs frequency figure'] , 'fig')
    %end
    
 
    
    
   
   
    results{I+1,1} = filename;
    results{I+1,2} = output_RMP;
    results{I+1,3} = Ri_Output_MegaOhms;
    results{I+1,4} = output_AP_amplitude;
    results{I+1,5} = output_ten_to_ninety_rise_time;
    results{I+1,6} = output_ninety_to_ten_decay_time;
    results{I+1,7} = output_half_width;
    results{I+1,8} = maxrise;
    results{I+1,9} = maxdecay;
    results{I+1,10} = output_FiringThreshold;
    results{I+1,11} = ThresholdFrame;
    
    %max firing freq
    if isempty(I_Freq_plot)==0  
    diffs = diff(AP.APsforframe17);
    freqs = 1/diffs;
    results{I+1,12} = max(freqs);
    end
    
    
   Completion = (I/length(directory))*100
    
end
