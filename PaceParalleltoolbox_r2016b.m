%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This script opens a local pool of Matlab workers (not Distributed     %
%%% Computing Server) on a desktop computer (Windows, MacOS) or a cluster %
%%% with a torque scheduler for Matlab versions r2013b, r2014a, r2015a,   %
%%% r2015a                                                                %
%%% and torque versions 4.2.5.h3 and 5.0.1.h2.                            %
%%%                                                                       %
%%% Blake Fleischer                                                       %
%%% Georgia Institute of Technology                                       %
%%% Partnership for an Advanced Computing Environment (PACE)              %
%%% http://www.pace.gatech.edu                                            %
%%% April 29, 2015                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PaceParallelToolbox_r2015a opens of matlab pool for local computers and
% Partnership for an Advanced Computing Environment (PACE) clusters.
%
% Functionality of the previous version is maintained - if you used an
% older version simply replace the function name (you can keep the
% previous inputs). The default behavior (i.e. running
% Pace_ParallelToolbox_r2015a) is Pace_ParallelToolbox_r2015a(false)
% for a desktop, and PaceParallelToolbox_r2015a(true) for the PACE
% cluster.
%
%
% INPUTS (all are optional):
% maxout : Applies to local computer only. True/False for setting the
% matlab pool equal to the number of cores on your local computer.
% If false, poolsize = total cores - 1. Default is false, so you don't
% totally max out your computer. Simply put true or false in front of all
% other arguments.
% For example: Pace_ParallelToolbox_r2015a(true)
%
% Parameter-specified inputs:
% cores : Manually specify the number of cores. If more are specified
% than are on the system, the max available cores will be used. This
% function will use the value of cores over the value of maxout if both
% are specified.
% For example: Pace_ParallelToolbox_r2015a('cores',3)
%
% job_storage : Applies to cluster only. Manually specify the matlab
% environment variable cluster.JobStorageLocation, i.e the folder where
% matlab's parallel processing files are kept. These files are only useful
% during a running job, and are used for matlab specific taks information
% is stored. This variable is intended for automating the cleaning up the
% files left behind many parallel job runs. Specify a full path.
% For example:
% Pace_ParallelToolbox_r2015a('job_storage','~/data/JOBID_58885-coolbeans')
%
% debug : Applies to cluster only. Optional parameter for running a series
% of tests to debug potential problems in the PACE cluster.
% For example: 'debug',true
%
%
% OUTPUTS:
% currCores : the size of the pool actually opened. Should be equal to or
% less than the number of cores.
%
% memRemain : the amount of ram on the system availble for computation (GB).
% Useful for determining the size of variables that are poolsize dependent.
%
%
% Changelog:
% PaceParallelToolbox_r2016b (only affects PACE jobs)
% -Queries job environment starting with torque (qstat -f), then
% moab (checkjob), and finally the pbs nodefile.
% -Opens matlabpools for jobs without memory specification, but throws warning.
% -Query environment for interactive job vs torque (qstat -f).
% -General tidying up source.
% -Prints command to close matlab pool for interactive jobs.

% PaceParallelToolbox_r2015a
% -MatLab r2015a compatible (and hopefully future versions).
% -Finds all logical cores on desktop computers - not just physical cores.
% -Mac Compatibility! Linux (64bit) might work too - it's not fully
% tested though.
% -Outputs available memory.
% -Lets you know if you're on a pace headnode (please don't run calculations on
% headnodes).
% -Deletion of the JobStorageLocation for PACE cluster runs is now
% possible because the folder can be set externally (and deleted
% afterwards via bash). This could prevent the filling of home folders on
% the cluster (preventing many other problems).

function [currCores, memRemain] = PaceParalleltoolbox_r2016b(varargin)

%% Check current version of matlab:
currVer = textscan(version,'%s','delimiter',' ');
currVer = textscan(currVer{1,1}{1,1},'%s','delimiter','.');
currVer = str2double([currVer{1,1}{1,1}, '.', currVer{1,1}{2,1}]);
% Some previous versions:
% 8.2 : R2013b
% 8.3 : R2014a
% 8.4 : r2015a
% 8.5 : R2015a

%% Matlab version specific variables settings:
if currVer == 8.2
    maxCores=12; % max allowable cores in r2013b,etc
    currCores = 0; % number of cores currently used in matlab pool
elseif currVer >= 8.3 % if matlab r2014a or later
    maxCores=512; % max allowable cores in r2014a
    currCores = 0;
    currClust = gcp('nocreate');
    if ~isempty(currClust)
        currCores = currClust.NumWorkers;
    end
else
    fprintf(['It looks like you''re using matlab (%s), which is older than ',...
        'r2013b. %s doesn''t currently support versions of matlab older than r2013b.\n',...
        'Sorry!\n'],version,mfilename)
    return
end

tic % Determine the amount of time it took to open the matlab pool.

%% Input parser - use matlab built-in input parser to simplify optional arguments
p = inputParser;
% Set defaults/conditions
defaultMaxOut = false;
defaultCores = inf;
defaultDebug = 0;
defaultJob_storage = 'default';

% Create variable parameters
addOptional(p,'maxout',defaultMaxOut,@islogical)
addParameter(p,'debug',defaultDebug,@islogical);
addParameter(p,'cores',defaultCores,@isnumeric)
addParameter(p,'job_storage',defaultJob_storage,@ischar)
% parse
parse(p,varargin{:})
% Set results:
maxout = p.Results.maxout;
debug = p.Results.debug;
desiredCores = p.Results.cores;
job_storage = p.Results.job_storage;


%% Check running environment - type of desktop or cluster?
currEnv = 'desktop';
if strfind(getenv('HOSTNAME'),'.pace.gatech.edu')
    currEnv = 'PACEcluster';
end

if strcmp(currEnv,'desktop') % Check for type of desktop environment
    currEnv = computer;
end

%% Determine resources available for different platforms:
memAvail = 0;
memRemain = memAvail;
if sum( strcmp(currEnv,{'PCWIN','PCWIN64'}) )>0
    import java.lang.*;
    r = Runtime.getRuntime;
    totalCores=r.availableProcessors;
    % totalCores = str2double(totalCores);
    % List free available memory (GB):
    [~, memInfo] = memory; %Determine total RAM available
    memAvail = memInfo.PhysicalMemory.Available;
    memAvail = memAvail/1024^3;
    
elseif sum( strcmp(currEnv,{'MACI64','GLNXA64'}) )>0
    import java.lang.*; % Determine desktop parameters via java, not matlab
    %function: matlab functions only detect physical procs.
    r=Runtime.getRuntime;
    totalCores=r.availableProcessors;
    % List free available memory (GB), parsing from top command
    [~, memAvail]=system('top -l 1 | head -n 10 | grep PhysMem');
    memAvail = strsplit(memAvail,' ');
    memAvail = str2double(memAvail{6}(1:end-1))/1024;
    
elseif strcmp(currEnv,'PACEcluster')
	% Obtain job info, set to: totalCores, memAvail, interActiv and vncServer.
	if strcmp(getenv('PBS_JOBID'),'') %if no pbs PBS_JOBID env variable, then must be on headnode.
        fprintf(['Can''t find an active job id for this job, and it looks like you''re running on a headnode.\n',...
            'Please run this script by submitting a job (and ',...
            'specifying memory>0.33*nprocs).\n'])
        return
    end
    
	% Determine how to get info from environment:
	qstat_working = true;
	checkjob_working = true;
	pbsndfle_working = true;
	defaultMem = 1; % Default on pace systems is 1 gb
	
	try % Try torque's qstat -f 
		[~, jobConds] = system(['qstat -f ',getenv('PBS_JOBID')]);
		%fprintf('Current running environment: \n%s',jobConds)
		jobCondsArr = strsplit(jobConds)';
		% Check number of cores
		totalCores=sscanf(jobCondsArr{find(strcmp('Resource_List.nodes',jobCondsArr))+2},'%*6c%u',[1, Inf]);
		% Check memory
		if isempty(find(strcmp('Resource_List.mem',jobCondsArr), 1))
			fprintf(['Please specify ram when making job submissions - use ''[qsub or msub] -l mem=XXgb.'' ',...
			'The parallel pool script needs to know if you have reserved enough ram for the ',...
			'job. Otherwise the job may unexpectedly quit. No memory specified - assuming default ',...
			'memory allocation of 1gb\n'])
			memAvail = defaultMem;
		else
			memAvail=sscanf(jobCondsArr{find(strcmp('Resource_List.mem',jobCondsArr))+2},'%u%*2c',[1, Inf]);
		end
		% Check for vnc server script (i.e. job is interactive, but running as vnc server)
		vncServer = any(~cellfun('isempty',strfind(jobCondsArr,'vncserverscript')));
    catch
		fprintf('\nTorque qstat not working!\n')
		fprintf('qstat error: %s\n', jobConds);
		qstat_working = false;
	end
	if ~qstat_working % Try moab's checkjob
		try
			fprintf('Attempting Moab checkjob...')
			[~, jobConds] = system(['checkjob ',getenv('PBS_JOBID')]);
			%fprintf('Current running environment: \n%s',jobConds)
			jobCondsArr = strsplit(jobConds)';
			% Check number of cores
			totalCores=sscanf(jobCondsArr{find(strcmp('TaskCount:',jobCondsArr))+1},'%u',[1,Inf]);
			% Check memory			
			memPcore=sscanf(jobCondsArr{find(strcmp('MEM:',jobCondsArr))+1},'%d%*1c',[1, Inf]);
			memAvail = [memPcore/1024*totalCores];
			fprintf(' success.\n')
		catch
			fprintf(' Moab checkjob not working!\n')
			fprintf('checkjob error: %s\n', jobConds);
			checkjob_working = false;
		end
		% Check for vnc server script (i.e. job is interactive, but running as vnc server)
		vncServer = any(~cellfun('isempty',strfind(jobCondsArr,'vncserverscript'))); %CHECK THIS!
	end	
	if ~checkjob_working %if no qstat or checkjob command
		try
			fprintf('Checking pbs nodefile...')
			% Check number of cores
			[~, totalCores] = system('cat ${PBS_NODEFILE} | wc -l');
			totalCores = str2num(totalCores);
			% Check Memory
			fprintf(' success.\n')
			fprintf(['!Important note! The job memory can''t be obtained.\n',...
			'Please use caution when using this script and be sure to specify enough memory\n',...
			'keeping in mind that matlab occupies approx 333 mb/core. Defaulting to 0.\n']);
		catch
			fprintf(' PBS NODEFILE not available: \n%s\n', getenv('PBS_NODEFILE'))
			system('cat ${PBS_NODEFILE}');
			fprintf('Hmm... we''re out of ideas at the moment. Looks like it''s time to contact pace support (if you''re at Georgia Tech. Exiting without creating MATLAB parpool.\n');
			return
		end
		% Check for vnc server script (i.e. job is interactive, but running as vnc server)
		vncServer = false; %can't really find out in this case, so err on side of caution
	end
	% Check for interactive jobs
	interActiv = strcmp(getenv('PBS_ENVIRONMENT'),'PBS_INTERACTIVE');
else
    fprintf(['The current achetecture of this computer (%s) isn''t ',...
        'currently supported. Sorry! \n'],currEnv)
    return
end

%% DebugFunction
if debug
    debugMatlab(currEnv);
end

%% Determine size of pool to create
if totalCores < maxCores
    maxCores = totalCores;
end

if sum( strcmp(currEnv,{'PCWIN','PCWIN64','MACI64','GLNXA64'}) )>0
    if desiredCores==inf %maxcores only used if desiredCores is set higher
        if maxout
            desiredCores = maxCores;
        else
            desiredCores = maxCores-1;
        end
    elseif desiredCores>maxCores
        desiredCores = maxCores;
    end
elseif strcmp(currEnv,'PACEcluster')
    if desiredCores>maxCores
        desiredCores = maxCores;
    end
end

if desiredCores < 2; %Don't open pools < 2 cores, but return memory information.
    fprintf('Poolsize is equal to or less than 1, not opening matlab pool. \n')
    return
end

%% Check if pool is already created that is of desired size
if currCores == desiredCores;
    fprintf('Matlabpool already created.\n')
    return
end

%% Close any open matlab pools that are open and not of the right size
if currCores >1 && currCores ~= desiredCores
    if currVer == 8.2
        matlabpool close
    elseif currVer >= 8.3 %if matlab r2014a or later
        delete(gcp('nocreate'));
    end
end
%% Open the matlab pool (start the headless versions of matlab)
if sum( strcmp(currEnv,{'PCWIN','PCWIN64','MACI64','GLNXA64'}) )>0
    if desiredCores==totalCores;
        fprintf(['%u cores available, using the max of %u cores ',...
            '(no multitask).\n'], totalCores, desiredCores)
    else desiredCores<totalCores;
        fprintf(['%u cores available, using only %u core(s), ('...
            'useful for multitasking).\n'],totalCores, desiredCores)
    end
    if currVer == 8.2;
        matlabpool(desiredCores)
    elseif currVer >= 8.3 % if matlab r2014a or later
        currClust = parcluster('local');
        if ~isempty(currClust.Jobs) % Delete any current jobs - not sure where from
            delete(currClust.Jobs)
        end
        fprintf('opening workers...\n')
        set(currClust,'NumWorkers',desiredCores);
        parpool('local',desiredCores);
        currPool = gcp('nocreate');
    end
    
elseif strcmp(currEnv,'PACEcluster')
    % Change location of job storage location so multiple jobs can run
    %simultaneously. If the folders don't get deleted, just remove the
    %jobs folders for completed or cancelled jobs.
    fprintf('\nInitializing local scheduler data.\n')
    currClust = parcluster('local');
    distcomp.feature( 'LocalUseMpiexec', false); %Unfortunately, this problem crops up about 50% of the time.
    % Can check this is a problem by running a GUI, validating the 'local' cluster profile, and checking the output.
    
    % Check user specified directory
    default_clust_data = fullfile( currClust.JobStorageLocation,...
        ['jobID_', strtok(getenv('PBS_JOBID'),'.')]);
    useUserJSL = false;
    if ~strcmp(job_storage,'default')
        directoryLegit = false;
        
        if ~any(regexp(job_storage,'[\/]'));
            fprintf(['File path (''%s'') must contain directory ',...
                'slashes (i.e. ''/'').\n',...
                '-----Using default instead-----\n'],job_storage)
        else
            jobStorageSpl = strsplit(job_storage,{'\','/'});
            
            % Recursively check to determine root directory (allows
            %subdirectories to be specified)
            for root = 1:length(jobStorageSpl)-1
                jobStorageRoot = fullfile(jobStorageSpl{1:end-root});
                if job_storage(1)=='/';
                    jobStorageRoot=['/',jobStorageRoot];
                end
                if exist(jobStorageRoot,'dir')
                    directoryLegit = true;
                    break
                end
            end
            
            if directoryLegit
                [~, jobStorageAttrib,~]=fileattrib(jobStorageRoot); %Check file attributes
                if jobStorageAttrib.UserWrite && jobStorageAttrib.UserRead
                    useUserJSL = true;
                else
                    fprintf(['\nMatlab doesn''t have read/write access to ''%s'' ',...
                        'for specified directory ''%s''.\n',...
                        '-----Using default instead-----\n'],jobStorageRoot,job_storage)
                end
            else
                fprintf(['\nCan''t find directory ''%s'' for specified '...
                    'directory ''%s'' (does it exist? need any slashes?).\n',...
                    '-----Using default instead-----.\n'], jobStorageRoot, job_storage)
            end
        end
    end
	
    % On PACE clusters, there is a wierd bug where if a user specifies /nv/
    %matlab thinks it has write access, and then can't open a folder
    %there, crashing. It's a problem with the matlab fileattrib function,
    %not this script.
    if useUserJSL
        local_clust_data = job_storage;
    else
        local_clust_data = default_clust_data;
    end
    
    % Create/set the job storage location
    mkdir(local_clust_data);
    currClust.JobStorageLocation = local_clust_data;
    fprintf('Matlabpool job_storage location is: %s\n',...
        currClust.JobStorageLocation);
    
    % Use a random time wait to open pools if not interactive jobs. Multiple
    %pools on the same node opened at the same time can fail.    
    if ~interActiv && ~vncServer % User probably won't open multiple interactive jobs on the
		% same node at the same time, but it is a problem on batch jobs
        fprintf('Job is not interactive.\n')
        rng(str2double(strtok(getenv('PBS_JOBID'),'.'))); % Uniquely seed the random
        %number generator
        randWaitTime = 1+60*rand();% Multiple parallel pools starting on the same
        %machine and time can cause problems...
        fprintf('Waiting %2.3f seconds to prevent parpool opening errors...',...
            randWaitTime)
        pause(randWaitTime);
        fprintf('done.\n')
    else
        fprintf('Job is interactive.\n')
    end

    % Open the matlab pool
	if isnan(memAvail)
		fprintf('%u/%u workers about to be opened, but job memory is unknown.\n',desiredCores,maxCores)
	else
		fprintf('%u/%u workers about to be opened for %g gb ram.\n',desiredCores,maxCores,memAvail)
	end
	
    if currVer == 8.2;
        matlabpool(desiredCores)
    elseif currVer >= 8.3 %if matlab r2014a or later
        currClust.NumWorkers=desiredCores;
        try % If parpool fails (has a history of intermittely failing on cluster, perhaps ulimit?)
            parpool(currClust, desiredCores);
        catch err
            % Get environment variables to resubmit job if this one fails
            fprintf('\n\nOriginal error Information:\n');
            disp(getReport(err,'extended'));
            % Run some other commands to help diagnose the problem...
            debugMatlab(currEnv)
        end % end of if parpool fails debugging info
        
        % Set pool properties
        currPool = gcp('nocreate');
    end
	else
		fprintf(['The current achetecture of this computer (%s) isn''t ',...
			'currently supported. Sorry! \n'],currEnv)
		return
	end

	% Set pool to not time out, check if pool is created first
	if currVer >= 8.3 && ~isempty(currPool) %if matlab r2014a or later and pool is created
		currPool.IdleTimeout=inf;
		if currPool.IdleTimeout==inf;
			fprintf('Pool will not time out.\n')
		else
			fprintf('Pool IdleTimeout is: %s minutes.\n',currPool.IdleTimeout);
		end
	end

	%% Report on opening of pool
	% Time Convert anonymous function
	timeConvert=@(secIn) (sprintf('%u:%u:%u:%2.4f (D:H:M:S)', [floor(secIn/(24*60*60)),...
		floor((secIn-24*60*60*floor(secIn/(24*60*60)))/(60*60)),...
		floor((secIn-(24*60*60)*floor(secIn/(24*60*60))-(60*60)*floor((secIn-24*60*60*floor(secIn/(24*60*60)))/(60*60)))/60),...
		secIn-(24*60*60*floor(secIn/(24*60*60)))-(60*60*floor((secIn-24*60*60*floor(secIn/(24*60*60)))/(60*60)))-(60*floor((secIn-(24*60*60)*floor(secIn/(24*60*60))-(60*60)*floor((secIn-24*60*60*floor(secIn/(24*60*60)))/(60*60)))/60))]));
	chckPool = gcp('nocreate');
	currCores = chckPool.NumWorkers;

	fprintf('Matlabpool opened with %u workers in: %s\n', currCores, timeConvert(toc));
	if ~isnan(memAvail)
		memRemain = memAvail - 0.330 * desiredCores; % Update memRemain to remaining memory
		memRemain(memRemain < 0) = 0;
		fprintf('Approx free memory is %g GB.\n',memRemain);

		if memRemain < (0.330 * currCores)
			fprintf(['\n\n!!Warning - ram allocation of %g gb is probably too ',...
				'low.\nConsider running 0.330*ppn gb of ram plus the size ',...
				'of any variables * ppn (%u gb + %u*variables).\nJob may '...
				'terminate or MatLab may crash unexpectedly as a result.\n'],...
				memAvail,ceil(0.330*currCores),currCores);
		end
	end

	% Print final messages about pool creation
	if currVer < 8.3
		close_command = 'matlabpool close;';
	else %currVer >=8.3
		close_command = 'delete(gcp(''nocreate''));';
	end
	if interActiv
		fprintf('To close this matlab pool, run: "%s".\n',close_command)
	end
	
	function debugMatlab(currEnv) %Error Checking Function
        if strcmp(currEnv,'PACEcluster')
            [~, jobCondsErr] = system(['qstat -f ',getenv('PBS_JOBID')]);
            fprintf('Current running environment: \n%s',jobCondsErr)
            fprintf('\ndone.\n\n')
        end
        
        if sum( strcmp(currEnv,{'GLNXA64','PACEcluster'}) )>0
            fprintf('\nRunning top to check other processes:\n');
            !top -bn 1 | head -n 70
            fprintf('done.\n\n')
        end
        
        fprintf(['\n\nRunning checks on the current matlab status on the cluster. ',...
            'These checks are from:\n',...
            'http://www.mathworks.com/matlabcentral/answers/92124-why-am-i-unable-to-use-matlabpool-or-parpool-with-the-local-scheduler-or-validate-my-local-configura',...
            '\n'])
        
        fprintf('Checking the parallel computing toolbox...\n')
        license checkout Distrib_Computing_Toolbox
        fprintf('done.\n\n')
        
        fprintf('Checking version of PCT...\n')
        versionM = ver;
        if sum( strcmp(currEnv,{'MACI64'}) )>0
            versionMname ={versionM.Name};
            versionMversion = {versionM.Version};
            fprintf([versionMname{9}, ' : ', versionMversion{9},'\n\n']);
        else
            versionM(39)
        end
        fprintf('Checking version of MATLAB...\n')
        version
        fprintf('done.\n\n')
        
        fprintf('Disabling LocalUseMpiexec...\n')
        distcomp.feature( 'LocalUseMpiexec', false )
        fprintf('done.\n\n')
        
        if strcmp(currEnv,'PACEcluster')
            localSchedFolder=strsplit(prefdir,'R201');
            localSchedFolder = [localSchedFolder{1}, 'local_scheduler_data'];
            fprintf(['One could delete the local_scheduler data, ',...
                'it is located in %s.\n\n'],localSchedFolder)
        end
        
        if sum( strcmp(currEnv,{'MACI64','GLNXA64','PACEcluster'}) )>0
            fprintf('Checking hostname and ping...\n')
            !echo $HOSTNAME
            !ping -c 8 $HOSTNAME
            fprintf('done.\n\n')
            
            fprintf('Checking ulimit settings, this has been problematic in the past:\n')
            !ulimit -a
            fprintf('done.\n\n')
        end
        
        fprintf('Checking for a startup.m file ...\n')
        which startup.m
        
        fprintf('Checking Java memory settings:\n')
        fprintf('java.lang.Runtime.getRuntime.maxMemory = %u\n',java.lang.Runtime.getRuntime.maxMemory)
        fprintf('java.lang.Runtime.getRuntime.totalMemory = %u\n',java.lang.Runtime.getRuntime.totalMemory)
        fprintf('java.lang.Runtime.getRuntime.freeMemory = %u\n',java.lang.Runtime.getRuntime.freeMemory)
        fprintf('done.\n')
        
        fprintf('End of MatLab checks.\n\n')
        fprintf(['If you see this output (privided you''re a are a pace ',...
            'user) and want to ',...
            'submit a help request, please add all of this output to your '...
            'help request.\n\n'])      
    end

end



% %To close pool:
% if r2014a
%     delete(gcp('nocreate'));
% else
%     matlabpool close;
% end