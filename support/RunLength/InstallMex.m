function Ok = InstallMex(SourceFile, varargin)
% Compile a C-Mex file
% This function calls MEX() to compile a C/C++/F-Mex source file.
% Advanced users can call "mex -O SourceFile.c ..." instead, but beginners are
% sometimes overwhelmed by compiling instructions.
%
% Ok = InstallMex(SourceFile, ...)
% INPUT:
%   SourceFile: Name of the source file, with or without absolute or partial
%               path. The default extension '.c' is appended on demand.
%   Optional arguments:
%   - Function name: Started after the compilation, e.g. a unit-test.
%   - Cell string:   Additional arguments for the compilation, e.g. libraries.
%   - '-debug':      Enabled debug mode.
%   - '-force32':    Use compatibleArrayDims under 64 bit Matlab.
%   - '-replace':    Overwrite existing mex file.
%
% OUTPUT:
%   Ok: Logical flag, TRUE if compilation was successful.
%
% COMPATIBILITY:
% - A compiler must be installed and setup before: mex -setup
% - For Linux and MacOS the C99 style is enabled for C-files.
% - The optimization flag -O is set.
% - The compiler directive -DMATLABVER<XYZ> is added to support pre-processor
%   switches, where <XYZ> is the current version, e.g. 708 for v7.8.
%
% EXAMPLES:
% Compile func1.c with LAPACK libraries:
%   InstallMex('func1', {'libmwlapack.lib', 'libmwblas.lib'})
% Compile func2.cpp, enable debugging and call a test function:
%   InstallMex('func2.cpp', '-debug', 'Test_func2');
% These lines can be appended after the help section of an M-file, when the
% compilation should be started automatically, if the compiled MEX is not found.
%
% Suggestions for improvements and comments are welcome!
% Feel free to add this function to your FEX submissions, when you change the
% "Precompiled" URL accordingly.
%
% Tested: Matlab 6.5, 7.7, 7.8, 7.13, WinXP/32, Win7/64
% Author: Jan Simon, Heidelberg, (C) 2012-2013 j@n-simon.de

% $JRev: R-t V:019 Sum:P+3HMiYKQ2d+ Date:17-Sep-2013 00:46:17 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $File: Tools\GLSource\InstallMex.m $
% History:
% 001: 27-Jul-2012 09:06, First version.
% 005: 29-Jul-2012 17:11, Run the unit-test instead of showing a link only.
% 006: 11-Aug-2012 23:59, Inputs are accepted in free order.

% Initialize: ==================================================================
% Global Interface: ------------------------------------------------------------
% URL to folder containing pre-compiled files, or the empty string if
% pre-compiled filesare not offered:
Precompiled = 'http://www.n-simon.de/mex';

% Initial values: --------------------------------------------------------------
bakCD   = cd;
matlabV = [100, 1] * sscanf(version, '%d.%d', 2);  % Numerical Matlab version

% Program Interface: -----------------------------------------------------------
% Parse inputs:
Param       = {};
UnitTestFcn = '';
doDebug     = false;
debugFlag   = {};
force32     = false;
replace     = false;

% First input is the name of the source file:
if ~ischar(SourceFile)
   error_L('BadTypeInput1', '1st input must be a string.');
end

% Additional inputs are identified by their type:
% String:      unit-test function or the flag to enable debugging.
% Cell string: additional parameters for the MEX command
for iArg = 1:numel(varargin)
   Arg = varargin{iArg};
   if ischar(Arg)
      if strcmpi(Arg, '-debug')
         doDebug     = true;
         debugFlag   = {'-v'};
      elseif strcmpi(Arg, '-force32')
         force32     = true;
      elseif strcmpi(Arg, '-replace')
         replace     = true;
      elseif exist(Arg, 'file') == 2
         UnitTestFcn = Arg;
      else
         error_L('MissFile', 'Unknown string or missing file: %s', Arg);
      end
   elseif iscellstr(Arg)  % As row cell:
      Param = Arg(:).';
   else
      error_L('BadInputType', 'Bad type of input.');
   end
end

% User Interface: --------------------------------------------------------------
hasHRef = usejava('jvm');   % Hyper-links in the command window?

% Do the work: =================================================================
% Search the source file, solve partial or relative path, get the real
% upper/lower case:
[dummy, dummy, Ext] = fileparts(SourceFile);  %#ok<ASGLU>
if isempty(Ext)
   SourceFile = [SourceFile, '.c'];
end

whichSource = which(SourceFile);
if isempty(whichSource)
   error_L('NoSource', 'Cannot find the source file: %s', SourceFile);
end
[SourcePath, SourceName, Ext] = fileparts(whichSource);
Source                        = [SourceName, Ext];
mexName                       = [SourceName, '.', mexext];

fprintf('== Compile: %s\n', fullfile(SourcePath, Source));

% Check if the compiled file is existing already:
whichMex = which([SourceName, '.', mexext]);
if ~isempty(whichMex)
   fprintf('::: Compiled file is existing already:\n    %s\n', whichMex);
   if replace
      fprintf('  Recompiled in: %s\n', SourcePath);
   else
      if nargout
         Ok = false;
      end
      return;
   end
end

if ~ispc && strcmpi(Ext, '.c')
   % C99 for the GCC and XCode compilers.
   % Note: 'CFLAGS="\$CFLAGS -std=c99"' must be separated to 2 strings!!!
   Flags = {'-O', 'CFLAGS="\$CFLAGS', '-std=c99"', Source};
else
   Flags = {'-O', Source};
end

% Large array dimensions under 64 bit, possible since R2007b:
matlabVDef = {sprintf('-DMATLABVER=%d', matlabV)};
if matlabV >= 705
   if any(strfind(computer, '64')) && ~force32
      Flags = cat(2, {'-largeArrayDims'}, Flags);
   else
      Flags = cat(2, {'-compatibleArrayDims'}, Flags);
   end
end

% Compile: ---------------------------------------------------------------------
% Display the compilation command:
Flags = cat(1, debugFlag, matlabVDef, Flags(:), Param(:));
cmd   = ['mex', sprintf(' %s', Flags{:})];
fprintf('  %s\n', cmd);

try    % Start the compilation:
   cd(SourcePath);
   mex(Flags{:});
   fprintf('Success:\n  %s\n', which(mexName));
   compiled = true;
   
   % Show other instances of the MEX and M, if the replace flag was set:
   if replace
      allWhich = which(SourceName, '-all');
      if length(allWhich) > 1
         fprintf('\n::: Multiple instances found:\n');
         fprintf('  %s\n', allWhich{:});
      end
   end
   
catch  % Compilation failed - no MException to support Matlab 6.5:
   compiled = false;
   err      = lasterror;
   fprintf(2, '\n*** Compilation failed:\n%s\n', err.message);
   if ~doDebug  % Compile again in debug mode if not done already:
      try
         mex(Flags{:}, '-v');
      catch  % Empty
      end
   end
   
   % Show commands for manual compilation and download pre-compiled files:
   fprintf('\n== The compilation failed! Possible solutions:\n');
   fprintf('  * Install and set up a compiler on demand:\n');
   if hasHRef
      fprintf('    <a href="matlab:mex -setup">mex -setup</a>\n');
      fprintf('  * Try to compile manually:\n    cd(''%s'')\n    %s -v\n', ...
         SourcePath, cmd);
      if ~isempty(Precompiled)
         fprintf('  * Or download the pre-compiled file %s:\n', mexName);
         fprintf('    <a href="matlab:web(''%s#%s'',''-browser'')">%s</a>\n', ...
            Precompiled, mexName, Precompiled);
      end
   else  % No hyper-references in the command window without Java:
      fprintf('  * mex -setup\n');
      fprintf('  * Try to compile manually:\n  cd(''%s'')\n  %s -v\n', ...
         SourcePath, cmd);
      if ~isempty(Precompiled)
         fprintf('  * Or download the pre-compiled file %s:\n  %s\n', ...
            mexName, Precompiled);
      end
   end
   fprintf('  * Or send this report to the author\n');
end

% Restore original directory:
cd(bakCD);

% Run the unit-test: -----------------------------------------------------------
if ~isempty(UnitTestFcn) && compiled
   fprintf('\n\n========== Post processing: ==========\n');
   [dum, UnitTestName] = fileparts(UnitTestFcn);  %#ok<ASGLU> % Remove extension
   if ~isempty(which(UnitTestName))
      fprintf('  Call: %s\n\n', UnitTestName);
      feval(UnitTestName);
   else
      fprintf(2, '??? Cannot find unit-test: %s\n', UnitTestFcn);
   end
end

% Return success of compilation: -----------------------------------------------
if nargout
   Ok = compiled;
end
fprintf('\n%s ready.\n', mfilename);

% end

% ******************************************************************************
function error_L(ID, Msg, varargin)
% Automatic error ID and mfilename in the message:
error(['JSimon:', mfilename, ':', ID], ...
   ['*** %s: ', Msg], mfilename, varargin{:});

% end
