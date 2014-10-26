function varargout = getOutput(funcName,outputNo,varargin)
    %get the nth(outputNo) output from a function
    %call syntax var = getOutput(@funcname, outputNo, vararin)
    %donot forget @ before the funcname when calling the function
    varargout = cell(max(outputNo),1);
    [varargout{:}] = funcName(varargin{:});
    varargout = varargout(outputNo);
end