function fprintAmplParamCLSU(fid, strName, dataInput, varargin)
% This m-file translates data in Matlab format into AMPL data format in a
% new AMPL data file.
% 
% source: Su and Judd (2011), Constrained Optimization Approaches to
% Estimation of Structural Models.
% Code Revised: Che-Lin Su, May 2010.

% Initialization
vecFirstIndex = [ 1 1 1 ];
intIncrement = 1;

% First indices
intArgLen = length(varargin);
for intI = 3:-1:1
	if intArgLen >= intI
        varargin{intI};
		vecFirstIndex(intI) = varargin{intI};	
	end
end

% Analyse input data

% intType = 0;
% if iscell(dataInput)
%	intType = 4;
% elseif isscalar(dataInput)
%	intType = 1;
% elseif isvector(dataInput)
% 	intType = 2;
% else
%	intType = 3; % Admittedly, this is a bit dodgy
% end

intType = 0;
if isscalar(dataInput)
	intType = 1;
elseif isvector(dataInput)
	intType = 2;
elseif length(size(dataInput))==2
	intType = 3; 
elseif length(size(dataInput))==3
	intType = 4; 
end

% Type 1: Scalar

if intType == 1
	% fprintf(fid, 'param %s := %s;\n', strName, num2str(dataInput));
    fprintf(fid, 'param %s := %16.12f;\n', strName, dataInput);
end

% Type 2: Vector

if intType == 2
	if ( size(dataInput, 1) > size(dataInput, 2) )
		dataInput = dataInput';
	end

	fprintf(fid, 'param %s := \n', strName);
        
	intCount = vecFirstIndex(1);
	for currentValue = dataInput
		% fprintf(fid, '\t%d\t%s\n', intCount, num2str(currentValue));
        fprintf(fid, '\t%d\t%16.12f\n', intCount, currentValue);
		intCount = intCount + intIncrement;
	end

	fprintf(fid, ';\n');
end

% Type 3: Matrix

if intType == 3
	intSizeX = size(dataInput, 1);
	intSizeY = size(dataInput, 2);
        
	fprintf(fid, 'param %s : \n', strName);
        
	for intJ = vecFirstIndex(2):intIncrement:vecFirstIndex(2)+(intSizeY-1)*intIncrement
		fprintf(fid, '\t%d\t', intJ);
	end
	fprintf(fid, ' := \n');
        
	for intI = vecFirstIndex(1):intIncrement:vecFirstIndex(1)+(intSizeX-1)*intIncrement
		fprintf(fid, '%d\t', intI);
		for intJ = vecFirstIndex(2):intIncrement:vecFirstIndex(2)+(intSizeY-1)*intIncrement
            fprintf(fid, '%6.4f\t', dataInput(intI+(1-vecFirstIndex(1)), intJ+(1-vecFirstIndex(2))));
		end
		fprintf(fid, '\n');
	end
	fprintf(fid, ';\n');
end

% type 4: 3-dimensional matrix

if intType == 4
	intSizeX = size(dataInput, 1);
	intSizeY = size(dataInput, 2);
    intSizeZ = size(dataInput, 3);
    
	fprintf(fid, 'param %s := \n', strName);
    for intZ = vecFirstIndex(3):intIncrement:vecFirstIndex(3)+(intSizeZ-1)*intIncrement   
        fprintf(fid, '[*,*,%d]: \n', intZ);     
        for intJ = vecFirstIndex(2):intIncrement:vecFirstIndex(2)+(intSizeY-1)*intIncrement
            fprintf(fid, '\t%d\t', intJ);
        end
        fprintf(fid, ' := \n');
        
        for intI = vecFirstIndex(1):intIncrement:vecFirstIndex(1)+(intSizeX-1)*intIncrement
            fprintf(fid, '%d\t', intI);
            for intJ = 1:intSizeY
                fprintf(fid, '%16.4f\t', dataInput(intI+(1-vecFirstIndex(1)), intJ+(1-vecFirstIndex(2)), intZ+(1-vecFirstIndex(3))));
            end
            fprintf(fid, '\n');
        end
    end
    fprintf(fid, ';\n');
end
