function dtype = matlab2FileArrayType(dtype)
% FORMAT fatype = matlab2fileArrayType(mtype)
% mtype  - Matlab datatype (ex: 'single')
% fatype - file_array data type name (ex: 'FLOAT32') or code (ex: 16)
    switch dtype
        case 'logical'
            dtype = 'BINARY';
        case 'int8'
            dtype =  'INT8';
        case 'uint8'
            dtype =  'UINT8';
        case 'int16'
            dtype =  'INT16';
        case 'uint16'
            dtype =  'UINT16';
        case 'int32'
            dtype =  'INT32';
        case 'uint32'
            dtype =  'UINT32';
        case 'int64'
            dtype = 'INT64';
        case 'uint64'
            dtype =  'UINT64';
        case 'single'
            dtype =  'FLOAT32';
        case 'double'
            dtype =  'FLOAT64';
        otherwise
            dtype = 'UNKNOWN';
    end
end