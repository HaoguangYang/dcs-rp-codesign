function data = network_stream_read(conn, dataSize, bufferSize, type)
    ptr = 1;
    data = uint8(zeros(dataSize, 1));
    while (dataSize > 0)
        count = min(int32(dataSize), int32(bufferSize));
        ptr1 = ptr + count;
        data(ptr : (ptr1 - 1)) = fread(conn, double(count), 'uint8');
        ptr = ptr1;
        dataSize = int32(dataSize) - int32(count);
    end
    fwrite(conn, uint32(ptr - 1), 'uint32');
    data = typecast(data, type);
end
