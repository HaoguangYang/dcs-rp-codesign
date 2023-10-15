function network_stream_write(conn, data, bufferSize)
    ptr = 1;
    data = typecast(data, 'uint8');
    dataSize = size(data, 1) * size(data, 2);
    while (dataSize > 0)
        count = min(dataSize, bufferSize);
        ptr1 = ptr + count;
        fwrite(conn, data(ptr : (ptr1 - 1)), 'uint8');
        ptr = ptr1;
        dataSize = dataSize - count;
    end
    count1 = fread(conn, 1, 'uint32');
    if (count1 ~= ptr - 1)
        ptr = ptr - 1;
        disp(ptr)
        disp(count1)
        error('Data Corrupted!');
    end
end
