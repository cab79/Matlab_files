function [] = printList(fid, list, perLine, indent)
% Output the values in list a specified number of values per line
count = 0;
fprintf(fid, '%s[ ', indent);
for k = 1:length(list)
    fprintf(fid, '%g ', list(k));
    count = count + 1;
    if mod(count, perLine) == 0
        fprintf(fid, '\n%s', indent);
        count = 0;
    end
end
fprintf(']\n');