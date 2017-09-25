function mymatr = returnnum(mystr)

j=0;
mystr=strtrim(mystr);
while j < size(mystr,2)
    j=j+1;
if isspace(mystr(j))&& isspace(mystr(j+1))
    temp1 = mystr(1:j);
    temp2 = mystr(j+2:size(mystr,2));
    mystr = cat (2,temp1, temp2);
j=0;     
end
end


mymatr= zeros(sum(isspace(mystr))+1,1);
to=0;
count=0;
frm=0;
for i = 1:size(mystr,2) 
    if ~isspace(mystr(i))
        frm=to+1;
        to=to+1;
        count= count+1;
        if to < size(mystr,2)
            while ~isspace(mystr(to)) && to<size(mystr,2)
                to = to+1;
            end
        else
            to = size(mystr,2);
        end
        if count <= sum(isspace(mystr))+1
        mymatr(count)= eval( mystr(frm:to));
        end
    end   
end

