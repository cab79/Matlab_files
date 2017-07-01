% excluding trials with eye movement

files=dir('E*.mat');


segmentlength = 200:850;  % segment for excluding eye movement from t= 0 to t= 1700ms

for k=1:length(files)
    load(files(k).name);

    accept=ones(1,size(el,2));
    
    for i=1:size(el,2)
        
       corr = abs(min(min(corrcoef(squeeze(el(35,i,segmentlength)), squeeze(el(62,i,segmentlength))  ))));
       
       maxAF7= max(abs(squeeze(el(35,i,:))));  maxAF8= max(abs(squeeze(el(36,i,:))));

       if corr >= 0.75 | maxAF7>= 70   |  maxAF8>= 70
           accept(i)=0;
       end
    end
    
    save(files(k).name, 'accept', '-append')
end
