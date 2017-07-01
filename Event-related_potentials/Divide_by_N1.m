load amplitudes
AMP = reshape(AMP2,length(subjects),2,5);
for i = 1:size(AMP,3)
    AMP_n1(:,:,i) = AMP(:,:,i)./AMP(:,:,3);
end
AMP(:,:,size(AMP,3)+1:size(AMP,3)*2) = AMP_n1;
AMP2 = reshape(AMP,length(subjects),20);
save 'amplitudes.mat' AMP2;