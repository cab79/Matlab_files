function MNI=vox2mni(M,VOX)

%function MNI=vox2mni(M,VOX)

%Voxel-space to MNI-space

 T=M(1:3,4);

 M=M(1:3,1:3);

 MNI=M*VOX;

 for i=1:3

MNI(i,:)=MNI(i,:)+T(i);

 end

return