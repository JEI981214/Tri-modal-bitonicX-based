% read in image
A = imread('01.jpg','jpg');

% Run filters and show results
for i = 1:3
  
  % Run appropriate processing
  switch (i)
    case 1
      title = 'Noisy image';
      B = A;
      
    case 2
      title = 'Bitonic X';
      fprintf(1,'Calculating Bitonic X ... ');
      tic;
      %mex xrankopen2.cpp
      mex anisotropic2_mex.cpp
      B = xbitonic2(A);
     
      t = toc;
      fprintf(1,'Completed in %.3f secs.\n', t);
      
    case 3
      title = 'Bitonic MX';
      fprintf(1,'Calculating Bitonic MX ... ');
      tic;
      B = mxbitonic2(A);
      t = toc;
      fprintf(1,'Completed in %.3f secs.\n', t);

  end
  
  % Display image with increased contrast to see residual errors
  figure(i);
  imagesc(3*double(B)/255-1); 
  colormap(gray);
  pos = get(i, 'Position');
  pos(3:4) = [512 512];
  set(i,'Position',pos);
  text(8,20,title,'Color',[1 1 0]);
  axis off;
  set(gca,'Position',[0 0 1 1]);
  set(gca,'Box','off');

end