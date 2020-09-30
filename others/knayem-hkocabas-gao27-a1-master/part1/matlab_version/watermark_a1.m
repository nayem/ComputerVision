image = imread('1.jpg');
figure
imshow(image);

F = fft2(image);
R = abs(F);
I = angle(F);
figure
imshow(R);
figure
imshow(I);

%F2 = fftshift(F);
%imshow(log(abs(F)),[-1 5]);
% colormap(jet); colorbar


