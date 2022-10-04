# Tracking Software

This matlab software tracks the position of A431 clusters migrating on E-cadherin coated polyacrylamide gels. Clusters are imaged with phase contrast using a 10x objective. The software generates a file named `Tracks.mat`. Inside this file there is an matlab array. Each element of the array represents the position of a cluster. Inside each cluster there are 4 columns (Frame, x, y and Area). x, y and Area are given in pixels. 

Besides the tracking data, the software generates a folder named `Tracking` with the images of the tracked clusters labeled with a number and its segmentation. These images will be used to discard those clusters with segmentation errors. 

This software need the following libraries:

- [Gray image thresholding using the Triangle Method - File Exchange - MATLAB Central (mathworks.com)](https://es.mathworks.com/matlabcentral/fileexchange/28047-gray-image-thresholding-using-the-triangle-method)

- [simpletracker - File Exchange - MATLAB Central (mathworks.com)](https://es.mathworks.com/matlabcentral/fileexchange/34040-simpletracker)

- [export_fig - File Exchange - MATLAB Central (mathworks.com)](https://es.mathworks.com/matlabcentral/fileexchange/23629-export_fig)


