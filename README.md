# resolventchannel
The following matlab codes are used for the resolvent analysis of turbulent channel flow at Re_tau = 182. The scripts are run using chebfun suites; therefore, chebfun must be downloaded from https://www.chebfun.org before running any code These matlab scripts should be run to derive the following plots.


- final4.png The 20 largest singular values of the resolvent operator for $(\lambda_x, \lambda_z, c) = (1000,100,10)$
- final5.png Color contour plot of streamwise velocity fluctuations from the most energetic resolvent mode at $(\lambda^+_x,\lambda^+_z,c) = (1000 ,100, 10)$
- final7.png Principal velocity response for $(\lambda^+_x,\lambda^+_z,c) = (1000 ,1000, 10)$ colored with the streamwise fluctuating velocity at 20% of its maximum

--- 

The matlab codes in this repo are 

- Re180.txt The mean velocity profile of the turbulent channel flow is obtained from https://turbulence.oden.utexas.edu
- singval.m code to plot final4.png
- velfluc.m code to plot final5.png
- isovel.m code to plot final7.png
- svdfr.m the svd function 

The codes are heavily influenced by ones found in https://gokulhari.github.io/svdfr/index.html, and a more detailed account of the behind-the-scene mechanisms of the code can be found in the URL. The svdfr.m was developed by Gokul Hariharan, Binh Lieu and Mihailo R. Jovanovic to accurately find the SVD of a resolvent operator. 
