%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
This code was developed by Li Yuanman 
yuanmanx.li@gmail.com
19/3/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Please cite these papers if you use this code:
[1] Y. M. Li and J. T. Zhou, "Fast and Effective Image Copy-Move Forgery Detection via Hierarchical Feature Point Matching",  To appear in IEEE Trans. on Inf. Forensics and Security (T-IFS), 2018.
[2] Y. M. Li and J. T. Zhou, "Image Copy-Move Forgery Detection Using Hierarchical Feature Point Matching", ASC, 2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
1. We have recomplied the vl_feat libary for adaption, which is also included our source code.  
2. The result of each run could be slightly different due to the randomness involved. 
3. all the pictures shown in our (journal) paper can be repoduced by running the script 'runDemo.m'

We test the results based on the current code (windows7 && matlab2017a): 
              TPR          FPR       F1
F220:         100%,      2.72%       98.65
F600:         97.50%     5.68%       91.50
GRIP:         100%,      0%          100
FAU:          100%,      2.08%       98.97
CMH+GRIPori:  96.3%,     0%          98.11
COVERAGE:     80.22%     41.76%      72.28
