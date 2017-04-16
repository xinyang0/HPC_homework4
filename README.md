# HPC_homework4

2. Weak scaling study:  
Fixed lN=100 and max_iter=10000.
tasks p   total size N   timings
 p=1         N=100       0.138146
 p=4         N=200       0.162171
 p=16        N=400       0.263917
 p=64        N=800       0.618267
 p=256       N=1600      1.530780
 p=1024      N=3200      3.297554
 
 Strong scaling study:
 Fixed N=3200 and max_iter=10000.
 tasks p          timings
   p=1            180.602354
   p=4             68.047975
   p=16            33.142351
   p=64            4.942223
   p=256           2.457551
   p=1024          
   
 
 
3. Fixed processor p=64.
 the number of elements per processor     timings
              N=100                      0.005930
              N=200                      0.006920
              N=400                      0.018462 
              N=800                      0.019282
              N=1600
  The timing may not increases with N sometimes. I noticed that this may caused by the unevenly distributed data among the tasks. If one processor has too much data, sorting this data in one processor takes long time. This happens more when N is small.
