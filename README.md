# How to Use Test Frame

First we obtain the codes using command ```wget```.
```
$ cd ~/
$ wget https://www.cs.ucr.edu/~yzhai015/lab1_frame.tar.gz
$ tar -xvf lab1_frame.tar.gz
$ cd ~/lab1
```

Then you should 
* Implement all functions in ```mygemm.c```.
* Modify optimal block size in ```line 39``` of ```cache_part4.c``` based on your own experiments in part3.
* If you managed to take care of the boundary condition, you can change the matrix size to {64, 128, 256, 512, 1024, 2048 and etc}. Some students may have difficulties in the boundary condition, that's OK. If you don't know how to deal with boundary condition issue (if you have), you may just use my test code.

Now you have your codes ready to run. Cong! Now you can run the code by
```
$ cd ~/lab1/build
$ export CXX=/act/gcc-4.7.2/bin/g++
$ export CC=/act/gcc-4.7.2/bin/gcc
$ cmake ..
$ make
```

Then you want to submit your code with ```SLURM```. You just need to do
```
$ cd ~/lab1
$ sbatch submit.sh
```

Finally your ouput can be found in ```~/lab1/data/```. You can check your output by type in
```
$ vim ~/lab1/data/mydata.txt
```