rsync -av *.cu *.cuh Makefile ../src/*.c ../src/*.h hasindu@kepler:/storage/hasindu/f5c/ && ssh kepler 'cd /storage/hasindu/f5c/ && make cuda=1'
rsync -av scripts/*.sh hasindu@kepler:/storage/hasindu/f5c/scripts/ 

rsync -av *.cu *.cuh Makefile ../src/*.c ../src/*.h hasindu@jetson:~/f5c/ && ssh jetson 'cd ~/f5c/ && make cuda=1'
rsync -av scripts/*.sh hasindu@jetson:~/f5c/scripts/
