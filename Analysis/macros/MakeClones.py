import os

os.system("rm -fr  ../PlotterSteer/*.steer")
os.system("ln -s /nfs/dust/cms/user/"+os.environ["USER"]+"/WorkingArea/SFramePlotter/Analysis/* ../PlotterSteer/")
